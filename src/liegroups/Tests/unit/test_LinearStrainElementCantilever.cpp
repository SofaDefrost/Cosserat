// Static rods and rod networks solved with Linear Strain Elements (liegroups/LinearStrainElement.h).
//
// A small Gauss-Newton solver (retraction g <- g Exp(δζ), tangent of Eq. 30) built on the kernel only,
// used to check the claims of Xun & Tamadazte, Computers & Structures 325 (2026) 108181, §6-7:
//  - a 4-element cantilever reproduces the large-deflection tip position of a fine mesh within 1 %,
//  - the strain-field energy error converges with order 4 for the linear strain element (LSE)
//    and 2 for the constant strain element (CSE, β = 0),
//  - a closed loop (thin ring under diametral compression) matches the analytic small-deflection solution.
// The same solver provides the reference values of the SOFA scenes of the mapping.
//
// Elements are node pairs, so the solver handles any rod graph; a chain is the special case of the cantilever.

#include <gtest/gtest.h>
#include <liegroups/LinearStrainElement.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <utility>
#include <vector>

using namespace sofa::component::cosserat::liegroups;

namespace {
using LSE = LinearStrainElementd;
using Vec6 = LSE::Vec6;
using Mat6 = LSE::Mat6;
using Vec3 = Eigen::Vector3d;
using VecX = Eigen::VectorXd;
using MatX = Eigen::MatrixXd;

/// Table 1 of the paper: L = 1 m, bending stiffness 0.2 N·m², shear and axial stiffness 1000 N.
Mat6 paperStiffness() {
	Vec6 k;
	k << 0.2, 0.2, 0.2, 1000.0, 1000.0, 1000.0; // GJ (unused in plane), EIy, EIz, EA, GA, GA
	return k.asDiagonal();
}

/// Strain of a straight unstretched rod: ε = (1, 0, 0).
Vec6 straightRestStrain() {
	Vec6 xi0 = Vec6::Zero();
	xi0[3] = 1.0;
	return xi0;
}

struct Solution {
	std::vector<SE3d> g;    ///< nodal poses
	std::vector<Vec6> beta; ///< slopes (zero for the CSE)
	bool converged = false;
	double residual = 0.0;
	double h = 0.0; ///< element length of a uniform chain (used by energyError)
};

/// A rod graph: initial nodal poses, node pairs, element lengths, clamped nodes and dead tip loads.
struct Network {
	std::vector<SE3d> g0;
	std::vector<std::array<int, 2>> elements;
	std::vector<double> lengths;
	std::vector<int> fixedNodes;
	std::vector<std::pair<int, Vec3>> loads; ///< world-frame force on a node
};

/// Gauss-Newton solve of the static equilibrium of a rod network. The rest strain of every element is the
/// mean strain of its initial poses with a zero slope (Remark 4 of the paper). With `withSlopes` false the
/// slopes are frozen to 0: constant strain element.
Solution solveNetwork(const Network &net, const Mat6 &K, bool withSlopes, int loadSteps = 10) {
	const int nNodes = static_cast<int>(net.g0.size());
	const int nElem = static_cast<int>(net.elements.size());

	std::vector<int> dofOfNode(nNodes, -1);
	int nFree = 0;
	for (int i = 0; i < nNodes; ++i)
		if (std::find(net.fixedNodes.begin(), net.fixedNodes.end(), i) == net.fixedNodes.end())
			dofOfNode[i] = 6 * nFree++;
	const int nPose = 6 * nFree;
	const int nU = nPose + (withSlopes ? 6 * nElem : 0);

	Solution s;
	s.g = net.g0;
	s.beta.assign(nElem, Vec6::Zero());
	s.h = nElem > 0 ? net.lengths[0] : 0.0;

	std::vector<Vec6> xi0(nElem);
	for (int e = 0; e < nElem; ++e)
		xi0[e] = LSE::meanStrain(net.g0[net.elements[e][0]], net.g0[net.elements[e][1]], Vec6::Zero(), net.lengths[e]);

	auto assemble = [&](double loadScale, VecX &r, MatX *T) {
		r = VecX::Zero(nU);
		if (T)
			*T = MatX::Zero(nU, nU);
		for (int e = 0; e < nElem; ++e) {
			const int a = net.elements[e][0], b = net.elements[e][1];
			const double h = net.lengths[e];
			const Vec6 beta = withSlopes ? s.beta[e] : Vec6::Zero();
			const LSE::Kinematics kin = LSE::evaluate(s.g[a], s.g[b], beta, h);
			const LSE::Residual re = LSE::residual(kin, xi0[e], beta, K, h);
			const int start[3] = {dofOfNode[a], dofOfNode[b], withSlopes ? nPose + 6 * e : -1};
			const Vec6 *part[3] = {&re.R1, &re.R2, &re.R3};
			for (int i = 0; i < 3; ++i)
				if (start[i] >= 0)
					r.segment<6>(start[i]) += *part[i];
			if (T) {
				const LSE::Mat18 Tg = LSE::gaussNewtonTangent(kin, K, h);
				for (int i = 0; i < 3; ++i)
					for (int j = 0; j < 3; ++j)
						if (start[i] >= 0 && start[j] >= 0)
							T->block<6, 6>(start[i], start[j]) += Tg.block<6, 6>(6 * i, 6 * j);
			}
		}
		// dead loads, expressed in the body frame of the loaded node (Eq. 27-28)
		for (const auto &[node, F] : net.loads) {
			if (dofOfNode[node] < 0)
				continue;
			Vec6 Fb = Vec6::Zero();
			Fb.tail<3>() = s.g[node].rotation().matrix().transpose() * (loadScale * F);
			r.segment<6>(dofOfNode[node]) -= Fb;
		}
	};

	auto apply = [&](const VecX &du, double alpha) {
		for (int i = 0; i < nNodes; ++i)
			if (dofOfNode[i] >= 0)
				s.g[i] = s.g[i].computeCompose(SE3d::computeExp(Vec6(alpha * du.segment<6>(dofOfNode[i]))));
		if (withSlopes)
			for (int e = 0; e < nElem; ++e)
				s.beta[e] += alpha * du.segment<6>(nPose + 6 * e);
	};

	VecX r, rTrial;
	MatX T;
	for (int step = 1; step <= loadSteps; ++step) {
		const double scale = double(step) / loadSteps;
		for (int it = 0; it < 300; ++it) {
			assemble(scale, r, &T);
			s.residual = r.norm();
			if (s.residual < 1e-10)
				break;
			const VecX du = T.ldlt().solve(-r);
			// backtracking on the residual norm
			double alpha = 1.0;
			const std::vector<SE3d> gSave = s.g;
			const std::vector<Vec6> bSave = s.beta;
			for (int k = 0; k < 12; ++k) {
				apply(du, alpha);
				assemble(scale, rTrial, nullptr);
				if (rTrial.norm() < s.residual)
					break;
				s.g = gSave;
				s.beta = bSave;
				alpha *= 0.5;
			}
		}
		assemble(scale, r, nullptr);
		s.residual = r.norm();
	}
	s.converged = s.residual < 1e-8;
	return s;
}

/// Clamped-free chain of N elements loaded by a dead (world-frame) force `fWorld` at the tip.
Solution solveCantilever(int N, double L, const Mat6 &K, const Vec3 &fWorld, bool withSlopes, int loadSteps = 10) {
	const double h = L / N;
	Network net;
	for (int i = 0; i <= N; ++i)
		net.g0.push_back(SE3d(SO3d(), Vec3(i * h, 0.0, 0.0)));
	for (int e = 0; e < N; ++e) {
		net.elements.push_back({e, e + 1});
		net.lengths.push_back(h);
	}
	net.fixedNodes = {0};
	net.loads = {{N, fWorld}};
	return solveNetwork(net, K, withSlopes, loadSteps);
}

/// Closed ring of radius R discretised with N elements (N multiple of 4), in the x-y plane. Node k sits at angle
/// -π/2 + 2πk/N; its frame has x along the counter-clockwise tangent, y towards the centre, z out of the plane.
/// Node 0 (bottom) is clamped and the force `-P e_y` acts on the top node N/2.
Network makeRing(int N, double R, double P) {
	const double pi = 3.14159265358979323846;
	Network net;
	for (int k = 0; k < N; ++k) {
		const double th = -0.5 * pi + 2.0 * pi * k / N;
		Eigen::Matrix3d Rot;
		Rot.col(0) = Vec3(-std::sin(th), std::cos(th), 0.0);
		Rot.col(1) = Vec3(-std::cos(th), -std::sin(th), 0.0);
		Rot.col(2) = Vec3(0.0, 0.0, 1.0);
		net.g0.push_back(SE3d(SO3d(Rot), Vec3(R * std::cos(th), R * std::sin(th), 0.0)));
	}
	for (int k = 0; k < N; ++k) {
		net.elements.push_back({k, (k + 1) % N}); // the last element closes the loop
		net.lengths.push_back(2.0 * pi * R / N);
	}
	net.fixedNodes = {0};
	net.loads = {{N / 2, Vec3(0.0, -P, 0.0)}};
	return net;
}

/// Strain energy error e_E = ∫ (ξ − ξ_r)ᵀ K (ξ − ξ_r) ds between a rod and a finer reference (Fig. 6).
double energyError(const Solution &coarse, const Solution &ref, bool coarseHasSlopes, const Mat6 &K) {
	const int Nc = static_cast<int>(coarse.beta.size());
	const int Nr = static_cast<int>(ref.beta.size());
	const double hc = coarse.h, hr = ref.h;
	// 4-point Gauss-Legendre on [0, 1]
	const double gp[4] = {0.0694318442029737, 0.3300094782075719, 0.6699905217924281, 0.9305681557970263};
	const double gw[4] = {0.1739274225687269, 0.3260725774312731, 0.3260725774312731, 0.1739274225687269};

	auto field = [&](const Solution &sol, bool slopes, int e, double sLocal) {
		const Vec6 beta = slopes ? sol.beta[e] : Vec6::Zero();
		const Vec6 xiBar = LSE::meanStrain(sol.g[e], sol.g[e + 1], beta, sol.h);
		return Vec6(xiBar + (sLocal - 0.5 * sol.h) * beta);
	};

	double err = 0.0;
	for (int e = 0; e < Nc; ++e) {
		for (int q = 0; q < 4; ++q) {
			const double s = (e + gp[q]) * hc;
			const int er = std::min(Nr - 1, static_cast<int>(std::floor(s / hr)));
			const Vec6 d = field(coarse, coarseHasSlopes, e, gp[q] * hc) - field(ref, true, er, s - er * hr);
			err += gw[q] * hc * d.dot(K * d);
		}
	}
	return err;
}

double slope(double e1, double e2, double n1, double n2) {
	return (std::log(e2) - std::log(e1)) / (std::log(n2) - std::log(n1));
}
} // namespace

TEST(LinearStrainElementCantilever, FourLinearElementsMatchFineMeshWithinOnePercent) {
	// Planar cantilever, tip force along z, F = 0.25 ... 3 N (paper §6.1, Fig. 3): relative tip error < 1 %.
	// Measured here with a 32-element reference: LSE 2e-5 ... 7e-4, CSE 0.5 % ... 1.7 %.
	const Mat6 K = paperStiffness();
	for (double F : {0.25, 0.5, 1.0, 2.0, 3.0}) {
		const Vec3 f(0.0, 0.0, F);
		const Solution ref = solveCantilever(32, 1.0, K, f, true, 4);
		const Solution lse = solveCantilever(4, 1.0, K, f, true);
		const Solution cse = solveCantilever(4, 1.0, K, f, false);
		ASSERT_TRUE(ref.converged) << "F=" << F << " residual " << ref.residual;
		ASSERT_TRUE(lse.converged) << "F=" << F << " residual " << lse.residual;
		ASSERT_TRUE(cse.converged) << "F=" << F << " residual " << cse.residual;

		const Vec3 tipRef = ref.g.back().translation();
		const double errLse = (lse.g.back().translation() - tipRef).norm() / tipRef.norm();
		const double errCse = (cse.g.back().translation() - tipRef).norm() / tipRef.norm();
		EXPECT_LT(errLse, 0.01) << "F=" << F;
		EXPECT_LT(errLse, errCse) << "LSE must be more accurate than CSE for the same number of elements, F=" << F;
		EXPECT_NEAR(tipRef.y(), 0.0, 1e-9) << "planar problem: no out-of-plane motion";
	}
}

TEST(LinearStrainElementCantilever, EnergyErrorConvergenceOrders) {
	// Strain-field error versus number of elements: order 4 (LSE) against order 2 (CSE), Fig. 6.
	const Mat6 K = paperStiffness();
	const Vec3 f(0.0, 0.0, 1.0);
	const Solution ref = solveCantilever(64, 1.0, K, f, true, 4);
	ASSERT_TRUE(ref.converged) << ref.residual;

	const int counts[] = {2, 4, 8};
	double eLse[3], eCse[3];
	for (int i = 0; i < 3; ++i) {
		const Solution lse = solveCantilever(counts[i], 1.0, K, f, true);
		const Solution cse = solveCantilever(counts[i], 1.0, K, f, false);
		ASSERT_TRUE(lse.converged) << "LSE N=" << counts[i] << " residual " << lse.residual;
		ASSERT_TRUE(cse.converged) << "CSE N=" << counts[i] << " residual " << cse.residual;
		eLse[i] = energyError(lse, ref, true, K);
		eCse[i] = energyError(cse, ref, false, K);
	}
	// Measured local orders: LSE 4.7, 4.2, 4.07 and CSE 1.98, 2.00, 2.00 for N = 2->4, 4->8, 8->16.
	const double pLse = -slope(eLse[0], eLse[2], counts[0], counts[2]);
	const double pCse = -slope(eCse[0], eCse[2], counts[0], counts[2]);
	EXPECT_GT(pLse, 3.5) << "LSE order " << pLse << " (expected ~4)";
	EXPECT_NEAR(pCse, 2.0, 0.4) << "CSE order " << pCse << " (expected ~2)";
	for (int i = 0; i < 3; ++i)
		EXPECT_LT(eLse[i], eCse[i]) << "N=" << counts[i];
}

TEST(LinearStrainElementCantilever, PureBendingGivesConstantMoment) {
	// A pure end moment M about y bends the rod into a circle: ξ_y = M/EIy for every element, with zero slope.
	const Mat6 K = paperStiffness();
	const int N = 4;
	const double L = 1.0, h = L / N, M = 0.3; // M/EI = 1.5 rad/m
	Vec6 xi;
	xi << 0.0, M / 0.2, 0.0, 1.0, 0.0, 0.0;
	SE3d g = SE3d();
	std::vector<SE3d> nodes{g};
	for (int e = 0; e < N; ++e) {
		g = LSE::endPose(g, xi, Vec6::Zero(), h);
		nodes.push_back(g);
	}
	// every element sees the same strain and a vanishing slope: the residual in the interior nodes is null
	for (int e = 0; e < N; ++e) {
		const LSE::Kinematics kin = LSE::evaluate(nodes[e], nodes[e + 1], Vec6::Zero(), h);
		EXPECT_LT((kin.xiBar - xi).norm(), 1e-12) << "element " << e;
	}
	const Vec6 r0 = LSE::residual(LSE::evaluate(nodes[0], nodes[1], Vec6::Zero(), h), straightRestStrain(), Vec6::Zero(), K, h).R2;
	const Vec6 r1 = LSE::residual(LSE::evaluate(nodes[1], nodes[2], Vec6::Zero(), h), straightRestStrain(), Vec6::Zero(), K, h).R1;
	// the force on the shared node cancels between the two elements: R2(e) + R1(e+1) = 0 (moment equilibrium)
	EXPECT_LT((r0 + r1).norm(), 1e-10);
}

TEST(LinearStrainElementCantilever, RingUnderDiametralCompressionMatchesAnalyticSolution) {
	// Closed loop (the last element joins node N-1 to node 0), thin ring, two opposite forces P (clamped bottom, load on
	// top). Small-deflection theory (Timoshenko): the vertical diameter shrinks by (π/4 − 2/π) P R³/EI = 0.1488 P R³/EI
	// and the horizontal one grows by (2/π − 1/2) P R³/EI = 0.1366 P R³/EI. P R²/EI = 0.005 keeps the problem linear
	// to ~0.1 %; the ring is thin (EI/(EA R²) = 2e-4) so shear and extension are negligible.
	const double pi = 3.14159265358979323846;
	const double R = 1.0, EI = 0.2, P = 0.001;
	const Mat6 K = paperStiffness();
	const double dVertical = (pi / 4.0 - 2.0 / pi) * P * R * R * R / EI;
	const double dHorizontal = (2.0 / pi - 0.5) * P * R * R * R / EI;

	for (int N : {16, 32}) {
		const Network net = makeRing(N, R, P);
		const Solution s = solveNetwork(net, K, true, 2);
		ASSERT_TRUE(s.converged) << "N=" << N << " residual " << s.residual;

		const double dv = net.g0[N / 2].translation().y() - s.g[N / 2].translation().y();
		const double dh = (s.g[N / 4].translation().x() - s.g[3 * N / 4].translation().x()) -
						  (net.g0[N / 4].translation().x() - net.g0[3 * N / 4].translation().x());
		EXPECT_NEAR(dv / dVertical, 1.0, 0.005) << "vertical diameter, N=" << N;
		EXPECT_NEAR(dh / dHorizontal, 1.0, 0.005) << "horizontal diameter, N=" << N;
		// the loop stays planar and symmetric about the vertical axis
		for (int k = 0; k < N; ++k)
			EXPECT_NEAR(s.g[k].translation().z(), 0.0, 1e-9) << "node " << k;
		EXPECT_NEAR(s.g[N / 4].translation().y(), s.g[3 * N / 4].translation().y(), 1e-9);
	}
}

TEST(LinearStrainElementCantilever, RingLinearElementsBeatConstantElements) {
	// Same coarse ring, with and without the slope: the linear strain element is the more accurate one.
	const double pi = 3.14159265358979323846;
	const double R = 1.0, EI = 0.2, P = 0.001;
	const Mat6 K = paperStiffness();
	const double dVertical = (pi / 4.0 - 2.0 / pi) * P * R * R * R / EI;
	const Network net = makeRing(8, R, P);
	const Solution lse = solveNetwork(net, K, true, 2);
	const Solution cse = solveNetwork(net, K, false, 2);
	ASSERT_TRUE(lse.converged && cse.converged) << lse.residual << " " << cse.residual;
	const double errLse = std::abs((net.g0[4].translation().y() - lse.g[4].translation().y()) / dVertical - 1.0);
	const double errCse = std::abs((net.g0[4].translation().y() - cse.g[4].translation().y()) / dVertical - 1.0);
	EXPECT_LT(errLse, errCse) << "LSE " << errLse << " vs CSE " << errCse;
}
