// Tests of the Linear Strain Element kernel (liegroups/LinearStrainElement.h):
//  - the closed-form fourth-order Magnus end pose against a fine RK4 integration of g' = g ξ̂(s),
//  - the explicit mean-strain reconstruction (round trip),
//  - the Jacobians J1, J2, J3, the energy gradient and the Gauss-Newton tangent against central
//    finite differences taken with right (body) perturbations g <- g Exp(ζ),
//  - the constant-strain limit (β = 0) and the determinant of A.
//
// Reference: Xun & Tamadazte, Computers & Structures 325 (2026) 108181.

#include <gtest/gtest.h>
#include <liegroups/LinearStrainElement.h>

#include <cmath>
#include <random>

using namespace sofa::component::cosserat::liegroups;

namespace {
using LSE = LinearStrainElementd;
using Vec6 = LSE::Vec6;
using Mat6 = LSE::Mat6;
using Vec18 = LSE::Vec18;
using Mat18 = LSE::Mat18;
using Mat4 = Eigen::Matrix4d;

struct Config {
	SE3d ga;
	SE3d gb;
	Vec6 beta;
	Vec6 xiBar;
	double h;
};

Vec6 randomVec6(std::mt19937_64 &rng, double angularScale, double linearScale) {
	std::uniform_real_distribution<double> u(-1.0, 1.0);
	Vec6 v;
	for (int i = 0; i < 3; ++i) {
		v[i] = angularScale * u(rng);
		v[3 + i] = linearScale * u(rng);
	}
	return v;
}

/// A generic element: random node a, random mean strain around a unit-stretch rod, random slope.
Config randomConfig(std::uint64_t seed, double h = 0.3) {
	std::mt19937_64 rng(seed);
	Config c;
	c.h = h;
	c.ga = SE3d::computeExp(randomVec6(rng, 1.0, 1.0));
	c.xiBar = randomVec6(rng, 1.5, 0.1);
	c.xiBar[3] += 1.0; // straight rod: ε = (1, 0, 0)
	c.beta = randomVec6(rng, 2.0, 0.5);
	c.gb = LSE::endPose(c.ga, c.xiBar, c.beta, h);
	return c;
}

Mat6 diagK() {
	Vec6 d;
	d << 0.5, 0.3, 0.3, 100.0, 20.0, 20.0; // GJ, EIy, EIz, EA, GA, GA
	return d.asDiagonal();
}

/// Right perturbation g <- g Exp(ζ_i e_i * eps).
SE3d perturbed(const SE3d &g, int i, double eps) {
	Vec6 z = Vec6::Zero();
	z[i] = eps;
	return g.computeCompose(SE3d::computeExp(z));
}

/// End pose by RK4 on the matrix ODE Y' = Y Ξ(s), Ξ(s) = hat(ξ̄ + (s - h/2) β), Y(0) = I.
SE3d referenceEndPose(const Vec6 &xiBar, const Vec6 &beta, double h, int steps = 4000) {
	auto Xi = [&](double s) { return SE3d::computeHat(Vec6(xiBar + (s - 0.5 * h) * beta)); };
	Mat4 Y = Mat4::Identity();
	const double dt = h / steps;
	for (int k = 0; k < steps; ++k) {
		const double s = k * dt;
		const Mat4 X0 = Xi(s), Xm = Xi(s + 0.5 * dt), X1 = Xi(s + dt);
		const Mat4 k1 = Y * X0;
		const Mat4 k2 = (Y + 0.5 * dt * k1) * Xm;
		const Mat4 k3 = (Y + 0.5 * dt * k2) * Xm;
		const Mat4 k4 = (Y + dt * k3) * X1;
		Y += dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
	}
	return SE3d(Y);
}

double poseDistance(const SE3d &a, const SE3d &b) {
	return a.computeInverse().computeCompose(b).computeLog().norm();
}

/// Element energy as a function of the 18 right perturbations (ζ_a, ζ_b, β) around a configuration.
double energyAt(const Config &c, const Vec6 &xi0, const Mat6 &K, const Vec18 &dq) {
	SE3d ga = c.ga.computeCompose(SE3d::computeExp(Vec6(dq.segment<6>(0))));
	SE3d gb = c.gb.computeCompose(SE3d::computeExp(Vec6(dq.segment<6>(6))));
	const Vec6 beta = c.beta + dq.segment<6>(12);
	const Vec6 xiBar = LSE::meanStrain(ga, gb, beta, c.h);
	return LSE::energy(xiBar, xi0, beta, K, c.h);
}
} // namespace

// ----------------------------------------------------------------------------------------------
// Magnus end pose
// ----------------------------------------------------------------------------------------------

TEST(LinearStrainElement, MagnusMatchesReferenceIntegration) {
	// Strong strains (‖κ‖ up to 1.5, ‖β‖ up to 2): the truncation error is C h^5 with a large C, about
	// 6e-4 at h = 0.4 and 6e-7 at h = 0.1 (the order itself is checked in MagnusIsFourthOrder).
	const Config c = randomConfig(7);
	const double h = 0.1;
	const SE3d ref = referenceEndPose(c.xiBar, c.beta, h);
	const SE3d magnus = LSE::endPose(SE3d(), c.xiBar, c.beta, h);
	EXPECT_LT(poseDistance(ref, magnus), 1e-5);
}

TEST(LinearStrainElement, MagnusIsFourthOrder) {
	// The closed form is a 4th-order Magnus truncation: error ~ h^5 for a fixed strain field.
	const Config c = randomConfig(11);
	double previous = 0.0;
	for (double h : {0.4, 0.2, 0.1}) {
		const SE3d ref = referenceEndPose(c.xiBar, c.beta, h);
		const double err = poseDistance(ref, LSE::endPose(SE3d(), c.xiBar, c.beta, h));
		if (previous > 0.0) {
			EXPECT_GT(previous / err, 24.0) << "expected a ratio close to 2^5 = 32 at h = " << h;
		}
		previous = err;
	}
	EXPECT_LT(previous, 1e-6);
}

TEST(LinearStrainElement, MagnusSignOfSlopeTermMatchesRightInvariantKinematics) {
	// g' = g ξ̂ is right-invariant: Ω = hξ̄ − (h³/12) ad_β ξ̄. The opposite sign is only 2nd order,
	// so it is much farther from the reference; this pins the sign down.
	const Config c = randomConfig(13);
	const double h = 0.2;
	const SE3d ref = referenceEndPose(c.xiBar, c.beta, h);
	const Vec6 omegaGood = LSE::integratedTwist(c.xiBar, c.beta, h);
	const Vec6 omegaBad = h * c.xiBar + (h * h * h / 12.0) * (SE3d::computeAd(c.beta) * c.xiBar);
	const double errGood = poseDistance(ref, SE3d::computeExp(omegaGood));
	const double errBad = poseDistance(ref, SE3d::computeExp(omegaBad));
	EXPECT_GT(errBad / errGood, 50.0);
}

// ----------------------------------------------------------------------------------------------
// Explicit reconstruction
// ----------------------------------------------------------------------------------------------

TEST(LinearStrainElement, MeanStrainRoundTrip) {
	for (std::uint64_t seed = 1; seed <= 20; ++seed) {
		const Config c = randomConfig(seed);
		const Vec6 xiBar = LSE::meanStrain(c.ga, c.gb, c.beta, c.h);
		EXPECT_LT((xiBar - c.xiBar).norm(), 1e-11) << "seed " << seed;
	}
}

TEST(LinearStrainElement, EvaluateAgreesWithMeanStrain) {
	const Config c = randomConfig(5);
	const LSE::Kinematics k = LSE::evaluate(c.ga, c.gb, c.beta, c.h);
	EXPECT_LT((k.xiBar - c.xiBar).norm(), 1e-11);
	EXPECT_LT((k.A * k.xiBar - k.omega).norm(), 1e-11);
}

TEST(LinearStrainElement, ConstantStrainLimit) {
	// β = 0: A = h·I, ξ̄ = Log(g_a⁻¹ g_b)/h, J2 = Jr⁻¹/h, J1 = −Jr⁻¹ Ad_{Exp(−Ω)}/h, J3 = −(h²/12) ad_ξ̄.
	Config c = randomConfig(17);
	c.beta.setZero();
	c.gb = LSE::endPose(c.ga, c.xiBar, c.beta, c.h);
	const LSE::Kinematics k = LSE::evaluate(c.ga, c.gb, c.beta, c.h);

	const Vec6 omega = c.ga.computeInverse().computeCompose(c.gb).computeLog();
	EXPECT_LT((k.xiBar - omega / c.h).norm(), 1e-12);

	const Mat6 JrInv = SE3d::computeRightJacobianInverse(omega);
	EXPECT_LT((k.J2 - JrInv / c.h).cwiseAbs().maxCoeff(), 1e-12);
	EXPECT_LT((k.J1 + JrInv * SE3d::computeExp(-omega).computeAdjoint() / c.h).cwiseAbs().maxCoeff(), 1e-12);
	EXPECT_LT((k.J3 + c.h * c.h / 12.0 * SE3d::computeAd(k.xiBar)).cwiseAbs().maxCoeff(), 1e-12);
}

TEST(LinearStrainElement, DeterminantOfAIsPositive) {
	// ad_β has the eigenvalues of [κ_β]× (0 and ±i‖κ_β‖), so A is invertible for every β.
	for (std::uint64_t seed = 1; seed <= 10; ++seed) {
		const Config c = randomConfig(seed, 0.5);
		for (double scale : {1.0, 10.0, 100.0}) {
			const Vec6 beta = scale * c.beta;
			const double expected =
					std::pow(c.h, 6) * std::pow(1.0 + std::pow(c.h * c.h * beta.head<3>().norm() / 12.0, 2), 2);
			const double det = LSE::computeA(beta, c.h).determinant();
			EXPECT_NEAR(det / expected, 1.0, 1e-9) << "seed " << seed << " scale " << scale;
		}
	}
}

// ----------------------------------------------------------------------------------------------
// Jacobians, energy gradient and tangent against finite differences
// ----------------------------------------------------------------------------------------------

TEST(LinearStrainElement, JacobiansMatchFiniteDifferences) {
	const double eps = 1e-5;
	for (std::uint64_t seed = 1; seed <= 10; ++seed) {
		const Config c = randomConfig(seed);
		const LSE::Kinematics k = LSE::evaluate(c.ga, c.gb, c.beta, c.h);
		for (int i = 0; i < 6; ++i) {
			Vec6 db = Vec6::Zero();
			db[i] = eps;
			const Vec6 colA = (LSE::meanStrain(perturbed(c.ga, i, eps), c.gb, c.beta, c.h) -
							   LSE::meanStrain(perturbed(c.ga, i, -eps), c.gb, c.beta, c.h)) / (2 * eps);
			const Vec6 colB = (LSE::meanStrain(c.ga, perturbed(c.gb, i, eps), c.beta, c.h) -
							   LSE::meanStrain(c.ga, perturbed(c.gb, i, -eps), c.beta, c.h)) / (2 * eps);
			const Vec6 colBeta = (LSE::meanStrain(c.ga, c.gb, c.beta + db, c.h) -
								  LSE::meanStrain(c.ga, c.gb, c.beta - db, c.h)) / (2 * eps);
			const double scale = 1.0 + k.J1.cwiseAbs().maxCoeff();
			EXPECT_LT((colA - k.J1.col(i)).norm(), 1e-7 * scale) << "J1 col " << i << " seed " << seed;
			EXPECT_LT((colB - k.J2.col(i)).norm(), 1e-7 * scale) << "J2 col " << i << " seed " << seed;
			EXPECT_LT((colBeta - k.J3.col(i)).norm(), 1e-7 * scale) << "J3 col " << i << " seed " << seed;
		}
	}
}

TEST(LinearStrainElement, ResidualIsEnergyGradient) {
	const double eps = 1e-6;
	const Mat6 K = diagK();
	for (std::uint64_t seed = 1; seed <= 10; ++seed) {
		const Config c = randomConfig(seed);
		// A rest strain different from the current strain so that the residual is not trivial.
		std::mt19937_64 rng(seed + 100);
		const Vec6 xi0 = c.xiBar + randomVec6(rng, 0.3, 0.05);
		const LSE::Kinematics k = LSE::evaluate(c.ga, c.gb, c.beta, c.h);
		const Vec18 r = LSE::residual(k, xi0, c.beta, K, c.h).stacked();
		for (int i = 0; i < 18; ++i) {
			Vec18 dq = Vec18::Zero();
			dq[i] = eps;
			const double fd = (energyAt(c, xi0, K, dq) - energyAt(c, xi0, K, -dq)) / (2 * eps);
			EXPECT_NEAR(fd, r[i], 1e-6 * (1.0 + std::abs(r[i]))) << "dof " << i << " seed " << seed;
		}
	}
}

TEST(LinearStrainElement, GaussNewtonTangentIsExactHessianAtRestStrain) {
	// When ξ̄ = ξ₀ the residual term of the Hessian vanishes: the Gauss-Newton tangent is exact.
	const double eps = 1e-4;
	const Mat6 K = diagK();
	for (std::uint64_t seed = 1; seed <= 5; ++seed) {
		const Config c = randomConfig(seed);
		const LSE::Kinematics k = LSE::evaluate(c.ga, c.gb, c.beta, c.h);
		const Mat18 T = LSE::gaussNewtonTangent(k, K, c.h);
		EXPECT_LT((T - T.transpose()).cwiseAbs().maxCoeff(), 1e-9 * (1.0 + T.cwiseAbs().maxCoeff()));

		Mat18 H;
		for (int i = 0; i < 18; ++i) {
			for (int j = 0; j < 18; ++j) {
				Vec18 ei = Vec18::Zero(), ej = Vec18::Zero();
				ei[i] = eps;
				ej[j] = eps;
				H(i, j) = (energyAt(c, k.xiBar, K, ei + ej) - energyAt(c, k.xiBar, K, ei - ej) -
						   energyAt(c, k.xiBar, K, -ei + ej) + energyAt(c, k.xiBar, K, -ei - ej)) /
						  (4 * eps * eps);
			}
		}
		const double scale = 1.0 + T.cwiseAbs().maxCoeff();
		EXPECT_LT((H - T).cwiseAbs().maxCoeff(), 5e-5 * scale) << "seed " << seed;
	}
}

TEST(LinearStrainElement, ElasticEnergyOfSlopeOnly) {
	// ξ̄ = ξ₀: only the slope term (h³/24) βᵀKβ remains, Eq. (21).
	const Mat6 K = diagK();
	const Config c = randomConfig(3);
	const double e = LSE::energy(c.xiBar, c.xiBar, c.beta, K, c.h);
	EXPECT_NEAR(e, std::pow(c.h, 3) / 24.0 * c.beta.dot(K * c.beta), 1e-14);
}
