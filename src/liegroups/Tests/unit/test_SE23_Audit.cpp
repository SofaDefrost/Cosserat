// Audit of SE_2(3) against long double matrix references.
// Tangent layout (documented in SE23.h): xi = [v ; w ; a] (linear, angular, acceleration),
// element = [[R, p, u], [0, 1, 0], [0, 0, 1]] (extendedMatrix()).

#include <gtest/gtest.h>
#include <liegroups/SE23.h>

#include <random>

#include "LieGroupReference.h"

using namespace sofa::component::cosserat::liegroups;
using namespace liegroup_ref;

namespace {
using SE23d = SE23<double>;
using V9 = Eigen::Matrix<double, 9, 1>;
using M9 = Eigen::Matrix<double, 9, 9>;
const double kPi = 3.14159265358979323846;
const double kAngles[] = {1e-9, 1e-6, 1e-4, 2e-4, 1e-3, 1e-2, 0.1, 0.3, 1.0, 2.0, 3.0, kPi - 1e-6, kPi - 1e-9};

V9 randomXi(std::mt19937_64 &rng, double angle) {
	std::uniform_real_distribution<double> u(-1.0, 1.0);
	Eigen::Vector3d axis(u(rng), u(rng), u(rng));
	axis.normalize();
	V9 x;
	x.segment<3>(0) = Eigen::Vector3d(u(rng), u(rng), u(rng)) * 2;
	x.segment<3>(3) = axis * angle;
	x.segment<3>(6) = Eigen::Vector3d(u(rng), u(rng), u(rng)) * 2;
	return x;
}
MatL<5> hat5(const VecL<9> &x) {
	MatL<5> M = MatL<5>::Zero();
	M.block<3, 3>(0, 0) = hat3(x.segment<3>(3));
	M.block<3, 1>(0, 3) = x.segment<3>(0);
	M.block<3, 1>(0, 4) = x.segment<3>(6);
	return M;
}
MatL<9> ad9(const VecL<9> &x) {
	const MatL<3> W = hat3(x.segment<3>(3)), Vh = hat3(x.segment<3>(0)), Ah = hat3(x.segment<3>(6));
	MatL<9> A = MatL<9>::Zero();
	A.block<3, 3>(0, 0) = W;
	A.block<3, 3>(0, 3) = Vh;
	A.block<3, 3>(3, 3) = W;
	A.block<3, 3>(6, 3) = Ah;
	A.block<3, 3>(6, 6) = W;
	return A;
}
/// Group element from the long double exponential (quaternion from the exact axis-angle).
SE23d referenceElement(const V9 &xi) {
	const MatL<5> T = expm<5>(hat5(xi.cast<LD>()));
	const VecL<3> phi = xi.segment<3>(3).cast<LD>();
	const LD th = phi.norm();
	Eigen::Quaterniond q = Eigen::Quaterniond::Identity();
	if (th > 0) {
		const VecL<3> ax = phi / th;
		const LD s = std::sin(th / 2);
		q = Eigen::Quaterniond(static_cast<double>(std::cos(th / 2)), static_cast<double>(ax[0] * s),
							   static_cast<double>(ax[1] * s), static_cast<double>(ax[2] * s));
	}
	return SE23d(SO3d(q), T.block<3, 1>(0, 3).cast<double>(), T.block<3, 1>(0, 4).cast<double>());
}
} // namespace

TEST(SE23Audit, ExpMatchesMatrixExponential) {
	std::mt19937_64 rng(51);
	for (double angle : kAngles) {
		for (int i = 0; i < 3; ++i) {
			const V9 xi = randomXi(rng, angle);
			EXPECT_LT(maxAbsDiff(SE23d::exp(xi).extendedMatrix(), expm<5>(hat5(xi.cast<LD>()))), 1e-14) << "angle " << angle;
		}
	}
}

TEST(SE23Audit, LogInvertsExp) {
	std::mt19937_64 rng(52);
	for (double angle : kAngles) {
		for (int i = 0; i < 3; ++i) {
			const V9 xi = randomXi(rng, angle);
			EXPECT_LT((referenceElement(xi).log() - xi).cwiseAbs().maxCoeff(), 2e-14) << "angle " << angle;
		}
	}
}

TEST(SE23Audit, HatVeeAndAdAreTheLieAlgebra) {
	std::mt19937_64 rng(53);
	for (double angle : {0.0, 0.4, 2.0}) {
		const V9 x = randomXi(rng, angle), y = randomXi(rng, angle / 2 + 0.1);
		EXPECT_LT((SE23d::vee(SE23d::hat(x)) - x).norm(), 1e-15);
		const M9 comm = SE23d::hat(x) * SE23d::hat(y) - SE23d::hat(y) * SE23d::hat(x);
		EXPECT_LT((SE23d::computeAd(x) * y - SE23d::vee(comm)).cwiseAbs().maxCoeff(), 1e-14) << angle;
		EXPECT_LT(maxAbsDiff(SE23d::computeAd(x), ad9(x.cast<LD>())), 1e-15) << angle;
	}
}

TEST(SE23Audit, AdjointIsExpOfAd) {
	std::mt19937_64 rng(54);
	for (double angle : kAngles) {
		const V9 xi = randomXi(rng, angle);
		EXPECT_LT(maxAbsDiff(SE23d::exp(xi).computeAdjoint(), expm<9>(ad9(xi.cast<LD>()))), 8e-15) << "angle " << angle;
	}
}

TEST(SE23Audit, ComposeAndInverseFollowTheMatrixProduct) {
	std::mt19937_64 rng(55);
	const SE23d g = SE23d::exp(randomXi(rng, 1.1)), h = SE23d::exp(randomXi(rng, 2.0));
	EXPECT_LT((g.compose(h).extendedMatrix() - g.extendedMatrix() * h.extendedMatrix()).cwiseAbs().maxCoeff(), 1e-14);
	EXPECT_LT((g.computeInverse().extendedMatrix() - g.extendedMatrix().inverse()).cwiseAbs().maxCoeff(), 1e-13);
}

TEST(SE23Audit, JacobiansMatchSeriesAndBchConverges) {
	std::mt19937_64 rng(56);
	for (double angle : kAngles) {
		const V9 xi = randomXi(rng, angle);
		const MatL<9> ad = ad9(xi.cast<LD>());
		EXPECT_LT(maxAbsDiff(SE23d::leftJacobian(xi), jacobianSeries<9>(ad)), 4e-15) << "angle " << angle;
		EXPECT_LT(maxAbsDiff(SE23d::rightJacobian(xi), jacobianSeries<9>(-ad)), 4e-15) << "angle " << angle;
	}
	std::uniform_real_distribution<double> u(-1.0, 1.0);
	V9 v0, w0;
	for (int i = 0; i < 9; ++i) {
		v0[i] = u(rng);
		w0[i] = u(rng);
	}
	auto err = [&](double s, int order) {
		const V9 v = v0 * s, w = w0 * s;
		return (SE23d::BCH(v, w, order) - (SE23d::exp(v) * SE23d::exp(w)).log()).norm();
	};
	for (int order = 2; order <= 4; ++order) {
		const double ratio = err(0.04, order) / err(0.02, order);
		EXPECT_GT(ratio, 0.7 * std::pow(2.0, order + 1)) << order;
		EXPECT_LT(ratio, 1.5 * std::pow(2.0, order + 1)) << order;
	}
}
