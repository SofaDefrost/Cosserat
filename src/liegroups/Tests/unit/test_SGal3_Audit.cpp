// Audit of SGal(3) against long double matrix references.
// Element = [[R, v, p], [0, 1, t], [0, 0, 1]] (extendedMatrix()),
// tangent xi = [rho ; w ; beta ; tau], hat(xi) = [[w^, beta, rho], [0, 0, tau], [0, 0, 0]].

#include <gtest/gtest.h>
#include <liegroups/SGal3.h>

#include <random>

#include "LieGroupReference.h"

using namespace sofa::component::cosserat::liegroups;
using namespace liegroup_ref;

namespace {
using V10 = Eigen::Matrix<double, 10, 1>;
using M10 = Eigen::Matrix<double, 10, 10>;
const double kPi = 3.14159265358979323846;
const double kAngles[] = {1e-9, 1e-6, 1e-4, 2e-4, 1e-3, 1e-2, 0.1, 0.3, 1.0, 2.0, 3.0, kPi - 1e-6, kPi - 1e-9};

V10 randomXi(std::mt19937_64 &rng, double angle) {
	std::uniform_real_distribution<double> u(-1.0, 1.0);
	Eigen::Vector3d axis(u(rng), u(rng), u(rng));
	axis.normalize();
	V10 x;
	x.segment<3>(0) = Eigen::Vector3d(u(rng), u(rng), u(rng)) * 2;
	x.segment<3>(3) = axis * angle;
	x.segment<3>(6) = Eigen::Vector3d(u(rng), u(rng), u(rng)) * 2;
	x[9] = 1.5 * u(rng);
	return x;
}
MatL<5> hat5(const VecL<10> &x) {
	MatL<5> M = MatL<5>::Zero();
	M.block<3, 3>(0, 0) = hat3(x.segment<3>(3));
	M.block<3, 1>(0, 3) = x.segment<3>(6);
	M.block<3, 1>(0, 4) = x.segment<3>(0);
	M(3, 4) = x[9];
	return M;
}
MatL<10> ad10(const VecL<10> &x) {
	const MatL<3> W = hat3(x.segment<3>(3));
	MatL<10> A = MatL<10>::Zero();
	A.block<3, 3>(0, 0) = W;
	A.block<3, 3>(0, 3) = hat3(x.segment<3>(0));
	A.block<3, 3>(0, 6) = -x[9] * MatL<3>::Identity();
	A.block<3, 1>(0, 9) = x.segment<3>(6);
	A.block<3, 3>(3, 3) = W;
	A.block<3, 3>(6, 3) = hat3(x.segment<3>(6));
	A.block<3, 3>(6, 6) = W;
	return A;
}
SGal3d referenceElement(const V10 &xi) {
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
	return SGal3d(SO3d(q), T.block<3, 1>(0, 4).cast<double>(), T.block<3, 1>(0, 3).cast<double>(),
				  static_cast<double>(T(3, 4)));
}
} // namespace

TEST(SGal3Audit, ClassCanBeInstantiated) {
	// Before: compose() did not exist and the CRTP hooks were private, so SGal3 could not be used at all.
	const SGal3d id = SGal3d::Identity();
	EXPECT_TRUE(id.isApprox(id * id));
	EXPECT_TRUE(id.inverse().isApprox(id));
	EXPECT_NEAR(id.log().norm(), 0.0, 1e-15);
}

TEST(SGal3Audit, ExpMatchesMatrixExponential) {
	std::mt19937_64 rng(61);
	for (double angle : kAngles) {
		for (int i = 0; i < 3; ++i) {
			const V10 xi = randomXi(rng, angle);
			EXPECT_LT(maxAbsDiff(SGal3d::exp(xi).extendedMatrix(), expm<5>(hat5(xi.cast<LD>()))), 1e-14) << "angle " << angle;
		}
	}
}

TEST(SGal3Audit, LogInvertsExp) {
	std::mt19937_64 rng(62);
	for (double angle : kAngles) {
		for (int i = 0; i < 3; ++i) {
			const V10 xi = randomXi(rng, angle);
			EXPECT_LT((referenceElement(xi).log() - xi).cwiseAbs().maxCoeff(), 5e-14) << "angle " << angle;
		}
	}
}

TEST(SGal3Audit, ComposeInverseAndActionFollowTheGroupLaw) {
	std::mt19937_64 rng(63);
	const SGal3d g = SGal3d::exp(randomXi(rng, 1.1)), h = SGal3d::exp(randomXi(rng, 2.0));
	EXPECT_LT((g.compose(h).extendedMatrix() - g.extendedMatrix() * h.extendedMatrix()).cwiseAbs().maxCoeff(), 1e-13);
	EXPECT_LT((g.inverse().extendedMatrix() - g.extendedMatrix().inverse()).cwiseAbs().maxCoeff(), 1e-12);
	EXPECT_TRUE((g * g.inverse()).isApprox(SGal3d::Identity(), 1e-12));
	SGal3d::ActionVector x;
	x << 0.3, -1.2, 0.8, 0.5, 0.1, -0.7, 0.9;
	EXPECT_LT((g.act(h.act(x)) - (g * h).act(x)).cwiseAbs().maxCoeff(), 1e-13);
}

TEST(SGal3Audit, HatVeeAdAndAdjointAreTheLieAlgebra) {
	std::mt19937_64 rng(64);
	for (double angle : {0.0, 0.4, 2.0}) {
		const V10 x = randomXi(rng, angle), y = randomXi(rng, angle / 2 + 0.1);
		EXPECT_LT((SGal3d::vee(SGal3d::hat(x)) - x).norm(), 1e-15);
		const M10 comm = SGal3d::hat(x) * SGal3d::hat(y) - SGal3d::hat(y) * SGal3d::hat(x);
		EXPECT_LT((SGal3d::computeAd(x) * y - SGal3d::vee(comm)).cwiseAbs().maxCoeff(), 1e-14) << angle;
		EXPECT_LT(maxAbsDiff(SGal3d::computeAd(x), ad10(x.cast<LD>())), 1e-15) << angle;
	}
	for (double angle : kAngles) {
		const V10 xi = randomXi(rng, angle);
		EXPECT_LT(maxAbsDiff(SGal3d::exp(xi).computeAdjoint(), expm<10>(ad10(xi.cast<LD>()))), 2e-14) << "angle " << angle;
	}
}

TEST(SGal3Audit, JacobiansMatchSeriesAndBchConverges) {
	std::mt19937_64 rng(65);
	for (double angle : kAngles) {
		const V10 xi = randomXi(rng, angle);
		const MatL<10> ad = ad10(xi.cast<LD>());
		EXPECT_LT(maxAbsDiff(SGal3d::leftJacobian(xi), jacobianSeries<10>(ad)), 1e-14) << "angle " << angle;
		EXPECT_LT(maxAbsDiff(SGal3d::rightJacobian(xi), jacobianSeries<10>(-ad)), 1e-14) << "angle " << angle;
	}
	std::uniform_real_distribution<double> u(-1.0, 1.0);
	V10 v0, w0;
	for (int i = 0; i < 10; ++i) {
		v0[i] = u(rng);
		w0[i] = u(rng);
	}
	auto err = [&](double s, int order) {
		const V10 v = v0 * s, w = w0 * s;
		return (SGal3d::BCH(v, w, order) - (SGal3d::exp(v) * SGal3d::exp(w)).log()).norm();
	};
	for (int order = 2; order <= 4; ++order) {
		const double ratio = err(0.04, order) / err(0.02, order);
		EXPECT_GT(ratio, 0.7 * std::pow(2.0, order + 1)) << order;
		EXPECT_LT(ratio, 1.5 * std::pow(2.0, order + 1)) << order;
	}
}
