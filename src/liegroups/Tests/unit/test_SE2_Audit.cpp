// Audit of SE2 (and the SO2 / RealSpace Jacobian hooks) against long double references.
// SE2 tangent layout: xi = [rho_x ; rho_y ; theta].

#include <gtest/gtest.h>
#include <liegroups/RealSpace.h>
#include <liegroups/SE2.h>
#include <liegroups/SO2.h>

#include <random>

#include "LieGroupReference.h"

using namespace sofa::component::cosserat::liegroups;
using namespace liegroup_ref;

namespace {
using V3 = Eigen::Vector3d;
using M3 = Eigen::Matrix3d;
const double kPi = 3.14159265358979323846;
// theta may be negative; |theta| < pi because SO2 stores an angle in [-pi, pi).
const double kThetas[] = {1e-12, 1e-9, 1e-7, 1e-5, 1e-4, 2e-4, 1e-3, 1e-2, 0.1, 0.3, 1.0, 2.0, 3.0,
						  3.14, kPi - 1e-6, kPi - 1e-9, -1e-9, -1e-3, -1.0, -3.0, -kPi + 1e-6};

V3 randomXi(std::mt19937_64 &rng, double theta) {
	std::uniform_real_distribution<double> u(-1.0, 1.0);
	return V3(2 * u(rng), 2 * u(rng), theta);
}
/// SE2 element built from the long double exponential of xi.
SE2d referenceElement(const V3 &xi) {
	const MatL<3> T = expm<3>(hatSE2(xi.cast<LD>()));
	return SE2d(xi[2], Eigen::Vector2d(static_cast<double>(T(0, 2)), static_cast<double>(T(1, 2))));
}
MatL<3> adSE2(const V3 &x) {
	// [x, y] = hat^-1(hat x hat y - hat y hat x): rho part theta_x J rho_y - theta_y J rho_x
	MatL<3> A = MatL<3>::Zero();
	const LD th = x[2];
	A(0, 1) = -th;
	A(1, 0) = th;
	A(0, 2) = x[1];
	A(1, 2) = -static_cast<LD>(x[0]);
	return A;
}
} // namespace

TEST(SE2Audit, ExpMatchesMatrixExponential) {
	std::mt19937_64 rng(41);
	for (double th : kThetas) {
		for (int i = 0; i < 5; ++i) {
			const V3 xi = randomXi(rng, th);
			EXPECT_LT(maxAbsDiff(SE2d::exp(xi).matrix(), expm<3>(hatSE2(xi.cast<LD>()))), 4e-15) << "theta " << th;
		}
	}
}

TEST(SE2Audit, LogInvertsExp) {
	std::mt19937_64 rng(42);
	for (double th : kThetas) {
		for (int i = 0; i < 5; ++i) {
			const V3 xi = randomXi(rng, th);
			EXPECT_LT((referenceElement(xi).log() - xi).cwiseAbs().maxCoeff(), 2e-14) << "theta " << th;
		}
	}
}

TEST(SE2Audit, AdIsTheLieBracketAndAdjointIsExpOfAd) {
	std::mt19937_64 rng(43);
	for (double th : {0.0, 0.3, -1.0, 2.5}) {
		const V3 x = randomXi(rng, th), y = randomXi(rng, th / 2 + 0.1);
		const M3 br = SE2d::computeHat(x) * SE2d::computeHat(y) - SE2d::computeHat(y) * SE2d::computeHat(x);
		EXPECT_LT((SE2d::computeAd(x) * y - SE2d::computeVee(br)).cwiseAbs().maxCoeff(), 1e-14) << th;
		EXPECT_LT((SE2d::ad(x) - SE2d::computeAd(x)).cwiseAbs().maxCoeff(), 1e-15) << th;
		// Ad_{exp(x)} = exp(ad_x)
		EXPECT_LT(maxAbsDiff(SE2d::exp(x).computeAdjoint(), expm<3>(adSE2(x))), 4e-15) << th;
	}
}

TEST(SE2Audit, MatrixAndTangentConstructorsCompileAndAgree) {
	std::mt19937_64 rng(44);
	const V3 xi = randomXi(rng, 0.8);
	const SE2d g = SE2d::exp(xi);
	const SE2d fromMatrix(g.matrix());
	const SE2d fromTangent(xi);
	EXPECT_TRUE(fromMatrix.isApprox(g, 1e-12));
	EXPECT_TRUE(fromTangent.isApprox(g, 1e-12));
}

