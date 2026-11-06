// Numerical audit of SO3/SE3 exp, log and Jacobians against long double references
// built from the matrix definition (LieGroupReference.h), over angles from 1e-9 to pi - 1e-9.

#include <gtest/gtest.h>
#include <liegroups/SE3.h>
#include <liegroups/SO3.h>

#include <random>

#include "LieGroupReference.h"

using namespace sofa::component::cosserat::liegroups;
using namespace liegroup_ref;

namespace {
using V6 = Eigen::Matrix<double, 6, 1>;
using V3 = Eigen::Vector3d;

const double kPi = 3.14159265358979323846;
// Includes both sides of the 1e-4 small-angle switch and of the thresholds of the series.
const double kAngles[] = {1e-9, 1e-7, 1e-5, 5e-5, 1e-4, 2e-4, 1e-3, 1e-2, 0.1, 0.2, 0.3, 1.0, 2.0, 3.0,
						  3.1,  3.14, kPi - 1e-6, kPi - 1e-9};

/// SE3 element whose quaternion and translation come from the long double exponential, so that
/// log() is tested without the rounding of a rotation matrix.
SE3d referenceElement(const V6 &xi) {
	const VecL<6> xl = xi.cast<LD>();
	const MatL<4> T = expm<4>(hatSE3(xl));
	const VecL<3> phi = xl.head<3>();
	const LD th = phi.norm();
	Eigen::Quaterniond q = Eigen::Quaterniond::Identity();
	if (th > 0) {
		const VecL<3> ax = phi / th;
		const LD s = std::sin(th / 2);
		q = Eigen::Quaterniond(static_cast<double>(std::cos(th / 2)), static_cast<double>(ax[0] * s),
							   static_cast<double>(ax[1] * s), static_cast<double>(ax[2] * s));
	}
	return SE3d(SO3d(q), T.block<3, 1>(0, 3).cast<double>());
}
} // namespace

TEST(SE3Numerics, SO3JacobiansMatchSeries) {
	std::mt19937_64 rng(31);
	for (double angle : kAngles) {
		const V3 w = randomSE3Tangent(rng, angle).head<3>();
		const MatL<3> Jl = jacobianSeries<3>(hat3(w.cast<LD>()));
		const MatL<3> Jr = jacobianSeries<3>(-hat3(w.cast<LD>()));
		EXPECT_LT(maxAbsDiff(SO3d::computeLeftJacobian(w), Jl), 3e-15) << "Jl angle " << angle;
		EXPECT_LT(maxAbsDiff(SO3d::computeRightJacobian(w), Jr), 3e-15) << "Jr angle " << angle;
		EXPECT_LT(maxAbsDiff(SO3d::computeLeftJacobianInverse(w), Jl.inverse()), 3e-15) << "Jl^-1 angle " << angle;
		EXPECT_LT(maxAbsDiff(SO3d::computeRightJacobianInverse(w), Jr.inverse()), 3e-15) << "Jr^-1 angle " << angle;
	}
}

TEST(SE3Numerics, ExpMatchesMatrixExponential) {
	std::mt19937_64 rng(32);
	for (double angle : kAngles) {
		for (int i = 0; i < 5; ++i) {
			const V6 xi = randomSE3Tangent(rng, angle);
			const MatL<4> T = expm<4>(hatSE3(xi.cast<LD>()));
			EXPECT_LT(maxAbsDiff(SE3d::exp(xi).matrix(), T), 4e-15) << "angle " << angle;
		}
	}
}

TEST(SE3Numerics, LogInvertsExp) {
	std::mt19937_64 rng(33);
	for (double angle : kAngles) {
		for (int i = 0; i < 5; ++i) {
			const V6 xi = randomSE3Tangent(rng, angle);
			const V6 l = referenceElement(xi).log();
			EXPECT_LT((l - xi).cwiseAbs().maxCoeff(), 2e-14) << "angle " << angle;
		}
	}
}

TEST(SE3Numerics, ExpOfLogIsIdentityOnTheGroup) {
	std::mt19937_64 rng(34);
	for (double angle : kAngles) {
		const SE3d g = referenceElement(randomSE3Tangent(rng, angle));
		const SE3d back = SE3d::exp(g.log());
		EXPECT_LT((g.inverse() * back).log().cwiseAbs().maxCoeff(), 3e-14) << "angle " << angle;
	}
}

TEST(SE3Numerics, JacobiansMatchSeries) {
	std::mt19937_64 rng(35);
	for (double angle : kAngles) {
		for (int i = 0; i < 3; ++i) {
			const V6 xi = randomSE3Tangent(rng, angle);
			const MatL<6> ad = adSE3(xi.cast<LD>());
			const MatL<6> Jl = jacobianSeries<6>(ad), Jr = jacobianSeries<6>(-ad);
			EXPECT_LT(maxAbsDiff(SE3d::computeLeftJacobian(xi), Jl), 4e-15) << "Jl angle " << angle;
			EXPECT_LT(maxAbsDiff(SE3d::computeRightJacobian(xi), Jr), 4e-15) << "Jr angle " << angle;
			EXPECT_LT(maxAbsDiff(SE3d::computeLeftJacobianInverse(xi), Jl.inverse()), 4e-15) << "Jl^-1 angle " << angle;
			EXPECT_LT(maxAbsDiff(SE3d::computeRightJacobianInverse(xi), Jr.inverse()), 4e-15) << "Jr^-1 angle " << angle;
		}
	}
}

TEST(SE3Numerics, AdjointIsExpOfAd) {
	// Ad_{exp(xi)} = exp(ad_xi)
	std::mt19937_64 rng(36);
	for (double angle : kAngles) {
		const V6 xi = randomSE3Tangent(rng, angle);
		EXPECT_LT(maxAbsDiff(SE3d::exp(xi).computeAdjoint(), expm<6>(adSE3(xi.cast<LD>()))), 8e-15) << "angle " << angle;
	}
}

TEST(SE3Numerics, LeftAndRightJacobianAreRelatedByAdjoint) {
	// Jl(xi) = Ad_exp(xi) Jr(xi)  and  Jl(xi) = Jr(-xi)
	std::mt19937_64 rng(37);
	for (double angle : kAngles) {
		const V6 xi = randomSE3Tangent(rng, angle);
		const auto Jl = SE3d::computeLeftJacobian(xi);
		EXPECT_LT((Jl - SE3d::exp(xi).computeAdjoint() * SE3d::computeRightJacobian(xi)).cwiseAbs().maxCoeff(), 1e-13)
				<< "angle " << angle;
		EXPECT_LT((Jl - SE3d::computeRightJacobian(-xi)).cwiseAbs().maxCoeff(), 1e-15) << "angle " << angle;
	}
}
