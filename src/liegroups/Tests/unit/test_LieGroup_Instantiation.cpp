// A class template member only fails to compile when it is instantiated, so
// -fsyntax-only on the headers misses members that cannot work. This file calls the
// public members that LieGroupBase promises to every group, for each group, and checks
// the result is meaningful (not only that it compiles).

#include <gtest/gtest.h>
#include <liegroups/SE23.h>
#include <liegroups/SE3.h>
#include <liegroups/SGal3.h>
#include <liegroups/SO3.h>

#include <random>

#include "LieGroupReference.h"

using namespace sofa::component::cosserat::liegroups;
using namespace liegroup_ref;

namespace {
using V6 = Eigen::Matrix<double, 6, 1>;
} // namespace

TEST(Instantiation, SE3HatVeeAndBracket) {
	// hat/vee are the CRTP hooks of LieGroupBase::hat/vee/ad (generic)/BCH/matrixLog.
	std::mt19937_64 rng(21);
	const V6 x = randomSE3Tangent(rng, 0.9), y = randomSE3Tangent(rng, 0.4);
	EXPECT_LT((SE3d::vee(SE3d::hat(x)) - x).norm(), 1e-15);
	EXPECT_LT((SE3d::hat(x) - hatSE3(x.cast<LD>()).cast<double>()).cwiseAbs().maxCoeff(), 1e-15);
	const Eigen::Matrix4d comm = SE3d::hat(x) * SE3d::hat(y) - SE3d::hat(y) * SE3d::hat(x);
	EXPECT_LT((SE3d::vee(comm) - SE3d::lieBracket(x, y)).norm(), 1e-14);
	const SE3d g = SE3d::exp(x);
	EXPECT_LT((g.matrixLog() - SE3d::hat(g.log())).cwiseAbs().maxCoeff(), 1e-15);
}

TEST(Instantiation, SE3BaseDistanceSquaredDistanceAndNearIdentity) {
	std::mt19937_64 rng(22);
	const SE3d g = SE3d::exp(randomSE3Tangent(rng, 0.9)), h = SE3d::exp(randomSE3Tangent(rng, 0.4));
	// The generic distance is the norm of log(g^-1 h); SE3::distance (weighted sum) hides it.
	using Base = LieGroupBase<SE3d, double, 4, 6, 3>;
	EXPECT_NEAR(g.Base::distance(h), (g.inverse() * h).log().norm(), 1e-12);
	EXPECT_NEAR(g.Base::squaredDistance(h), (g.inverse() * h).log().squaredNorm(), 1e-12);
	// A group element within rounding of the identity (translation 1e-17, not exactly 0) is near identity.
	const SE3d tiny = SE3d::exp(randomSE3Tangent(rng, 1e-12, 1e-12));
	EXPECT_TRUE(tiny.isNearIdentity(1e-9));
	EXPECT_TRUE(SE3d::Identity().isApprox(tiny, 1e-9));
	EXPECT_FALSE(SE3d::Identity().isApprox(g, 1e-9));
}

// ---- the free functions of SE3.inl / SE23.inl / SGal3.inl (these files used to be included by nobody,
// and parts of them did not compile) ----

TEST(Instantiation, SE3InlFreeFunctions) {
	std::mt19937_64 rng(23);
	const double roll = 0.3, pitch = -0.5, yaw = 1.2;
	const Eigen::Vector3d pos(1.0, -2.0, 0.5);
	const SE3d g = fromPositionEulerZYX<double>(pos, roll, pitch, yaw);
	// ZYX convention: R = Rz(yaw) Ry(pitch) Rx(roll)
	const Eigen::Matrix3d R = (Eigen::AngleAxisd(yaw, Eigen::Vector3d::UnitZ()) *
							   Eigen::AngleAxisd(pitch, Eigen::Vector3d::UnitY()) *
							   Eigen::AngleAxisd(roll, Eigen::Vector3d::UnitX()))
									  .toRotationMatrix();
	EXPECT_LT((g.rotation().matrix() - R).cwiseAbs().maxCoeff(), 1e-14);
	EXPECT_LT((g.translation() - pos).norm(), 1e-15);
	const V6 rpy = toPositionEulerZYX(g);
	EXPECT_LT((rpy.head<3>() - pos).norm(), 1e-15);
	EXPECT_NEAR(rpy[3], roll, 1e-13);
	EXPECT_NEAR(rpy[4], pitch, 1e-13);
	EXPECT_NEAR(rpy[5], yaw, 1e-13);
	const SE3d h = SE3d::exp(randomSE3Tangent(rng, 0.9));
	EXPECT_TRUE(interpolate(g, h, 0.0).isApprox(g, 1e-12));
	EXPECT_TRUE(interpolate(g, h, 1.0).isApprox(h, 1e-12));
	const Eigen::Vector3d p(0.1, 0.2, 0.3);
	EXPECT_LT((p / g - g.inverse().act(p)).norm(), 1e-15);
}

TEST(Instantiation, SE23InlFreeFunctions) {
	const Eigen::Vector3d pos(1.0, -2.0, 0.5), vel(0.3, 0.1, -0.4);
	const SO3d R = SO3d::exp(Eigen::Vector3d(0.2, -0.1, 0.4));
	const SE23<double> g = fromComponents<double>(pos, R, vel);
	EXPECT_LT((g.velocity() - vel).norm(), 1e-15);
	const SE23<double> e = fromPositionEulerVelocity<double>(pos, 0.3, -0.5, 1.2, vel);
	const Eigen::Matrix<double, 9, 1> v = toPositionEulerVelocity(e);
	EXPECT_NEAR(v[3], 0.3, 1e-13);
	EXPECT_NEAR(v[5], 1.2, 1e-13);
	EXPECT_LT((v.tail<3>() - vel).norm(), 1e-15);
	EXPECT_TRUE(interpolate(g, e, 0.0).isApprox(g, 1e-12));
	EXPECT_TRUE(interpolate(g, e, 1.0).isApprox(e, 1e-12));
}

TEST(Instantiation, SGal3InlFreeFunctions) {
	const Eigen::Vector3d pos(1.0, -2.0, 0.5), vel(0.3, 0.1, -0.4);
	const SO3d R = SO3d::exp(Eigen::Vector3d(0.2, -0.1, 0.4));
	const SGal3d g = fromComponents<double>(pos, R, vel, 2.0);
	EXPECT_NEAR(g.time(), 2.0, 1e-15);
	const SGal3d e = fromPositionEulerVelocityTime<double>(pos, 0.3, -0.5, 1.2, vel, 0.7);
	const Eigen::Matrix<double, 10, 1> v = toPositionEulerVelocityTime(e);
	EXPECT_NEAR(v[4], -0.5, 1e-13);
	EXPECT_NEAR(v[9], 0.7, 1e-15);
}
