// A class template member only fails to compile when it is instantiated, so
// -fsyntax-only on the headers misses members that cannot work. This file calls the
// public members that LieGroupBase promises to every group, for each group, and checks
// the result is meaningful (not only that it compiles).

#include <gtest/gtest.h>
#include <liegroups/SE3.h>
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
