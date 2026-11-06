// Audit tests of the generic helpers of LieGroupBase (LieGroupBase.inl) that are
// reached through SO3 / SE3: dexp_alternative, dexpInv_alternative, dlog, BCH and slerp.
// The references are the series of LieGroupReference.h (long double), not closed forms.

#include <gtest/gtest.h>
#include <liegroups/SE3.h>
#include <liegroups/SO3.h>

#include <random>

#include "LieGroupReference.h"

using namespace sofa::component::cosserat::liegroups;
using namespace liegroup_ref;

namespace {
using V6 = Eigen::Matrix<double, 6, 1>;
using M6 = Eigen::Matrix<double, 6, 6>;
const double kAngles[] = {1e-9, 1e-6, 1e-3, 0.1, 1.0, 2.0, 3.0};
} // namespace

TEST(LieGroupBaseGeneric, SlerpIsTheGeodesicAtParameterT) {
	const SO3d a = SO3d::exp(Eigen::Vector3d(0.3, 0.5, -0.2));
	const SO3d b = SO3d::exp(Eigen::Vector3d(-1.0, 0.4, 1.2));
	for (double t : {0.1, 0.25, 0.5, 0.9})
		EXPECT_LT((a.slerp(b, t).log() - a.interpolate(b, t).log()).norm(), 1e-12) << t;
	std::mt19937_64 rng(16);
	const SE3d g = SE3d::exp(randomSE3Tangent(rng, 0.7)), h = SE3d::exp(randomSE3Tangent(rng, 1.1));
	for (double t : {0.1, 0.25, 0.5, 0.9}) {
		// the geodesic g exp(t log(g^-1 h)) is at relative distance t from g
		const V6 d = (g.inverse() * g.slerp(h, t)).log();
		const V6 full = (g.inverse() * h).log();
		EXPECT_LT((d - t * full).norm(), 1e-10) << t;
	}
}
