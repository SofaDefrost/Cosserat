// Bundle<G1, G2, ...> is a product group: every operation must act component-wise, with the
// algebra / action offsets taken from each group's AlgebraDim / ActionDim (not Dim: SE3 has
// Dim 4 but 6 algebra coordinates, SO2 has Dim 2 but 1, SE2 acts on 2 coordinates).
// The class used to fail to compile for every combination of groups.

#include <gtest/gtest.h>
#include <liegroups/Bundle.h>

#include <random>

#include "LieGroupReference.h"

using namespace sofa::component::cosserat::liegroups;
using namespace liegroup_ref;

namespace {
using R2 = RealSpace<double, 2>;
using PoseAndJoints = Bundle<SE3d, R2>;           // algebra 6 + 2, action 3 + 2
using PlanarChain = Bundle<SE2d, SO2d, SO2d>;     // algebra 3 + 1 + 1, action 2 + 2 + 2
using Mixed = Bundle<SO3d, SE2d, R2>;             // algebra 3 + 3 + 2, action 3 + 2 + 2
} // namespace

TEST(BundleAudit, DimensionsAreTheSumsOfAlgebraAndActionDimensions) {
	static_assert(PoseAndJoints::AlgebraDim == 8);
	static_assert(PoseAndJoints::ActionDim == 5);
	static_assert(PlanarChain::AlgebraDim == 5);
	static_assert(PlanarChain::ActionDim == 6);
	static_assert(Mixed::AlgebraDim == 8);
	static_assert(Mixed::ActionDim == 7);
	EXPECT_EQ(PoseAndJoints().algebraDimension(), 8);
	EXPECT_EQ(PoseAndJoints().actionDimension(), 5);
}

TEST(BundleAudit, ExpLogAndCompositionActComponentWise) {
	std::mt19937_64 rng(71);
	const Eigen::Matrix<double, 6, 1> xa = randomSE3Tangent(rng, 1.3), xb = randomSE3Tangent(rng, 0.4);
	PoseAndJoints::TangentVector a, b;
	a << xa, 0.3, -0.7;
	b << xb, -1.1, 0.2;
	const PoseAndJoints ga = PoseAndJoints::exp(a), gb = PoseAndJoints::exp(b);
	// component-wise exp
	EXPECT_TRUE(ga.get<0>().isApprox(SE3d::exp(xa), 1e-14));
	EXPECT_NEAR(ga.get<1>().log()[0], 0.3, 1e-15);
	EXPECT_LT((ga.log() - a).cwiseAbs().maxCoeff(), 1e-13);
	// composition and inverse act per component
	const PoseAndJoints prod = ga * gb;
	EXPECT_TRUE(prod.get<0>().isApprox(SE3d::exp(xa) * SE3d::exp(xb), 1e-14));
	EXPECT_LT((prod.get<1>().log() - (R2::TangentVector(0.3, -0.7) + R2::TangentVector(-1.1, 0.2))).norm(), 1e-15);
	EXPECT_TRUE((ga * ga.inverse()).isApprox(PoseAndJoints(), 1e-13));
	// generic plus / minus of the base class work through the hooks
	EXPECT_LT((ga.minus(prod) - gb.log()).cwiseAbs().maxCoeff(), 1e-12);
}

TEST(BundleAudit, AdjointIsBlockDiagonalAndActOffsetsFollowActionDim) {
	std::mt19937_64 rng(72);
	const Eigen::Matrix<double, 6, 1> xa = randomSE3Tangent(rng, 0.8);
	PoseAndJoints::TangentVector a;
	a << xa, 0.5, 0.25;
	const PoseAndJoints g = PoseAndJoints::exp(a);
	const auto Ad = g.adjoint();
	EXPECT_LT((Ad.block<6, 6>(0, 0) - SE3d::exp(xa).computeAdjoint()).cwiseAbs().maxCoeff(), 1e-14);
	EXPECT_LT((Ad.block<2, 2>(6, 6) - Eigen::Matrix2d::Identity()).cwiseAbs().maxCoeff(), 1e-15);
	const double offBlocks = Ad.block<6, 2>(0, 6).cwiseAbs().maxCoeff() + Ad.block<2, 6>(6, 0).cwiseAbs().maxCoeff();
	EXPECT_LT(offBlocks, 1e-15);
	// action: SE3 on the first 3 coordinates, translation of R^2 on the last 2
	PoseAndJoints::ActionVector x;
	x << 0.1, 0.2, 0.3, 1.0, -1.0;
	const auto y = g.act(x);
	EXPECT_LT((y.head<3>() - SE3d::exp(xa).act(x.head<3>())).norm(), 1e-14);
	EXPECT_LT((y.tail<2>() - (x.tail<2>() + Eigen::Vector2d(0.5, 0.25))).norm(), 1e-15);
}

TEST(BundleAudit, GroupsWithDimDifferentFromAlgebraDim) {
	// SO2 (Dim 2, algebra 1) and SE2 (action 2) inside a bundle
	PlanarChain::TangentVector a;
	a << 0.4, -0.2, 0.9, 0.3, -1.2;
	const PlanarChain g = PlanarChain::exp(a);
	EXPECT_LT((g.log() - a).cwiseAbs().maxCoeff(), 1e-13);
	PlanarChain::ActionVector x;
	x << 1, 2, 3, 4, 5, 6;
	const auto y = g.act(x);
	EXPECT_LT((y.segment<2>(0) - g.get<0>().act(x.segment<2>(0))).norm(), 1e-14);
	EXPECT_LT((y.segment<2>(2) - g.get<1>().act(x.segment<2>(2))).norm(), 1e-14);
	EXPECT_LT((y.segment<2>(4) - g.get<2>().act(x.segment<2>(4))).norm(), 1e-14);
}

TEST(BundleAudit, JacobiansAndBchUseTheBlockDiagonalAlgebra) {
	std::mt19937_64 rng(73);
	const Eigen::Matrix<double, 6, 1> xa = randomSE3Tangent(rng, 0.9);
	PoseAndJoints::TangentVector a;
	a << xa, 0.5, 0.25;
	const auto Jl = PoseAndJoints::leftJacobian(a);
	EXPECT_LT((Jl.block<6, 6>(0, 0) - SE3d::computeLeftJacobian(xa)).cwiseAbs().maxCoeff(), 1e-13);
	EXPECT_LT((Jl.block<2, 2>(6, 6) - Eigen::Matrix2d::Identity()).cwiseAbs().maxCoeff(), 1e-15);
	PoseAndJoints::TangentVector b = a * 0.02, c = a.reverse() * 0.02;
	const auto ref = (PoseAndJoints::exp(b) * PoseAndJoints::exp(c)).log();
	EXPECT_LT((PoseAndJoints::BCH(b, c, 3) - ref).norm(), 1e-6);
}

TEST(BundleAudit, RandomIdentityAndStreaming) {
	std::mt19937_64 rng(74);
	const Mixed g = Mixed::random(rng);
	EXPECT_TRUE((g * g.inverse()).isApprox(Mixed::identity(), 1e-12));
	EXPECT_TRUE(Mixed::Identity().isApprox(Mixed(), 1e-15));
	std::ostringstream os;
	os << g;
	EXPECT_FALSE(os.str().empty());
}
