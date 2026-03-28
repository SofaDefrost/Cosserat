#include <Cosserat/mapping/BeamShapeInterpolation.h>
#include <gtest/gtest.h>
#include <liegroups/SE3.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/VecTypes.h>

using namespace sofa::component::cosserat::liegroups;
using namespace Cosserat::mapping;

TEST(AdvancedFeaturesTest, BCHApproximation) {
	using SE3Type = SE3<double>;
	using TangentVector = SE3Type::TangentVector;

	TangentVector x = TangentVector::Zero();
	x[0] = 0.1; // small translation

	TangentVector y = TangentVector::Zero();
	y[1] = 0.1; // small translation

	// For pure translations, [x,y] = 0, so BCH(x,y) = x+y
	TangentVector z = SE3Type::computeBCH(x, y);
	TangentVector expected = x + y;

	EXPECT_TRUE(z.isApprox(expected));

	// For rotations
	x = TangentVector::Zero();
	x[3] = 0.1; // rotation X
	y = TangentVector::Zero();
	y[4] = 0.1; // rotation Y

	// [x,y] should be non-zero (rotation Z)
	z = SE3Type::computeBCH(x, y);
	EXPECT_FALSE(z.isApprox(x + y));
}

TEST(AdvancedFeaturesTest, ParallelTransport) {
	using SE3Type = SE3<double>;
	using TangentVector = SE3Type::TangentVector;

	SE3Type pose = SE3Type::computeIdentity();
	TangentVector v = TangentVector::Ones();

	// Transport along 0-length geodesic (at identity) should be identity
	TangentVector v_trans = pose.parallelTransport(v);
	EXPECT_TRUE(v_trans.isApprox(v));

	// Transport along a path
	TangentVector u = TangentVector::Zero();
	u[3] = M_PI / 2.0; // 90 deg rotation X
	pose = SE3Type::computeExp(u);

	v = TangentVector::Zero();
	v[4] = 1.0; // Y axis vector

	// Transporting Y vector along X-axis rotation by 90 degrees
	// Should result in Z vector?
	// Parallel transport preserves the angle with the tangent of the curve?
	// This depends on the definition.
	// Our implementation: v - 0.5 * [u, v]
	v_trans = pose.parallelTransport(v);

	// [u, v]: u=[0,0,0, 1,0,0], v=[0,0,0, 0,1,0]
	// phi_u=[1,0,0], phi_v=[0,1,0]
	// phi_bracket = [1,0,0]x[0,1,0] = [0,0,1]
	// bracket = [0,0,0, 0,0,1]
	// v_trans = v - 0.5 * [0,0,1] = [0,1,-0.5]

	EXPECT_DOUBLE_EQ(v_trans[5], -0.5 * M_PI / 2.0); // Wait, u has magnitude PI/2
}

TEST(AdvancedFeaturesTest, ShapeInterpolation) {
	using SE3Type = SE3<double>;
	using TangentVector = SE3Type::TangentVector;

	std::vector<SE3Type> shape1(2);
	shape1[0] = SE3Type::computeIdentity();
	shape1[1] = SE3Type::computeIdentity();

	std::vector<SE3Type> shape2(2);
	TangentVector u = TangentVector::Zero();
	u[0] = 1.0;
	shape2[0] = SE3Type::computeExp(u);
	shape2[1] = SE3Type::computeExp(u);

	// Interpolate at 0.5
	auto interpolated = BeamShapeInterpolation::interpolateShapes(shape1, shape2, 0.5);

	EXPECT_EQ(interpolated.size(), 2);

	TangentVector expected_log = u * 0.5;
	EXPECT_TRUE(interpolated[0].log().isApprox(expected_log));
}
