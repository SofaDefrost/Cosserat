/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                  *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/\>.        *
******************************************************************************/

#include <sofa/testing/BaseTest.h>
#include <Cosserat/liegroups/SO3.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>

namespace sofa::component::cosserat::liegroups::testing {

using namespace sofa::testing;

/**
 * Test suite for SO3 Lie group implementation
 */
class SO3Test : public BaseTest
{
protected:
    using SO3d = SO3<double>;
    using Vector3 = Eigen::Vector3d;
    using Matrix3 = Eigen::Matrix3d;
    using Quaternion = Eigen::Quaterniond;

    const double pi = M_PI;
    const double eps = 1e-10;

    // Helper function to create rotation vector
    Vector3 rotationVector(double angle, const Vector3& axis) {
        return axis.normalized() * angle;
    }

    void SetUp() override {}
    void TearDown() override {}
};

/**
 * Test identity element properties
 */
TEST_F(SO3Test, Identity)
{
    SO3d id = SO3d::identity();
    
    // Test identity quaternion properties
    EXPECT_NEAR(id.quaternion().w(), 1.0, eps);
    EXPECT_NEAR(id.quaternion().x(), 0.0, eps);
    EXPECT_NEAR(id.quaternion().y(), 0.0, eps);
    EXPECT_NEAR(id.quaternion().z(), 0.0, eps);
    
    // Test identity matrix
    EXPECT_TRUE(id.matrix().isApprox(Matrix3::Identity()));
    
    // Test composition with identity
    Vector3 axis = Vector3(1, 1, 1).normalized();
    SO3d rot(pi/4, axis);  // 45-degree rotation around (1,1,1)
    EXPECT_TRUE((rot * id).isApprox(rot));
    EXPECT_TRUE((id * rot).isApprox(rot));
}

/**
 * Test group operation (rotation composition)
 */
TEST_F(SO3Test, GroupOperation)
{
    // Create two rotations with different axes
    Vector3 axis1 = Vector3(1, 0, 0);  // X-axis
    Vector3 axis2 = Vector3(0, 1, 0);  // Y-axis
    SO3d rx(pi/2, axis1);  // 90° around X
    SO3d ry(pi/2, axis2);  // 90° around Y
    
    // Test composition
    SO3d rxy = rx * ry;
    
    // Verify using matrix multiplication
    Matrix3 Rx = rx.matrix();
    Matrix3 Ry = ry.matrix();
    Matrix3 Rxy = rxy.matrix();
    EXPECT_TRUE((Rx * Ry).isApprox(Rxy));
    
    // Test non-commutativity
    SO3d ryx = ry * rx;
    EXPECT_FALSE(rxy.isApprox(ryx));
}

/**
 * Test inverse operation
 */
TEST_F(SO3Test, Inverse)
{
    Vector3 axis = Vector3(1, 1, 0).normalized();
    SO3d rot(pi/3, axis);  // 60° rotation
    SO3d inv = rot.inverse();
    
    // Test inverse properties
    EXPECT_TRUE((rot * inv).isApprox(SO3d::identity()));
    EXPECT_TRUE((inv * rot).isApprox(SO3d::identity()));
    
    // Test matrix inverse
    EXPECT_TRUE(inv.matrix().isApprox(rot.matrix().inverse()));
    
    // Test quaternion conjugate
    EXPECT_TRUE(inv.quaternion().coeffs().isApprox(rot.quaternion().conjugate().coeffs()));
}

/**
 * Test exponential and logarithm maps
 */
TEST_F(SO3Test, ExpLog)
{
    // Test exp(log(g)) = g
    Vector3 axis = Vector3(1, 2, 3).normalized();
    SO3d rot(pi/4, axis);  // 45° rotation
    Vector3 omega = rot.log();
    SO3d rot2 = SO3d().exp(omega);
    EXPECT_TRUE(rot.isApprox(rot2));
    
    // Test small rotation approximation
    Vector3 small_omega = Vector3(0.001, 0.001, 0.001);
    SO3d small_rot = SO3d().exp(small_omega);
    Matrix3 expected = Matrix3::Identity() + SO3d::hat(small_omega);
    EXPECT_TRUE(small_rot.matrix().isApprox(expected, 1e-6));
    
    // Test rotation vector recovery
    Vector3 rot_vec = rotationVector(pi/3, Vector3(1,0,0));
    SO3d g = SO3d().exp(rot_vec);
    EXPECT_TRUE(g.log().isApprox(rot_vec));
}

/**
 * Test adjoint representation
 */
TEST_F(SO3Test, Adjoint)
{
    Vector3 axis = Vector3(0, 0, 1);  // Z-axis
    SO3d rot(pi/2, axis);  // 90° around Z
    Matrix3 Ad = rot.adjoint();
    
    // Adjoint should be the rotation matrix itself for SO(3)
    EXPECT_TRUE(Ad.isApprox(rot.matrix()));
    
    // Test adjoint action on vectors
    Vector3 v(1, 0, 0);
    EXPECT_TRUE((Ad * v).isApprox(rot.act(v)));
}

/**
 * Test group action on points
 */
TEST_F(SO3Test, Action)
{
    // 90° rotation around Z-axis
    SO3d rot(pi/2, Vector3(0, 0, 1));
    Vector3 p(1, 0, 0);  // Point on x-axis
    
    // Should map (1,0,0) to (0,1,0)
    Vector3 q = rot.act(p);
    EXPECT_NEAR(q[0], 0.0, eps);
    EXPECT_NEAR(q[1], 1.0, eps);
    EXPECT_NEAR(q[2], 0.0, eps);
    
    // Test that action preserves length
    EXPECT_NEAR(p.norm(), q.norm(), eps);
    
    // Test that action matches matrix multiplication
    EXPECT_TRUE(q.isApprox(rot.matrix() * p));
}

/**
 * Test hat and vee operators
 */
TEST_F(SO3Test, HatVee)
{
    Vector3 omega(1, 2, 3);
    Matrix3 Omega = SO3d::hat(omega);
    
    // Test skew-symmetry
    EXPECT_TRUE(Omega.transpose().isApprox(-Omega));
    
    // Test vee operator (inverse of hat)
    EXPECT_TRUE(SO3d::vee(Omega).isApprox(omega));
    
    // Test that hat(omega) * v = omega × v
    Vector3 v(4, 5, 6);
    EXPECT_TRUE((Omega * v).isApprox(omega.cross(v)));
}

/**
 * Test interpolation
 */
TEST_F(SO3Test, Interpolation)
{
    // Create start and end rotations
    Vector3 axis = Vector3(1, 1, 1).normalized();
    SO3d start = SO3d::identity();
    SO3d end(pi/2, axis);  // 90° rotation
    
    // Test midpoint interpolation
    SO3d mid = interpolate(start, end, 0.5);
    
    // Midpoint should be 45° rotation around same axis
    SO3d expected(pi/4, axis);
    EXPECT_TRUE(mid.isApprox(expected));
    
    // Test boundary conditions
    EXPECT_TRUE(interpolate(start, end, 0.0).isApprox(start));
    EXPECT_TRUE(interpolate(start, end, 1.0).isApprox(end));
}

/**
 * Test conversion between different rotation representations
 */
TEST_F(SO3Test, Conversions)
{
    // Create rotation from angle-axis
    Vector3 axis = Vector3(1, 0, 0);  // X-axis
    double angle = pi/3;  // 60°
    SO3d rot(angle, axis);
    
    // Test conversion to/from quaternion
    Quaternion quat(Eigen::AngleAxisd(angle, axis));
    EXPECT_TRUE(rot.quaternion().coeffs().isApprox(quat.coeffs()));
    
    // Test conversion to/from rotation matrix
    Matrix3 R = quat.toRotationMatrix();
    EXPECT_TRUE(rot.matrix().isApprox(R));
    
    // Test conversion to/from angle-axis
    Eigen::AngleAxisd aa = rot.angleAxis();
    EXPECT_NEAR(aa.angle(), angle, eps);
    EXPECT_TRUE(aa.axis().isApprox(axis));
}

} // namespace sofa::component::cosserat::liegroups::testing
