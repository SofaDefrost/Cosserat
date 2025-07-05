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
#include <Cosserat/liegroups/SE3.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>

namespace sofa::component::cosserat::liegroups::testing {

using namespace sofa::testing;

/**
 * Test suite for SE3 Lie group implementation
 */
class SE3Test : public BaseTest
{
protected:
    using SE3d = SE3<double>;
    using SO3d = SO3<double>;
    using Vector3 = Eigen::Vector3d;
    using Vector6 = Eigen::Matrix<double, 6, 1>;
    using Matrix3 = Eigen::Matrix3d;
    using Matrix4 = Eigen::Matrix4d;
    using Matrix6 = Eigen::Matrix<double, 6, 6>;
    using Quaternion = Eigen::Quaterniond;

    const double pi = M_PI;
    const double eps = 1e-10;

    // Helper function to create twist vector
    Vector6 twist(const Vector3& v, const Vector3& omega) {
        Vector6 xi;
        xi << v, omega;
        return xi;
    }

    void SetUp() override {}
    void TearDown() override {}
};

/**
 * Test identity element properties
 */
TEST_F(SE3Test, Identity)
{
    SE3d id = SE3d::identity();
    
    // Test identity components
    EXPECT_TRUE(id.rotation().quaternion().isApprox(Quaternion::Identity()));
    EXPECT_TRUE(id.translation().isZero());
    
    // Test identity matrix
    EXPECT_TRUE(id.matrix().isApprox(Matrix4::Identity()));
    
    // Test composition with identity
    Vector3 axis = Vector3(1, 1, 1).normalized();
    Vector3 trans(1, 2, 3);
    SE3d g(SO3d(pi/4, axis), trans);  // 45° rotation + translation
    EXPECT_TRUE((g * id).isApprox(g));
    EXPECT_TRUE((id * g).isApprox(g));
}

/**
 * Test group operation (rigid transformation composition)
 */
TEST_F(SE3Test, GroupOperation)
{
    // Create transformations with different rotations and translations
    Vector3 axis1 = Vector3(1, 0, 0);  // X-axis
    Vector3 axis2 = Vector3(0, 1, 0);  // Y-axis
    Vector3 t1(1, 0, 0);  // X-translation
    Vector3 t2(0, 1, 0);  // Y-translation
    
    SE3d g1(SO3d(pi/2, axis1), t1);  // 90° around X + X-translation
    SE3d g2(SO3d(pi/2, axis2), t2);  // 90° around Y + Y-translation
    
    // Test composition
    SE3d g12 = g1 * g2;
    
    // Verify using homogeneous matrices
    Matrix4 T1 = g1.matrix();
    Matrix4 T2 = g2.matrix();
    Matrix4 T12 = g12.matrix();
    EXPECT_TRUE((T1 * T2).isApprox(T12));
    
    // Test non-commutativity
    SE3d g21 = g2 * g1;
    EXPECT_FALSE(g12.isApprox(g21));
}

/**
 * Test inverse operation
 */
TEST_F(SE3Test, Inverse)
{
    Vector3 axis = Vector3(1, 1, 0).normalized();
    Vector3 trans(1, 2, 3);
    SE3d g(SO3d(pi/3, axis), trans);  // 60° rotation + translation
    SE3d inv = g.inverse();
    
    // Test inverse properties
    EXPECT_TRUE((g * inv).isApprox(SE3d::identity()));
    EXPECT_TRUE((inv * g).isApprox(SE3d::identity()));
    
    // Test matrix inverse
    EXPECT_TRUE(inv.matrix().isApprox(g.matrix().inverse()));
    
    // Test inverse translation
    Vector3 expected_trans = -(g.rotation().inverse().act(trans));
    EXPECT_TRUE(inv.translation().isApprox(expected_trans));
}

/**
 * Test exponential and logarithm maps
 */
TEST_F(SE3Test, ExpLog)
{
    // Test exp(log(g)) = g
    Vector3 axis = Vector3(1, 2, 3).normalized();
    Vector3 trans(4, 5, 6);
    SE3d g(SO3d(pi/4, axis), trans);
    Vector6 xi = g.log();
    SE3d g2 = SE3d().exp(xi);
    EXPECT_TRUE(g.isApprox(g2));
    
    // Test pure translation
    Vector6 xi_trans = twist(Vector3(1, 2, 3), Vector3::Zero());
    SE3d g_trans = SE3d().exp(xi_trans);
    EXPECT_TRUE(g_trans.rotation().isApprox(SO3d::identity()));
    EXPECT_TRUE(g_trans.translation().isApprox(Vector3(1, 2, 3)));
    
    // Test pure rotation
    Vector6 xi_rot = twist(Vector3::Zero(), Vector3(0, 0, pi/2));
    SE3d g_rot = SE3d().exp(xi_rot);
    EXPECT_TRUE(g_rot.translation().isZero());
    EXPECT_TRUE(g_rot.rotation().isApprox(SO3d(pi/2, Vector3(0, 0, 1))));
    
    // Test small motion approximation
    Vector6 xi_small = twist(Vector3(0.001, 0.001, 0.001),
                           Vector3(0.001, 0.001, 0.001));
    SE3d g_small = SE3d().exp(xi_small);
    Matrix4 expected = Matrix4::Identity();
    expected.block<3,3>(0,0) = Matrix3::Identity() + SO3d::hat(xi_small.tail<3>());
    expected.block<3,1>(0,3) = xi_small.head<3>();
    EXPECT_TRUE(g_small.matrix().isApprox(expected, 1e-6));
}

/**
 * Test adjoint representation
 */
TEST_F(SE3Test, Adjoint)
{
    Vector3 axis = Vector3(0, 0, 1);  // Z-axis
    Vector3 trans(1, 2, 0);
    SE3d g(SO3d(pi/2, axis), trans);  // 90° around Z + translation
    Matrix6 Ad = g.adjoint();
    
    // Test adjoint structure
    Matrix3 R = g.rotation().matrix();
    Matrix3 t_hat = SO3d::hat(trans);
    
    EXPECT_TRUE(Ad.block<3,3>(0,0).isApprox(R));
    EXPECT_TRUE(Ad.block<3,3>(0,3).isApprox(t_hat * R));
    EXPECT_TRUE(Ad.block<3,3>(3,0).isZero());
    EXPECT_TRUE(Ad.block<3,3>(3,3).isApprox(R));
    
    // Test adjoint action on twists
    Vector6 xi = twist(Vector3(1, 0, 0), Vector3(0, 0, 1));
    Vector6 xi_transformed = Ad * xi;
    
    // Verify transformation matches conjugation
    SE3d h = SE3d().exp(xi);
    SE3d h_transformed = g * h * g.inverse();
    EXPECT_TRUE(h_transformed.isApprox(SE3d().exp(xi_transformed)));
}

/**
 * Test group action on points
 */
TEST_F(SE3Test, Action)
{
    // 90° rotation around Z-axis + translation in X
    SE3d g(SO3d(pi/2, Vector3(0, 0, 1)), Vector3(1, 0, 0));
    Vector3 p(1, 0, 0);  // Point on x-axis
    
    // Should map (1,0,0) to (1,1,0)
    Vector3 q = g.act(p);
    EXPECT_NEAR(q[0], 1.0, eps);
    EXPECT_NEAR(q[1], 1.0, eps);
    EXPECT_NEAR(q[2], 0.0, eps);
    
    // Test that action matches homogeneous transformation
    Vector4 p_hom;
    p_hom << p, 1.0;
    Vector4 q_hom = g.matrix() * p_hom;
    EXPECT_TRUE(q.isApprox(q_hom.head<3>()));
}

/**
 * Test interpolation
 */
TEST_F(SE3Test, Interpolation)
{
    // Create start and end transformations
    Vector3 axis = Vector3(1, 1, 1).normalized();
    SE3d start = SE3d::identity();
    SE3d end(SO3d(pi/2, axis), Vector3(2, 0, 0));
    
    // Test midpoint interpolation
    SE3d mid = interpolate(start, end, 0.5);
    
    // Verify midpoint properties
    EXPECT_TRUE(mid.rotation().isApprox(SO3d(pi/4, axis)));
    EXPECT_TRUE(mid.translation().isApprox(Vector3(1, 0, 0)));
    
    // Test boundary conditions
    EXPECT_TRUE(interpolate(start, end, 0.0).isApprox(start));
    EXPECT_TRUE(interpolate(start, end, 1.0).isApprox(end));
}

/**
 * Test conversion between different representations
 */
TEST_F(SE3Test, Conversions)
{
    // Create transformation from rotation and translation
    Vector3 axis = Vector3(1, 0, 0);  // X-axis
    double angle = pi/3;  // 60°
    Vector3 trans(1, 2, 3);
    SE3d g(SO3d(angle, axis), trans);
    
    // Test conversion to/from homogeneous matrix
    Matrix4 T = Matrix4::Identity();
    T.block<3,3>(0,0) = Eigen::AngleAxisd(angle, axis).toRotationMatrix();
    T.block<3,1>(0,3) = trans;
    
    SE3d g2(T);
    EXPECT_TRUE(g.isApprox(g2));
    EXPECT_TRUE(g.matrix().isApprox(T));
}

} // namespace sofa::component::cosserat::liegroups::testing
