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
#include <Cosserat/liegroups/SE23.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>

namespace sofa::component::cosserat::liegroups::testing {

using namespace sofa::testing;

/**
 * Test suite for SE_2(3) Lie group implementation
 */
class SE23Test : public BaseTest
{
protected:
    using SE23d = SE23<double>;
    using SE3d = SE3<double>;
    using SO3d = SO3<double>;
    using Vector3 = Eigen::Vector3d;
    using Vector9 = Eigen::Matrix<double, 9, 1>;
    using Matrix3 = Eigen::Matrix3d;
    using Matrix4 = Eigen::Matrix4d;
    using Matrix9 = Eigen::Matrix<double, 9, 9>;
    using Quaternion = Eigen::Quaterniond;

    const double pi = M_PI;
    const double eps = 1e-10;

    // Helper function to create twist vector with acceleration
    Vector9 twist(const Vector3& v, const Vector3& omega, const Vector3& a) {
        Vector9 xi;
        xi << v, omega, a;
        return xi;
    }

    void SetUp() override {}
    void TearDown() override {}
};

/**
 * Test identity element properties
 */
TEST_F(SE23Test, Identity)
{
    SE23d id = SE23d::identity();
    
    // Test identity components
    EXPECT_TRUE(id.pose().rotation().quaternion().isApprox(Quaternion::Identity()));
    EXPECT_TRUE(id.pose().translation().isZero());
    EXPECT_TRUE(id.velocity().isZero());
    
    // Test identity matrix
    EXPECT_TRUE(id.matrix().isIdentity());
    
    // Test composition with identity
    Vector3 axis = Vector3(1, 1, 1).normalized();
    Vector3 trans(1, 2, 3);
    Vector3 vel(0.1, 0.2, 0.3);
    SE23d g(SE3d(SO3d(pi/4, axis), trans), vel);
    EXPECT_TRUE((g * id).isApprox(g));
    EXPECT_TRUE((id * g).isApprox(g));
}

/**
 * Test group operation (extended pose composition)
 */
TEST_F(SE23Test, GroupOperation)
{
    // Create transformations with different rotations, translations, and velocities
    Vector3 axis1 = Vector3(1, 0, 0);  // X-axis
    Vector3 axis2 = Vector3(0, 1, 0);  // Y-axis
    Vector3 t1(1, 0, 0);  // X-translation
    Vector3 t2(0, 1, 0);  // Y-translation
    Vector3 v1(0.1, 0, 0);  // X-velocity
    Vector3 v2(0, 0.2, 0);  // Y-velocity
    
    SE23d g1(SE3d(SO3d(pi/2, axis1), t1), v1);
    SE23d g2(SE3d(SO3d(pi/2, axis2), t2), v2);
    
    // Test composition
    SE23d g12 = g1 * g2;
    
    // Verify pose composition
    EXPECT_TRUE(g12.pose().isApprox(g1.pose() * g2.pose()));
    
    // Verify velocity transformation
    Vector3 expected_vel = v1 + g1.pose().rotation().act(v2);
    EXPECT_TRUE(g12.velocity().isApprox(expected_vel));
    
    // Test non-commutativity
    SE23d g21 = g2 * g1;
    EXPECT_FALSE(g12.isApprox(g21));
}

/**
 * Test inverse operation
 */
TEST_F(SE23Test, Inverse)
{
    Vector3 axis = Vector3(1, 1, 0).normalized();
    Vector3 trans(1, 2, 3);
    Vector3 vel(0.1, 0.2, 0.3);
    SE23d g(SE3d(SO3d(pi/3, axis), trans), vel);
    SE23d inv = g.inverse();
    
    // Test inverse properties
    EXPECT_TRUE((g * inv).isApprox(SE23d::identity()));
    EXPECT_TRUE((inv * g).isApprox(SE23d::identity()));
    
    // Test inverse pose
    EXPECT_TRUE(inv.pose().isApprox(g.pose().inverse()));
    
    // Test inverse velocity
    Vector3 expected_vel = -g.pose().rotation().inverse().act(vel);
    EXPECT_TRUE(inv.velocity().isApprox(expected_vel));
}

/**
 * Test exponential and logarithm maps
 */
TEST_F(SE23Test, ExpLog)
{
    // Test exp(log(g)) = g
    Vector3 axis = Vector3(1, 2, 3).normalized();
    Vector3 trans(4, 5, 6);
    Vector3 vel(0.1, 0.2, 0.3);
    SE23d g(SE3d(SO3d(pi/4, axis), trans), vel);
    Vector9 xi = g.log();
    SE23d g2 = SE23d().exp(xi);
    EXPECT_TRUE(g.isApprox(g2));
    
    // Test pure translation
    Vector9 xi_trans = twist(Vector3(1, 2, 3), Vector3::Zero(), Vector3::Zero());
    SE23d g_trans = SE23d().exp(xi_trans);
    EXPECT_TRUE(g_trans.pose().rotation().isApprox(SO3d::identity()));
    EXPECT_TRUE(g_trans.pose().translation().isApprox(Vector3(1, 2, 3)));
    EXPECT_TRUE(g_trans.velocity().isZero());
    
    // Test pure velocity
    Vector9 xi_vel = twist(Vector3::Zero(), Vector3::Zero(), Vector3(0.1, 0.2, 0.3));
    SE23d g_vel = SE23d().exp(xi_vel);
    EXPECT_TRUE(g_vel.pose().isApprox(SE3d::identity()));
    EXPECT_TRUE(g_vel.velocity().isApprox(Vector3(0.1, 0.2, 0.3)));
    
    // Test small motion approximation
    Vector9 xi_small = twist(Vector3(0.001, 0.001, 0.001),
                           Vector3(0.001, 0.001, 0.001),
                           Vector3(0.001, 0.001, 0.001));
    SE23d g_small = SE23d().exp(xi_small);
    Matrix4 expected_pose = Matrix4::Identity();
    expected_pose.block<3,3>(0,0) = Matrix3::Identity() + SO3d::hat(xi_small.segment<3>(3));
    expected_pose.block<3,1>(0,3) = xi_small.head<3>();
    EXPECT_TRUE(g_small.pose().matrix().isApprox(expected_pose, 1e-6));
    EXPECT_TRUE(g_small.velocity().isApprox(xi_small.tail<3>(), 1e-6));
}

/**
 * Test adjoint representation
 */
TEST_F(SE23Test, Adjoint)
{
    Vector3 axis = Vector3(0, 0, 1);  // Z-axis
    Vector3 trans(1, 2, 0);
    Vector3 vel(0.1, 0.2, 0);
    SE23d g(SE3d(SO3d(pi/2, axis), trans), vel);
    Matrix9 Ad = g.adjoint();
    
    // Test adjoint structure
    Matrix3 R = g.pose().rotation().matrix();
    Matrix3 t_hat = SO3d::hat(trans);
    Matrix3 v_hat = SO3d::hat(vel);
    
    EXPECT_TRUE(Ad.block<3,3>(0,0).isApprox(R));
    EXPECT_TRUE(Ad.block<3,3>(0,3).isApprox(t_hat * R));
    EXPECT_TRUE(Ad.block<3,3>(0,6).isApprox(v_hat * R));
    EXPECT_TRUE(Ad.block<3,3>(3,0).isZero());
    EXPECT_TRUE(Ad.block<3,3>(3,3).isApprox(R));
    EXPECT_TRUE(Ad.block<3,3>(6,6).isApprox(R));
    
    // Test adjoint action
    Vector9 xi = twist(Vector3(1, 0, 0), Vector3(0, 0, 1), Vector3(0.1, 0, 0));
    Vector9 xi_transformed = Ad * xi;
    
    // Verify transformation matches conjugation
    SE23d h = SE23d().exp(xi);
    SE23d h_transformed = g * h * g.inverse();
    EXPECT_TRUE(h_transformed.isApprox(SE23d().exp(xi_transformed)));
}

/**
 * Test group action on points with velocity
 */
TEST_F(SE23Test, Action)
{
    // 90° rotation around Z-axis + translation in X + velocity
    SE23d g(SE3d(SO3d(pi/2, Vector3(0, 0, 1)), Vector3(1, 0, 0)),
            Vector3(0.1, 0.2, 0));
            
    // Point with velocity
    Vector3 p(1, 0, 0);
    Vector3 v(0.1, 0, 0);
    Vector9 state;
    state << p, v, Vector3::Zero();
    
    // Test transformation
    Vector9 transformed = g.act(state);
    Vector3 p_new = transformed.head<3>();
    Vector3 v_new = transformed.segment<3>(3);
    
    // Check position transformation
    EXPECT_NEAR(p_new[0], 1.0, eps);  // Original x-translation + point x
    EXPECT_NEAR(p_new[1], 1.0, eps);  // Rotated point x
    EXPECT_NEAR(p_new[2], 0.0, eps);
    
    // Check velocity transformation
    Vector3 expected_vel = g.pose().rotation().act(v) + g.velocity();
    EXPECT_TRUE(v_new.isApprox(expected_vel));
}

/**
 * Test interpolation
 */
TEST_F(SE23Test, Interpolation)
{
    // Create start and end states
    Vector3 axis = Vector3(1, 1, 1).normalized();
    SE23d start = SE23d::identity();
    SE23d end(SE3d(SO3d(pi/2, axis), Vector3(2, 0, 0)),
              Vector3(0.2, 0, 0));
    
    // Test midpoint interpolation
    SE23d mid = interpolate(start, end, 0.5);
    
    // Verify midpoint properties
    EXPECT_TRUE(mid.pose().rotation().isApprox(SO3d(pi/4, axis)));
    EXPECT_TRUE(mid.pose().translation().isApprox(Vector3(1, 0, 0)));
    EXPECT_TRUE(mid.velocity().isApprox(Vector3(0.1, 0, 0)));
    
    // Test boundary conditions
    EXPECT_TRUE(interpolate(start, end, 0.0).isApprox(start));
    EXPECT_TRUE(interpolate(start, end, 1.0).isApprox(end));
}

} // namespace sofa::component::cosserat::liegroups::testing
