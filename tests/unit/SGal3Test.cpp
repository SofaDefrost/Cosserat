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
#include <Cosserat/liegroups/SGal3.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>

namespace sofa::component::cosserat::liegroups::testing {

using namespace sofa::testing;

/**
 * Test suite for SGal(3) Lie group implementation
 */
class SGal3Test : public BaseTest
{
protected:
    using SGal3d = SGal3<double>;
    using SE3d = SE3<double>;
    using SO3d = SO3<double>;
    using Vector3 = Eigen::Vector3d;
    using Vector10 = Eigen::Matrix<double, 10, 1>;
    using Matrix3 = Eigen::Matrix3d;
    using Matrix4 = Eigen::Matrix4d;
    using Matrix10 = Eigen::Matrix<double, 10, 10>;
    using Quaternion = Eigen::Quaterniond;

    const double pi = M_PI;
    const double eps = 1e-10;

    // Helper function to create Galilean twist vector
    Vector10 twist(const Vector3& v, const Vector3& omega, const Vector3& beta, double tau) {
        Vector10 xi;
        xi << v, omega, beta, tau;
        return xi;
    }

    void SetUp() override {}
    void TearDown() override {}
};

/**
 * Test identity element properties
 */
TEST_F(SGal3Test, Identity)
{
    SGal3d id = SGal3d::identity();
    
    // Test identity components
    EXPECT_TRUE(id.pose().rotation().quaternion().isApprox(Quaternion::Identity()));
    EXPECT_TRUE(id.pose().translation().isZero());
    EXPECT_TRUE(id.velocity().isZero());
    EXPECT_NEAR(id.time(), 0.0, eps);
    
    // Test composition with identity
    Vector3 axis = Vector3(1, 1, 1).normalized();
    Vector3 trans(1, 2, 3);
    Vector3 vel(0.1, 0.2, 0.3);
    double t = 1.0;
    SGal3d g(SE3d(SO3d(pi/4, axis), trans), vel, t);
    EXPECT_TRUE((g * id).isApprox(g));
    EXPECT_TRUE((id * g).isApprox(g));
}

/**
 * Test group operation (Galilean transformation composition)
 */
TEST_F(SGal3Test, GroupOperation)
{
    // Create transformations with different components
    Vector3 axis1 = Vector3(1, 0, 0);  // X-axis
    Vector3 axis2 = Vector3(0, 1, 0);  // Y-axis
    Vector3 t1(1, 0, 0);  // X-translation
    Vector3 t2(0, 1, 0);  // Y-translation
    Vector3 v1(0.1, 0, 0);  // X-velocity
    Vector3 v2(0, 0.2, 0);  // Y-velocity
    double time1 = 1.0;
    double time2 = 2.0;
    
    SGal3d g1(SE3d(SO3d(pi/2, axis1), t1), v1, time1);
    SGal3d g2(SE3d(SO3d(pi/2, axis2), t2), v2, time2);
    
    // Test composition
    SGal3d g12 = g1 * g2;
    
    // Verify pose composition
    EXPECT_TRUE(g12.pose().isApprox(g1.pose() * g2.pose()));
    
    // Verify velocity transformation
    Vector3 expected_vel = v1 + g1.pose().rotation().act(v2);
    EXPECT_TRUE(g12.velocity().isApprox(expected_vel));
    
    // Verify time addition
    EXPECT_NEAR(g12.time(), time1 + time2, eps);
    
    // Test non-commutativity
    SGal3d g21 = g2 * g1;
    EXPECT_FALSE(g12.isApprox(g21));
}

/**
 * Test inverse operation
 */
TEST_F(SGal3Test, Inverse)
{
    Vector3 axis = Vector3(1, 1, 0).normalized();
    Vector3 trans(1, 2, 3);
    Vector3 vel(0.1, 0.2, 0.3);
    double time = 1.5;
    SGal3d g(SE3d(SO3d(pi/3, axis), trans), vel, time);
    SGal3d inv = g.inverse();
    
    // Test inverse properties
    EXPECT_TRUE((g * inv).isApprox(SGal3d::identity()));
    EXPECT_TRUE((inv * g).isApprox(SGal3d::identity()));
    
    // Test inverse pose
    EXPECT_TRUE(inv.pose().isApprox(g.pose().inverse()));
    
    // Test inverse velocity
    Vector3 expected_vel = -g.pose().rotation().inverse().act(vel);
    EXPECT_TRUE(inv.velocity().isApprox(expected_vel));
    
    // Test inverse time
    EXPECT_NEAR(inv.time(), -time, eps);
}

/**
 * Test exponential and logarithm maps
 */
TEST_F(SGal3Test, ExpLog)
{
    // Test exp(log(g)) = g
    Vector3 axis = Vector3(1, 2, 3).normalized();
    Vector3 trans(4, 5, 6);
    Vector3 vel(0.1, 0.2, 0.3);
    double time = 1.0;
    SGal3d g(SE3d(SO3d(pi/4, axis), trans), vel, time);
    Vector10 xi = g.log();
    SGal3d g2 = SGal3d().exp(xi);
    EXPECT_TRUE(g.isApprox(g2));
    
    // Test pure translation
    Vector10 xi_trans = twist(Vector3(1, 2, 3), Vector3::Zero(), Vector3::Zero(), 0);
    SGal3d g_trans = SGal3d().exp(xi_trans);
    EXPECT_TRUE(g_trans.pose().rotation().isApprox(SO3d::identity()));
    EXPECT_TRUE(g_trans.pose().translation().isApprox(Vector3(1, 2, 3)));
    EXPECT_TRUE(g_trans.velocity().isZero());
    EXPECT_NEAR(g_trans.time(), 0.0, eps);
    
    // Test pure velocity and time
    Vector10 xi_vel = twist(Vector3::Zero(), Vector3::Zero(), Vector3(0.1, 0.2, 0.3), 1.0);
    SGal3d g_vel = SGal3d().exp(xi_vel);
    EXPECT_TRUE(g_vel.pose().isApprox(SE3d::identity()));
    EXPECT_TRUE(g_vel.velocity().isApprox(Vector3(0.1, 0.2, 0.3)));
    EXPECT_NEAR(g_vel.time(), 1.0, eps);
}

/**
 * Test adjoint representation
 */
TEST_F(SGal3Test, Adjoint)
{
    Vector3 axis = Vector3(0, 0, 1);  // Z-axis
    Vector3 trans(1, 2, 0);
    Vector3 vel(0.1, 0.2, 0);
    double time = 1.0;
    SGal3d g(SE3d(SO3d(pi/2, axis), trans), vel, time);
    Matrix10 Ad = g.adjoint();
    
    // Test adjoint structure
    Matrix3 R = g.pose().rotation().matrix();
    Matrix3 t_hat = SO3d::hat(trans);
    Matrix3 v_hat = SO3d::hat(vel);
    
    // Verify block structure
    EXPECT_TRUE(Ad.block<3,3>(0,0).isApprox(R));
    EXPECT_TRUE(Ad.block<3,3>(0,3).isApprox(t_hat * R));
    EXPECT_TRUE(Ad.block<3,3>(0,6).isApprox(v_hat * R));
    EXPECT_TRUE(Ad.block<3,3>(3,0).isZero());
    EXPECT_TRUE(Ad.block<3,3>(3,3).isApprox(R));
    EXPECT_TRUE(Ad.block<3,3>(6,6).isApprox(R));
    EXPECT_NEAR(Ad(9,9), 1.0, eps);
}

/**
 * Test group action on points with velocity and time
 */
TEST_F(SGal3Test, Action)
{
    // Create a Galilean transformation
    SGal3d g(SE3d(SO3d(pi/2, Vector3(0, 0, 1)), Vector3(1, 0, 0)),
             Vector3(0.1, 0.2, 0), 1.0);
            
    // Point with velocity and time
    Vector3 p(1, 0, 0);
    Vector3 v(0.1, 0, 0);
    Vector3 boost(0.01, 0, 0);
    double t = 0.5;
    Vector10 state;
    state << p, v, boost, t;
    
    // Test transformation
    Vector10 transformed = g.act(state);
    Vector3 p_new = transformed.head<3>();
    Vector3 v_new = transformed.segment<3>(3);
    Vector3 boost_new = transformed.segment<3>(6);
    double t_new = transformed[9];
    
    // Check position transformation with time evolution
    Vector3 expected_pos = g.pose().act(p) + g.velocity() * t;
    EXPECT_TRUE(p_new.isApprox(expected_pos));
    
    // Check velocity transformation
    Vector3 expected_vel = g.pose().rotation().act(v) + g.velocity();
    EXPECT_TRUE(v_new.isApprox(expected_vel));
    
    // Check boost transformation
    Vector3 expected_boost = g.pose().rotation().act(boost);
    EXPECT_TRUE(boost_new.isApprox(expected_boost));
    
    // Check time transformation
    EXPECT_NEAR(t_new, t + g.time(), eps);
}

/**
 * Test interpolation
 */
TEST_F(SGal3Test, Interpolation)
{
    // Create start and end states
    Vector3 axis = Vector3(1, 1, 1).normalized();
    SGal3d start = SGal3d::identity();
    SGal3d end(SE3d(SO3d(pi/2, axis), Vector3(2, 0, 0)),
               Vector3(0.2, 0, 0), 2.0);
    
    // Test midpoint interpolation
    SGal3d mid = interpolate(start, end, 0.5);
    
    // Verify midpoint properties
    EXPECT_TRUE(mid.pose().rotation().isApprox(SO3d(pi/4, axis)));
    EXPECT_TRUE(mid.pose().translation().isApprox(Vector3(1, 0, 0)));
    EXPECT_TRUE(mid.velocity().isApprox(Vector3(0.1, 0, 0)));
    EXPECT_NEAR(mid.time(), 1.0, eps);
    
    // Test boundary conditions
    EXPECT_TRUE(interpolate(start, end, 0.0).isApprox(start));
    EXPECT_TRUE(interpolate(start, end, 1.0).isApprox(end));
}

} // namespace sofa::component::cosserat::liegroups::testing
