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
#include <Cosserat/liegroups/Bundle.h>
#include <Cosserat/liegroups/RealSpace.h>
#include <Cosserat/liegroups/SO3.h>
#include <Cosserat/liegroups/SE3.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>

namespace sofa::component::cosserat::liegroups::testing {

using namespace sofa::testing;

/**
 * Test suite for Bundle Lie group implementation
 */
class BundleTest : public BaseTest
{
protected:
    // Define common types
    using RealSpace3d = RealSpace<double, 3>;
    using SO3d = SO3<double>;
    using SE3d = SE3<double>;
    using Vector3 = Eigen::Vector3d;
    using Matrix3 = Eigen::Matrix3d;
    using Quaternion = Eigen::Quaterniond;

    // Define test bundle types
    using PoseVel = Bundle<SE3d, RealSpace3d>;  // Pose + velocity
    using MultiBody = Bundle<SE3d, SE3d, SE3d>;  // Three rigid bodies
    using ComplexSystem = Bundle<SE3d, SO3d, RealSpace3d>;  // Mixed system

    const double pi = M_PI;
    const double eps = 1e-10;

    void SetUp() override {}
    void TearDown() override {}
};

/**
 * Test identity element properties
 */
TEST_F(BundleTest, Identity)
{
    // Test PoseVel identity
    PoseVel id_pv = PoseVel::identity();
    EXPECT_TRUE(id_pv.get<0>().isApprox(SE3d::identity()));
    EXPECT_TRUE(id_pv.get<1>().isApprox(RealSpace3d::identity()));
    
    // Test MultiBody identity
    MultiBody id_mb = MultiBody::identity();
    EXPECT_TRUE(id_mb.get<0>().isApprox(SE3d::identity()));
    EXPECT_TRUE(id_mb.get<1>().isApprox(SE3d::identity()));
    EXPECT_TRUE(id_mb.get<2>().isApprox(SE3d::identity()));
    
    // Test right and left identity
    Vector3 axis = Vector3(1, 1, 1).normalized();
    PoseVel g(SE3d(SO3d(pi/4, axis), Vector3(1, 2, 3)),
              RealSpace3d(Vector3(0.1, 0.2, 0.3)));
    
    EXPECT_TRUE((g * id_pv).isApprox(g));
    EXPECT_TRUE((id_pv * g).isApprox(g));
}

/**
 * Test group operation (component-wise composition)
 */
TEST_F(BundleTest, GroupOperation)
{
    // Test PoseVel composition
    Vector3 axis1(1, 0, 0), axis2(0, 1, 0);
    Vector3 trans1(1, 0, 0), trans2(0, 1, 0);
    Vector3 vel1(0.1, 0, 0), vel2(0, 0.2, 0);
    
    PoseVel g1(SE3d(SO3d(pi/2, axis1), trans1), RealSpace3d(vel1));
    PoseVel g2(SE3d(SO3d(pi/2, axis2), trans2), RealSpace3d(vel2));
    
    PoseVel g12 = g1 * g2;
    
    // Verify component-wise composition
    EXPECT_TRUE(g12.get<0>().isApprox(g1.get<0>() * g2.get<0>()));
    EXPECT_TRUE(g12.get<1>().isApprox(g1.get<1>() * g2.get<1>()));
    
    // Test MultiBody composition
    MultiBody mb1(SE3d(SO3d(pi/2, axis1), trans1),
                 SE3d(SO3d(pi/3, axis1), trans1),
                 SE3d(SO3d(pi/4, axis1), trans1));
    
    MultiBody mb2(SE3d(SO3d(pi/2, axis2), trans2),
                 SE3d(SO3d(pi/3, axis2), trans2),
                 SE3d(SO3d(pi/4, axis2), trans2));
    
    MultiBody mb12 = mb1 * mb2;
    
    // Verify all bodies composed correctly
    EXPECT_TRUE(mb12.get<0>().isApprox(mb1.get<0>() * mb2.get<0>()));
    EXPECT_TRUE(mb12.get<1>().isApprox(mb1.get<1>() * mb2.get<1>()));
    EXPECT_TRUE(mb12.get<2>().isApprox(mb1.get<2>() * mb2.get<2>()));
}

/**
 * Test inverse operation
 */
TEST_F(BundleTest, Inverse)
{
    // Test PoseVel inverse
    Vector3 axis = Vector3(1, 1, 0).normalized();
    Vector3 trans(1, 2, 3);
    Vector3 vel(0.1, 0.2, 0.3);
    
    PoseVel g(SE3d(SO3d(pi/3, axis), trans), RealSpace3d(vel));
    PoseVel inv = g.inverse();
    
    // Test inverse properties
    EXPECT_TRUE((g * inv).isApprox(PoseVel::identity()));
    EXPECT_TRUE((inv * g).isApprox(PoseVel::identity()));
    
    // Verify component-wise inverse
    EXPECT_TRUE(inv.get<0>().isApprox(g.get<0>().inverse()));
    EXPECT_TRUE(inv.get<1>().isApprox(g.get<1>().inverse()));
}

/**
 * Test exponential and logarithm maps
 */
TEST_F(BundleTest, ExpLog)
{
    // Test PoseVel exp/log
    Vector3 axis = Vector3(1, 2, 3).normalized();
    Vector3 trans(4, 5, 6);
    Vector3 vel(0.1, 0.2, 0.3);
    
    PoseVel g(SE3d(SO3d(pi/4, axis), trans), RealSpace3d(vel));
    auto xi = g.log();
    PoseVel g2 = PoseVel().exp(xi);
    
    EXPECT_TRUE(g.isApprox(g2));
    
    // Test MultiBody exp/log
    MultiBody mb(SE3d(SO3d(pi/3, axis), trans),
                SE3d(SO3d(pi/4, axis), trans),
                SE3d(SO3d(pi/6, axis), trans));
    
    auto xi_mb = mb.log();
    MultiBody mb2 = MultiBody().exp(xi_mb);
    
    EXPECT_TRUE(mb.isApprox(mb2));
}

/**
 * Test adjoint representation
 */
TEST_F(BundleTest, Adjoint)
{
    // Test PoseVel adjoint (should be block diagonal)
    Vector3 axis = Vector3(0, 0, 1);
    Vector3 trans(1, 2, 0);
    Vector3 vel(0.1, 0.2, 0);
    
    PoseVel g(SE3d(SO3d(pi/2, axis), trans), RealSpace3d(vel));
    auto Ad = g.adjoint();
    
    // Verify block diagonal structure
    auto Ad1 = g.get<0>().adjoint();
    auto Ad2 = g.get<1>().adjoint();
    
    EXPECT_TRUE(Ad.block(0, 0, Ad1.rows(), Ad1.cols()).isApprox(Ad1));
    EXPECT_TRUE(Ad.block(Ad1.rows(), Ad1.cols(), Ad2.rows(), Ad2.cols()).isApprox(Ad2));
}

/**
 * Test group action
 */
TEST_F(BundleTest, Action)
{
    // Test PoseVel action
    Vector3 axis = Vector3(0, 0, 1);
    Vector3 trans(1, 0, 0);
    Vector3 vel(0.1, 0, 0);
    
    PoseVel g(SE3d(SO3d(pi/2, axis), trans), RealSpace3d(vel));
    
    // Create point and velocity
    Vector3 p(1, 0, 0);
    Vector3 v(0.1, 0, 0);
    Eigen::VectorXd state(6);
    state << p, v;
    
    // Test action
    auto result = g.act(state);
    
    // Verify components transformed correctly
    Vector3 p_new = result.head<3>();
    Vector3 v_new = result.tail<3>();
    
    EXPECT_TRUE(p_new.isApprox(g.get<0>().act(p)));
    EXPECT_TRUE(v_new.isApprox(v + vel));
}

/**
 * Test interpolation
 */
TEST_F(BundleTest, Interpolation)
{
    // Test PoseVel interpolation
    Vector3 axis = Vector3(1, 1, 1).normalized();
    PoseVel start = PoseVel::identity();
    PoseVel end(SE3d(SO3d(pi/2, axis), Vector3(2, 0, 0)),
                RealSpace3d(Vector3(0.2, 0, 0)));
    
    // Test midpoint
    PoseVel mid = interpolate(start, end, 0.5);
    
    // Verify component-wise interpolation
    EXPECT_TRUE(mid.get<0>().isApprox(interpolate(start.get<0>(), end.get<0>(), 0.5)));
    EXPECT_TRUE(mid.get<1>().isApprox(interpolate(start.get<1>(), end.get<1>(), 0.5)));
    
    // Test boundary conditions
    EXPECT_TRUE(interpolate(start, end, 0.0).isApprox(start));
    EXPECT_TRUE(interpolate(start, end, 1.0).isApprox(end));
}

/**
 * Test mixed bundle operations
 */
TEST_F(BundleTest, ComplexSystem)
{
    // Create a complex system with different types
    Vector3 axis = Vector3(1, 1, 1).normalized();
    ComplexSystem sys(
        SE3d(SO3d(pi/4, axis), Vector3(1, 2, 3)),    // Rigid body
        SO3d(pi/3, Vector3(0, 0, 1)),                // Pure rotation
        RealSpace3d(Vector3(0.1, 0.2, 0.3))          // Vector space
    );
    
    // Test identity composition
    EXPECT_TRUE((sys * ComplexSystem::identity()).isApprox(sys));
    
    // Test inverse
    auto inv = sys.inverse();
    EXPECT_TRUE((sys * inv).isApprox(ComplexSystem::identity()));
    
    // Test exp/log
    auto xi = sys.log();
    auto sys2 = ComplexSystem().exp(xi);
    EXPECT_TRUE(sys.isApprox(sys2));
}

} // namespace sofa::component::cosserat::liegroups::testing
