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
#include <Cosserat/liegroups/SE23.h>
#include <Cosserat/liegroups/SGal3.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
#include <cmath>

namespace sofa::component::cosserat::liegroups::testing {

using namespace sofa::testing;

/**
 * Integration tests demonstrating real-world applications of Lie groups
 */
class LieGroupIntegrationTest : public BaseTest
{
protected:
    // Define common types
    using Vector3 = Eigen::Vector3d;
    using Vector6 = Eigen::Matrix<double, 6, 1>;
    using Matrix3 = Eigen::Matrix3d;
    using Matrix4 = Eigen::Matrix4d;
    using Quaternion = Eigen::Quaterniond;

    // Basic Lie groups
    using SO3d = SO3<double>;
    using SE3d = SE3<double>;
    using SE23d = SE23<double>;
    using SGal3d = SGal3<double>;
    using RealSpace3d = RealSpace<double, 3>;

    // Specialized bundle types for Cosserat mechanics
    using CosseratSection = Bundle<SE3d, RealSpace3d>;  // Position, strain
    using CosseratNode = Bundle<SE23d, RealSpace3d>;    // Configuration with velocity
    using CosseratRod = Bundle<SE3d, RealSpace3d, RealSpace3d>; // Position, strain, stress

    const double pi = M_PI;
    const double eps = 1e-10;

    /**
     * Helper to create a rotation from axis-angle
     */
    SO3d makeRotation(double angle, const Vector3& axis) {
        return SO3d(angle, axis.normalized());
    }

    /**
     * Helper to create a material frame
     */
    SE3d makeMaterialFrame(const Vector3& position, const Vector3& tangent,
                          const Vector3& normal) {
        // Create rotation that aligns z-axis with tangent
        Vector3 z_axis(0, 0, 1);
        Vector3 rot_axis = z_axis.cross(tangent);
        double rot_angle = std::acos(z_axis.dot(tangent));
        SO3d R(rot_angle, rot_axis.normalized());

        // Additional rotation around tangent to align normal
        Vector3 current_normal = R.act(Vector3(1, 0, 0));
        Vector3 target_normal = normal - normal.dot(tangent) * tangent;
        target_normal.normalize();
        
        double twist_angle = std::acos(current_normal.dot(target_normal));
        if (std::abs(twist_angle) > eps) {
            SO3d twist(twist_angle, tangent);
            R = twist * R;
        }

        return SE3d(R, position);
    }

    void SetUp() override {}
    void TearDown() override {}
};

/**
 * Test Cosserat rod centerline representation
 */
TEST_F(LieGroupIntegrationTest, CosseratCenterline)
{
    // Create a circular arc centerline
    const int num_points = 50;
    const double radius = 1.0;
    const double arc_angle = pi;  // Half circle
    
    std::vector<SE3d> frames;
    frames.reserve(num_points);

    for (int i = 0; i < num_points; ++i) {
        double t = static_cast<double>(i) / (num_points - 1);
        double theta = t * arc_angle;
        
        // Position on circle
        Vector3 position(radius * std::cos(theta), radius * std::sin(theta), 0);
        
        // Tangent vector (normalized derivative)
        Vector3 tangent(-std::sin(theta), std::cos(theta), 0);
        
        // Normal vector (points inward)
        Vector3 normal(-std::cos(theta), -std::sin(theta), 0);
        
        frames.push_back(makeMaterialFrame(position, tangent, normal));
    }

    // Verify geometric properties
    for (size_t i = 1; i < frames.size(); ++i) {
        // Get relative transform between consecutive frames
        SE3d rel = frames[i-1].inverse() * frames[i];
        
        // Extract curvature from relative rotation
        Vector3 omega = rel.rotation().log();
        double curvature = omega.norm() * (num_points - 1) / arc_angle;
        
        // Should be constant and equal to 1/radius
        EXPECT_NEAR(curvature, 1.0/radius, 0.1);
    }
}

/**
 * Test Cosserat rod dynamics
 */
TEST_F(LieGroupIntegrationTest, CosseratDynamics)
{
    // Create a rod configuration with velocity
    Vector3 position(0, 0, 0);
    Vector3 tangent(0, 0, 1);
    Vector3 normal(1, 0, 0);
    
    // Initial material frame
    SE3d frame = makeMaterialFrame(position, tangent, normal);
    
    // Initial velocity state (linear and angular velocity)
    Vector3 linear_vel(0.1, 0, 0);    // Moving in x direction
    Vector3 angular_vel(0, 0.5, 0);    // Rotating around y axis
    
    // Create extended pose with velocity
    SE23d state(frame, linear_vel);
    
    // Initial strain and stress
    Vector3 strain(0, 0, 1);   // Unit stretch
    Vector3 stress(0, 0, 0);   // No initial stress
    
    // Create full Cosserat state
    CosseratRod rod(frame, strain, stress);
    
    // Simulate simple time evolution
    const double dt = 0.01;
    const int steps = 100;
    
    for (int i = 0; i < steps; ++i) {
        // Update position and orientation using current velocity
        Vector6 twist;
        twist << linear_vel, angular_vel;
        frame = frame * SE3d().exp(dt * twist);
        
        // Update strain based on deformation
        strain = frame.rotation().inverse().act(tangent);
        
        // Simple linear stress response
        stress = 100.0 * (strain - Vector3(0, 0, 1));  // Hook's law
        
        // Store new state
        rod = CosseratRod(frame, strain, stress);
        
        // Verify physical constraints
        EXPECT_NEAR(strain.norm(), 1.0, 0.1);  // Length preservation
        EXPECT_TRUE(frame.rotation().matrix().determinant() > 0);  // Proper rotation
    }
}

/**
 * Test multi-body articulation
 */
TEST_F(LieGroupIntegrationTest, ArticulatedSystem)
{
    // Create a simple 3-link articulated system
    using ThreeLink = Bundle<SE3d, SE3d, SE3d>;
    
    // Initial configuration (vertical stack)
    std::vector<SE3d> links;
    const double link_length = 1.0;
    
    for (int i = 0; i < 3; ++i) {
        Vector3 position(0, 0, i * link_length);
        SO3d orientation = SO3d::identity();
        links.push_back(SE3d(orientation, position));
    }
    
    ThreeLink system(links[0], links[1], links[2]);
    
    // Apply joint motions
    const double joint_angle = pi/4;  // 45 degrees
    
    // Rotate first joint around y-axis
    SO3d R1(joint_angle, Vector3(0, 1, 0));
    links[0] = SE3d(R1, links[0].translation());
    
    // Rotate second joint around x-axis
    SO3d R2(joint_angle, Vector3(1, 0, 0));
    SE3d T1 = links[0];  // Transform from first link
    links[1] = T1 * SE3d(R2, Vector3(0, 0, link_length));
    
    // Rotate third joint around y-axis
    SO3d R3(-joint_angle, Vector3(0, 1, 0));
    SE3d T2 = links[1];  // Transform from second link
    links[2] = T2 * SE3d(R3, Vector3(0, 0, link_length));
    
    // Update system state
    system = ThreeLink(links[0], links[1], links[2]);
    
    // Verify kinematic chain properties
    for (int i = 1; i < 3; ++i) {
        // Check link connections
        Vector3 parent_end = links[i-1].translation() + 
                           links[i-1].rotation().act(Vector3(0, 0, link_length));
        Vector3 child_start = links[i].translation();
        
        EXPECT_TRUE((parent_end - child_start).norm() < eps);
    }
}

/**
 * Test time-varying trajectories
 */
TEST_F(LieGroupIntegrationTest, TimeVaryingTrajectory)
{
    // Create a helix trajectory using SGal(3)
    const double radius = 1.0;
    const double pitch = 0.5;
    const double angular_vel = 1.0;
    const double vertical_vel = pitch * angular_vel;
    
    std::vector<SGal3d> trajectory;
    const int num_points = 50;
    
    for (int i = 0; i < num_points; ++i) {
        double t = static_cast<double>(i) / (num_points - 1);
        double theta = t * 2 * pi;
        
        // Position on helix
        Vector3 position(radius * std::cos(theta),
                        radius * std::sin(theta),
                        pitch * theta);
                        
        // Tangent vector
        Vector3 tangent(-radius * std::sin(theta),
                        radius * std::cos(theta),
                        pitch);
        tangent.normalize();
        
        // Create frame from position and tangent
        SE3d frame = makeMaterialFrame(position, tangent, Vector3(0, 0, 1));
        
        // Velocity components
        Vector3 linear_vel(-radius * angular_vel * std::sin(theta),
                          radius * angular_vel * std::cos(theta),
                          vertical_vel);
                          
        // Create Galilean transformation
        trajectory.push_back(SGal3d(frame, linear_vel, t));
    }
    
    // Verify trajectory properties
    for (size_t i = 1; i < trajectory.size(); ++i) {
        const auto& g1 = trajectory[i-1];
        const auto& g2 = trajectory[i];
        
        // Time should increase monotonically
        EXPECT_GT(g2.time(), g1.time());
        
        // Velocity should be constant in magnitude
        double vel_mag1 = g1.velocity().norm();
        double vel_mag2 = g2.velocity().norm();
        EXPECT_NEAR(vel_mag1, vel_mag2, eps);
        
        // Position should follow helix equation
        Vector3 pos = g1.pose().translation();
        double theta = std::atan2(pos.y(), pos.x());
        double expected_z = pitch * theta;
        EXPECT_NEAR(pos.z(), expected_z, 0.1);
    }
}

} // namespace sofa::component::cosserat::liegroups::testing
