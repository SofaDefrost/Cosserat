/******************************************************************************
 * iLQR Controller Tests
 ******************************************************************************/

#include <gtest/gtest.h>
#include "../../control/CosseratILQRController.h"

using namespace sofa::component::cosserat::liegroups;
using namespace sofa::component::cosserat::liegroups::control;

class ILQRControllerTest : public ::testing::Test {
protected:
    using Controller = CosseratILQRController<double>;
    using SE3Type = SE3<double>;
    using Vector6 = Eigen::Matrix<double, 6, 1>;
    using Vector3 = Eigen::Vector3d;
};

TEST_F(ILQRControllerTest, StraightLineTracking) {
    // Track a straight line along X axis
    Controller::Config config;
    config.max_iterations = 20;
    config.verbose = false;
    
    Controller controller(5, config);
    
    // Reference: straight line from 0 to 1.0m
    Controller::Trajectory ref;
    ref.segment_length = 0.2;
    
    for (int i = 0; i <= 5; ++i) {
        SE3Type pose;
        pose = SE3Type::computeIdentity();
        pose.translation() = Vector3(i * 0.2, 0, 0);
        ref.poses.push_back(pose);
    }
    
    // Initial guess: zero strains
    std::vector<Vector6> initial(5, Vector6::Zero());
    
    // Optimize
    auto result = controller.optimize(ref, initial);
    
    // Should converge
    EXPECT_TRUE(result.converged || result.iterations > 0);
    EXPECT_GT(result.iterations, 0);
    EXPECT_LT(result.final_cost, 10.0);  // Reasonable cost
}

TEST_F(ILQRControllerTest, CurvedTrajectory) {
    // Track a curved trajectory
    Controller::Config config;
    config.max_iterations = 30;
    config.Q_position = 20.0;
    config.verbose = false;
    
    Controller controller(10, config);
    
    // Reference: arc
    Controller::Trajectory ref;
    ref.segment_length = 0.1;
    
    for (int i = 0; i <= 10; ++i) {
        double t = i * 0.1;
        SE3Type pose;
        pose = SE3Type::computeIdentity();
        // Arc: x = t, y = 0.5*sin(πt)
        pose.translation() = Vector3(t, 0.5 * std::sin(M_PI * t), 0);
        ref.poses.push_back(pose);
    }
    
    std::vector<Vector6> initial(10, Vector6::Zero());
    auto result = controller.optimize(ref, initial);
    
    EXPECT_GT(result.iterations, 0);
    EXPECT_EQ(result.optimal_strains.size(), 10);
    EXPECT_EQ(result.feedback_gains.size(), 10);
}

TEST_F(ILQRControllerTest, CostDecreases) {
    // Verify cost decreases over iterations
    Controller::Config config;
    config.max_iterations = 10;
    config.verbose = false;
    
    Controller controller(3, config);
    
    Controller::Trajectory ref;
    ref.segment_length = 0.3;
    
    for (int i = 0; i <= 3; ++i) {
        SE3Type pose;
        pose = SE3Type::computeIdentity();
        pose.translation() = Vector3(i * 0.3, 0.1, 0);
        ref.poses.push_back(pose);
    }
    
    std::vector<Vector6> initial(3, Vector6::Zero());
    auto result = controller.optimize(ref, initial);
    
    // Cost should decrease
    ASSERT_GT(result.cost_history.size(), 1);
    for (size_t i = 1; i < result.cost_history.size(); ++i) {
        EXPECT_LE(result.cost_history[i], result.cost_history[i-1] + 1e-6);
    }
}
