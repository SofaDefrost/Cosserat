/**
 * @file test_analytical_jacobians.cpp
 * @brief Test suite for validating newly implemented analytical Jacobians
 * 
 * This file tests the analytical Jacobians for:
 * - Group composition (compose)
 * - Group inverse (inverse)
 * - Group action (act)
 * 
 * All tests compare analytical jacobians with numerical finite differences.
 */

#include <gtest/gtest.h>
#include "DifferentiationTestUtils.h"
#include "../../SO3.h"
#include "../../SE3.h"

using namespace sofa::component::cosserat::liegroups;
using namespace sofa::component::cosserat::liegroups::testing;

class AnalyticalJacobianTest : public ::testing::Test {
protected:
    using Scalar = double;
    using TestUtils = DifferentiationTestUtils<Scalar>;
    
    static constexpr Scalar tolerance = 1e-5;
    static constexpr Scalar relaxed_tolerance = 1e-4;
    
    void SetUp() override {
        std::srand(42);
    }
};

// =============================================================================
// SO3 Jacobian Tests
// =============================================================================

TEST_F(AnalyticalJacobianTest, SO3_ComposeJacobians) {
    using SO3d = SO3<Scalar>;
    using Vector3 = typename SO3d::TangentVector;
    
    // Create two random rotations
    Vector3 omega1(0.5, 0.3, -0.2);
    Vector3 omega2(0.2, -0.4, 0.3);
    SO3d R = SO3d::exp(omega1);
    SO3d S = SO3d::exp(omega2);
    
    // Get analytical jacobians
    auto [J_left_analytical, J_right_analytical] = R.composeJacobians(S);
    
    // Test left Jacobian: ∂(R*S)/∂R
    auto compose_wrt_R = [&S](const Vector3& delta) -> Vector3 {
        SO3d R_perturbed = SO3d::exp(delta);
        return (R_perturbed * S).log();
    };
    
    Vector3 zero = Vector3::Zero();
    auto J_left_numerical = TestUtils::centralDifferenceJacobian<3, 3>(
        compose_wrt_R, zero
    );
    
    bool left_passed = TestUtils::compareMatrices<3, 3>(
        J_left_analytical, J_left_numerical, tolerance, false
    );
    
    EXPECT_TRUE(left_passed) << "Left Jacobian test failed\n"
        << "Analytical:\n" << J_left_analytical << "\n"
        << "Numerical:\n" << J_left_numerical;
    
    // Test right Jacobian: ∂(R*S)/∂S
    auto compose_wrt_S = [&R](const Vector3& delta) -> Vector3 {
        SO3d S_perturbed = SO3d::exp(delta);
        return (R * S_perturbed).log();
    };
    
    auto J_right_numerical = TestUtils::centralDifferenceJacobian<3, 3>(
        compose_wrt_S, zero
    );
    
    bool right_passed = TestUtils::compareMatrices<3, 3>(
        J_right_analytical, J_right_numerical, tolerance, false
    );
    
    EXPECT_TRUE(right_passed) << "Right Jacobian test failed";
}

TEST_F(AnalyticalJacobianTest, SO3_InverseJacobian) {
    using SO3d = SO3<Scalar>;
    using Vector3 = typename SO3d::TangentVector;
    
    Vector3 omega(0.5, 0.3, -0.2);
    SO3d R = SO3d::exp(omega);
    
    // Analytical Jacobian
    auto J_analytical = R.inverseJacobian();
    
    // Numerical Jacobian: ∂(R^{-1})/∂R
    auto inverse_func = [&R](const Vector3& delta) -> Vector3 {
        SO3d R_perturbed = SO3d::exp(delta) * R;
        return R_perturbed.inverse().log();
    };
    
    Vector3 zero = Vector3::Zero();
    auto J_numerical = TestUtils::centralDifferenceJacobian<3, 3>(
        inverse_func, zero
    );
    
    bool passed = TestUtils::compareMatrices<3, 3>(
        J_analytical, J_numerical, tolerance, false
    );
    
    EXPECT_TRUE(passed) << "Inverse Jacobian test failed\n"
        << "Analytical:\n" << J_analytical << "\n"
        << "Numerical:\n" << J_numerical;
}

TEST_F(AnalyticalJacobianTest, SO3_ActionJacobians) {
    using SO3d = SO3<Scalar>;
    using Vector3 = typename SO3d::TangentVector;
    
    Vector3 omega(0.3, 0.2, -0.1);
    SO3d R = SO3d::exp(omega);
    Vector3 p(1.0, 2.0, 3.0);
    
    // Get analytical jacobians
    auto [J_rot_analytical, J_point_analytical] = R.actionJacobians(p);
    
    // Test Jacobian w.r.t. rotation
    auto action_wrt_R = [&R, &p](const Vector3& delta) -> Vector3 {
        return (SO3d::exp(delta) * R).act(p);
    };
    
    Vector3 zero = Vector3::Zero();
    auto J_rot_numerical = TestUtils::centralDifferenceJacobian<3, 3>(
        action_wrt_R, zero
    );
    
    bool rot_passed = TestUtils::compareMatrices<3, 3>(
        J_rot_analytical, J_rot_numerical, relaxed_tolerance, false
    );
    
    EXPECT_TRUE(rot_passed) << "Action Jacobian w.r.t. rotation failed\n"
        << "Analytical:\n" << J_rot_analytical << "\n"
        << "Numerical:\n" << J_rot_numerical;
    
    // Test Jacobian w.r.t. point (should be rotation matrix)
    auto action_wrt_p = [&R](const Vector3& delta_p) -> Vector3 {
        return R.act(delta_p);
    };
    
    auto J_point_numerical = TestUtils::centralDifferenceJacobian<3, 3>(
        action_wrt_p, p
    );
    
    bool point_passed = TestUtils::compareMatrices<3, 3>(
        J_point_analytical, J_point_numerical, tolerance, false
    );
    
    EXPECT_TRUE(point_passed) << "Action Jacobian w.r.t. point failed";
}

// =============================================================================
// SE3 Jacobian Tests
// =============================================================================

TEST_F(AnalyticalJacobianTest, SE3_ComposeJacobians) {
    using SE3d = SE3<Scalar>;
    using Vector6 = typename SE3d::TangentVector;
    using Vector3 = Eigen::Matrix<Scalar, 3, 1>;
    
    // Create two random transformations
    Vector6 xi1;
    xi1 << 0.1, 0.2, 0.15, 0.3, -0.2, 0.1;
    Vector6 xi2;
    xi2 << -0.05, 0.15, 0.1, 0.2, 0.1, -0.15;
    
    SE3d g = SE3d::exp(xi1);
    SE3d h = SE3d::exp(xi2);
    
    // Get analytical jacobians
    auto [J_left_analytical, J_right_analytical] = g.composeJacobians(h);
    
    // Test left Jacobian: ∂(g*h)/∂g
    auto compose_wrt_g = [&h](const Vector6& delta) -> Vector6 {
        SE3d g_perturbed = SE3d::exp(delta);
        return (g_perturbed * h).log();
    };
    
    Vector6 zero = Vector6::Zero();
    auto J_left_numerical = TestUtils::centralDifferenceJacobian<6, 6>(
        compose_wrt_g, zero
    );
    
    bool left_passed = TestUtils::compareMatrices<6, 6>(
        J_left_analytical, J_left_numerical, tolerance, false
    );
    
    EXPECT_TRUE(left_passed) << "SE3 Left Jacobian test failed\n"
        << "Max error: " << (J_left_analytical - J_left_numerical).cwiseAbs().maxCoeff();
    
    // Test right Jacobian: ∂(g*h)/∂h
    auto compose_wrt_h = [&g](const Vector6& delta) -> Vector6 {
        SE3d h_perturbed = SE3d::exp(delta);
        return (g * h_perturbed).log();
    };
    
    auto J_right_numerical = TestUtils::centralDifferenceJacobian<6, 6>(
        compose_wrt_h, zero
    );
    
    bool right_passed = TestUtils::compareMatrices<6, 6>(
        J_right_analytical, J_right_numerical, tolerance, false
    );
    
    EXPECT_TRUE(right_passed) << "SE3 Right Jacobian test failed";
}

TEST_F(AnalyticalJacobianTest, SE3_InverseJacobian) {
    using SE3d = SE3<Scalar>;
    using Vector6 = typename SE3d::TangentVector;
    
    Vector6 xi;
    xi << 0.5, 1.0, 0.3, 0.2, -0.1, 0.15;
    SE3d g = SE3d::exp(xi);
    
    // Analytical Jacobian
    auto J_analytical = g.inverseJacobian();
    
    // Numerical Jacobian
    auto inverse_func = [&g](const Vector6& delta) -> Vector6 {
        SE3d g_perturbed = SE3d::exp(delta) * g;
        return g_perturbed.inverse().log();
    };
    
    Vector6 zero = Vector6::Zero();
    auto J_numerical = TestUtils::centralDifferenceJacobian<6, 6>(
        inverse_func, zero
    );
    
    bool passed = TestUtils::compareMatrices<6, 6>(
        J_analytical, J_numerical, tolerance, false
    );
    
    EXPECT_TRUE(passed) << "SE3 Inverse Jacobian test failed\n"
        << "Max error: " << (J_analytical - J_numerical).cwiseAbs().maxCoeff();
}

TEST_F(AnalyticalJacobianTest, SE3_ActionJacobians) {
    using SE3d = SE3<Scalar>;
    using Vector6 = typename SE3d::TangentVector;
    using Vector3 = typename SE3d::ActionVector;
    using Matrix3 = Eigen::Matrix<Scalar, 3, 3>;
    using Matrix3x6 = Eigen::Matrix<Scalar, 3, 6>;
    
    Vector6 xi;
    xi << 0.5, 1.0, 0.3, 0.2, -0.1, 0.15;
    SE3d g = SE3d::exp(xi);
    Vector3 p(1.0, 2.0, 3.0);
    
    // Get analytical jacobians
    auto [J_group_analytical, J_point_analytical] = g.actionJacobians(p);
    
    // Test Jacobian w.r.t. group element
    auto action_wrt_g = [&g, &p](const Vector6& delta) -> Vector3 {
        return (SE3d::exp(delta) * g).act(p);
    };
    
    Vector6 zero = Vector6::Zero();
    auto J_group_numerical = TestUtils::centralDifferenceJacobian<3, 6>(
        action_wrt_g, zero
    );
    
    bool group_passed = TestUtils::compareMatrices<3, 6>(
        J_group_analytical, J_group_numerical, relaxed_tolerance, false
    );
    
    EXPECT_TRUE(group_passed) << "SE3 Action Jacobian w.r.t. group failed\n"
        << "Max error: " << (J_group_analytical - J_group_numerical).cwiseAbs().maxCoeff();
    
    // Test Jacobian w.r.t. point
    auto action_wrt_p = [&g](const Vector3& delta_p) -> Vector3 {
        return g.act(delta_p);
    };
    
    auto J_point_numerical = TestUtils::centralDifferenceJacobian<3, 3>(
        action_wrt_p, p
    );
    
    bool point_passed = TestUtils::compareMatrices<3, 3>(
        J_point_analytical, J_point_numerical, tolerance, false
    );
    
    EXPECT_TRUE(point_passed) << "SE3 Action Jacobian w.r.t. point failed";
}

// =============================================================================
// Consistency Tests
// =============================================================================

TEST_F(AnalyticalJacobianTest, SO3_JacobiansConsistency) {
    using SO3d = SO3<Scalar>;
    using Vector3 = typename SO3d::TangentVector;
    
    // Test: (R^{-1})^{-1} = R, so jacobians should be inverses
    Vector3 omega(0.3, 0.2, -0.1);
    SO3d R = SO3d::exp(omega);
    
    auto J_inv = R.inverseJacobian();
    auto J_inv_inv = R.inverse().inverseJacobian();
    
    // These should be negatives: J(R^{-1}) = -J(R)
    auto product = J_inv * J_inv_inv;
    auto expected = -SO3d::AdjointMatrix::Identity();
    
    bool consistent = product.isApprox(expected, 0.01);
    EXPECT_TRUE(consistent) << "Inverse Jacobian consistency failed";
}

TEST_F(AnalyticalJacobianTest, SE3_ActionWithIdentity) {
    using SE3d = SE3<Scalar>;
    using Vector3 = typename SE3d::ActionVector;
    
    // Action by identity should have specific structure
    SE3d identity = SE3d::Identity();
    Vector3 p(1.0, 2.0, 3.0);
    
    auto [J_group, J_point] = identity.actionJacobians(p);
    
    // For identity transformation:
    // J_point should be identity matrix
    auto I3 = Eigen::Matrix<Scalar, 3, 3>::Identity();
    EXPECT_TRUE(J_point.isApprox(I3, tolerance)) 
        << "Identity action Jacobian w.r.t. point should be identity";
}
