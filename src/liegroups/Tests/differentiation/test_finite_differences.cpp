/**
 * @file test_finite_differences.cpp
 * @brief Test suite for validating analytical Jacobians using finite differences
 * 
 * This file contains tests that compare analytically computed Jacobians
 * with numerically computed ones using finite differences.
 */

#include <gtest/gtest.h>
#include "DifferentiationTestUtils.h"
#include "../../SO3.h"
#include "../../SE3.h"
#include "../../SO2.h"
#include "../../SE2.h"

using namespace sofa::component::cosserat::liegroups;
using namespace sofa::component::cosserat::liegroups::testing;

/**
 * @brief Test fixture for differentiation tests
 */
class DifferentiationTest : public ::testing::Test {
protected:
    using Scalar = double;
    using TestUtils = DifferentiationTestUtils<Scalar>;
    
    // Tolerance for numerical comparisons
    static constexpr Scalar tolerance = 1e-5;
    
    void SetUp() override {
        // Seed random number generator for reproducibility
        std::srand(42);
    }
};

// =============================================================================
// SO3 Tests
// =============================================================================

TEST_F(DifferentiationTest, SO3_ExpJacobian) {
    using SO3d = SO3<Scalar>;
    using Vector3 = typename SO3d::TangentVector;
    
    // Test at several points in the tangent space
    std::vector<Vector3> test_points = {
        Vector3(0.1, 0.2, 0.3),
        Vector3(1.0, 0.5, -0.3),
        Vector3(-0.5, 0.8, 0.2),
        Vector3(0.0, 0.0, 0.1),  // Near identity
    };
    
    for (const auto& omega : test_points) {
        // Function: exp map from tangent space to SO3
        auto exp_func = [](const Vector3& w) -> Vector3 {
            return SO3d::exp(w).log();  // Round trip should be identity
        };
        
        // Analytical Jacobian at omega (should be close to identity for small omega)
        // For now, we test if the numerical derivative matches a known pattern
        auto numerical_jacobian = TestUtils::centralDifferenceJacobian<3, 3>(
            exp_func, omega
        );
        
        // The Jacobian of exp->log composition should be close to identity
        Eigen::Matrix<Scalar, 3, 3> expected = Eigen::Matrix<Scalar, 3, 3>::Identity();
        
        bool passed = TestUtils::compareMatrices<3, 3>(
            expected, numerical_jacobian, tolerance, false
        );
        
        if (!passed) {
            std::cout << "Failed at omega = " << omega.transpose() << std::endl;
        }
        
        EXPECT_TRUE(passed);
    }
}

TEST_F(DifferentiationTest, SO3_LogJacobian) {
    using SO3d = SO3<Scalar>;
    using Vector3 = typename SO3d::TangentVector;
    
    // Test log Jacobian
    std::vector<Vector3> test_angles = {
        Vector3(0.1, 0.2, 0.3),
        Vector3(0.5, -0.3, 0.4),
    };
    
    for (const auto& omega : test_angles) {
        SO3d R = SO3d::exp(omega);
        
        // Function: Variation in log around R
        auto log_func = [&R](const Vector3& delta) -> Vector3 {
            SO3d R_delta = SO3d::exp(delta) * R;
            return R_delta.log();
        };
        
        Vector3 zero = Vector3::Zero();
        auto numerical_jacobian = TestUtils::centralDifferenceJacobian<3, 3>(
            log_func, zero
        );
        
        // Near identity perturbation, Jacobian should be close to identity
        Eigen::Matrix<Scalar, 3, 3> expected = Eigen::Matrix<Scalar, 3, 3>::Identity();
        
        bool passed = TestUtils::compareMatrices<3, 3>(
            expected, numerical_jacobian, tolerance, false
        );
        
        EXPECT_TRUE(passed);
    }
}

TEST_F(DifferentiationTest, SO3_ActionJacobian) {
    using SO3d = SO3<Scalar>;
    using Vector3 = typename SO3d::TangentVector;
    
    // Test action Jacobian: R * p
    Vector3 omega(0.5, 0.3, -0.2);
    SO3d R = SO3d::exp(omega);
    Vector3 p(1.0, 2.0, 3.0);
    
    // Function: How action changes with perturbation in tangent space
    auto action_func = [&R, &p](const Vector3& delta) -> Vector3 {
        SO3d R_delta = SO3d::exp(delta) * R;
        return R_delta.act(p);
    };
    
    Vector3 zero = Vector3::Zero();
    auto numerical_jacobian = TestUtils::centralDifferenceJacobian<3, 3>(
        action_func, zero
    );
    
    // Expected: [omega]× * R * p (cross product form)
    Vector3 Rp = R.act(p);
    Eigen::Matrix<Scalar, 3, 3> expected = -SO3d::buildAntisymmetric(Rp);
    
    bool passed = TestUtils::compareMatrices<3, 3>(
        expected, numerical_jacobian, tolerance * 10, false  // Slightly relaxed tolerance
    );
    
    EXPECT_TRUE(passed);
}

// =============================================================================
// SO2 Tests
// =============================================================================

TEST_F(DifferentiationTest, SO2_ExpLogConsistency) {
    using SO2d = SO2<Scalar>;
    using TangentVector = typename SO2d::TangentVector;
    
    std::vector<Scalar> test_angles = {0.1, 0.5, 1.0, M_PI/2, -M_PI/4};
    
    for (Scalar theta : test_angles) {
        TangentVector tangent;
        tangent[0] = theta;
        
        // Function: exp->log should be identity
        auto round_trip = [](const TangentVector& v) -> TangentVector {
            return SO2d::exp(v).log();
        };
        
        auto numerical_jacobian = TestUtils::centralDifferenceJacobian<1, 1>(
            round_trip, tangent
        );
        
        Eigen::Matrix<Scalar, 1, 1> expected;
        expected(0, 0) = Scalar(1);
        
        bool passed = TestUtils::compareMatrices<1, 1>(
            expected, numerical_jacobian, tolerance, false
        );
        
        EXPECT_TRUE(passed);
    }
}

// =============================================================================
// SE3 Tests
// =============================================================================

TEST_F(DifferentiationTest, SE3_ExpLogConsistency) {
    using SE3d = SE3<Scalar>;
    using TangentVector = typename SE3d::TangentVector;
    
    TangentVector xi;
    xi << 0.1, 0.2, 0.15,  // translation
          0.3, -0.2, 0.1;   // rotation
    
    // Function: exp->log round trip
    auto round_trip = [](const TangentVector& v) -> TangentVector {
        return SE3d::exp(v).log();
    };
    
    auto numerical_jacobian = TestUtils::centralDifferenceJacobian<6, 6>(
        round_trip, xi
    );
    
    // Should be close to identity
    Eigen::Matrix<Scalar, 6, 6> expected = Eigen::Matrix<Scalar, 6, 6>::Identity();
    
    bool passed = TestUtils::compareMatrices<6, 6>(
        expected, numerical_jacobian, tolerance, false
    );
    
    EXPECT_TRUE(passed);
}

TEST_F(DifferentiationTest, SE3_TranslationExtraction) {
    using SE3d = SE3<Scalar>;
    using TangentVector = typename SE3d::TangentVector;
    using Vector3 = Eigen::Matrix<Scalar, 3, 1>;
    
    TangentVector xi;
    xi << 0.5, 1.0, 0.3,   // translation
          0.2, -0.1, 0.15;  // rotation
    
    // Function: Extract translation component
    auto translation_func = [](const TangentVector& v) -> Vector3 {
        SE3d g = SE3d::exp(v);
        return g.translation();
    };
    
    auto numerical_jacobian = TestUtils::centralDifferenceJacobian<3, 6>(
        translation_func, xi
    );
    
    // We expect specific structure but this is a sanity check
    EXPECT_GT(numerical_jacobian.norm(), 0.0);
    EXPECT_LT(numerical_jacobian.norm(), 100.0);  // Reasonable magnitude
}

// =============================================================================
// Utility Tests
// =============================================================================

TEST_F(DifferentiationTest, FiniteDifferenceAccuracy) {
    // Test that central differences are more accurate than forward differences
    
    auto test_func = [](const Eigen::Matrix<Scalar, 1, 1>& x) -> Scalar {
        return std::sin(x(0)) * std::exp(-x(0) * x(0));
    };
    
    Eigen::Matrix<Scalar, 1, 1> x;
    x(0) = 0.5;
    
    // Analytical derivative: cos(x)*exp(-x^2) - 2*x*sin(x)*exp(-x^2)
    Scalar analytical = std::cos(x(0)) * std::exp(-x(0) * x(0)) 
                       - 2 * x(0) * std::sin(x(0)) * std::exp(-x(0) * x(0));
    
    auto grad_forward = TestUtils::finiteDifferenceGradient<1>(test_func, x);
    auto grad_central = TestUtils::centralDifferenceGradient<1>(test_func, x);
    
    Scalar error_forward = std::abs(grad_forward(0) - analytical);
    Scalar error_central = std::abs(grad_central(0) - analytical);
    
    // Central differences should be significantly more accurate
    EXPECT_LT(error_central, error_forward * 0.01);
    EXPECT_LT(error_central, 1e-8);
}

// =============================================================================
// Main
// =============================================================================

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
