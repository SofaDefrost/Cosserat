/**
 * @file test_jacobian_verification.cpp
 * @brief Step 1: Verify Jacobian Implementation
 *
 * This test verifies that computeTangExpImplementation() produces correct results
 * by comparing with numerical finite differences and checking mathematical properties.
 */

#include <Cosserat/mapping/CosseratGeometryMapping.h>
#include <gtest/gtest.h>
#include <liegroups/SE3.h>
#include <liegroups/SO3.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/VecTypes.h>

#include <Cosserat/mapping/CosseratGeometryMapping.inl>
#include <Eigen/Dense>
#include <cmath>
#include <iostream>

using namespace Cosserat::mapping;
using namespace sofa::component::cosserat::liegroups;

// Test fixture for Jacobian verification
class JacobianVerificationTest : public ::testing::Test {
   protected:
    using SE3Type = SE3<double>;
    using SO3Type = SO3<double>;
    using TangentVector = typename SE3Type::TangentVector;
    using AdjointMatrix = typename SE3Type::AdjointMatrix;
    using Vector3 = typename SE3Type::Vector3;
    using Matrix3 = typename SE3Type::Matrix3;

    // Use actual SOFA types for the mapping
    using MappingType =
        CosseratGeometryMapping<sofa::defaulttype::Vec3Types, sofa::defaulttype::Rigid3Types,
                                sofa::defaulttype::Rigid3Types>;

    void SetUp() override { tolerance_ = 1e-6; }

    // Helper to compute numerical Jacobian using finite differences
    AdjointMatrix computeNumericalJacobian(double curv_abs, const TangentVector &strain,
                                           double eps = 1e-7) {
        AdjointMatrix numerical_jac = AdjointMatrix::Zero();

        // For each dimension of the tangent space
        for (int i = 0; i < 6; ++i) {
            TangentVector strain_plus = strain;
            TangentVector strain_minus = strain;

            strain_plus[i] += eps;
            strain_minus[i] -= eps;

            // Compute SE3 transformations
            SE3Type g_plus = SE3Type::computeExp(strain_plus * curv_abs);
            SE3Type g_minus = SE3Type::computeExp(strain_minus * curv_abs);

            // Compute difference
            SE3Type g_diff = g_plus * g_minus.computeInverse();
            TangentVector diff = g_diff.computeLog();

            numerical_jac.col(i) = diff / (2.0 * eps);
        }

        return numerical_jac;
    }

    double tolerance_;
};

// Test 1: Zero strain - should be exactly curv_abs * I
TEST_F(JacobianVerificationTest, ZeroStrain) {
    double curv_abs = 1.0;
    TangentVector strain = TangentVector::Zero();

    AdjointMatrix tang_adjoint;
    MappingType::computeTangExpImplementation(curv_abs, strain, tang_adjoint);

    AdjointMatrix expected = curv_abs * AdjointMatrix::Identity();

    EXPECT_TRUE(tang_adjoint.isApprox(expected, tolerance_))
        << "Zero strain Jacobian should be curv_abs * I\n"
        << "Expected:\n"
        << expected << "\n"
        << "Got:\n"
        << tang_adjoint;
}

// Test 2: Small strain (near zero) - should use first-order approximation
TEST_F(JacobianVerificationTest, SmallStrain) {
    double curv_abs = 0.1;
    TangentVector strain;
    strain << 1e-7, 1e-7, 1e-7, 0.0, 0.0, 0.0;  // Very small angular strain

    AdjointMatrix tang_adjoint;
    MappingType::computeTangExpImplementation(curv_abs, strain, tang_adjoint);

    // For small strains, should be approximately curv_abs * I
    AdjointMatrix expected = curv_abs * AdjointMatrix::Identity();

    EXPECT_TRUE(tang_adjoint.isApprox(expected, 1e-4)) << "Small strain Jacobian failed\n"
                                                       << "Expected:\n"
                                                       << expected << "\n"
                                                       << "Got:\n"
                                                       << tang_adjoint << "\n"
                                                       << "Difference:\n"
                                                       << (tang_adjoint - expected);
}

// Test 3: Moderate strain - typical use case
TEST_F(JacobianVerificationTest, ModerateStrain) {
    double curv_abs = 0.5;
    TangentVector strain;
    strain << 0.1, 0.05, 0.02, 0.01, 0.01, 0.01;  // Moderate strain

    AdjointMatrix tang_adjoint;

    MappingType::computeTangExpImplementation(curv_abs, strain, tang_adjoint);

    // Verify it's not zero or identity
    EXPECT_FALSE(tang_adjoint.isZero(tolerance_));
    EXPECT_FALSE(tang_adjoint.isApprox(AdjointMatrix::Identity(), tolerance_));

    // Verify all elements are finite
    EXPECT_TRUE(tang_adjoint.allFinite()) << "Jacobian contains non-finite values";
}

// Test 4: Large strain - stress test
TEST_F(JacobianVerificationTest, LargeStrain) {
    double curv_abs = 1.0;
    TangentVector strain;
    strain << 1.0, 0.5, 0.3, 0.1, 0.1, 0.1;  // Large angular strain

    AdjointMatrix tang_adjoint;
    MappingType::computeTangExpImplementation(curv_abs, strain, tang_adjoint);

    // Should still be finite and non-zero
    EXPECT_TRUE(tang_adjoint.allFinite()) << "Large strain Jacobian contains non-finite values";
    EXPECT_FALSE(tang_adjoint.isZero(tolerance_));
}

// Test 5: Numerical accuracy vs finite differences
TEST_F(JacobianVerificationTest, NumericalAccuracy) {
    double curv_abs = 0.3;
    TangentVector strain;
    strain << 0.2, 0.1, 0.05, 0.02, 0.02, 0.02;

    // Analytical Jacobian
    AdjointMatrix analytical_jac;
    MappingType::computeTangExpImplementation(curv_abs, strain, analytical_jac);

    // Numerical Jacobian (finite differences)
    AdjointMatrix numerical_jac = computeNumericalJacobian(curv_abs, strain, 1e-7);

    // Compare
    double max_error = (analytical_jac - numerical_jac).cwiseAbs().maxCoeff();

    EXPECT_LT(max_error, 1e-4) << "Jacobian numerical accuracy test failed\n"
                               << "Max error: " << max_error << "\n"
                               << "Analytical:\n"
                               << analytical_jac << "\n"
                               << "Numerical:\n"
                               << numerical_jac << "\n"
                               << "Difference:\n"
                               << (analytical_jac - numerical_jac);
}

// Test 6: Consistency across multiple calls
TEST_F(JacobianVerificationTest, Consistency) {
    double curv_abs = 0.5;
    TangentVector strain;
    strain << 0.1, 0.05, 0.02, 0.01, 0.01, 0.01;

    AdjointMatrix tang_adjoint1, tang_adjoint2;

    // Call twice with same inputs
    MappingType::computeTangExpImplementation(curv_abs, strain, tang_adjoint1);

    MappingType::computeTangExpImplementation(curv_abs, strain, tang_adjoint2);

    // Results should be identical
    EXPECT_TRUE(tang_adjoint1.isApprox(tang_adjoint2, 1e-15))
        << "Jacobian computation is not deterministic";
}

// Test 7: Scaling properties
TEST_F(JacobianVerificationTest, ScalingProperties) {
    TangentVector strain;
    strain << 0.1, 0.05, 0.02, 0.01, 0.01, 0.01;

    double curv_abs1 = 0.5;
    double curv_abs2 = 1.0;

    AdjointMatrix tang_adjoint1, tang_adjoint2;

    MappingType::computeTangExpImplementation(curv_abs1, strain, tang_adjoint1);

    MappingType::computeTangExpImplementation(curv_abs2, strain, tang_adjoint2);

    // Both should be finite
    EXPECT_TRUE(tang_adjoint1.allFinite());
    EXPECT_TRUE(tang_adjoint2.allFinite());

    // They should be different (not a simple scaling)
    EXPECT_FALSE(tang_adjoint1.isApprox(tang_adjoint2, tolerance_));
}

// Test 8: Performance benchmark
TEST_F(JacobianVerificationTest, PerformanceBenchmark) {
    double curv_abs = 0.5;
    TangentVector strain;
    strain << 0.1, 0.05, 0.02, 0.01, 0.01, 0.01;

    AdjointMatrix tang_adjoint;

    const int iterations = 10000;
    auto start = std::chrono::high_resolution_clock::now();

    for (int i = 0; i < iterations; ++i) {
        MappingType::computeTangExpImplementation(curv_abs, strain, tang_adjoint);
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    double avg_time = duration.count() / static_cast<double>(iterations);

    std::cout << "Average Jacobian computation time: " << avg_time << " microseconds\n";

    // Expect it to be reasonably fast (< 10 microseconds per call)
    EXPECT_LT(avg_time, 10.0) << "Jacobian computation is too slow: " << avg_time << " μs";
}
