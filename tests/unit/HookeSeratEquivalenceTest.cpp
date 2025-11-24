#include <gtest/gtest.h>
#include <Cosserat/mapping/HookeSeratBaseMapping.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <vector>
#include <chrono>
#include <random>

// Test fixture for HookeSerat equivalence testing
class HookeSeratEquivalenceTest : public ::testing::Test {
protected:
    using Mapping = Cosserat::mapping::HookeSeratBaseMapping<
        sofa::defaulttype::Vec3Types,
        sofa::defaulttype::Rigid3Types,
        sofa::defaulttype::Rigid3Types>;

    using SE3Type = sofa::component::cosserat::liegroups::SE3<double>;
    using TangentVector = SE3Type::TangentVector;
    using AdjointMatrix = SE3Type::AdjointMatrix;

    void SetUp() override {
        // Generate test data
        generateTestCases();
    }

    void generateTestCases() {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dist(-1.0, 1.0);

        // Generate various test strains
        test_strains = {
            TangentVector::Zero(), // Zero strain
            TangentVector::Random(), // Random strain
            TangentVector::UnitX(), // Unit X rotation
            TangentVector::UnitY(), // Unit Y rotation
            TangentVector::UnitZ(), // Unit Z rotation
        };

        // Add small strains
        for (int i = 0; i < 5; ++i) {
            TangentVector small_strain;
            for (int j = 0; j < 6; ++j) {
                small_strain[j] = dist(gen) * 0.1; // Small values
            }
            test_strains.push_back(small_strain);
        }

        // Generate test curvatures
        test_curvatures = {0.0, 0.1, 0.5, 1.0, 2.0, 5.0};

        // Generate test adjoint matrices
        for (int i = 0; i < 10; ++i) {
            AdjointMatrix adj = AdjointMatrix::Random();
            // Ensure it's a valid adjoint matrix (skew-symmetric rotation part)
            adj.block<3,3>(0,0) = (adj.block<3,3>(0,0) - adj.block<3,3>(0,0).transpose()) * 0.5;
            test_adjoint_matrices.push_back(adj);
        }
    }

    std::vector<TangentVector> test_strains;
    std::vector<double> test_curvatures;
    std::vector<AdjointMatrix> test_adjoint_matrices;
};

// Test equivalence between new and legacy tangent exponential implementations
TEST_F(HookeSeratEquivalenceTest, TangentExponentialEquivalence) {
    const double tolerance = 1e-6;

    for (const auto& strain : test_strains) {
        for (const auto& curvature : test_curvatures) {
            for (const auto& adjoint : test_adjoint_matrices) {
                EXPECT_TRUE(Mapping::testTangExpImplementationEquivalence(
                    curvature, strain, adjoint, tolerance))
                    << "Failed for strain: " << strain.transpose()
                    << ", curvature: " << curvature;
            }
        }
    }
}

// Test performance comparison between implementations
TEST_F(HookeSeratEquivalenceTest, PerformanceComparison) {
    const int num_iterations = 1000;
    const TangentVector test_strain = TangentVector::Random();
    const double test_curvature = 1.0;
    const AdjointMatrix test_adjoint = AdjointMatrix::Random();

    // Time new implementation
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < num_iterations; ++i) {
        AdjointMatrix result;
        Mapping::computeTangExpImplementation(test_curvature, test_strain, test_adjoint, result);
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto new_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    // Time legacy implementation
    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < num_iterations; ++i) {
        AdjointMatrix result;
        Mapping::computeTangExpImplementationLegacy(test_curvature, test_strain, test_adjoint, result);
    }
    end = std::chrono::high_resolution_clock::now();
    auto legacy_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    std::cout << "New implementation time: " << new_time.count() << " microseconds" << std::endl;
    std::cout << "Legacy implementation time: " << legacy_time.count() << " microseconds" << std::endl;
    std::cout << "Performance ratio (legacy/new): " << static_cast<double>(legacy_time.count()) / new_time.count() << std::endl;

    // New implementation should be faster
    EXPECT_LT(new_time.count(), legacy_time.count());
}

// Test numerical stability at edge cases
TEST_F(HookeSeratEquivalenceTest, NumericalStability) {
    const double tolerance = 1e-4; // Relaxed tolerance for edge cases

    // Test with very small strains (near identity)
    TangentVector tiny_strain = TangentVector::Zero();
    tiny_strain[0] = 1e-12;

    // Test with large curvatures
    std::vector<double> large_curvatures = {10.0, 50.0, 100.0};

    AdjointMatrix identity_adjoint = AdjointMatrix::Identity();

    for (double curvature : large_curvatures) {
        EXPECT_TRUE(Mapping::testTangExpImplementationEquivalence(
            curvature, tiny_strain, identity_adjoint, tolerance))
            << "Failed numerical stability test for curvature: " << curvature;
    }
}

// Test SE3 right Jacobian implementation directly
TEST_F(HookeSeratEquivalenceTest, SE3RightJacobian) {
    // Test that SE3::rightJacobian produces reasonable results
    TangentVector test_strain = TangentVector::UnitX() * 0.1; // Small rotation around X

    AdjointMatrix jacobian = SE3Type::rightJacobian(test_strain);

    // Jacobian should be close to identity for small strains
    AdjointMatrix identity = AdjointMatrix::Identity();
    AdjointMatrix error = jacobian - identity;
    double max_error = error.cwiseAbs().maxCoeff();

    EXPECT_LT(max_error, 0.01) << "Jacobian not close to identity for small strain";
}

// Test adjoint matrix properties
TEST_F(HookeSeratEquivalenceTest, AdjointMatrixProperties) {
    // Test that computed adjoint matrices have expected properties
    SE3Type test_transform = SE3Type::computeExp(TangentVector::Random() * 0.1);

    AdjointMatrix adjoint = test_transform.computeAdjoint();

    // Adjoint should be invertible
    EXPECT_TRUE(adjoint.determinant() != 0.0) << "Adjoint matrix should be invertible";

    // Adjoint of inverse should be inverse of adjoint
    SE3Type inverse_transform = test_transform.computeInverse();
    AdjointMatrix adjoint_inverse = inverse_transform.computeAdjoint();
    AdjointMatrix expected_inverse = adjoint.inverse();

    AdjointMatrix error = adjoint_inverse - expected_inverse;
    double max_error = error.cwiseAbs().maxCoeff();

    EXPECT_LT(max_error, 1e-10) << "Adjoint of inverse should equal inverse of adjoint";
}

// Test strain state bundle operations
TEST_F(HookeSeratEquivalenceTest, StrainStateBundleOperations) {
    using StrainState = Cosserat::mapping::StrainState;
    using SO3Type = sofa::component::cosserat::liegroups::SO3<double>;
    using RealSpace = sofa::component::cosserat::liegroups::RealSpace<double, 3>;
    using Vector3 = sofa::component::cosserat::liegroups::SE3<double>::Vector3;

    // Create test components
    Vector3 angular_strain(0.1, 0.2, 0.3);
    Vector3 linear_strain(0.0, 0.0, 0.1);

    SO3Type so3_component(angular_strain);
    RealSpace real_component(linear_strain);

    // Create bundle
    StrainState bundle(so3_component, real_component);

    // Test that we can extract components
    const auto& extracted_so3 = bundle.template get<0>();
    const auto& extracted_real = bundle.template get<1>();

    // Verify components
    TangentVector extracted_angular = extracted_so3.log();
    TangentVector extracted_linear = extracted_real.log();

    EXPECT_TRUE(extracted_angular.isApprox(angular_strain, 1e-10));
    EXPECT_TRUE(extracted_linear.isApprox(linear_strain, 1e-10));
}

// Test SectionInfo strain state operations
TEST_F(HookeSeratEquivalenceTest, SectionInfoStrainState) {
    using SectionInfo = Cosserat::mapping::SectionInfo;
    using StrainState = Cosserat::mapping::StrainState;
    using TangentVector = sofa::component::cosserat::liegroups::SE3<double>::TangentVector;

    // Create section with known strain
    TangentVector strain = TangentVector::Zero();
    strain << 0.1, 0.2, 0.3, 0.0, 0.0, 0.1;

    SectionInfo section(1.0, strain, 0);

    // Get strain state
    StrainState strain_state = section.getStrainState();

    // Create new section and set strain state
    SectionInfo new_section(1.0, TangentVector::Zero(), 0);
    new_section.setStrainState(strain_state);

    // Verify that strains match
    EXPECT_TRUE(section.getStrainsVec().isApprox(new_section.getStrainsVec(), 1e-10));
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}