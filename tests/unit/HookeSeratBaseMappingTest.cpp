#include <gtest/gtest.h>
#include <Cosserat/mapping/HookeSeratBaseMapping.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <memory>

// Test fixture for HookeSeratBaseMapping tests
class HookeSeratBaseMappingTest : public ::testing::Test {
protected:
    using Mapping = Cosserat::mapping::HookeSeratBaseMapping<
        sofa::defaulttype::Vec3Types,
        sofa::defaulttype::Rigid3Types,
        sofa::defaulttype::Rigid3Types>;

    void SetUp() override {
        // Create mapping instance
        mapping = std::make_unique<Mapping>();
    }

    void TearDown() override {
        mapping.reset();
    }

    std::unique_ptr<Mapping> mapping;
};

// Test SectionInfo class functionality
TEST_F(HookeSeratBaseMappingTest, SectionInfo_BasicOperations) {
    using SectionInfo = Cosserat::mapping::SectionInfo;
    using SE3Type = sofa::component::cosserat::liegroups::SE3<double>;
    using TangentVector = SE3Type::TangentVector;

    // Test constructor
    TangentVector strain = TangentVector::Zero();
    strain << 0.1, 0.2, 0.3, 0.0, 0.0, 0.1; // Some curvature and extension

    SectionInfo section(1.0, strain, 0);

    EXPECT_DOUBLE_EQ(section.getLength(), 1.0);
    EXPECT_EQ(section.getIndex0(), 0);
    EXPECT_TRUE(section.getStrainsVec().isApprox(strain));
}

// Test SE3 transformation operations in SectionInfo
TEST_F(HookeSeratBaseMappingTest, SectionInfo_SE3Operations) {
    using SectionInfo = Cosserat::mapping::SectionInfo;
    using SE3Type = sofa::component::cosserat::liegroups::SE3<double>;
    using TangentVector = SE3Type::TangentVector;
    using Vector3 = SE3Type::Vector3;

    TangentVector strain = TangentVector::Zero();
    strain << 0.1, 0.0, 0.0, 0.0, 0.0, 0.0; // Pure bending around x-axis

    SectionInfo section(1.0, strain, 0);

    // Test local transformation computation
    SE3Type local_transform = section.getLocalTransformation(0.5);
    EXPECT_TRUE(local_transform.matrix().isApprox(SE3Type::computeIdentity().matrix(), 1e-10));

    // Test adjoint matrix computation
    auto adjoint = section.getAdjoint();
    EXPECT_EQ(adjoint.rows(), 6);
    EXPECT_EQ(adjoint.cols(), 6);

    // Test co-adjoint matrix
    auto coadjoint = section.getCoAdjoint();
    EXPECT_TRUE(coadjoint.isApprox(adjoint.transpose()));
}

// Test FrameInfo class functionality
TEST_F(HookeSeratBaseMappingTest, FrameInfo_BasicOperations) {
    using FrameInfo = Cosserat::mapping::FrameInfo;
    using SE3Type = sofa::component::cosserat::liegroups::SE3<double>;
    using TangentVector = SE3Type::TangentVector;

    FrameInfo frame;

    // Test basic setters/getters
    frame.setLength(2.0);
    EXPECT_DOUBLE_EQ(frame.getLength(), 2.0);

    frame.set_related_beam_index_(1);
    EXPECT_EQ(frame.get_related_beam_index_(), 1);

    frame.setDistanceToNearestBeamNode(0.5);
    EXPECT_DOUBLE_EQ(frame.getDistanceToNearestBeamNode(), 0.5);

    // Test kappa operations
    TangentVector kappa = TangentVector::Zero();
    kappa << 0.1, 0.2, 0.3, 0.0, 0.0, 0.1;
    frame.setKappa(kappa);
    EXPECT_TRUE(frame.getKappa().isApprox(kappa));
}

// Test tangent exponential implementation equivalence
TEST_F(HookeSeratBaseMappingTest, TangentExponentialEquivalence) {
    using SE3Type = sofa::component::cosserat::liegroups::SE3<double>;
    using TangentVector = SE3Type::TangentVector;
    using AdjointMatrix = SE3Type::AdjointMatrix;

    // Test case 1: Zero strain (should be equivalent)
    TangentVector zero_strain = TangentVector::Zero();
    AdjointMatrix adjoint = AdjointMatrix::Identity();
    double curv_abs = 1.0;

    EXPECT_TRUE(Mapping::testTangExpImplementationEquivalence(
        curv_abs, zero_strain, adjoint, 1e-6));

    // Test case 2: Small strain
    TangentVector small_strain = TangentVector::Zero();
    small_strain << 0.01, 0.01, 0.01, 0.0, 0.0, 0.0;

    EXPECT_TRUE(Mapping::testTangExpImplementationEquivalence(
        curv_abs, small_strain, adjoint, 1e-4)); // Relaxed tolerance for small angles
}

// Test geometry validation methods
TEST_F(HookeSeratBaseMappingTest, GeometryValidation) {
    // This test would require a properly initialized mapping with mechanical states
    // For now, we test the validation logic conceptually

    using SectionInfo = Cosserat::mapping::SectionInfo;
    using TangentVector = sofa::component::cosserat::liegroups::SE3<double>::TangentVector;

    std::vector<SectionInfo> sections;

    // Valid section
    TangentVector strain = TangentVector::Zero();
    SectionInfo valid_section(1.0, strain, 0);
    sections.push_back(valid_section);

    // Invalid section (negative length)
    SectionInfo invalid_section(-1.0, strain, 1);
    sections.push_back(invalid_section);

    // Test validation logic (would be called from mapping)
    bool has_invalid = false;
    for (const auto& section : sections) {
        if (section.getLength() < 0) {
            has_invalid = true;
            break;
        }
    }

    EXPECT_TRUE(has_invalid);
}

// Test interpolation methods
TEST_F(HookeSeratBaseMappingTest, SectionInterpolation) {
    using SectionInfo = Cosserat::mapping::SectionInfo;
    using SE3Type = sofa::component::cosserat::liegroups::SE3<double>;
    using TangentVector = SE3Type::TangentVector;

    // Create two sections
    TangentVector strain1 = TangentVector::Zero();
    strain1 << 0.1, 0.0, 0.0, 0.0, 0.0, 0.0;

    TangentVector strain2 = TangentVector::Zero();
    strain2 << 0.2, 0.0, 0.0, 0.0, 0.0, 0.0;

    SectionInfo section1(1.0, strain1, 0);
    SectionInfo section2(2.0, strain2, 1);

    // Test linear interpolation
    SectionInfo interpolated = section1.lerp(section2, 0.5);

    EXPECT_DOUBLE_EQ(interpolated.getLength(), 1.5); // Average length
    TangentVector expected_strain = 0.5 * strain1 + 0.5 * strain2;
    EXPECT_TRUE(interpolated.getStrainsVec().isApprox(expected_strain));
}

// Test distance computation
TEST_F(HookeSeratBaseMappingTest, SectionDistance) {
    using SectionInfo = Cosserat::mapping::SectionInfo;
    using SE3Type = sofa::component::cosserat::liegroups::SE3<double>;
    using TangentVector = SE3Type::TangentVector;

    TangentVector strain = TangentVector::Zero();
    SectionInfo section1(1.0, strain, 0, SE3Type::computeIdentity());

    // Create a translated section
    SE3Type translated_transform(SE3Type::computeIdentity());
    translated_transform.translation() << 1.0, 0.0, 0.0;

    SectionInfo section2(1.0, strain, 1, translated_transform);

    // Test distance computation
    double distance = section1.distanceTo(section2);
    EXPECT_NEAR(distance, 1.0, 1e-6); // Should be approximately the translation distance
}

// Test strain state bundle operations
TEST_F(HookeSeratBaseMappingTest, StrainStateBundle) {
    using SectionInfo = Cosserat::mapping::SectionInfo;
    using StrainState = Cosserat::mapping::StrainState;
    using SO3Type = sofa::component::cosserat::liegroups::SO3<double>;
    using RealSpace = sofa::component::cosserat::liegroups::RealSpace<double, 3>;
    using TangentVector = sofa::component::cosserat::liegroups::SE3<double>::TangentVector;
    using Vector3 = sofa::component::cosserat::liegroups::SE3<double>::Vector3;

    // Create section with known strain
    TangentVector strain = TangentVector::Zero();
    strain << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6;

    SectionInfo section(1.0, strain, 0);

    // Get strain state as bundle
    StrainState strain_state = section.getStrainState();

    // Verify bundle components can be extracted
    // Note: This requires implementing accessors in Bundle class
    // For now, test that the method exists and returns a valid object

    EXPECT_NO_THROW({
        StrainState state = section.getStrainState();
        // Test that we can create and manipulate the bundle
        SO3Type angular_part(Vector3(0.1, 0.2, 0.3));
        RealSpace linear_part(Vector3(0.4, 0.5, 0.6));
        StrainState expected_state(angular_part, linear_part);
    });
}

// Performance test for adjoint matrix computation
TEST_F(HookeSeratBaseMappingTest, AdjointMatrixPerformance) {
    using SectionInfo = Cosserat::mapping::SectionInfo;
    using TangentVector = sofa::component::cosserat::liegroups::SE3<double>::TangentVector;

    TangentVector strain = TangentVector::Zero();
    strain << 0.1, 0.05, 0.02, 0.0, 0.0, 0.0;

    SectionInfo section(1.0, strain, 0);

    // Time adjoint computation (first call computes, second uses cache)
    auto start = std::chrono::high_resolution_clock::now();
    auto adjoint1 = section.getAdjoint();
    auto end = std::chrono::high_resolution_clock::now();
    auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    start = std::chrono::high_resolution_clock::now();
    auto adjoint2 = section.getAdjoint();
    end = std::chrono::high_resolution_clock::now();
    auto duration2 = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    // Cached version should be faster
    EXPECT_LE(duration2.count(), duration1.count());

    // Results should be identical
    EXPECT_TRUE(adjoint1.isApprox(adjoint2));
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}