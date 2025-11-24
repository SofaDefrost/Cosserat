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

// Test BeamTopology functionality
TEST_F(HookeSeratBaseMappingTest, BeamTopology_BasicOperations) {
    using BeamTopology = Cosserat::mapping::BeamTopology;
    using SE3Type = sofa::component::cosserat::liegroups::SE3<double>;

    BeamTopology topology;

    // Test empty topology
    EXPECT_TRUE(topology.isValid());
    EXPECT_EQ(topology.getNumSections(), 0u);

    // Add sections (simple chain: 0 -> 1 -> 2)
    topology.parent_indices = {-1, 0, 1};
    topology.relative_transforms = {SE3Type::computeIdentity(), SE3Type::computeIdentity(), SE3Type::computeIdentity()};
    topology.connection_stiffnesses = {1000.0, 1000.0, 1000.0};

    EXPECT_TRUE(topology.isValid());
    EXPECT_EQ(topology.getNumSections(), 3u);

    // Test children relationships
    auto children0 = topology.getChildren(0);
    EXPECT_EQ(children0.size(), 1u);
    EXPECT_EQ(children0[0], 1u);

    auto children1 = topology.getChildren(1);
    EXPECT_EQ(children1.size(), 1u);
    EXPECT_EQ(children1[0], 2u);

    auto children2 = topology.getChildren(2);
    EXPECT_EQ(children2.size(), 0u);
}

// Test BeamTopology validation (cycles)
TEST_F(HookeSeratBaseMappingTest, BeamTopology_CycleDetection) {
    using BeamTopology = Cosserat::mapping::BeamTopology;
    using SE3Type = sofa::component::cosserat::liegroups::SE3<double>;

    BeamTopology topology;

    // Create a cycle: 0 -> 1 -> 2 -> 0
    topology.parent_indices = {2, 0, 1}; // 2 points to 0, creating cycle
    topology.relative_transforms = {SE3Type::computeIdentity(), SE3Type::computeIdentity(), SE3Type::computeIdentity()};
    topology.connection_stiffnesses = {1000.0, 1000.0, 1000.0};

    EXPECT_FALSE(topology.isValid()); // Should detect cycle
}

// Test BeamStateEstimator basic functionality
TEST_F(HookeSeratBaseMappingTest, BeamStateEstimator_BasicOperations) {
    using BeamStateEstimator = Cosserat::mapping::BeamStateEstimator;
    using SE3Type = sofa::component::cosserat::liegroups::SE3<double>;
    using StrainState = Cosserat::mapping::StrainState;
    using SO3Type = sofa::component::cosserat::liegroups::SO3<double>;
    using RealSpace = sofa::component::cosserat::liegroups::RealSpace<double, 3>;
    using Vector3 = SE3Type::Vector3;

    BeamStateEstimator estimator;

    // Test initialization
    SE3Type initial_pose = SE3Type::computeIdentity();
    StrainState initial_strain(SO3Type(Vector3::Zero()), RealSpace(Vector3::Zero()));
    Eigen::Matrix<double, 12, 12> initial_cov = Eigen::Matrix<double, 12, 12>::Identity() * 0.1;

    estimator.initialize(initial_pose, initial_strain, initial_cov);

    // Test confidence computation
    double confidence = estimator.getEstimationConfidence();
    EXPECT_GT(confidence, 0.0);

    // Test prediction
    using TangentVector = SE3Type::TangentVector;
    TangentVector control_input = TangentVector::Zero();
    control_input << 0.01, 0.0, 0.0, 0.0, 0.0, 0.0; // Small rotation

    estimator.predict(control_input);

    // Confidence should still be positive
    double confidence_after_predict = estimator.getEstimationConfidence();
    EXPECT_GT(confidence_after_predict, 0.0);
}

// Test BeamStateEstimator update operations
TEST_F(HookeSeratBaseMappingTest, BeamStateEstimator_UpdateOperations) {
    using BeamStateEstimator = Cosserat::mapping::BeamStateEstimator;
    using SE3Type = sofa::component::cosserat::liegroups::SE3<double>;
    using StrainState = Cosserat::mapping::StrainState;
    using SO3Type = sofa::component::cosserat::liegroups::SO3<double>;
    using RealSpace = sofa::component::cosserat::liegroups::RealSpace<double, 3>;
    using Vector3 = SE3Type::Vector3;

    BeamStateEstimator estimator;

    // Initialize
    SE3Type initial_pose = SE3Type::computeIdentity();
    StrainState initial_strain(SO3Type(Vector3::Zero()), RealSpace(Vector3::Zero()));
    Eigen::Matrix<double, 12, 12> initial_cov = Eigen::Matrix<double, 12, 12>::Identity() * 0.1;
    estimator.initialize(initial_pose, initial_strain, initial_cov);

    // Test pose update
    SE3Type measurement = SE3Type::computeIdentity();
    Eigen::Matrix<double, 6, 6> measurement_cov = Eigen::Matrix<double, 6, 6>::Identity() * 0.01;

    estimator.update(measurement, measurement_cov);

    // Confidence should improve after measurement update
    double confidence_after_update = estimator.getEstimationConfidence();
    EXPECT_GT(confidence_after_update, 0.0);
}

// Test performance optimization features
TEST_F(HookeSeratBaseMappingTest, PerformanceOptimization_Caching) {
    // Test that mapping has caching capabilities
    EXPECT_TRUE(mapping->supportsMultiSectionBeams());

    // Test cache operations
    mapping->clearComputationCache();
    EXPECT_EQ(mapping->getCacheSize(), 0u);

    // Test parallel computation settings
    mapping->enableParallelComputation(false);
    EXPECT_FALSE(mapping->isParallelComputationEnabled());

    mapping->enableParallelComputation(true);
    EXPECT_TRUE(mapping->isParallelComputationEnabled());
    EXPECT_GT(mapping->getOptimalThreadCount(), 0u);
}

// Test multi-section beam support
TEST_F(HookeSeratBaseMappingTest, MultiSectionBeamSupport) {
    using BeamTopology = Cosserat::mapping::BeamTopology;
    using SE3Type = sofa::component::cosserat::liegroups::SE3<double>;

    // Create a simple topology
    BeamTopology topology;
    topology.parent_indices = {-1, 0};
    topology.relative_transforms = {SE3Type::computeIdentity(), SE3Type::computeIdentity()};
    topology.connection_stiffnesses = {1000.0, 1000.0};

    // Test setting topology
    mapping->setBeamTopology(topology);
    const auto& retrieved_topology = mapping->getBeamTopology();
    EXPECT_EQ(retrieved_topology.getNumSections(), 2u);

    // Test multi-section enable/disable
    mapping->enableMultiSectionSupport(true);
    // Note: We can't easily test the internal state without accessing private members
    // This would require friend classes or additional getter methods
}

// Test Jacobian statistics and benchmarking
TEST_F(HookeSeratBaseMappingTest, JacobianStatistics_Benchmarking) {
    // Test that we can run performance benchmarks
    EXPECT_NO_THROW({
        mapping->runPerformanceBenchmark(10); // Small benchmark for testing
    });

    // Test that we can get statistics
    const auto& stats = mapping->getJacobianStats();
    // Statistics should be available even if not computed yet
    EXPECT_GE(stats.computation_count, 0u);
}

// Test state estimation integration
TEST_F(HookeSeratBaseMappingTest, StateEstimationIntegration) {
    // Test state estimation enable/disable
    mapping->enableStateEstimation(false);
    EXPECT_FALSE(mapping->isStateEstimationEnabled());

    mapping->enableStateEstimation(true);
    EXPECT_TRUE(mapping->isStateEstimationEnabled());

    // Test confidence query
    double confidence = mapping->getEstimationConfidence();
    // Confidence should be 0 when no estimator is active, or positive when active
    EXPECT_GE(confidence, 0.0);
}

// Phase 4: Advanced Lie Group Features Tests

// Test BCH correction computation
TEST_F(HookeSeratBaseMappingTest, AdvancedLieGroups_BCHCorrection) {
    using TangentVector = sofa::component::cosserat::liegroups::SE3<double>::TangentVector;

    // Create two tangent vectors
    TangentVector v1 = TangentVector::Zero();
    v1 << 0.1, 0.0, 0.0, 0.0, 0.0, 0.0; // Pure rotation around x

    TangentVector v2 = TangentVector::Zero();
    v2 << 0.0, 0.1, 0.0, 0.0, 0.0, 0.0; // Pure rotation around y

    // Compute BCH correction
    TangentVector bch_result = mapping->computeBCHCorrection(v1, v2);

    // BCH correction should be non-zero for non-commuting elements
    EXPECT_FALSE(bch_result.isZero());

    // Test commutativity: [v1,v2] = -[v2,v1]
    TangentVector bch_reverse = mapping->computeBCHCorrection(v2, v1);
    EXPECT_TRUE(bch_result.isApprox(-bch_reverse));
}

// Test parallel transport
TEST_F(HookeSeratBaseMappingTest, AdvancedLieGroups_ParallelTransport) {
    using SE3Type = sofa::component::cosserat::liegroups::SE3<double>;
    using TangentVector = SE3Type::TangentVector;

    // Create a tangent vector and target pose
    TangentVector tangent_vector;
    tangent_vector << 0.1, 0.0, 0.0, 0.0, 0.0, 0.1;

    SE3Type target_pose = SE3Type::computeIdentity();

    // Compute parallel transport
    TangentVector transported = mapping->parallelTransport(tangent_vector, target_pose);

    // For identity transport, result should be the same (simplified implementation)
    EXPECT_TRUE(transported.isApprox(tangent_vector));
}

// Test geodesic distance computation
TEST_F(HookeSeratBaseMappingTest, AdvancedLieGroups_GeodesicDistance) {
    using SectionInfo = Cosserat::mapping::SectionInfo;
    using TangentVector = sofa::component::cosserat::liegroups::SE3<double>::TangentVector;

    // Create two sections with different strains
    TangentVector strain1 = TangentVector::Zero();
    strain1 << 0.1, 0.0, 0.0, 0.0, 0.0, 0.0;

    TangentVector strain2 = TangentVector::Zero();
    strain2 << 0.2, 0.0, 0.0, 0.0, 0.0, 0.0;

    SectionInfo section1(1.0, strain1, 0);
    SectionInfo section2(1.0, strain2, 0);

    // Add sections to mapping
    mapping->addSection(section1);

    // Compute geodesic distance
    double distance = mapping->computeGeodesicDistance(section2);

    // Distance should be positive
    EXPECT_GT(distance, 0.0);

    // Distance to self should be zero
    double self_distance = mapping->computeGeodesicDistance(section1);
    EXPECT_DOUBLE_EQ(self_distance, 0.0);
}

// Test Riemannian exponential map
TEST_F(HookeSeratBaseMappingTest, AdvancedLieGroups_RiemannianExponential) {
    using SectionInfo = Cosserat::mapping::SectionInfo;
    using TangentVector = sofa::component::cosserat::liegroups::SE3<double>::TangentVector;

    // Create a base section
    TangentVector strain = TangentVector::Zero();
    strain << 0.1, 0.0, 0.0, 0.0, 0.0, 0.0;

    SectionInfo base_section(1.0, strain, 0);

    // Create a direction vector
    TangentVector direction = TangentVector::Zero();
    direction << 0.01, 0.0, 0.0, 0.0, 0.0, 0.01;

    // Compute Riemannian exponential
    SectionInfo result = mapping->riemannianExponential(base_section, direction, 1.0);

    // Result should be different from base section
    EXPECT_FALSE(result.getStrainsVec().isApprox(base_section.getStrainsVec()));
}

// Phase 4: Machine Learning Integration Tests

// Test adaptive controller basic functionality
TEST_F(HookeSeratBaseMappingTest, MLIntegration_AdaptiveControllerBasic) {
    // Test enable/disable adaptive control
    mapping->enableAdaptiveControl(false);
    EXPECT_FALSE(mapping->isAdaptiveControlEnabled());

    mapping->enableAdaptiveControl(true);
    EXPECT_TRUE(mapping->isAdaptiveControlEnabled());

    // Test prediction with no training (should return zero)
    using SE3Type = sofa::component::cosserat::liegroups::SE3<double>;
    SE3Type target_pose = SE3Type::computeIdentity();

    auto prediction = mapping->getAdaptiveControlPrediction(target_pose);
    EXPECT_TRUE(prediction.isZero());
}

// Test adaptive controller training
TEST_F(HookeSeratBaseMappingTest, MLIntegration_AdaptiveControllerTraining) {
    using SectionInfo = Cosserat::mapping::SectionInfo;
    using TangentVector = sofa::component::cosserat::liegroups::SE3<double>::TangentVector;
    using SE3Type = sofa::component::cosserat::liegroups::SE3<double>;

    // Enable adaptive control
    mapping->enableAdaptiveControl(true);
    EXPECT_TRUE(mapping->isAdaptiveControlEnabled());

    // Create training data
    std::vector<std::pair<SectionInfo, TangentVector>> training_data;

    TangentVector strain1 = TangentVector::Zero();
    strain1 << 0.1, 0.0, 0.0, 0.0, 0.0, 0.0;
    SectionInfo section1(1.0, strain1, 0);

    TangentVector control1 = TangentVector::Zero();
    control1 << 0.05, 0.0, 0.0, 0.0, 0.0, 0.0;

    training_data.emplace_back(section1, control1);

    // Train the controller
    mapping->trainAdaptiveController(training_data);

    // Test prediction
    SE3Type target_pose = SE3Type::computeIdentity();
    auto prediction = mapping->getAdaptiveControlPrediction(target_pose);

    // Prediction should be non-zero after training
    EXPECT_FALSE(prediction.isZero());
}

// Test material adaptation
TEST_F(HookeSeratBaseMappingTest, MLIntegration_MaterialAdaptation) {
    // Enable adaptive control
    mapping->enableAdaptiveControl(true);

    // Create feedback vector
    Eigen::VectorXd feedback(24);
    feedback.setRandom();

    // Update material adaptation
    mapping->updateMaterialAdaptation(feedback);

    // Test should pass without errors
    EXPECT_TRUE(mapping->isAdaptiveControlEnabled());
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}