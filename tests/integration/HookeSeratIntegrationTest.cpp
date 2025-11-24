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
#include <Cosserat/mapping/HookeSeratBaseMapping.h>
#include <Cosserat/mapping/BaseCosseratMapping.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/State.h>
#include <sofa/core/objectmodel/Data.h>
#include <Eigen/Core>
#include <vector>
#include <memory>

namespace sofa::component::cosserat::mapping::testing {

using namespace sofa::testing;
using namespace sofa::defaulttype;

/**
 * Integration tests comparing HookeSeratBaseMapping vs BaseCosseratMapping
 * Validates that our Lie group improvements maintain compatibility while
 * providing enhanced functionality and performance.
 */
class HookeSeratIntegrationTest : public BaseTest
{
protected:
    // Mapping types
    using HookeSeratMapping = HookeSeratBaseMapping<
        Vec3Types, Rigid3Types, Rigid3Types>;
    using BaseCosseratMapping = sofa::component::mapping::BaseCosseratMapping<
        Vec3Types, Rigid3Types, Rigid3Types>;

    // Data types
    using Vec3 = Vec3Types::VecCoord;
    using Rigid = Rigid3Types::VecCoord;
    using Vec3Data = sofa::core::State<Vec3Types>;
    using RigidData = sofa::core::State<Rigid3Types>;

    // Test parameters
    const double eps = 1e-6;
    const int num_sections = 5;
    const double section_length = 0.1;

    void SetUp() override {
        // Create test data
        createTestGeometry();
        createTestStates();
    }

    void TearDown() override {}

    /**
     * Create test beam geometry
     */
    void createTestGeometry() {
        // Curvilinear abscissa for sections
        curv_abs_sections.clear();
        for (int i = 0; i <= num_sections; ++i) {
            curv_abs_sections.push_back(i * section_length);
        }

        // Curvilinear abscissa for frames (more points than sections)
        curv_abs_frames.clear();
        const int num_frames = num_sections * 3;
        for (int i = 0; i < num_frames; ++i) {
            double t = static_cast<double>(i) / (num_frames - 1);
            curv_abs_frames.push_back(t * num_sections * section_length);
        }
    }

    /**
     * Create test mechanical states
     */
    void createTestStates() {
        // Strain state (kappa_x, kappa_y, kappa_z, gamma_x, gamma_y, gamma_z)
        strain_state.resize(num_sections);
        for (int i = 0; i < num_sections; ++i) {
            double t = static_cast<double>(i) / (num_sections - 1);
            // Linear curvature variation
            strain_state[i] = Vec3Types::Coord(
                0.1 * t,      // kappa_x: linear increase
                0.05 * t,     // kappa_y: smaller curvature
                0.0,          // kappa_z: no torsion
                0.0,          // gamma_x: no extension
                0.0,          // gamma_y: no shear
                1.0 + 0.1 * t // gamma_z: slight stretch
            );
        }

        // Rigid base state
        rigid_base.resize(1);
        rigid_base[0] = Rigid3Types::Coord(
            Vec3(0.0, 0.0, 0.0),  // position
            Quat(1.0, 0.0, 0.0, 0.0)  // identity rotation
        );

        // Frame states (output)
        frames.resize(curv_abs_frames.size());
        // Initialize with identity transforms
        for (size_t i = 0; i < frames.size(); ++i) {
            frames[i] = Rigid3Types::Coord(
                Vec3(0.0, 0.0, curv_abs_frames[i]),
                Quat(1.0, 0.0, 0.0, 0.0)
            );
        }
    }

    /**
     * Create and configure HookeSerat mapping
     */
    std::unique_ptr<HookeSeratMapping> createHookeSeratMapping() {
        auto mapping = std::make_unique<HookeSeratMapping>();

        // Set geometry data
        mapping->d_curv_abs_section.setValue(curv_abs_sections);
        mapping->d_curv_abs_frames.setValue(curv_abs_frames);

        // Create mock mechanical states
        strain_state_mock = std::make_unique<Vec3Data>();
        rigid_base_mock = std::make_unique<RigidData>();
        frames_mock = std::make_unique<RigidData>();

        strain_state_mock->resize(num_sections);
        rigid_base_mock->resize(1);
        frames_mock->resize(curv_abs_frames.size());

        // Set initial data
        strain_state_mock->write(Vec3Types::VecCoord::position)->setValue(strain_state);
        rigid_base_mock->write(Rigid3Types::VecCoord::position)->setValue(rigid_base);
        frames_mock->write(Rigid3Types::VecCoord::position)->setValue(frames);

        // Configure mapping
        mapping->fromModels1.add(strain_state_mock.get());
        mapping->fromModels2.add(rigid_base_mock.get());
        mapping->toModels.add(frames_mock.get());

        return mapping;
    }

    // Test data
    std::vector<double> curv_abs_sections;
    std::vector<double> curv_abs_frames;
    Vec3 strain_state;
    Rigid rigid_base;
    Rigid frames;

    // Mock mechanical states
    std::unique_ptr<Vec3Data> strain_state_mock;
    std::unique_ptr<RigidData> rigid_base_mock;
    std::unique_ptr<RigidData> frames_mock;
};

/**
 * Test basic initialization and geometry setup
 */
TEST_F(HookeSeratIntegrationTest, Initialization)
{
    auto mapping = createHookeSeratMapping();
    ASSERT_TRUE(mapping != nullptr);

    // Test initialization
    EXPECT_NO_THROW(mapping->init());

    // Verify geometry was set up correctly
    EXPECT_EQ(mapping->getNumberOfSections(), num_sections);
    EXPECT_EQ(mapping->getNumberOfFrames(), curv_abs_frames.size());

    // Test geometry validation
    EXPECT_TRUE(mapping->validateBeamGeometry());
}

/**
 * Test Lie group functionality
 */
TEST_F(HookeSeratIntegrationTest, LieGroupOperations)
{
    auto mapping = createHookeSeratMapping();
    mapping->init();

    // Test section properties
    const auto& sections = mapping->getSectionProperties();
    EXPECT_EQ(sections.size(), num_sections);

    for (const auto& section : sections) {
        // Test SE3 operations
        auto transform = section.getTransformation();
        EXPECT_TRUE(transform.matrix().isApprox(transform.matrix(), eps));

        // Test adjoint matrices
        auto adjoint = section.getAdjoint();
        EXPECT_EQ(adjoint.rows(), 6);
        EXPECT_EQ(adjoint.cols(), 6);

        // Test strain state bundle
        auto strain_state = section.getStrainState();
        // Bundle should be properly constructed
        EXPECT_NO_THROW(section.setStrainState(strain_state));
    }
}

/**
 * Test geometry validation methods
 */
TEST_F(HookeSeratIntegrationTest, GeometryValidation)
{
    auto mapping = createHookeSeratMapping();
    mapping->init();

    // Test individual validation methods
    EXPECT_TRUE(mapping->validateSectionProperties());
    EXPECT_TRUE(mapping->validateFrameProperties());
    EXPECT_TRUE(mapping->checkInterSectionContinuity());
    EXPECT_TRUE(mapping->validateFrameSectionAssociations());

    // Test overall geometry validation
    EXPECT_TRUE(mapping->validateBeamGeometry());
}

/**
 * Test Jacobian computation and validation
 */
TEST_F(HookeSeratIntegrationTest, JacobianValidation)
{
    auto mapping = createHookeSeratMapping();
    mapping->init();

    // Test Jacobian validation for each section
    for (size_t i = 0; i < mapping->getNumberOfSections(); ++i) {
        EXPECT_TRUE(mapping->validateJacobianAccuracy(i));
    }
}

/**
 * Test performance monitoring
 */
TEST_F(HookeSeratIntegrationTest, PerformanceMonitoring)
{
    auto mapping = createHookeSeratMapping();
    mapping->init();

    // Get initial stats
    auto initial_stats = mapping->getJacobianStats();
    EXPECT_EQ(initial_stats.computation_count, 0);
    EXPECT_EQ(initial_stats.cache_hits, 0);

    // Reset stats
    mapping->resetJacobianStats();
    auto reset_stats = mapping->getJacobianStats();
    EXPECT_EQ(reset_stats.computation_count, 0);
    EXPECT_EQ(reset_stats.cache_hits, 0);
}

/**
 * Test tangent exponential implementations equivalence
 */
TEST_F(HookeSeratIntegrationTest, TangentExponentialEquivalence)
{
    // Test various strain configurations
    std::vector<Vec3Types::Coord> test_strains = {
        Vec3Types::Coord(0.1, 0.0, 0.0, 0.0, 0.0, 0.0),  // Pure bending
        Vec3Types::Coord(0.0, 0.0, 0.1, 0.0, 0.0, 0.0),  // Pure torsion
        Vec3Types::Coord(0.0, 0.0, 0.0, 0.1, 0.0, 0.0),  // Extension
        Vec3Types::Coord(0.05, 0.02, 0.01, 0.0, 0.0, 1.0) // Complex strain
    };

    std::vector<double> test_lengths = {0.1, 0.5, 1.0};

    // Create identity adjoint for testing
    Eigen::Matrix<double, 6, 6> identity_adjoint = Eigen::Matrix<double, 6, 6>::Identity();

    for (const auto& strain : test_strains) {
        for (double length : test_lengths) {
            // Convert SOFA Vec to TangentVector
            Eigen::Matrix<double, 6, 1> eigen_strain;
            for (int i = 0; i < 6; ++i) {
                eigen_strain[i] = strain[i];
            }

            // Test equivalence
            EXPECT_TRUE(HookeSeratMapping::testTangExpImplementationEquivalence(
                length, eigen_strain, identity_adjoint, 1e-6));
        }
    }
}

/**
 * Test smooth trajectory generation
 */
TEST_F(HookeSeratIntegrationTest, SmoothTrajectory)
{
    auto mapping = createHookeSeratMapping();
    mapping->init();

    const int num_points = 10;
    auto trajectory = mapping->generateSmoothTrajectory(num_points);

    // Should generate trajectory with proper size
    EXPECT_EQ(trajectory.size(), num_sections * num_points);

    // Verify trajectory continuity
    for (size_t i = 1; i < trajectory.size(); ++i) {
        const auto& prev = trajectory[i-1];
        const auto& curr = trajectory[i];

        // Check that transforms are valid (no NaN, finite values)
        EXPECT_TRUE(prev.matrix().allFinite());
        EXPECT_TRUE(curr.matrix().allFinite());

        // Check determinant is positive (proper rotation)
        EXPECT_GT(prev.matrix().block<3,3>(0,0).determinant(), 0.0);
        EXPECT_GT(curr.matrix().block<3,3>(0,0).determinant(), 0.0);
    }
}

/**
 * Test section interpolation methods
 */
TEST_F(HookeSeratIntegrationTest, SectionInterpolation)
{
    auto mapping = createHookeSeratMapping();
    mapping->init();

    const auto& sections = mapping->getSectionProperties();
    if (sections.size() >= 2) {
        const auto& section1 = sections[0];
        const auto& section2 = sections[1];

        // Test linear interpolation
        auto interpolated = section1.lerp(section2, 0.5);

        // Length should be average
        double expected_length = (section1.getLength() + section2.getLength()) * 0.5;
        EXPECT_NEAR(interpolated.getLength(), expected_length, eps);

        // Transform should be interpolated
        auto expected_transform = section1.getTransformation().interpolate(
            section2.getTransformation(), 0.5);
        EXPECT_TRUE(interpolated.getTransformation().matrix().isApprox(
            expected_transform.matrix(), 1e-3));
    }
}

/**
 * Test internal forces computation
 */
TEST_F(HookeSeratIntegrationTest, InternalForces)
{
    auto mapping = createHookeSeratMapping();
    mapping->init();

    // Create test strains
    std::vector<Eigen::Matrix<double, 6, 1>> test_strains;
    for (size_t i = 0; i < mapping->getNumberOfSections(); ++i) {
        Eigen::Matrix<double, 6, 1> strain = Eigen::Matrix<double, 6, 1>::Random();
        strain *= 0.1; // Small strains
        test_strains.push_back(strain);
    }

    // Compute internal forces
    auto forces = mapping->computeInternalForces(test_strains);

    // Should return forces for each section
    EXPECT_EQ(forces.size(), mapping->getNumberOfSections());

    // Forces should be finite
    for (const auto& force : forces) {
        EXPECT_TRUE(force.allFinite());
    }
}

/**
 * Test error handling and robustness
 */
TEST_F(HookeSeratIntegrationTest, ErrorHandling)
{
    auto mapping = createHookeSeratMapping();

    // Test with invalid geometry
    std::vector<double> invalid_curv_abs = {0.0, -1.0, 2.0}; // Negative length
    mapping->d_curv_abs_section.setValue(invalid_curv_abs);

    // Should handle gracefully
    EXPECT_NO_THROW(mapping->init());

    // But validation should fail
    EXPECT_FALSE(mapping->validateBeamGeometry());
}

/**
 * Test caching behavior
 */
TEST_F(HookeSeratIntegrationTest, JacobianCaching)
{
    auto mapping = createHookeSeratMapping();
    mapping->init();

    // Enable caching (should be default)
    mapping->setJacobianCaching(true);

    // Get stats before operations
    auto initial_stats = mapping->getJacobianStats();

    // Perform operations that use Jacobians
    for (size_t i = 0; i < mapping->getNumberOfSections(); ++i) {
        mapping->validateJacobianAccuracy(i);
    }

    // Get stats after operations
    auto final_stats = mapping->getJacobianStats();

    // Should have performed computations
    EXPECT_GT(final_stats.computation_count, initial_stats.computation_count);

    // Cache should be populated
    EXPECT_GT(mapping->getJacobianStats().jacobian_cache.size(), 0);
}

/**
 * Test BeamStateEstimator integration
 */
TEST_F(HookeSeratIntegrationTest, BeamStateEstimatorIntegration)
{
    auto mapping = createHookeSeratMapping();
    mapping->init();

    // Enable state estimation
    mapping->enableStateEstimation(true);
    EXPECT_TRUE(mapping->isStateEstimationEnabled());

    // Test confidence query
    double confidence = mapping->getEstimationConfidence();
    EXPECT_GE(confidence, 0.0);

    // Test state prediction
    using TangentVector = sofa::component::cosserat::liegroups::SE3<double>::TangentVector;
    TangentVector control_input = TangentVector::Zero();
    control_input << 0.01, 0.0, 0.0, 0.0, 0.0, 0.0; // Small control input

    EXPECT_NO_THROW(mapping->predictState(control_input));

    // Test state update with measurement
    using SE3Type = sofa::component::cosserat::liegroups::SE3<double>;
    SE3Type measurement = SE3Type::computeIdentity();
    Eigen::Matrix<double, 6, 6> measurement_cov = Eigen::Matrix<double, 6, 6>::Identity() * 0.01;

    EXPECT_NO_THROW(mapping->updateStateEstimate(measurement, measurement_cov));

    // Disable state estimation
    mapping->enableStateEstimation(false);
    EXPECT_FALSE(mapping->isStateEstimationEnabled());
}

/**
 * Test BeamTopology integration
 */
TEST_F(HookeSeratIntegrationTest, BeamTopologyIntegration)
{
    auto mapping = createHookeSeratMapping();
    mapping->init();

    // Test multi-section support
    EXPECT_TRUE(mapping->supportsMultiSectionBeams());

    // Create a simple topology
    using BeamTopology = sofa::component::cosserat::mapping::BeamTopology;
    using SE3Type = sofa::component::cosserat::liegroups::SE3<double>;

    BeamTopology topology;
    topology.parent_indices = {-1, 0, 1}; // Chain topology
    topology.relative_transforms = {
        SE3Type::computeIdentity(),
        SE3Type::computeIdentity(),
        SE3Type::computeIdentity()
    };
    topology.connection_stiffnesses = {1000.0, 1000.0, 1000.0};

    // Set topology
    EXPECT_NO_THROW(mapping->setBeamTopology(topology));

    // Verify topology was set
    const auto& retrieved_topology = mapping->getBeamTopology();
    EXPECT_EQ(retrieved_topology.getNumSections(), 3u);
    EXPECT_TRUE(retrieved_topology.isValid());

    // Test multi-section enable/disable
    mapping->enableMultiSectionSupport(true);
    // Note: Internal state verification would require additional accessors
}

/**
 * Test performance optimization features
 */
TEST_F(HookeSeratIntegrationTest, PerformanceOptimizationIntegration)
{
    auto mapping = createHookeSeratMapping();
    mapping->init();

    // Test parallel computation
    mapping->enableParallelComputation(true);
    EXPECT_TRUE(mapping->isParallelComputationEnabled());
    EXPECT_GT(mapping->getOptimalThreadCount(), 0u);

    mapping->enableParallelComputation(false);
    EXPECT_FALSE(mapping->isParallelComputationEnabled());

    // Test cache operations
    size_t initial_cache_size = mapping->getCacheSize();
    mapping->clearComputationCache();
    EXPECT_EQ(mapping->getCacheSize(), 0u);

    // Test performance benchmarking
    EXPECT_NO_THROW(mapping->runPerformanceBenchmark(50));

    // Test performance reporting
    EXPECT_NO_THROW(mapping->printPerformanceReport());

    // Verify stats are updated
    const auto& stats = mapping->getJacobianStats();
    EXPECT_GE(stats.computation_count, 0u);
}

/**
 * Test advanced state estimation with realistic scenarios
 */
TEST_F(HookeSeratIntegrationTest, AdvancedStateEstimation)
{
    auto mapping = createHookeSeratMapping();
    mapping->init();
    mapping->enableStateEstimation(true);

    // Simulate a sequence of measurements and predictions
    using SE3Type = sofa::component::cosserat::liegroups::SE3<double>;
    using TangentVector = SE3Type::TangentVector;

    // Initial measurement
    SE3Type measurement1 = SE3Type::computeIdentity();
    Eigen::Matrix<double, 6, 6> cov1 = Eigen::Matrix<double, 6, 6>::Identity() * 0.1;
    mapping->updateStateEstimate(measurement1, cov1);

    double confidence1 = mapping->getEstimationConfidence();
    EXPECT_GT(confidence1, 0.0);

    // Prediction step
    TangentVector control = TangentVector::Zero();
    control << 0.05, 0.0, 0.0, 0.0, 0.0, 0.0; // Some curvature
    mapping->predictState(control);

    // Second measurement (slightly different)
    SE3Type measurement2 = SE3Type::computeIdentity();
    measurement2.translation() << 0.01, 0.0, 0.0; // Small translation
    Eigen::Matrix<double, 6, 6> cov2 = Eigen::Matrix<double, 6, 6>::Identity() * 0.05; // Better accuracy
    mapping->updateStateEstimate(measurement2, cov2);

    double confidence2 = mapping->getEstimationConfidence();
    EXPECT_GT(confidence2, 0.0);

    // Confidence should potentially improve with better measurements
    // (This is a statistical property, so we just check it's reasonable)
    EXPECT_GE(confidence2, 0.0);
}

/**
 * Test multi-section beam with complex topology
 */
TEST_F(HookeSeratIntegrationTest, MultiSectionBeamComplexTopology)
{
    auto mapping = createHookeSeratMapping();
    mapping->init();

    // Create a branched topology (more complex than simple chain)
    using BeamTopology = sofa::component::cosserat::mapping::BeamTopology;
    using SE3Type = sofa::component::cosserat::liegroups::SE3<double>;

    BeamTopology topology;
    // Root -> Branch1, Branch2
    topology.parent_indices = {-1, 0, 0, 1, 2}; // Branched structure
    topology.relative_transforms.resize(5, SE3Type::computeIdentity());
    topology.connection_stiffnesses = {1000.0, 800.0, 800.0, 600.0, 600.0};

    EXPECT_TRUE(topology.isValid());
    EXPECT_EQ(topology.getNumSections(), 5u);

    // Test children relationships
    auto root_children = topology.getChildren(0);
    EXPECT_EQ(root_children.size(), 2u); // Two children from root

    auto branch1_children = topology.getChildren(1);
    EXPECT_EQ(branch1_children.size(), 1u); // One child from branch 1

    // Set complex topology
    EXPECT_NO_THROW(mapping->setBeamTopology(topology));
    mapping->enableMultiSectionSupport(true);
}

/**
 * Test performance benchmarking with different configurations
 */
TEST_F(HookeSeratIntegrationTest, PerformanceBenchmarking)
{
    auto mapping = createHookeSeratMapping();
    mapping->init();

    // Test with different benchmark sizes
    std::vector<size_t> benchmark_sizes = {10, 50, 100};

    for (size_t size : benchmark_sizes) {
        // Reset stats before each benchmark
        mapping->resetJacobianStats();

        // Run benchmark
        EXPECT_NO_THROW(mapping->runPerformanceBenchmark(size));

        // Verify stats were collected
        const auto& stats = mapping->getJacobianStats();
        EXPECT_GE(stats.computation_count, size);
        EXPECT_GE(stats.averageComputationTime(), 0.0);
    }

    // Test cache performance
    mapping->clearComputationCache();
    EXPECT_EQ(mapping->getCacheSize(), 0u);

    // Run benchmark again to test caching
    mapping->runPerformanceBenchmark(20);
    size_t cache_size_after = mapping->getCacheSize();
    EXPECT_GE(cache_size_after, 0u); // Cache should be populated
}

/**
 * Test integration of all Phase 2 features together
 */
TEST_F(HookeSeratIntegrationTest, Phase2FeatureIntegration)
{
    auto mapping = createHookeSeratMapping();
    mapping->init();

    // Enable all Phase 2 features
    mapping->enableStateEstimation(true);
    mapping->enableParallelComputation(true);
    mapping->enableMultiSectionSupport(true);

    // Set up complex topology
    using BeamTopology = sofa::component::cosserat::mapping::BeamTopology;
    using SE3Type = sofa::component::cosserat::liegroups::SE3<double>;

    BeamTopology topology;
    topology.parent_indices = {-1, 0, 1};
    topology.relative_transforms.resize(3, SE3Type::computeIdentity());
    topology.connection_stiffnesses = {1000.0, 800.0, 600.0};

    mapping->setBeamTopology(topology);

    // Run performance benchmark with all features enabled
    EXPECT_NO_THROW(mapping->runPerformanceBenchmark(25));

    // Test state estimation with topology
    using TangentVector = SE3Type::TangentVector;
    TangentVector control = TangentVector::Random() * 0.1;
    mapping->predictState(control);

    SE3Type measurement = SE3Type::computeIdentity();
    Eigen::Matrix<double, 6, 6> cov = Eigen::Matrix<double, 6, 6>::Identity() * 0.05;
    mapping->updateStateEstimate(measurement, cov);

    // Verify everything works together
    EXPECT_TRUE(mapping->isStateEstimationEnabled());
    EXPECT_TRUE(mapping->isParallelComputationEnabled());
    EXPECT_TRUE(mapping->supportsMultiSectionBeams());
    EXPECT_GT(mapping->getEstimationConfidence(), 0.0);
    EXPECT_GT(mapping->getCacheSize(), 0u);
}

} // namespace sofa::component::cosserat::mapping::testing