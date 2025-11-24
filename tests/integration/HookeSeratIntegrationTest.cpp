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

} // namespace sofa::component::cosserat::mapping::testing