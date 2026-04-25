/******************************************************************************
 *       SOFA, Simulation Open-Framework Architecture, development version     *
 *                (c) 2006-2019 INRIA, USTL, UJF, CNRS, MGH                    *
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
 * along with this program. If not, see <http://www.gnu.org/licenses/>.        *
 ******************************************************************************/

#include <Cosserat/mapping/Strain2RigidCosseratMapping.h>
#include <Cosserat/mapping/Strain2RigidCosseratMapping.cpp>
#include <gtest/gtest.h>
#include <sofa/component/statecontainer/MechanicalObject.h>
#include <sofa/core/MechanicalParams.h>
#include <sofa/core/VecId.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/helper/accessor.h>
#include <sofa/helper/system/FileRepository.h>
#include <sofa/simulation/Node.h>
#include <sofa/simulation/Simulation.h>
#include <sofa/simulation/graph/DAGSimulation.h>

using namespace Cosserat::mapping;
using namespace sofa::defaulttype;
using namespace sofa::type;
using namespace sofa::simulation;
using namespace sofa::component::statecontainer;

/**
 * @brief Test fixture for Strain2RigidCosseratMapping
 */
class Strain2RigidCosseratMappingTest : public ::testing::Test {
   protected:
    using Mapping = Strain2RigidCosseratMapping<Vec3Types, Rigid3Types, Rigid3Types>;
    using StrainMO = MechanicalObject<Vec3Types>;
    using RigidMO = MechanicalObject<Rigid3Types>;

    sofa::simulation::Node::SPtr root;
    typename Mapping::SPtr mapping;
    typename StrainMO::SPtr strainState;
    typename RigidMO::SPtr rigidBase;
    typename RigidMO::SPtr outputFrames;

    void SetUp() override {
        // Create simulation root
        root = sofa::simulation::getSimulation()->createNewNode("root");

        // Create mechanical objects
        strainState = sofa::core::objectmodel::New<StrainMO>();
        rigidBase = sofa::core::objectmodel::New<RigidMO>();
        outputFrames = sofa::core::objectmodel::New<RigidMO>();

        // Add to scene graph
        root->addObject(strainState);
        root->addObject(rigidBase);
        root->addObject(outputFrames);

        // Create mapping
        mapping = sofa::core::objectmodel::New<Mapping>();
        root->addObject(mapping);

        // Link inputs and outputs
        mapping->setModels(strainState.get(), rigidBase.get(), outputFrames.get());
    }

    void TearDown() override {
        if (root) {
            sofa::simulation::node::unload(root);
        }
    }

    /**
     * @brief Setup a simple straight beam configuration
     */
    void setupStraightBeam(int numSections = 5) {
        // Setup curvilinear abscissas
        sofa::type::vector<double> curvAbsSection;
        sofa::type::vector<double> curvAbsFrames;

        double sectionLength = 1.0;
        for (int i = 0; i <= numSections; ++i) {
            curvAbsSection.push_back(i * sectionLength);
            curvAbsFrames.push_back(i * sectionLength);
        }

        mapping->d_curv_abs_section.setValue(curvAbsSection);
        mapping->d_curv_abs_frames.setValue(curvAbsFrames);

        // Initialize strain state (zero strain = straight beam)
        strainState->resize(numSections);
        {
            sofa::helper::WriteAccessor<sofa::Data<sofa::type::vector<Vec3Types::Coord>>> writer =
                *strainState->write(sofa::core::vec_id::write_access::position);
            for (int i = 0; i < numSections; ++i) {
                writer[i] = Vec3Types::Coord(0, 0, 0);
            }
        }

        // Initialize rigid base (identity)
        rigidBase->resize(1);
        {
            sofa::helper::WriteAccessor<sofa::Data<sofa::type::vector<Rigid3Types::Coord>>> writer =
                *rigidBase->write(sofa::core::vec_id::write_access::position);
            writer[0] = Rigid3Types::Coord(sofa::type::Vec3(0, 0, 0), Quat<SReal>(0, 0, 0, 1));
        }

        // Initialize output frames
        outputFrames->resize(numSections + 1);

        // Initialize mapping
        mapping->init();
    }
};

/**
 * @brief Test basic initialization
 */
TEST_F(Strain2RigidCosseratMappingTest, Initialization) {
    setupStraightBeam(5);

    EXPECT_NE(mapping, nullptr);
    EXPECT_EQ(mapping->getNumberOfSections(), 6);  // 5 sections + base
    EXPECT_EQ(mapping->getNumberOfFrames(), 6);
}

/**
 * @brief Test apply() with zero strain (straight beam)
 */
TEST_F(Strain2RigidCosseratMappingTest, ApplyZeroStrain) {
    setupStraightBeam(5);

    // Apply mapping
    sofa::core::MechanicalParams mparams;
    mapping->apply(&mparams, {outputFrames->write(sofa::core::vec_id::write_access::position)},
                   {strainState->read(sofa::core::vec_id::read_access::position)},
                   {rigidBase->read(sofa::core::vec_id::read_access::position)});

    const auto &frames = outputFrames->read(sofa::core::vec_id::read_access::position)->getValue();

    for (size_t i = 0; i < frames.size(); ++i) {
        const auto &frame = frames[i];

        // Check position is along x-axis
        EXPECT_NEAR(frame.getCenter()[0], i * 1.0, 1e-6) << "Frame " << i << " x position";
        EXPECT_NEAR(frame.getCenter()[1], 0.0, 1e-6) << "Frame " << i << " y position";
        EXPECT_NEAR(frame.getCenter()[2], 0.0, 1e-6) << "Frame " << i << " z position";

        // Check orientation is identity
        const auto &quat = frame.getOrientation();
        EXPECT_NEAR(quat[0], 0.0, 1e-6) << "Frame " << i << " quat x";
        EXPECT_NEAR(quat[1], 0.0, 1e-6) << "Frame " << i << " quat y";
        EXPECT_NEAR(quat[2], 0.0, 1e-6) << "Frame " << i << " quat z";
        EXPECT_NEAR(quat[3], 1.0, 1e-6) << "Frame " << i << " quat w";
    }
}

/**
 * @brief Test applyJ() with finite differences
 */
TEST_F(Strain2RigidCosseratMappingTest, JacobianFiniteDifference) {
    setupStraightBeam(3);

    const double epsilon = 1e-7;
    const double tolerance = 1e-5;

    // Get initial positions
    sofa::core::MechanicalParams mparams;
    mapping->apply(&mparams, {outputFrames->write(sofa::core::vec_id::write_access::position)},
                   {strainState->read(sofa::core::vec_id::read_access::position)},
                   {rigidBase->read(sofa::core::vec_id::read_access::position)});

    const auto &frames0 = outputFrames->read(sofa::core::vec_id::read_access::position)->getValue();

    // Test Jacobian for each strain component
    for (int strainIdx = 0; strainIdx < 3; ++strainIdx) {
        for (int component = 0; component < 3; ++component) {
            // Perturb strain
            {
                sofa::helper::WriteAccessor<sofa::Data<sofa::type::vector<Vec3Types::Coord>>>
                    writer = *strainState->write(sofa::core::vec_id::write_access::position);
                writer[strainIdx][component] += epsilon;
            }
            mapping->apply(&mparams,
                           {outputFrames->write(sofa::core::vec_id::write_access::position)},
                           {strainState->read(sofa::core::vec_id::read_access::position)},
                           {rigidBase->read(sofa::core::vec_id::read_access::position)});

            const auto &framesPerturbed =
                outputFrames->read(sofa::core::vec_id::read_access::position)->getValue();

            // Compute finite difference
            sofa::type::vector<Rigid3Types::Deriv> fdJacobian;
            fdJacobian.resize(framesPerturbed.size());

            for (size_t i = 0; i < framesPerturbed.size(); ++i) {
                // Approximate derivative
                for (int k = 0; k < 6; ++k) {
                    if (k < 3) {
                        fdJacobian[i][k] =
                            (framesPerturbed[i].getCenter()[k] - frames0[i].getCenter()[k]) /
                            epsilon;
                    } else {
                        // For orientation, use quaternion difference (simplified)
                        fdJacobian[i][k] = (framesPerturbed[i].getOrientation()[k - 3] -
                                            frames0[i].getOrientation()[k - 3]) /
                                           epsilon;
                    }
                }
            }

            // Reset strain
            {
                sofa::helper::WriteAccessor<sofa::Data<sofa::type::vector<Vec3Types::Coord>>>
                    writer = *strainState->write(sofa::core::vec_id::write_access::position);
                writer[strainIdx][component] -= epsilon;
            }
            // Compute analytical Jacobian using applyJ
            sofa::type::vector<Vec3Types::Deriv> strainVel;
            strainVel.resize(3);
            strainVel[strainIdx][component] = 1.0;

            sofa::type::vector<Rigid3Types::Deriv> baseVel;
            baseVel.resize(1);
            baseVel[0] = Rigid3Types::Deriv(sofa::type::Vec3(0, 0, 0), sofa::type::Vec3(0, 0, 0));

            sofa::type::vector<Rigid3Types::Deriv> frameVel;
            frameVel.resize(framesPerturbed.size());

            mapping->applyJ(&mparams,
                            {outputFrames->write(sofa::core::vec_id::write_access::velocity)},
                            {strainState->read(sofa::core::vec_id::read_access::velocity)},
                            {rigidBase->read(sofa::core::vec_id::read_access::velocity)});

            // Compare (simplified - full comparison would need proper SE(3) metrics)
            // This is a basic sanity check
            for (size_t i = 0; i < frameVel.size(); ++i) {
                for (int k = 0; k < 3; ++k) {
                    double diff = std::abs(frameVel[i][k] - fdJacobian[i][k]);
                    EXPECT_LT(diff, tolerance)
                        << "Jacobian mismatch at frame " << i << " component " << k
                        << " for strain " << strainIdx << "," << component;
                }
            }
        }
    }
}

/**
 * @brief Test applyJT() is transpose of applyJ()
 */
TEST_F(Strain2RigidCosseratMappingTest, JacobianTranspose) {
    setupStraightBeam(3);

    sofa::core::MechanicalParams mparams;

    // Create random velocities
    sofa::type::vector<Vec3Types::Deriv> strainVel;
    strainVel.resize(3);
    for (int i = 0; i < 3; ++i) {
        strainVel[i] = Vec3Types::Deriv(0.1 * i, 0.2 * i, 0.3 * i);
    }

    sofa::type::vector<Rigid3Types::Deriv> baseVel;
    baseVel.resize(1);
    baseVel[0] =
        Rigid3Types::Deriv(sofa::type::Vec3(0.1, 0.2, 0.3), sofa::type::Vec3(0.01, 0.02, 0.03));

    // Apply J
    sofa::type::vector<Rigid3Types::Deriv> frameVel;
    frameVel.resize(4);

    // TODO: Complete this test when applyJ is fully functional
    // Test: <J*v, f> = <v, J^T*f>
}

/**
 * @brief Test with curved beam (non-zero strain)
 */
TEST_F(Strain2RigidCosseratMappingTest, CurvedBeam) {
    setupStraightBeam(5);

    // Set constant curvature (bending in y-direction)
    {
        sofa::helper::WriteAccessor<sofa::Data<sofa::type::vector<Vec3Types::Coord>>> writer =
            *strainState->write(sofa::core::vec_id::write_access::position);
        for (int i = 0; i < 5; ++i) {
            writer[i] = Vec3Types::Coord(0, 0.1, 0);  // Curvature around z-axis
        }
    }

    // Apply mapping
    sofa::core::MechanicalParams mparams;
    mapping->apply(&mparams, {outputFrames->write(sofa::core::vec_id::write_access::position)},
                   {strainState->read(sofa::core::vec_id::read_access::position)},
                   {rigidBase->read(sofa::core::vec_id::read_access::position)});

    const auto &frames = outputFrames->read(sofa::core::vec_id::read_access::position)->getValue();

    // Verify beam is curved (not straight)
    bool isCurved = false;
    for (size_t i = 1; i < frames.size() - 1; ++i) {
        // Check if middle frames deviate from straight line
        double expectedY = 0.0;  // Straight beam would have y=0
        if (std::abs(frames[i].getCenter()[1] - expectedY) > 0.01) {
            isCurved = true;
            break;
        }
    }

    EXPECT_TRUE(isCurved) << "Beam should be curved with non-zero strain";
}

/**
 * @brief Test validateJacobianAccuracy method
 */
TEST_F(Strain2RigidCosseratMappingTest, ValidateJacobianAccuracy) {
    setupStraightBeam(3);

    // This test verifies the built-in numerical validation
    bool isValid = mapping->validateJacobianAccuracy(1e-5);

    EXPECT_TRUE(isValid) << "Jacobian accuracy validation should pass";
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
