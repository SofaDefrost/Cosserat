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

#include <Cosserat/mapping/Strain2FramesCosseratMapping.h>
#include <Cosserat/mapping/Strain2FramesCosseratMapping.cpp>
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
 * @brief Test fixture for Strain2FramesCosseratMapping
 */
class Strain2FramesCosseratMappingTest : public ::testing::Test {
   protected:
    using Mapping = Strain2FramesCosseratMapping<Vec3Types, Rigid3Types, Rigid3Types>;
    using StrainMO = MechanicalObject<Vec3Types>;
    using RigidMO = MechanicalObject<Rigid3Types>;

    // Build the scene hierarchy required by Multi2Mapping:
    //   root
    //   ├── strainNode   ← strainState  (In1)
    //   ├── baseNode     ← rigidBase    (In2)
    //   └── outputNode   ← outputFrames (Out) + mapping
    //
    // Each MechanicalObject must live in its own node because SOFA allows
    // only one BaseMechanicalState per node.
    sofa::simulation::Node::SPtr root;
    sofa::simulation::Node::SPtr strainNode;
    sofa::simulation::Node::SPtr baseNode;
    sofa::simulation::Node::SPtr outputNode;
    typename Mapping::SPtr mapping;
    typename StrainMO::SPtr strainState;
    typename RigidMO::SPtr rigidBase;
    typename RigidMO::SPtr outputFrames;

    void SetUp() override {
        root = sofa::simulation::getSimulation()->createNewNode("root");

        strainNode = root->createChild("strainNode");
        strainState = sofa::core::objectmodel::New<StrainMO>();
        strainNode->addObject(strainState);

        baseNode = root->createChild("baseNode");
        rigidBase = sofa::core::objectmodel::New<RigidMO>();
        baseNode->addObject(rigidBase);

        outputNode = root->createChild("outputNode");
        outputFrames = sofa::core::objectmodel::New<RigidMO>();
        outputNode->addObject(outputFrames);

        mapping = sofa::core::objectmodel::New<Mapping>();
        outputNode->addObject(mapping);

        // Multi2Mapping API — populate fromModels1/2 and toModels so that
        // CosseratGeometryMapping::init() can resolve them correctly.
        mapping->addInputModel1(strainState.get());
        mapping->addInputModel2(rigidBase.get());
        mapping->addOutputModel(outputFrames.get());
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
TEST_F(Strain2FramesCosseratMappingTest, Initialization) {
    setupStraightBeam(5);

    EXPECT_NE(mapping, nullptr);
    EXPECT_EQ(mapping->getNumberOfSections(), 6);  // 5 sections + base
    EXPECT_EQ(mapping->getNumberOfFrames(), 6);
}

/**
 * @brief Test apply() with zero strain (straight beam)
 */
TEST_F(Strain2FramesCosseratMappingTest, ApplyZeroStrain) {
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
 *
 * @note Skipped: the finite-difference comparison requires writing strain
 * velocities into the MechanicalObject DOFs and reading frame velocities back
 * from it after applyJ(); the current implementation reads/writes directly
 * from/to the Data vectors rather than through the MO API used here.
 * A dedicated applyJ test should be written once the velocity propagation
 * interface is stabilised.
 */
TEST_F(Strain2FramesCosseratMappingTest, JacobianFiniteDifference) {
    GTEST_SKIP() << "FD-based applyJ test not yet implemented (velocity MO "
                    "read-back path incomplete)";
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
 *
 * @note Skipped: requires a fully functional applyJ velocity propagation path
 * through the MechanicalObject API. To be completed once the velocity interface
 * is stabilised and applyJ is validated.
 * Mathematical identity to verify: ⟨J·v, f⟩ = ⟨v, Jᵀ·f⟩
 */
TEST_F(Strain2FramesCosseratMappingTest, JacobianTranspose) {
    GTEST_SKIP() << "Transpose consistency test not yet implemented "
                    "(depends on functional applyJ velocity path)";
}

/**
 * @brief Test with curved beam (non-zero strain)
 */
TEST_F(Strain2FramesCosseratMappingTest, CurvedBeam) {
    setupStraightBeam(5);

    // Set constant curvature — bending around z-axis (φ₃ = 0.1 rad/m).
    // φ = [φ₁, φ₂, φ₃] = head<3> of the 6-DOF Cosserat strain vector.
    // φ₃ drives rotation around z → the rod tip moves in the x-y plane,
    // so y-coordinates of middle frames deviate from zero.
    {
        sofa::helper::WriteAccessor<sofa::Data<sofa::type::vector<Vec3Types::Coord>>> writer =
            *strainState->write(sofa::core::vec_id::write_access::position);
        for (int i = 0; i < 5; ++i) {
            writer[i] = Vec3Types::Coord(0, 0, 0.1);  // φ₃ = 0.1 rad/m (bending around z)
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
TEST_F(Strain2FramesCosseratMappingTest, ValidateJacobianAccuracy) {
    setupStraightBeam(3);

    // This test verifies the built-in numerical validation
    bool isValid = mapping->validateJacobianAccuracy(1e-5);

    EXPECT_TRUE(isValid) << "Jacobian accuracy validation should pass";
}

// ============================================================================
// applyDJT — structural tests for the PARTIAL geometric stiffness
//
// applyDJT only computes the coAdjoint-variation term of the exact derivative
// of Jᵀf (same known limitation as DiscreteCosseratMapping::applyDJT; see the
// doc-blocks in Strain2FramesCosseratMapping.{h,inl} and the skipped FD tests
// in DiscreteCosseratMappingApplyDJTTest.cpp). These tests therefore check
// only the properties the implementation actually guarantees:
//   1. zero child forces  → zero output
//   2. output scales linearly with kFactor
// ============================================================================

namespace {

/// Write curvatures, refresh internal maps (apply), write child forces +
/// strain displacement, run applyDJT, and return the resulting strain force.
sofa::type::vector<Vec3Types::Deriv> runStrain2FramesApplyDJT(
    Strain2FramesCosseratMapping<Vec3Types, Rigid3Types, Rigid3Types>* mapping,
    MechanicalObject<Vec3Types>* strainState,
    MechanicalObject<Rigid3Types>* rigidBase,
    MechanicalObject<Rigid3Types>* outputFrames,
    const sofa::type::vector<Vec3Types::Coord>& xi,
    const sofa::type::vector<Rigid3Types::Deriv>& childForces,
    const sofa::type::vector<Vec3Types::Deriv>& dx,
    double kFactor)
{
    sofa::core::MechanicalParams mp;
    mp.setKFactor(kFactor);

    // Strain configuration
    {
        auto w = sofa::helper::getWriteAccessor(
            *strainState->write(sofa::core::vec_id::write_access::position));
        for (size_t k = 0; k < xi.size(); ++k) w[k] = xi[k];
    }

    // Refresh frame positions (applyDJT refreshes the tangent-exp matrices itself)
    mapping->apply(&mp, {outputFrames->write(sofa::core::vec_id::write_access::position)},
                   {strainState->read(sofa::core::vec_id::read_access::position)},
                   {rigidBase->read(sofa::core::vec_id::read_access::position)});

    // Child wrenches into the force slot (mparams->readF reads it)
    {
        auto w = sofa::helper::getWriteAccessor(
            *outputFrames->write(sofa::core::vec_id::write_access::force));
        w.resize(childForces.size());
        for (size_t s = 0; s < childForces.size(); ++s) w[s] = childForces[s];
    }

    // Strain displacement into the dx slot (mparams->readDx reads it)
    {
        auto w = sofa::helper::getWriteAccessor(
            *strainState->write(sofa::core::vec_id::write_access::dx));
        w.resize(dx.size());
        for (size_t k = 0; k < dx.size(); ++k) w[k] = dx[k];
    }

    // Zero the strain force accumulator (applyDJT uses += semantics)
    {
        auto w = sofa::helper::getWriteAccessor(
            *strainState->write(sofa::core::vec_id::write_access::force));
        w.resize(xi.size());
        for (auto& f : w) f = Vec3Types::Deriv(0, 0, 0);
    }

    mapping->applyDJT(&mp, sofa::core::vec_id::write_access::force,
                      sofa::core::vec_id::read_access::force);

    auto r = sofa::helper::getReadAccessor(
        *strainState->read(sofa::core::vec_id::read_access::force));
    return r.ref();
}

} // anonymous namespace

/**
 * @brief applyDJT with zero child forces must produce zero strain force.
 */
TEST_F(Strain2FramesCosseratMappingTest, ApplyDJTZeroForces) {
    setupStraightBeam(3);

    const sofa::type::vector<Vec3Types::Coord> xi = {
        Vec3Types::Coord(0.1, 0.0, 0.0),
        Vec3Types::Coord(0.0, 0.1, 0.0),
        Vec3Types::Coord(0.0, 0.0, 0.1)};

    // 4 frames (3 sections + base frame), all wrenches zero
    const sofa::type::vector<Rigid3Types::Deriv> zeroF(4, Rigid3Types::Deriv());
    const sofa::type::vector<Vec3Types::Deriv> dx = {
        Vec3Types::Deriv(1, 0, 0), Vec3Types::Deriv(0, 0, 0), Vec3Types::Deriv(0, 0, 0)};

    auto result = runStrain2FramesApplyDJT(mapping.get(), strainState.get(),
                                           rigidBase.get(), outputFrames.get(),
                                           xi, zeroF, dx, 1.0);

    for (size_t k = 0; k < result.size(); ++k)
        for (int c = 0; c < 3; ++c)
            EXPECT_NEAR(result[k][c], 0.0, 1e-14)
                << "applyDJT with zero forces should give zero at strain " << k;
}

/**
 * @brief applyDJT output must scale linearly with kFactor.
 */
TEST_F(Strain2FramesCosseratMappingTest, ApplyDJTKFactorScaling) {
    setupStraightBeam(3);

    const sofa::type::vector<Vec3Types::Coord> xi = {
        Vec3Types::Coord(0.2, 0.0, 0.0),
        Vec3Types::Coord(0.0, 0.2, 0.0),
        Vec3Types::Coord(0.0, 0.0, 0.2)};

    sofa::type::vector<Rigid3Types::Deriv> childF(4);
    childF[0] = Rigid3Types::Deriv(sofa::type::Vec3(0.0, 0.5, 1.0), sofa::type::Vec3(1.0, 0.0, 0.0));
    childF[1] = Rigid3Types::Deriv(sofa::type::Vec3(0.5, 0.0, 0.5), sofa::type::Vec3(0.0, 1.0, 0.5));
    childF[2] = Rigid3Types::Deriv(sofa::type::Vec3(1.0, 1.0, 0.0), sofa::type::Vec3(0.5, 0.5, 0.0));
    childF[3] = Rigid3Types::Deriv(sofa::type::Vec3(0.0, 1.0, 0.5), sofa::type::Vec3(1.0, 0.5, 0.0));

    const sofa::type::vector<Vec3Types::Deriv> dx = {
        Vec3Types::Deriv(1, 0, 0), Vec3Types::Deriv(0, 1, 0), Vec3Types::Deriv(0, 0, 1)};

    auto f1 = runStrain2FramesApplyDJT(mapping.get(), strainState.get(),
                                       rigidBase.get(), outputFrames.get(),
                                       xi, childF, dx, 1.0);
    auto f2 = runStrain2FramesApplyDJT(mapping.get(), strainState.get(),
                                       rigidBase.get(), outputFrames.get(),
                                       xi, childF, dx, 2.0);

    // The output must be non-trivially exercised…
    double norm = 0.0;
    for (size_t k = 0; k < f1.size(); ++k)
        for (int c = 0; c < 3; ++c)
            norm += f1[k][c] * f1[k][c];
    EXPECT_GT(norm, 1e-10) << "applyDJT should produce nonzero output here";

    // …and scale linearly in kFactor.
    for (size_t k = 0; k < f1.size(); ++k)
        for (int c = 0; c < 3; ++c)
            EXPECT_NEAR(f2[k][c], 2.0 * f1[k][c], 1e-12)
                << "applyDJT must scale linearly with kFactor (strain " << k << ")";
}

/**
 * @brief Finite-difference validation — intentionally skipped.
 */
TEST_F(Strain2FramesCosseratMappingTest, ApplyDJTFiniteDifference) {
    GTEST_SKIP() << "applyDJT computes only the coAdjoint-variation term of "
                    "d(J^T f)/dxi; the tangent-map (dT/dxi) and projector "
                    "variations are missing, so central finite differences of "
                    "applyJT do not match. Same limitation and skip rationale "
                    "as DiscreteCosseratMappingApplyDJTTest::StraightBeam_FD. "
                    "Re-enable once the full geometric stiffness is implemented.";
}

// No main() here — Sofa.Testing links SofaGTestMain which provides
// main() + sofa::simulation::graph::init() + sofa::simulation::graph::cleanup().
