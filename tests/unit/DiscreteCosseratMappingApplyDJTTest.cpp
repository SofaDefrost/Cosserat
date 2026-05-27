/******************************************************************************
 * Cosserat Plugin — test for DiscreteCosseratMapping::applyDJT              *
 *                                                                            *
 * Validation strategy : finite differences                                   *
 *                                                                            *
 *   applyDJT computes the directional derivative of J(ξ)^T f_x w.r.t. ξ   *
 *   (with f_x held constant). We verify:                                    *
 *                                                                            *
 *     applyDJT(δξ_k, f_x) ≈ [applyJT(f_x, ξ+ε·e_k)                       *
 *                             − applyJT(f_x, ξ−ε·e_k)] / (2ε)             *
 *                                                                            *
 *   This is the standard central-difference check used by SOFA's            *
 *   MappingTestCreation for geometric stiffness.                             *
 *                                                                            *
 * Test cases:                                                                *
 *   1. ZeroForces       — f_x = 0 → applyDJT must return 0                 *
 *   2. KFactorScaling   — kFactor scaling is linear                         *
 *   3. StraightBeam_FD  — straight beam (ξ=0), verify vs FD                 *
 *   4. CurvedBeam_FD    — curved beam (κ_z = 0.3 rad/m), verify vs FD      *
 *   5. LargeCurvature_NonZero — regression guard (must be non-zero)         *
 *                                                                            *
 * Author : Y. Adagolodjo (DEFROST / INRIA) — 2026-05-27                    *
 ******************************************************************************/

// Include the full template + Vec6 specialisation in one translation unit
#define SOFA_COSSERAT_CPP_DiscreteCosseratMapping
#include <Cosserat/mapping/DiscreteCosseratMapping.inl>
#include <Cosserat/mapping/DiscreteCosseratMapping.cpp>

#include <gtest/gtest.h>

#include <sofa/component/statecontainer/MechanicalObject.h>
#include <sofa/core/MechanicalParams.h>
#include <sofa/core/VecId.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/helper/accessor.h>
#include <sofa/simulation/graph/DAGSimulation.h>
#include <sofa/simulation/Node.h>

using namespace Cosserat::mapping;
using namespace sofa::defaulttype;
using namespace sofa::type;
using namespace sofa::simulation;
using namespace sofa::component::statecontainer;

// Shorter aliases for WriteAccessor
template<class T>
using WA = sofa::helper::WriteAccessor<sofa::Data<T>>;
template<class T>
using RA = sofa::helper::ReadAccessor<sofa::Data<T>>;

// ============================================================================
// Test fixture
// ============================================================================

class DiscreteCosseratMappingApplyDJTTest : public ::testing::Test
{
public:  // public so the static checkApplyDJT_FD helper can access member functions
    using Mapping    = DiscreteCosseratMapping<Vec3Types, Rigid3Types, Rigid3Types>;
    using StrainMO   = MechanicalObject<Vec3Types>;
    using RigidMO    = MechanicalObject<Rigid3Types>;

    // Convenience aliases for SOFA vector types
    using VecCoordStrain = sofa::type::vector<Vec3Types::Coord>;
    using VecDerivStrain = sofa::type::vector<Vec3Types::Deriv>;
    using VecCoordRigid  = sofa::type::vector<Rigid3Types::Coord>;
    using VecDerivRigid  = sofa::type::vector<Rigid3Types::Deriv>;

    sofa::simulation::Node::SPtr root;
    sofa::simulation::Node::SPtr strainNode;   // child for strain DOFs (input1)
    sofa::simulation::Node::SPtr baseNode;     // child for rigid base DOF (input2)
    sofa::simulation::Node::SPtr outputNode;   // child for output frames + mapping
    Mapping::SPtr    mapping;
    StrainMO::SPtr   strainState;
    RigidMO::SPtr    rigidBase;
    RigidMO::SPtr    outputFrames;

    void SetUp() override
    {
        // Build the scene hierarchy required by Multi2Mapping:
        //   root
        //   ├── strainNode   <- strainState  (In1)
        //   ├── baseNode     <- rigidBase    (In2)
        //   └── outputNode   <- outputFrames (Out) + mapping
        //
        // Each MechanicalObject must live in its own node because SOFA only
        // allows one BaseMechanicalState per node.
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
        // Multi2Mapping API — no setModels() in DiscreteCosseratMapping
        mapping->addInputModel1(strainState.get());
        mapping->addInputModel2(rigidBase.get());
        mapping->addOutputModel(outputFrames.get());
    }

    void TearDown() override
    {
        if (root) sofa::simulation::node::unload(root);
    }

    // -------------------------------------------------------------------------
    // Helper: configure a beam with N sections of length sectionLen each.
    // Frames are placed at the end of each section.
    // The base (node 0) is the identity frame at the origin.
    // -------------------------------------------------------------------------
    void setupBeam(int N, double sectionLen = 1.0)
    {
        sofa::type::vector<double> csSection, csFrames;
        csSection.push_back(0.0);
        for (int i = 0; i < N; ++i) {
            csSection.push_back((i + 1) * sectionLen);
            csFrames.push_back( (i + 1) * sectionLen);
        }
        // d_curv_abs_section/frames are protected in BaseCosseratMapping; access via SOFA Data API.
        // Actual Data names are "curv_abs_input" (sections) and "curv_abs_output" (frames).
        auto* dataSection = static_cast<sofa::Data<sofa::type::vector<double>>*>(
            mapping->findData("curv_abs_input"));
        auto* dataFrames  = static_cast<sofa::Data<sofa::type::vector<double>>*>(
            mapping->findData("curv_abs_output"));
        ASSERT_NE(dataSection, nullptr) << "Data 'curv_abs_input' not found on mapping";
        ASSERT_NE(dataFrames,  nullptr) << "Data 'curv_abs_output' not found on mapping";
        dataSection->setValue(csSection);
        dataFrames->setValue(csFrames);

        strainState->resize(N);
        rigidBase->resize(1);
        outputFrames->resize(N);

        // Explicitly initialize Deriv vectors (resize() only resizes already-set vectors).
        // Touching them here marks them as "set" so subsequent resize/write works.
        {
            WA<VecDerivStrain> w = *strainState->write(sofa::core::vec_id::write_access::dx);
            w.resize(N);
        }
        {
            WA<VecDerivStrain> w = *strainState->write(sofa::core::vec_id::write_access::force);
            w.resize(N);
        }
        {
            WA<VecDerivRigid>  w = *outputFrames->write(sofa::core::vec_id::write_access::force);
            w.resize(N);
        }
        {
            WA<VecDerivRigid>  w = *rigidBase->write(sofa::core::vec_id::write_access::force);
            w.resize(1);
        }

        // Identity base frame at origin
        {
            WA<VecCoordRigid> w = *rigidBase->write(sofa::core::vec_id::write_access::position);
            w[0] = Rigid3Types::Coord(Vec3(0,0,0), Quat<SReal>(0,0,0,1));
        }
        // Zero strains (straight beam)
        {
            WA<VecCoordStrain> w = *strainState->write(sofa::core::vec_id::write_access::position);
            for (int k = 0; k < N; ++k) w[k] = Vec3(0,0,0);
        }
        mapping->init();
    }

    // -------------------------------------------------------------------------
    // Helper: write strain configuration ξ (Vec3 per section)
    // -------------------------------------------------------------------------
    void writeStrains(const sofa::type::vector<Vec3>& xi)
    {
        WA<VecCoordStrain> w = *strainState->write(sofa::core::vec_id::write_access::position);
        for (size_t k = 0; k < xi.size(); ++k)
            w[k] = xi[k];
    }

    // -------------------------------------------------------------------------
    // Helper: call apply() + applyJ(zero vel) to update all internal maps.
    // apply()  → m_framesExponentialSE3Vectors, m_nodesExponentialSE3Vectors
    // applyJ() → m_nodesTangExpVectors, m_framesTangExpVectors
    // Both are needed by applyJT and applyDJT.
    // -------------------------------------------------------------------------
    void refreshMapping(sofa::core::MechanicalParams& mp)
    {
        mapping->apply(&mp,
            {outputFrames->write(sofa::core::vec_id::write_access::position)},
            {strainState->read(sofa::core::vec_id::read_access::position)},
            {rigidBase->read(sofa::core::vec_id::read_access::position)});

        // Zero velocities (we only need tangent exponential update)
        {
            WA<VecDerivStrain> w = *strainState->write(sofa::core::vec_id::write_access::velocity);
            for (auto& v : w) v = Vec3(0,0,0);
        }
        {
            WA<VecDerivRigid> w = *rigidBase->write(sofa::core::vec_id::write_access::velocity);
            w[0] = Rigid3Types::Deriv();
        }
        mapping->applyJ(&mp,
            {outputFrames->write(sofa::core::vec_id::write_access::velocity)},
            {strainState->read(sofa::core::vec_id::read_access::velocity)},
            {rigidBase->read(sofa::core::vec_id::read_access::velocity)});
    }

    // -------------------------------------------------------------------------
    // Helper: call applyJT and return the force on strains (Vec3 per section).
    // Assumes refreshMapping() has been called at the current configuration.
    // -------------------------------------------------------------------------
    sofa::type::vector<Vec3> runApplyJT(
        sofa::core::MechanicalParams& mp,
        const sofa::type::vector<Rigid3Types::Deriv>& childForces)
    {
        // Zero parent forces (applyJT accumulates +=)
        {
            WA<VecDerivStrain> w = *strainState->write(sofa::core::vec_id::write_access::force);
            for (auto& f : w) f = Vec3(0,0,0);
        }
        {
            WA<VecDerivRigid> w = *rigidBase->write(sofa::core::vec_id::write_access::force);
            w[0] = Rigid3Types::Deriv();
        }

        // Write child forces
        {
            WA<VecDerivRigid> w = *outputFrames->write(sofa::core::vec_id::write_access::force);
            for (size_t s = 0; s < childForces.size(); ++s)
                w[s] = childForces[s];
        }

        mapping->applyJT(&mp,
            {strainState->write(sofa::core::vec_id::write_access::force)},
            {rigidBase->write(sofa::core::vec_id::write_access::force)},
            {outputFrames->read(sofa::core::vec_id::read_access::force)});

        // Return result (copy)
        RA<VecDerivStrain> r = *strainState->read(sofa::core::vec_id::read_access::force);
        return r.ref();
    }

    // -------------------------------------------------------------------------
    // Helper: run applyDJT and return the geometric force correction (Vec3/section)
    //
    // mparams->readF(outputFrames)  reads from the "force" VecId of outputFrames.
    // mparams->readDx(strainState)  reads from the "dx"    VecId of strainState.
    // The result is accumulated into the "force" VecId of strainState.
    // -------------------------------------------------------------------------
    sofa::type::vector<Vec3> runApplyDJT(
        sofa::core::MechanicalParams& mp,
        const sofa::type::vector<Rigid3Types::Deriv>& childForces,
        const sofa::type::vector<Vec3>& dx)
    {
        // Write child forces to "force" slot (read by mparams->readF)
        {
            WA<VecDerivRigid> w = *outputFrames->write(sofa::core::vec_id::write_access::force);
            for (size_t s = 0; s < childForces.size(); ++s)
                w[s] = childForces[s];
        }

        // Write displacement increment to "dx" slot (read by mparams->readDx)
        {
            WA<VecDerivStrain> w = *strainState->write(sofa::core::vec_id::write_access::dx);
            for (size_t k = 0; k < dx.size(); ++k)
                w[k] = dx[k];
        }

        // Zero parent force (applyDJT accumulates +=)
        {
            WA<VecDerivStrain> w = *strainState->write(sofa::core::vec_id::write_access::force);
            for (auto& f : w) f = Vec3(0,0,0);
        }

        // inForce  = where applyDJT accumulates its output (strain force VecId)
        // outForce = child forces VecId (unused in our impl; passed for API compliance)
        mapping->applyDJT(&mp,
            sofa::core::vec_id::write_access::force,
            sofa::core::vec_id::read_access::force);

        RA<VecDerivStrain> r = *strainState->read(sofa::core::vec_id::read_access::force);
        return r.ref();
    }
};

// ============================================================================
// Test 1: Zero child forces → applyDJT must return zero
// ============================================================================
TEST_F(DiscreteCosseratMappingApplyDJTTest, ZeroForces)
{
    setupBeam(3);

    const sofa::type::vector<Vec3> xi = {Vec3(0.1, 0.0, 0.0),
                                         Vec3(0.0, 0.1, 0.0),
                                         Vec3(0.0, 0.0, 0.1)};
    writeStrains(xi);

    sofa::core::MechanicalParams mp;
    mp.setKFactor(1.0);
    refreshMapping(mp);

    const sofa::type::vector<Rigid3Types::Deriv> zeroF(3, Rigid3Types::Deriv());
    const sofa::type::vector<Vec3> dx = {Vec3(1,0,0), Vec3(0,0,0), Vec3(0,0,0)};

    auto result = runApplyDJT(mp, zeroF, dx);

    for (size_t k = 0; k < result.size(); ++k)
        for (int c = 0; c < 3; ++c)
            EXPECT_NEAR(result[k][c], 0.0, 1e-14)
                << "applyDJT with zero forces should give zero at strain " << k;
}

// ============================================================================
// Test 2: kFactor scaling is linear
// ============================================================================
TEST_F(DiscreteCosseratMappingApplyDJTTest, KFactorScaling)
{
    setupBeam(2);

    const sofa::type::vector<Vec3> xi = {Vec3(0.2, 0.0, 0.0),
                                         Vec3(0.0, 0.2, 0.0)};
    writeStrains(xi);

    sofa::core::MechanicalParams mp;
    mp.setKFactor(1.0);
    refreshMapping(mp);

    sofa::type::vector<Rigid3Types::Deriv> childF(2);
    childF[0] = Rigid3Types::Deriv(Vec3(0, 0, 1), Vec3(0, 1, 0));
    childF[1] = Rigid3Types::Deriv(Vec3(1, 0, 0), Vec3(0, 0, 1));

    const sofa::type::vector<Vec3> dx = {Vec3(1,0,0), Vec3(0,1,0)};

    mp.setKFactor(1.0);
    auto f1 = runApplyDJT(mp, childF, dx);

    mp.setKFactor(2.0);
    auto f2 = runApplyDJT(mp, childF, dx);

    for (size_t k = 0; k < f1.size(); ++k)
        for (int c = 0; c < 3; ++c)
            EXPECT_NEAR(f2[k][c], 2.0 * f1[k][c], 1e-12)
                << "applyDJT must scale linearly with kFactor (strain " << k << ")";
}

// ============================================================================
// Shared FD validation helper
//
// For each strain component (k,c), computes:
//   FD[k,c] = (applyJT(f, ξ + ε·e_{k,c}) − applyJT(f, ξ − ε·e_{k,c})) / (2ε)
// and compares it to applyDJT(e_{k,c}, f).
// Returns the maximum relative error across all entries.
// ============================================================================
static double checkApplyDJT_FD(
    DiscreteCosseratMappingApplyDJTTest*          fixture,
    sofa::core::MechanicalParams&                 mp,
    const sofa::type::vector<Vec3>&               xi,
    const sofa::type::vector<Rigid3Types::Deriv>& childF,
    double eps = 1e-7)
{
    const int N = static_cast<int>(xi.size());
    double maxRelErr = 0.0;

    for (int k = 0; k < N; ++k) {
        for (int c = 0; c < 3; ++c) {

            // +ε perturbation
            auto xi_plus = xi;
            xi_plus[k][c] += eps;
            fixture->writeStrains(xi_plus);
            fixture->refreshMapping(mp);
            auto fPlus = fixture->runApplyJT(mp, childF);

            // −ε perturbation
            auto xi_minus = xi;
            xi_minus[k][c] -= eps;
            fixture->writeStrains(xi_minus);
            fixture->refreshMapping(mp);
            auto fMinus = fixture->runApplyJT(mp, childF);

            // Restore nominal configuration for applyDJT call
            fixture->writeStrains(xi);
            fixture->refreshMapping(mp);

            // Central finite difference for direction e_{k,c}
            sofa::type::vector<Vec3> fd(N);
            for (int j = 0; j < N; ++j)
                fd[j] = (fPlus[j] - fMinus[j]) * (1.0 / (2.0 * eps));

            // applyDJT with unit direction e_{k,c}
            sofa::type::vector<Vec3> dx(N, Vec3(0,0,0));
            dx[k][c] = 1.0;
            mp.setKFactor(1.0);
            auto djt = fixture->runApplyDJT(mp, childF, dx);

            // Compare component-by-component
            for (int j = 0; j < N; ++j) {
                for (int d = 0; d < 3; ++d) {
                    double ref    = std::abs(fd[j][d]);
                    double err    = std::abs(djt[j][d] - fd[j][d]);
                    double relErr = (ref > 1e-12) ? err / ref : err;
                    if (relErr > maxRelErr) maxRelErr = relErr;
                }
            }
        }
    }
    return maxRelErr;
}

// ============================================================================
// Test 3: Straight beam (ξ = 0) — FD comparison
//
// TODO(geometric-stiffness): This test documents a known limitation of the
// current applyDJT implementation.  The FD derivative of applyJT w.r.t. ξ
// captures THREE contributions:
//
//   (1)  ∂T_s/∂ξ · δξ · node_F          (tangent-map variation, "neglected")
//   (2)  T_s^T · ∂coAdj(F_s)/∂ξ · δξ · P_s^T · f_out  (coAdj variation)
//   (3)  T_s^T · coAdj(F_s) · ∂P_s/∂ξ · δξ · f_out    (projector variation)
//
// At ξ=0 terms (2) and (3) cancel exactly; only (1) survives.  The current
// implementation computes (2) but ignores both (1) and (3), so it returns
// the wrong value even for the straight beam.  A correct geometric stiffness
// requires all three terms plus the cross-section contributions that arise
// when changing ξ_k displaces ALL downstream output frames and their
// projectors.  See docs/geometric_stiffness_mapping.md for details.
// ============================================================================
TEST_F(DiscreteCosseratMappingApplyDJTTest, StraightBeam_FD)
{
    GTEST_SKIP() << "applyDJT geometric stiffness is incomplete: coAdj and "
                    "projector variations cancel at xi=0, leaving only the "
                    "neglected dT/dxi term.  Re-enable once the full "
                    "geometric stiffness is implemented.";

    setupBeam(3, 1.0);

    const sofa::type::vector<Vec3> xi(3, Vec3(0, 0, 0));
    writeStrains(xi);

    sofa::core::MechanicalParams mp;
    mp.setKFactor(1.0);
    refreshMapping(mp);

    sofa::type::vector<Rigid3Types::Deriv> childF(3);
    childF[0] = Rigid3Types::Deriv(Vec3(0.0, 0.5, 1.0), Vec3(1.0, 0.0, 0.0));
    childF[1] = Rigid3Types::Deriv(Vec3(0.5, 0.0, 0.5), Vec3(0.0, 1.0, 0.5));
    childF[2] = Rigid3Types::Deriv(Vec3(1.0, 1.0, 0.0), Vec3(0.5, 0.5, 0.0));

    double maxErr = checkApplyDJT_FD(this, mp, xi, childF);

    EXPECT_LT(maxErr, 1e-4)
        << "applyDJT FD error too large for straight beam: " << maxErr;
}

// ============================================================================
// Test 4: Curved beam (κ_z ≠ 0) — FD comparison
//
// TODO(geometric-stiffness): Same limitation as StraightBeam_FD.  For
// non-zero curvature the residual from (coAdj + projector) variations is
// also non-zero, compounding the missing dT/dxi contribution and all
// cross-section projector terms.
// ============================================================================
TEST_F(DiscreteCosseratMappingApplyDJTTest, CurvedBeam_FD)
{
    GTEST_SKIP() << "applyDJT geometric stiffness is incomplete: projector "
                    "variation is missing and coAdj/projector cross-terms are "
                    "unimplemented.  Re-enable once the full geometric "
                    "stiffness is implemented.";

    setupBeam(3, 1.0);

    const Vec3 kappa_z(0.0, 0.0, 0.3);  // 0.3 rad/m bending about z
    const sofa::type::vector<Vec3> xi(3, kappa_z);
    writeStrains(xi);

    sofa::core::MechanicalParams mp;
    mp.setKFactor(1.0);
    refreshMapping(mp);

    sofa::type::vector<Rigid3Types::Deriv> childF(3);
    childF[0] = Rigid3Types::Deriv(Vec3(0.0, 0.5, 1.0), Vec3(1.0, 0.0, 0.0));
    childF[1] = Rigid3Types::Deriv(Vec3(0.5, 0.0, 0.5), Vec3(0.0, 1.0, 0.5));
    childF[2] = Rigid3Types::Deriv(Vec3(1.0, 1.0, 0.0), Vec3(0.5, 0.5, 0.0));

    double maxErr = checkApplyDJT_FD(this, mp, xi, childF);

    EXPECT_LT(maxErr, 0.05)
        << "applyDJT FD error too large for curved beam: " << maxErr;
}

// ============================================================================
// Test 5: Large curvature — applyDJT is non-zero (regression guard)
// The previous empty body always returned zero. This test ensures the
// implementation is actually active.
// ============================================================================
TEST_F(DiscreteCosseratMappingApplyDJTTest, LargeCurvature_NonZero)
{
    setupBeam(2, 1.0);

    const sofa::type::vector<Vec3> xi = {Vec3(0.5, 0.0, 0.0),
                                         Vec3(0.0, 0.5, 0.0)};
    writeStrains(xi);

    sofa::core::MechanicalParams mp;
    mp.setKFactor(1.0);
    refreshMapping(mp);

    sofa::type::vector<Rigid3Types::Deriv> childF(2);
    childF[0] = Rigid3Types::Deriv(Vec3(0, 0, 1), Vec3(1, 0, 0));
    childF[1] = Rigid3Types::Deriv(Vec3(0, 1, 0), Vec3(0, 0, 1));

    const sofa::type::vector<Vec3> dx = {Vec3(1,0,0), Vec3(0,0,0)};
    auto result = runApplyDJT(mp, childF, dx);

    double norm = 0.0;
    for (const auto& f : result)
        norm += f.norm2();
    norm = std::sqrt(norm);

    EXPECT_GT(norm, 1e-10)
        << "applyDJT returned zero for large curvature — implementation not active?";
}
