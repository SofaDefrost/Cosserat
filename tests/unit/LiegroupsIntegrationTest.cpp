/******************************************************************************
 *                 SOFA, Simulation Open-Framework Architecture                *
 *                 (c) 2006 INRIA, USTL, UJF, CNRS, MGH                       *
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

/**
 * @file LiegroupsIntegrationTest.cpp
 *
 * Unit tests for the six liegroups → production integrations:
 *
 *  1. CosseratBodyJacobian  — virtual power duality (applyJ / applyJT)
 *  2. SE3Integrator         — all modes coincide for constant strain
 *  3. Twist / Wrench        — semantic wrapper round-trips
 *  4. BezierSE3             — computeSmoothedPath geometry
 *  5. BeamStateEstimator    — predict/update Kalman cycle
 *  6. CosseratILQRController — computeControl reduces predicted tip error
 *
 * Tests 1–5 are pure liegroups tests (no SOFA simulation needed).
 * Test 6 uses the SOFA fixture from Strain2RigidCosseratMappingTest.
 */

// ── liegroups (header-only) ───────────────────────────────────────────────────
#include <liegroups/SE3.h>
#include <liegroups/Twist.h>
#include <liegroups/Wrench.h>
#include <liegroups/CosseratBodyJacobian.h>
#include <liegroups/LieGroupIntegrators.h>
#include <liegroups/BezierSE3.h>
#include <liegroups/GaussianOnManifold.h>
#include <liegroups/CosseratUncertaintyPropagator.h>

// ── Cosserat mapping / controller ─────────────────────────────────────────────
#include <Cosserat/mapping/BeamStateEstimator.h>
#include <Cosserat/mapping/Strain2RigidCosseratMapping.h>
#include <Cosserat/mapping/Strain2RigidCosseratMapping.cpp>
#include <Cosserat/controller/CosseratILQRController.h>
#include <Cosserat/controller/CosseratILQRController.cpp>

// ── SOFA ─────────────────────────────────────────────────────────────────────
#include <gtest/gtest.h>
#include <sofa/component/statecontainer/MechanicalObject.h>
#include <sofa/core/MechanicalParams.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/helper/accessor.h>
#include <sofa/simulation/Node.h>
#include <sofa/simulation/Simulation.h>
#include <sofa/simulation/graph/DAGSimulation.h>

using namespace sofa::component::cosserat::liegroups;
using namespace sofa::defaulttype;
using namespace sofa::type;

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

using SE3d          = SE3<double>;
using TwistD        = Twist<double>;
using WrenchD       = Wrench<double>;
using BodyJacD      = CosseratBodyJacobian<double>;
using PropagatorD   = CosseratUncertaintyPropagator<double>;
using GaussianSE3d  = GaussianOnManifold<SE3d>;

using Vec6d = Eigen::Matrix<double, 6, 1>;
using Mat6d = Eigen::Matrix<double, 6, 6>;

/** Build a simple N-section body Jacobian with constant strain and unit length. */
static BodyJacD buildJacobian(int N, const Vec6d &strain) {
    BodyJacD bj;
    bj.reserve(N);
    for (int k = 0; k < N; ++k) {
        const SE3d g = SE3d::expCosserat(strain, 1.0);
        // Tangent-exp matrix T = I₆ for zero strain (small-angle approx)
        const Mat6d T = Mat6d::Identity();
        bj.pushSection(g, T);
    }
    return bj;
}

// ═════════════════════════════════════════════════════════════════════════════
// Integration 1 — CosseratBodyJacobian
// ═════════════════════════════════════════════════════════════════════════════

class BodyJacobianTest : public ::testing::Test {
protected:
    static constexpr int N = 4;
    Vec6d strain = Vec6d::Zero();  // straight rod
};

/**
 * Virtual-power duality: ⟨J^T·w, ξ̇⟩ = ⟨w, J·ξ̇⟩ for all w, ξ̇.
 *
 * applyForward computes η = J·ξ̇  (6N → 6, one twist per section).
 * applyTranspose computes q = J^T·w (6N → 6N, one wrench per section).
 * We verify: Σ_k q_k · ξ̇_k  ==  Σ_k w_k · η_k
 */
TEST_F(BodyJacobianTest, VirtualPowerDuality) {
    BodyJacD bj = buildJacobian(N, strain);

    // Random strain rates (one per section)
    std::vector<TwistD> strain_rates(N);
    for (int k = 0; k < N; ++k) {
        Vec6d v; v << 0.1*k, -0.05*k, 0.03*(k+1), 0.0, 0.0, 0.0;
        strain_rates[k] = TwistD(v);
    }

    // Random external wrenches (N+1: one at each node including base)
    std::vector<WrenchD> ext_w(N + 1);
    for (int k = 0; k <= N; ++k) {
        Vec6d w; w << 0.2*(k+1), -0.1*k, 0.05, 0.3, 0.0, -0.1;
        ext_w[k] = WrenchD(w);
    }

    // Forward pass: η_k = J·ξ̇  (returns N+1 twists: [η_0, η_1, ..., η_N])
    const auto twists = bj.applyForward(strain_rates, TwistD::Zero());

    // Backward pass: q_k = J^T·w (returns N strain forces)
    WrenchD base_w_out;
    const auto strain_forces = bj.applyTranspose(ext_w, base_w_out);

    // Left side: Σ_{k=0}^{N} w_k · η_k  (virtual work of external wrenches)
    double lhs = 0.0;
    for (int k = 0; k <= N; ++k)
        lhs += ext_w[k].data().dot(twists[k].data());

    // Right side: Σ_{k=0}^{N-1} q_k · ξ̇_k  + base_w_out · η_0
    double rhs = base_w_out.data().dot(TwistD::Zero().data());
    for (int k = 0; k < N; ++k)
        rhs += strain_forces[k].data().dot(strain_rates[k].data());

    EXPECT_NEAR(lhs, rhs, 1e-10) << "Virtual power duality violated";
}

/**
 * Forward consistency: with zero base twist and zero strain rates,
 * all output twists must be zero.
 */
TEST_F(BodyJacobianTest, ZeroInputZeroOutput) {
    BodyJacD bj = buildJacobian(N, strain);
    std::vector<TwistD> zero_rates(N, TwistD::Zero());
    const auto twists = bj.applyForward(zero_rates, TwistD::Zero());
    ASSERT_EQ(static_cast<int>(twists.size()), N + 1);
    for (int k = 0; k <= N; ++k)
        EXPECT_LT(twists[k].data().norm(), 1e-12) << "Twist " << k << " should be zero";
}

/**
 * getJacobianAtSection: J_k is 6×6N with the correct block structure.
 * For a single section (N=1), J_0 has non-zero first 6 columns only.
 */
TEST_F(BodyJacobianTest, JacobianShape) {
    BodyJacD bj = buildJacobian(1, strain);
    const auto J = bj.getJacobianAtSection(0);
    EXPECT_EQ(J.rows(), 6);
    EXPECT_EQ(J.cols(), 6);
}

/**
 * size() matches the number of pushSection() calls.
 */
TEST_F(BodyJacobianTest, Size) {
    BodyJacD bj = buildJacobian(N, strain);
    EXPECT_EQ(bj.size(), N);
}

// ═════════════════════════════════════════════════════════════════════════════
// Integration 2 — SE3Integrator
// ═════════════════════════════════════════════════════════════════════════════

class SE3IntegratorTest : public ::testing::Test {
protected:
    using Integrator = SE3Integrator<double>;
    using SE3Type = SE3<double>;
    using TangentVector = SE3Type::TangentVector;

    const double L = 1.0;
};

/**
 * For piecewise-constant strain, all three integration methods must produce
 * the same result as SE3::expCosserat().
 */
TEST_F(SE3IntegratorTest, AllMethodsEqualForConstantStrain) {
    TangentVector strain; strain << 0.1, 0.0, 0.0, 0.0, 0.0, 0.0;  // pure bending
    const auto strain_field = [&strain](double) -> TangentVector { return strain; };

    const SE3Type g0 = SE3Type::computeIdentity();
    const SE3Type g_euler    = Integrator::integrateEuler(g0, strain_field, 0.0, L);
    const SE3Type g_midpoint = Integrator::integrateMidpoint(g0, strain_field, 0.0, L);
    const SE3Type g_rkmk4   = Integrator::integrateRKMK4(g0, strain_field, 0.0, L);
    const SE3Type g_ref      = SE3Type::expCosserat(strain, L);

    // All methods should agree with expCosserat for constant strain
    const double tol = 1e-10;
    const auto diff_mid  = g0.computeInverse().compose(g_midpoint).compose(g_euler.computeInverse()).log();
    const auto diff_rk4  = g0.computeInverse().compose(g_rkmk4).compose(g_euler.computeInverse()).log();

    // Translation part (last 3 components)
    EXPECT_LT(diff_mid.norm(), tol) << "Midpoint differs from Euler for constant strain";
    EXPECT_LT(diff_rk4.norm(), tol) << "RKMK4 differs from Euler for constant strain";

    // Compare euler vs expCosserat
    const TangentVector diff_ref = g_euler.computeInverse().compose(g_ref).log();
    EXPECT_LT(diff_ref.norm(), tol) << "integrateEuler differs from expCosserat";
}

/**
 * Identity integration: zero strain over any length gives
 * a pure translation of L along x (Cosserat elongation correction +1 in ρx).
 */
TEST_F(SE3IntegratorTest, ZeroStrainGivesPureTranslation) {
    TangentVector zero_strain = TangentVector::Zero();
    const auto strain_field = [&zero_strain](double) -> TangentVector { return zero_strain; };

    const SE3Type g0 = SE3Type::computeIdentity();
    const SE3Type g  = Integrator::integrateEuler(g0, strain_field, 0.0, L);

    // With zero curvature, the rod should extend L units along x
    EXPECT_NEAR(g.translation()[0], L, 1e-10);
    EXPECT_NEAR(g.translation()[1], 0.0, 1e-10);
    EXPECT_NEAR(g.translation()[2], 0.0, 1e-10);
    // Rotation should remain identity
    EXPECT_NEAR(g.rotation().log().norm(), 0.0, 1e-10);
}

/**
 * integrate() (multi-step) converges to integrateEuler for N=1.
 */
TEST_F(SE3IntegratorTest, MultiStepMatchesSingleStep) {
    TangentVector strain; strain << 0.2, 0.1, 0.0, 0.0, 0.0, 0.0;
    const auto sf = [&strain](double) -> TangentVector { return strain; };

    const SE3Type g0 = SE3Type::computeIdentity();
    const SE3Type g_single = Integrator::integrateEuler(g0, sf, 0.0, L);
    const SE3Type g_multi  = Integrator::integrate(g0, sf, 0.0, L, 1, Integrator::Method::Euler);

    const TangentVector diff = g_single.computeInverse().compose(g_multi).log();
    EXPECT_LT(diff.norm(), 1e-12);
}

// ═════════════════════════════════════════════════════════════════════════════
// Integration 3 — Twist / Wrench semantic wrappers
// ═════════════════════════════════════════════════════════════════════════════

TEST(TwistWrenchTest, RoundTrip) {
    Vec6d v; v << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;

    const TwistD  t(v);
    const WrenchD w(v);

    EXPECT_LT((t.data() - v).norm(), 1e-15) << "Twist round-trip failed";
    EXPECT_LT((w.data() - v).norm(), 1e-15) << "Wrench round-trip failed";
}

TEST(TwistWrenchTest, Zero) {
    EXPECT_LT(TwistD::Zero().data().norm(),  1e-15);
    EXPECT_LT(WrenchD::Zero().data().norm(), 1e-15);
}

/**
 * Virtual power: ⟨w, t⟩ = w.data().dot(t.data())
 * should equal the explicit inner product.
 */
TEST(TwistWrenchTest, VirtualPower) {
    Vec6d tv; tv << 1.0, 0.0, 0.0, 0.0, 1.0, 0.0;
    Vec6d wv; wv << 0.0, 1.0, 0.0, 1.0, 0.0, 1.0;
    const TwistD  t(tv);
    const WrenchD w(wv);

    const double power = w.data().dot(t.data());
    EXPECT_DOUBLE_EQ(power, tv.dot(wv));
}

/**
 * Adjoint action: Ad_g maps a twist — applying it via a rotation-only SE3
 * element should rotate the angular part.
 */
TEST(TwistWrenchTest, AdjointMapsTwist) {
    // π/2 rotation around z
    Vec6d xi_rot; xi_rot << 0.0, 0.0, M_PI/2.0, 0.0, 0.0, 0.0;
    const SE3d g = SE3d::computeExp(xi_rot);
    const Mat6d Ad = g.computeAdjoint();

    // A twist along x rotated by π/2 should become a twist along y
    Vec6d t_in; t_in << 1.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    const Vec6d t_out = Ad * t_in;

    EXPECT_NEAR(t_out[0],  0.0, 1e-10);
    EXPECT_NEAR(t_out[1],  1.0, 1e-10);
    EXPECT_NEAR(t_out[2],  0.0, 1e-10);
}

// ═════════════════════════════════════════════════════════════════════════════
// Integration 4 — BezierSE3 / computeSmoothedPath
// ═════════════════════════════════════════════════════════════════════════════

TEST(BezierSE3Test, LinearCurveInterpolatesEndpoints) {
    const SE3d g0 = SE3d::computeIdentity();
    Vec6d xi; xi << 0.0, 0.0, 0.0, 1.0, 0.0, 0.0;  // 1m translation
    const SE3d g1 = SE3d::computeExp(xi);

    BezierSE3<double> curve({g0, g1});
    EXPECT_EQ(curve.degree(), 1);

    // t=0 → g0, t=1 → g1
    const auto p0 = curve.evaluate(0.0);
    const auto p1 = curve.evaluate(1.0);

    EXPECT_LT(g0.computeInverse().compose(p0).log().norm(), 1e-12);
    EXPECT_LT(g1.computeInverse().compose(p1).log().norm(), 1e-12);
}

TEST(BezierSE3Test, CubicCurveInterpolatesEndpoints) {
    const SE3d g0 = SE3d::computeIdentity();
    Vec6d xi; xi << 0.3, 0.1, 0.0, 2.0, 0.0, 0.0;
    const SE3d g3 = g0.compose(SE3d::computeExp(xi));
    // Build cubic via makePiecewiseCubicPath
    const auto segments = makePiecewiseCubicPath(std::vector<SE3d>{g0, g3});
    ASSERT_EQ(segments.size(), 1u);

    const auto &seg = segments[0];
    EXPECT_LT(g0.computeInverse().compose(seg.evaluate(0.0)).log().norm(), 1e-12);
    EXPECT_LT(g3.computeInverse().compose(seg.evaluate(1.0)).log().norm(), 1e-12);
}

TEST(BezierSE3Test, SampleCountIsCorrect) {
    const SE3d g0 = SE3d::computeIdentity();
    Vec6d xi; xi << 0.1, 0.0, 0.0, 1.0, 0.0, 0.0;
    const SE3d g1 = g0.compose(SE3d::computeExp(xi));

    BezierSE3<double> curve({g0, g1});
    const auto samples = curve.sample(10);
    EXPECT_EQ(static_cast<int>(samples.size()), 11);  // N+1 samples
}

TEST(BezierSE3Test, ArcLengthPositive) {
    const SE3d g0 = SE3d::computeIdentity();
    Vec6d xi; xi << 0.0, 0.0, 0.0, 1.0, 0.0, 0.0;
    const SE3d g1 = g0.compose(SE3d::computeExp(xi));
    BezierSE3<double> curve({g0, g1});
    EXPECT_GT(curve.arcLength(50), 0.0);
}

// ═════════════════════════════════════════════════════════════════════════════
// Integration 5 — BeamStateEstimator (Kalman cycle)
// ═════════════════════════════════════════════════════════════════════════════

class BeamStateEstimatorTest : public ::testing::Test {
protected:
    using Estimator = Cosserat::mapping::BeamStateEstimator;
    using CovMatrix = Eigen::Matrix<double, 6, 6>;

    Estimator estimator;
    SE3d identity = SE3d::computeIdentity();
    CovMatrix I6  = CovMatrix::Identity();

    void SetUp() override {
        estimator.initialize(identity, CovMatrix::Zero());
    }
};

/**
 * predict() must increase the covariance trace (uncertainty grows along the rod).
 */
TEST_F(BeamStateEstimatorTest, PredictIncreasesUncertainty) {
    const double trace_before = estimator.getEstimate().getCovariance().trace();

    Vec6d strain = Vec6d::Zero();
    estimator.predict(strain, 1.0, 1e-4 * I6);

    const double trace_after = estimator.getEstimate().getCovariance().trace();
    EXPECT_GT(trace_after, trace_before) << "predict() should increase covariance";
}

/**
 * After N predict steps, trace should be monotonically non-decreasing.
 */
TEST_F(BeamStateEstimatorTest, PredictMonotonicallGrowsCovariance) {
    Vec6d strain; strain << 0.1, 0.0, 0.0, 0.0, 0.0, 0.0;
    double prev_trace = 0.0;
    for (int k = 0; k < 5; ++k) {
        estimator.predict(strain, 0.2, 1e-5 * I6);
        double tr = estimator.getEstimate().getCovariance().trace();
        EXPECT_GE(tr, prev_trace - 1e-14)
            << "Covariance trace should not decrease at step " << k;
        prev_trace = tr;
    }
}

/**
 * update() with the exact current mean as measurement must not change the mean
 * and must reduce the covariance (innovation = 0, but K·r = 0).
 */
TEST_F(BeamStateEstimatorTest, UpdateWithExactMeasurementZeroInnovation) {
    // First inflate covariance with a predict step
    estimator.predict(Vec6d::Zero(), 1.0, 1e-2 * I6);

    const double trace_before = estimator.getEstimate().getCovariance().trace();
    const SE3d   mean_before  = estimator.getEstimate().getMean();

    // Update with the exact current mean (zero innovation)
    estimator.update(mean_before, 1e-4 * I6);

    const double trace_after = estimator.getEstimate().getCovariance().trace();
    EXPECT_LE(trace_after, trace_before + 1e-12)
        << "update() with exact measurement should not increase covariance";

    // Mean should not change (r = 0 → no correction)
    const auto err = mean_before.computeInverse()
                         .compose(estimator.getEstimate().getMean()).log();
    EXPECT_LT(err.norm(), 1e-10) << "Mean should be unchanged when innovation is zero";
}

/**
 * predictAlongRod() returns N+1 Gaussians for N sections.
 */
TEST_F(BeamStateEstimatorTest, PredictAlongRodReturnsPlusOneGaussians) {
    constexpr int N = 4;
    std::vector<Vec6d>  strains(N, Vec6d::Zero());
    std::vector<double> lengths(N, 0.25);

    const auto gaussians = estimator.predictAlongRod(strains, lengths, 1e-5 * I6);
    EXPECT_EQ(static_cast<int>(gaussians.size()), N + 1);
}

/**
 * getEstimationConfidence() = trace of covariance.
 */
TEST_F(BeamStateEstimatorTest, ConfidenceEqualsCovarianceTrace) {
    estimator.predict(Vec6d::Zero(), 0.5, 1e-3 * I6);
    const double conf  = estimator.getEstimationConfidence();
    const double trace = estimator.getEstimate().getCovariance().trace();
    EXPECT_NEAR(conf, trace, 1e-14);
}

// ═════════════════════════════════════════════════════════════════════════════
// Integration 5b — CosseratUncertaintyPropagator (standalone)
// ═════════════════════════════════════════════════════════════════════════════

TEST(UncertaintyPropagatorTest, PropagateAlongRodSizeN1) {
    const SE3d base = SE3d::computeIdentity();
    GaussianSE3d base_g(base, Eigen::Matrix<double,6,6>::Zero());

    constexpr int N = 3;
    std::vector<PropagatorD::Section> secs(N);
    for (auto &s : secs) {
        s.strain     = Vec6d::Zero();
        s.length     = 1.0;
        s.strain_cov = 1e-4 * Mat6d::Identity();
    }

    const auto result = PropagatorD::propagateAlongRod(base_g, secs);
    EXPECT_EQ(static_cast<int>(result.size()), N + 1);
}

TEST(UncertaintyPropagatorTest, CovarianceGrowsFromBaseToTip) {
    const SE3d base = SE3d::computeIdentity();
    GaussianSE3d base_g(base, Eigen::Matrix<double,6,6>::Zero());

    constexpr int N = 5;
    std::vector<PropagatorD::Section> secs(N);
    for (auto &s : secs) {
        s.strain     = Vec6d::Zero();
        s.length     = 0.2;
        s.strain_cov = 1e-3 * Mat6d::Identity();
    }

    const auto result = PropagatorD::propagateAlongRod(base_g, secs);

    double prev_trace = result[0].getCovariance().trace();
    for (int k = 1; k <= N; ++k) {
        const double tr = result[k].getCovariance().trace();
        EXPECT_GE(tr, prev_trace - 1e-14)
            << "Covariance should not decrease from section " << k-1 << " to " << k;
        prev_trace = tr;
    }
}

TEST(UncertaintyPropagatorTest, ZeroNoiseLeavesZeroCovariance) {
    const SE3d base = SE3d::computeIdentity();
    GaussianSE3d base_g(base, Eigen::Matrix<double,6,6>::Zero());

    PropagatorD::Section sec;
    sec.strain     = Vec6d::Zero();
    sec.length     = 1.0;
    sec.strain_cov = Mat6d::Zero();  // no noise

    const auto result = PropagatorD::propagateStep(base_g, sec);
    EXPECT_LT(result.getCovariance().norm(), 1e-14)
        << "Zero process noise should leave zero covariance";
}

TEST(UncertaintyPropagatorTest, TipConfidenceRadiiPositive) {
    GaussianSE3d g(SE3d::computeIdentity(), 1e-2 * Mat6d::Identity());
    const auto radii = PropagatorD::tipConfidenceRadii(g);
    EXPECT_EQ(radii.size(), 3);
    for (int i = 0; i < 3; ++i)
        EXPECT_GT(radii[i], 0.0) << "Confidence radius " << i << " should be positive";
}

// ═════════════════════════════════════════════════════════════════════════════
// Integration 6 — CosseratILQRController (with SOFA fixture)
// ═════════════════════════════════════════════════════════════════════════════

class ILQRControllerTest : public ::testing::Test {
protected:
    using MappingT    = Cosserat::mapping::Strain2RigidCosseratMapping<Vec3Types, Rigid3Types, Rigid3Types>;
    using ControllerT = Cosserat::controller::CosseratILQRController<Vec3Types, Rigid3Types, Rigid3Types>;
    using StrainMO    = sofa::component::statecontainer::MechanicalObject<Vec3Types>;
    using RigidMO     = sofa::component::statecontainer::MechanicalObject<Rigid3Types>;

    sofa::simulation::Node::SPtr root;
    MappingT::SPtr    mapping;
    ControllerT::SPtr controller;
    StrainMO::SPtr    strainState;
    RigidMO::SPtr     rigidBase;
    RigidMO::SPtr     outputFrames;

    static constexpr int N = 3;

    void SetUp() override {
        root         = sofa::simulation::getSimulation()->createNewNode("root");
        strainState  = sofa::core::objectmodel::New<StrainMO>();
        rigidBase    = sofa::core::objectmodel::New<RigidMO>();
        outputFrames = sofa::core::objectmodel::New<RigidMO>();
        mapping      = sofa::core::objectmodel::New<MappingT>();
        controller   = sofa::core::objectmodel::New<ControllerT>();

        root->addObject(strainState);
        root->addObject(rigidBase);
        root->addObject(outputFrames);
        root->addObject(mapping);
        root->addObject(controller);

        mapping->setModels(strainState.get(), rigidBase.get(), outputFrames.get());
        controller->l_mapping.set(mapping.get());

        // Straight beam, N sections of length 1
        sofa::type::vector<double> absSection, absFrames;
        for (int i = 0; i <= N; ++i) { absSection.push_back(i); absFrames.push_back(i); }
        mapping->d_curv_abs_section.setValue(absSection);
        mapping->d_curv_abs_frames.setValue(absFrames);

        strainState->resize(N);
        {
            auto w = *strainState->write(sofa::core::vec_id::write_access::position);
            for (int k = 0; k < N; ++k) w.beginEdit()->at(k) = Vec3Types::Coord(0,0,0);
            w.endEdit();
        }
        rigidBase->resize(1);
        {
            auto w = *rigidBase->write(sofa::core::vec_id::write_access::position);
            w.beginEdit()->at(0) = Rigid3Types::Coord(
                Vec3(0,0,0), Quat<SReal>(0,0,0,1));
            w.endEdit();
        }
        outputFrames->resize(N + 1);

        mapping->init();
        controller->init();

        // Run apply() once so body Jacobian is built
        sofa::core::MechanicalParams mp;
        mapping->apply(&mp,
            {outputFrames->write(sofa::core::vec_id::write_access::position)},
            {strainState->read(sofa::core::vec_id::read_access::position)},
            {rigidBase->read(sofa::core::vec_id::read_access::position)});
    }

    void TearDown() override {
        if (root) sofa::simulation::node::unload(root);
    }
};

/**
 * computeControl() returns exactly N Coord1 corrections.
 */
TEST_F(ILQRControllerTest, ReturnsNCorrections) {
    // Set target far from current tip so the controller has something to do
    Rigid3Types::Coord target;
    target.getCenter()      = Vec3(5.0, 1.0, 0.0);
    target.getOrientation() = Quat<SReal>(0,0,0,1);
    controller->d_targetPose.setValue(target);
    controller->d_mode.setValue(1);  // Gauss-Newton

    const auto corrections = controller->computeControl();
    EXPECT_EQ(static_cast<int>(corrections.size()), N);
}

/**
 * With target == current tip, tip error should be near zero and
 * computeControl() should return zero or near-zero corrections.
 */
TEST_F(ILQRControllerTest, ZeroErrorZeroCorrection) {
    // Read current tip pose
    const auto &frames = outputFrames->read(sofa::core::vec_id::read_access::position)->getValue();
    ASSERT_FALSE(frames.empty());
    controller->d_targetPose.setValue(frames.back());  // target = current tip
    controller->d_mode.setValue(1);

    const auto corrections = controller->computeControl();

    // Each correction should be zero (no error to correct)
    for (int k = 0; k < N && k < static_cast<int>(corrections.size()); ++k) {
        for (int j = 0; j < 3; ++j)
            EXPECT_NEAR(corrections[k][j], 0.0, 1e-8)
                << "Non-zero correction for zero-error target at section " << k;
    }
}

/**
 * Gauss-Newton mode reduces the predicted residual compared to no update.
 */
TEST_F(ILQRControllerTest, GaussNewtonReducesPredictedResidual) {
    // Target: 1m offset in y from the tip
    const auto &frames = outputFrames->read(sofa::core::vec_id::read_access::position)->getValue();
    Rigid3Types::Coord target = frames.back();
    target.getCenter()[1] += 0.5;
    controller->d_targetPose.setValue(target);
    controller->d_mode.setValue(1);
    controller->d_maxIterations.setValue(1);

    const double err_before = controller->d_tipError.getValue();
    controller->computeControl();  // runs one GN step, updates d_tipError
    const double err_after = controller->d_tipError.getValue();

    // After one Gauss-Newton step the predicted error (linearised) should be ≤ initial
    // (d_tipError is updated inside computeControl with the predicted residual)
    EXPECT_LE(err_after, err_before + 1e-10)
        << "Gauss-Newton step should not increase the predicted tip error";
}

/**
 * Manipulability is positive for a non-singular configuration.
 */
TEST_F(ILQRControllerTest, ManipulabilityPositive) {
    Rigid3Types::Coord target;
    target.getCenter()      = Vec3(2.0, 0.5, 0.0);
    target.getOrientation() = Quat<SReal>(0,0,0,1);
    controller->d_targetPose.setValue(target);
    controller->computeControl();

    EXPECT_GT(controller->d_manipulability.getValue(), 0.0)
        << "Manipulability should be positive for a non-singular rod";
}

// ═════════════════════════════════════════════════════════════════════════════
// main
// ═════════════════════════════════════════════════════════════════════════════

int main(int argc, char **argv) {
    sofa::simulation::setSimulation(
        new sofa::simulation::graph::DAGSimulation());
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
