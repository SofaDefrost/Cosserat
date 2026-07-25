#include <gtest/gtest.h>
#include <liegroups/LieGroupIntegrators.h>
#include <liegroups/SE3.h>
#include <cmath>
#include <vector>

using namespace sofa::component::cosserat::liegroups;
using SE3d = SE3<double>;
using Integrator = SE3Integrator<double>;
using TangentVector = SE3d::TangentVector;

// ── Helpers ───────────────────────────────────────────────────────────────────

// Pose error metric: norm of the log of the relative transformation
static double poseError(const SE3d& g_ref, const SE3d& g_approx) {
    return g_ref.inverse().compose(g_approx).log().norm();
}

// Constant strain field (piecewise-constant elements)
static Integrator::StrainField constantField(const TangentVector& xi) {
    return [xi](double) { return xi; };
}

// Linearly varying strain field: xi(s) = xi0 + s * dxi
static Integrator::StrainField linearField(const TangentVector& xi0, const TangentVector& dxi) {
    return [xi0, dxi](double s) { return xi0 + s * dxi; };
}

// ── Constant strain — all methods should give identical results ───────────────

TEST(IntegratorTest, ConstantStrainEulerEqualsExpCosserat) {
    TangentVector strain;
    strain << 0.1, 0.0, 0.2,   // angular
              0.0, 0.0, 0.0;   // linear (pure bending)
    double ds = 0.05;

    SE3d g0 = SE3d::Identity();
    SE3d g_euler = Integrator::integrateEuler(g0, constantField(strain), 0.0, ds);
    SE3d g_exp   = g0.compose(SE3d::expCosserat(strain, ds));

    EXPECT_NEAR(poseError(g_exp, g_euler), 0.0, 1e-14);
}

TEST(IntegratorTest, ConstantStrainAllMethodsAgree) {
    TangentVector strain;
    strain << 0.2, 0.1, 0.05,  // angular
              0.01, 0.0, 0.0;  // linear
    double ds = 0.1;

    SE3d g0 = SE3d::Identity();
    SE3d g_euler  = Integrator::integrateEuler(   g0, constantField(strain), 0.0, ds);
    SE3d g_mid    = Integrator::integrateMidpoint( g0, constantField(strain), 0.0, ds);
    SE3d g_rkmk4  = Integrator::integrateRKMK4(   g0, constantField(strain), 0.0, ds);

    // For constant strain [k1,k4]=0, commutator vanishes → all three agree
    EXPECT_NEAR(poseError(g_euler, g_mid),   0.0, 1e-14);
    EXPECT_NEAR(poseError(g_euler, g_rkmk4), 0.0, 1e-14);
}

// ── Order verification — varying strain ──────────────────────────────────────
//
// Strategy: compute a "reference" solution with RKMK4 and many sub-steps,
// then compare coarse single-step solutions. The error ratio between step
// sizes h and 2h should be ~2^p for a p-th order method.

class IntegratorOrderTest : public ::testing::Test {
protected:
    // Linearly varying strain over [0, L]
    TangentVector xi0, dxi;
    double L = 0.2;

    void SetUp() override {
        xi0 << 0.3, 0.1, 0.05,   // angular strain at s=0
               0.02, 0.0, 0.01;  // linear strain at s=0
        dxi << 0.2, -0.1, 0.1,   // variation per unit length
               0.01, 0.02, 0.0;
    }

    // High-accuracy reference: RKMK4 with 1000 sub-steps
    SE3d reference() {
        return Integrator::integrate(
            SE3d::Identity(), linearField(xi0, dxi), 0.0, L, 1000,
            Integrator::Method::RKMK4);
    }

    // Single-step error for a given method and step size h
    double singleStepError(Integrator::Method method, double h) {
        SE3d g_ref = Integrator::integrate(
            SE3d::Identity(), linearField(xi0, dxi), 0.0, h, 1000,
            Integrator::Method::RKMK4);
        SE3d g_approx = Integrator::integrate(
            SE3d::Identity(), linearField(xi0, dxi), 0.0, h, 1,
            method);
        return poseError(g_ref, g_approx);
    }
};

TEST_F(IntegratorOrderTest, EulerIsOrder1) {
    double h1 = 0.04, h2 = 0.02;
    double e1 = singleStepError(Integrator::Method::Euler, h1);
    double e2 = singleStepError(Integrator::Method::Euler, h2);
    double observed_order = std::log(e1 / e2) / std::log(h1 / h2);
    // Expect order ≈ 1, allow generous tolerance due to commutator effects
    EXPECT_GT(observed_order, 0.7);
    EXPECT_LT(observed_order, 2.0);
}

TEST_F(IntegratorOrderTest, MidpointIsOrder2) {
    double h1 = 0.08, h2 = 0.04;
    double e1 = singleStepError(Integrator::Method::Midpoint, h1);
    double e2 = singleStepError(Integrator::Method::Midpoint, h2);
    double observed_order = std::log(e1 / e2) / std::log(h1 / h2);
    // Expect order ≈ 2
    EXPECT_GT(observed_order, 1.5);
    EXPECT_LT(observed_order, 3.0);
}

TEST_F(IntegratorOrderTest, RKMK4IsOrder4) {
    double h1 = 0.1, h2 = 0.05;
    double e1 = singleStepError(Integrator::Method::RKMK4, h1);
    double e2 = singleStepError(Integrator::Method::RKMK4, h2);
    double observed_order = std::log(e1 / e2) / std::log(h1 / h2);
    // Expect order ≈ 4
    EXPECT_GT(observed_order, 3.0);
    EXPECT_LT(observed_order, 5.5);
}

TEST_F(IntegratorOrderTest, RKMK4MoreAccurateThanEuler) {
    // At a fixed step, RKMK4 must be strictly more accurate than Euler
    double h = 0.1;
    double err_euler = singleStepError(Integrator::Method::Euler,   h);
    double err_rkmk4 = singleStepError(Integrator::Method::RKMK4,  h);
    EXPECT_LT(err_rkmk4, err_euler);
}

// ── integrate() chains multiple steps ────────────────────────────────────────

TEST(IntegratorTest, MultiStepEquivalentToChainedSingleSteps) {
    TangentVector strain;
    strain << 0.1, 0.05, 0.0,
              0.0, 0.0,  0.0;
    double L = 0.5;
    int N = 5;

    auto field = constantField(strain);
    SE3d g0 = SE3d::Identity();

    SE3d g_multi = Integrator::integrate(g0, field, 0.0, L, N, Integrator::Method::RKMK4);

    // Chain manually
    SE3d g_manual = g0;
    double ds = L / N;
    for (int i = 0; i < N; ++i)
        g_manual = Integrator::integrateRKMK4(g_manual, field, i * ds, ds);

    EXPECT_NEAR(poseError(g_multi, g_manual), 0.0, 1e-13);
}

// ── integratePath() ───────────────────────────────────────────────────────────

TEST(IntegratorTest, IntegratePathSize) {
    TangentVector strain;
    strain << 0.1, 0.0, 0.0, 0.0, 0.0, 0.0;
    int N = 10;

    auto path = Integrator::integratePath(
        SE3d::Identity(), constantField(strain), 0.0, 1.0, N,
        Integrator::Method::RKMK4);

    EXPECT_EQ(static_cast<int>(path.size()), N + 1);
    EXPECT_TRUE(path.front().isApprox(SE3d::Identity()));
}

TEST(IntegratorTest, IntegratePathFinalMatchesIntegrate) {
    TangentVector strain;
    strain << 0.05, 0.1, 0.0,
              0.0,  0.0, 0.01;
    int N = 8;
    double L = 0.4;

    auto field = constantField(strain);
    SE3d g_final = Integrator::integrate(
        SE3d::Identity(), field, 0.0, L, N, Integrator::Method::RKMK4);

    auto path = Integrator::integratePath(
        SE3d::Identity(), field, 0.0, L, N, Integrator::Method::RKMK4);

    EXPECT_NEAR(poseError(g_final, path.back()), 0.0, 1e-14);
}

// ── Pure rotation — known analytical solution ─────────────────────────────────

TEST(IntegratorTest, PureConstantTorsionAnalyticalSolution) {
    // Constant torsion phi_x = kappa, no other strain
    // g(s) = (Rx(kappa*s), 0)  — rotation around X by kappa*s
    double kappa = 0.5;  // rad/m
    double L     = 0.4;  // m

    TangentVector strain = TangentVector::Zero();
    strain[0] = kappa;  // phi_x = torsion

    SE3d g_rkmk4 = Integrator::integrate(
        SE3d::Identity(), constantField(strain), 0.0, L, 1,
        Integrator::Method::RKMK4);

    // Analytical: pure rotation Rx(kappa*L)
    Eigen::AngleAxisd rot(kappa * L, Eigen::Vector3d::UnitX());
    SE3d g_analytic(SO3<double>(rot.toRotationMatrix()), Eigen::Vector3d::Zero());

    EXPECT_NEAR(poseError(g_analytic, g_rkmk4), 0.0, 1e-12);
}

// ── Backward compatibility ────────────────────────────────────────────────────

TEST(IntegratorTest, EulerMatchesExpCosseratWithNonzeroLinearStrain) {
    // Elongation along X
    TangentVector strain;
    strain << 0.0, 0.0, 0.0,
              0.05, 0.0, 0.0;  // rho_x = elongation deviation
    double ds = 0.1;

    SE3d g0 = SE3d::Identity();
    SE3d g_euler = Integrator::integrateEuler(g0, constantField(strain), 0.0, ds);
    SE3d g_exp   = g0.compose(SE3d::expCosserat(strain, ds));

    EXPECT_NEAR(poseError(g_exp, g_euler), 0.0, 1e-14);
}
