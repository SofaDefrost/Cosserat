/******************************************************************************
 *                 SOFA, Simulation Open-Framework Architecture                *
 *                 (c) 2006 INRIA, USTL, UJF, CNRS, MGH                       *
 ******************************************************************************/

/**
 * @file LogExpRoundTripTest.cpp
 *
 * Round-trip tests for the log/exp pair on SO(3) and SE(3).
 *
 *   Given a group element A, we form B = log(A) ∈ algebra, then A' = exp(B).
 *   We then verify  A' == A  to numerical tolerance.
 *
 * These tests are the canonical sanity check for any Lie group implementation
 * and would have caught the 2026-05-31 SE3 bug where V and V⁻¹ used K = [axis]×
 * (unit) with coefficients meant for K = [φ]× (non-unit) — making
 * log(exp(ξ)) ≠ ξ for every θ ≠ 1.
 *
 * Test fixtures cover:
 *   • Identity element
 *   • Pure rotations (no translation)
 *   • Pure translations (no rotation)
 *   • Mixed twists at small θ (Taylor branch)
 *   • Mixed twists at general θ (general branch)
 *   • Boundary of branches (θ ≈ SMALL_ANGLE_THRESHOLD)
 *   • Near singularities (θ → π)
 *   • Random elements (stress test, 100 samples)
 */

#include <gtest/gtest.h>
#include <liegroups/SE3.h>
#include <liegroups/SO3.h>

#include <Eigen/Core>
#include <cmath>
#include <random>

using SO3d = sofa::component::cosserat::liegroups::SO3<double>;
using SE3d = sofa::component::cosserat::liegroups::SE3<double>;

using Vec3d = Eigen::Matrix<double, 3, 1>;
using Vec6d = Eigen::Matrix<double, 6, 1>;
using Mat3d = Eigen::Matrix<double, 3, 3>;

namespace {

/** Build a se(3) twist from angular and linear parts (Cosserat convention: φ first). */
inline Vec6d makeXi(const Vec3d& phi, const Vec3d& rho) {
    Vec6d xi;
    xi.head<3>() = phi;
    xi.tail<3>() = rho;
    return xi;
}

/** Frobenius norm of the rotation discrepancy: ‖R₁·R₂ᵀ − I‖. */
inline double rotationError(const SO3d& a, const SO3d& b) {
    return (a.matrix() * b.matrix().transpose() - Mat3d::Identity()).norm();
}

/** ‖t₁ − t₂‖₂ on translation parts. */
inline double translationError(const Vec3d& t1, const Vec3d& t2) {
    return (t1 - t2).norm();
}

/**
 * Check round-trip on SE3 for a given twist:
 *   ξ → A = exp(ξ) → B = log(A) → A' = exp(B)
 * Verifies both A' ≈ A (group element) and B ≈ ξ (algebra element).
 */
void assertSE3RoundTrip(const Vec6d& xi_in, double tol_rot, double tol_trans,
                        double tol_xi, const std::string& label) {
    const SE3d A     = SE3d::computeExp(xi_in);
    const Vec6d xi_out = A.computeLog();
    const SE3d A_back = SE3d::computeExp(xi_out);

    // (1) Group element round-trip : A' == A
    const double err_rot = rotationError(A.rotation(), A_back.rotation());
    EXPECT_LT(err_rot, tol_rot)
        << label << " — rotation part: ‖A·A'ᵀ − I‖ = " << err_rot;

    // Translation part : extract from internal storage (not exposed publicly, so
    // we compare via action on origin).
    const Vec3d t_A      = A.computeAction(Vec3d::Zero());
    const Vec3d t_A_back = A_back.computeAction(Vec3d::Zero());
    const double err_trans = translationError(t_A, t_A_back);
    EXPECT_LT(err_trans, tol_trans)
        << label << " — translation part: ‖t_A − t_A'‖ = " << err_trans;

    // (2) Algebra round-trip : B == ξ
    const double err_xi = (xi_in - xi_out).norm();
    EXPECT_LT(err_xi, tol_xi)
        << label << " — algebra: ‖ξ_out − ξ_in‖ = " << err_xi
        << "\n    ξ_in  = (" << xi_in.transpose() << ")"
        << "\n    ξ_out = (" << xi_out.transpose() << ")";
}

/** Same idea for SO3 : ω → R = exp(ω) → ω' = log(R) → R' = exp(ω'). */
void assertSO3RoundTrip(const Vec3d& omega_in, double tol_rot, double tol_omega,
                        const std::string& label) {
    const SO3d R       = SO3d::exp(omega_in);
    const Vec3d omega_out = R.log();
    const SO3d R_back  = SO3d::exp(omega_out);

    const double err_rot = rotationError(R, R_back);
    EXPECT_LT(err_rot, tol_rot)
        << label << " — ‖R·R'ᵀ − I‖ = " << err_rot;

    const double err_omega = (omega_in - omega_out).norm();
    EXPECT_LT(err_omega, tol_omega)
        << label << " — ‖ω_out − ω_in‖ = " << err_omega
        << "\n    ω_in  = (" << omega_in.transpose() << ")"
        << "\n    ω_out = (" << omega_out.transpose() << ")";
}

}  // anonymous namespace

// ═════════════════════════════════════════════════════════════════════════════
// SO(3) round-trip tests
// ═════════════════════════════════════════════════════════════════════════════

TEST(SO3LogExpRoundTrip, Identity) {
    assertSO3RoundTrip(Vec3d::Zero(), 1e-14, 1e-14, "identity");
}

TEST(SO3LogExpRoundTrip, SmallRotationsTaylorBranch) {
    // |ω| below SMALL_ANGLE_THRESHOLD (1e-4) → first-order quaternion branch
    for (double theta : {1e-12, 1e-8, 1e-6, 1e-5, 1e-4}) {
        const Vec3d axis = Vec3d(0.3, 0.5, -0.4).normalized();
        const Vec3d omega = axis * theta;
        assertSO3RoundTrip(omega, 1e-13, std::max(1e-15, 10.0 * theta * 1e-15),
                           "small ω θ=" + std::to_string(theta));
    }
}

TEST(SO3LogExpRoundTrip, GeneralBranch) {
    // |ω| ∈ [1e-3, 2.5] → AngleAxis branch with comfortable margin from the
    // log singularity at θ = π. The near-π regime is tested separately in
    // NearPiSingularity with widened tolerances.
    const Vec3d axis = Vec3d(0.5, -0.7, 0.2).normalized();
    for (double theta : {1e-3, 0.1, 0.5, 1.0, 1.5, 2.0, 2.5}) {
        const Vec3d omega = axis * theta;
        assertSO3RoundTrip(omega, 1e-13, 1e-12,
                           "general ω θ=" + std::to_string(theta));
    }
}

TEST(SO3LogExpRoundTrip, NearPiSingularity) {
    // sin(θ) → 0 in the log formula → loss of precision near π. We verify
    // graceful degradation rather than full precision.
    const Vec3d axis = Vec3d(0.5, -0.7, 0.2).normalized();
    for (double theta : {3.0, M_PI - 1e-2, M_PI - 1e-4}) {
        // Tolerance scales roughly with √(π-θ) since the V⁻¹ coefficient
        // diverges like 1/(π-θ). A 1e-6 tolerance is generous.
        const Vec3d omega = axis * theta;
        assertSO3RoundTrip(omega, 1e-6, 1e-6,
                           "near π ω θ=" + std::to_string(theta));
    }
}

TEST(SO3LogExpRoundTrip, RandomStress100) {
    std::mt19937 rng(42);
    std::uniform_real_distribution<double> dist_angle(-M_PI + 0.01, M_PI - 0.01);
    std::uniform_real_distribution<double> dist_axis(-1.0, 1.0);

    for (int k = 0; k < 100; ++k) {
        Vec3d axis(dist_axis(rng), dist_axis(rng), dist_axis(rng));
        if (axis.norm() < 1e-10) continue;
        axis.normalize();
        const double theta = dist_angle(rng);
        const Vec3d omega = axis * theta;
        assertSO3RoundTrip(omega, 1e-12, 1e-11,
                           "random k=" + std::to_string(k));
    }
}

// ═════════════════════════════════════════════════════════════════════════════
// SE(3) round-trip tests — these would have caught the V/V_inv bug
// ═════════════════════════════════════════════════════════════════════════════

TEST(SE3LogExpRoundTrip, Identity) {
    assertSE3RoundTrip(Vec6d::Zero(), 1e-14, 1e-14, 1e-14, "identity");
}

TEST(SE3LogExpRoundTrip, PureRotationNoTranslation) {
    // Pure rotation: ρ = 0. Even with the old V bug, V·0 = 0, so this passed.
    // Included for completeness — should remain robust.
    // θ = 3.0 (close to π) is excluded — tested in NearPiSingularity below.
    const Vec3d axis = Vec3d(1.0, 1.0, 1.0).normalized();
    for (double theta : {0.1, 0.5, 1.0, 2.0}) {
        const Vec6d xi = makeXi(axis * theta, Vec3d::Zero());
        assertSE3RoundTrip(xi, 1e-13, 1e-13, 1e-12,
                           "pure rotation θ=" + std::to_string(theta));
    }
}

TEST(SE3LogExpRoundTrip, PureTranslationNoRotation) {
    // φ = 0 → V = I trivially. Should be exact.
    const Vec6d xi = makeXi(Vec3d::Zero(), Vec3d(1.5, -2.3, 0.7));
    assertSE3RoundTrip(xi, 1e-14, 1e-14, 1e-14, "pure translation");
}

/**
 * @brief THE BUG SENTINEL : mixed twist at θ = 0.5
 *
 * Before the 2026-05-31 fix, V_code · V_inv_code ≠ I for θ ≠ 1, so the
 * round-trip diverged with error ~3.6e-1. This test would have caught it
 * immediately.
 */
TEST(SE3LogExpRoundTrip, MixedTwistTheta05_BUG_SENTINEL) {
    const Vec6d xi = makeXi(Vec3d(0.5, 0.0, 0.0), Vec3d(1.0, 2.0, -0.5));
    assertSE3RoundTrip(xi, 1e-12, 1e-12, 1e-11, "mixed θ=0.5 (bug sentinel)");
}

TEST(SE3LogExpRoundTrip, MixedTwistTheta1) {
    // θ = 1 was the ONLY angle where the old buggy code happened to work.
    // After the fix, it should still work.
    const Vec6d xi = makeXi(Vec3d(1.0, 0.0, 0.0), Vec3d(2.0, 0.5, -1.0));
    assertSE3RoundTrip(xi, 1e-12, 1e-12, 1e-11, "mixed θ=1.0");
}

TEST(SE3LogExpRoundTrip, MixedTwistTaylorBranch) {
    // θ < SMALL_ANGLE_THRESHOLD → forces the Taylor branch we added.
    // ρ is non-trivial to also exercise V·ρ ≠ 0.
    const Vec3d axis = Vec3d(0.7, -0.4, 0.6).normalized();
    for (double theta : {1e-6, 1e-5, 1e-4, 0.99e-4}) {
        const Vec6d xi = makeXi(axis * theta, Vec3d(0.3, -0.2, 0.5));
        assertSE3RoundTrip(xi, 1e-12, 1e-12, 1e-11,
                           "Taylor branch θ=" + std::to_string(theta));
    }
}

TEST(SE3LogExpRoundTrip, BranchBoundaryAgreement) {
    // At θ ≈ SMALL_ANGLE_THRESHOLD, Taylor and general branches must agree.
    // We round-trip just above and just below the threshold and check the
    // result is continuous.
    const Vec3d axis = Vec3d(1.0, 0.0, 0.0);
    const Vec3d rho(0.5, -1.0, 0.2);

    const Vec6d xi_low  = makeXi(axis * 9.9e-5, rho);
    const Vec6d xi_high = makeXi(axis * 1.01e-4, rho);

    const Vec6d back_low  = SE3d::computeExp(xi_low).computeLog();
    const Vec6d back_high = SE3d::computeExp(xi_high).computeLog();

    EXPECT_LT((xi_low  - back_low).norm(),  1e-11) << "Taylor branch boundary";
    EXPECT_LT((xi_high - back_high).norm(), 1e-11) << "general branch boundary";
}

TEST(SE3LogExpRoundTrip, NearPiSingularity) {
    // θ close to π is a known difficult region : sin(θ) → 0 in V⁻¹ coefficient.
    // We verify graceful degradation only. θ = π is a genuine log singularity
    // (the rotation axis is ambiguous, sign flip possible).
    //
    // Tolerance follows the empirical loss-of-precision rate ~ 1/(π-θ).
    const Vec3d axis = Vec3d(0.3, 0.5, -0.4).normalized();
    struct Case { double theta; double tol; };
    for (const auto& c : {Case{3.0,         1e-10},
                          Case{M_PI - 1e-2, 1e-10},
                          Case{M_PI - 1e-4, 1e-6}}) {
        const Vec6d xi = makeXi(axis * c.theta, Vec3d(0.5, -0.3, 0.2));
        assertSE3RoundTrip(xi, c.tol, c.tol, c.tol,
                           "near π θ=" + std::to_string(c.theta));
    }
}

TEST(SE3LogExpRoundTrip, RandomStress100) {
    std::mt19937 rng(12345);
    std::uniform_real_distribution<double> dist_angle(-2.5, 2.5);
    std::uniform_real_distribution<double> dist_coord(-3.0, 3.0);

    int failed = 0;
    for (int k = 0; k < 100; ++k) {
        Vec3d axis(dist_coord(rng), dist_coord(rng), dist_coord(rng));
        if (axis.norm() < 1e-10) continue;
        axis.normalize();
        const Vec3d phi = axis * dist_angle(rng);
        const Vec3d rho(dist_coord(rng), dist_coord(rng), dist_coord(rng));

        const Vec6d xi = makeXi(phi, rho);
        const SE3d A      = SE3d::computeExp(xi);
        const Vec6d xi_out = A.computeLog();
        const SE3d A_back = SE3d::computeExp(xi_out);

        const double err_group   = rotationError(A.rotation(), A_back.rotation());
                                  + translationError(A.computeAction(Vec3d::Zero()),
                                                     A_back.computeAction(Vec3d::Zero()));
        const double err_algebra = (xi - xi_out).norm();
        if (err_group > 1e-10 || err_algebra > 1e-10) {
            ++failed;
            ADD_FAILURE() << "random k=" << k
                          << "  err_group=" << err_group
                          << "  err_algebra=" << err_algebra
                          << "  ξ=(" << xi.transpose() << ")";
        }
    }
    EXPECT_EQ(failed, 0) << "Number of failed random round-trips";
}

// ═════════════════════════════════════════════════════════════════════════════
// Direct matrix-level check : A = exp(log(A))  on a concrete matrix
// ═════════════════════════════════════════════════════════════════════════════

/**
 * @brief « Take a matrix A and show that, if B = log(A), then A = exp(B) ».
 *
 * Builds a specific SE(3) element A (rotation 30° around an arbitrary axis,
 * translation (1, 0.5, -0.2)) and verifies the matrix-level identity by
 * comparing the 4×4 homogeneous matrices.
 */
TEST(SE3LogExpRoundTrip, MatrixLevelIdentity_A_equals_exp_log_A) {
    // ── Build A by hand : rotation 30° around (1, 1, 0)/√2, translation t
    const double angle = M_PI / 6.0;                              // 30°
    const Vec3d axis(1.0, 1.0, 0.0);
    const Vec3d axis_n = axis.normalized();
    const SO3d  R(angle, axis_n);
    const Vec3d t(1.0, 0.5, -0.2);
    const SE3d  A(R, t);

    // ── A's 4×4 homogeneous matrix (reference)
    Eigen::Matrix4d M_A = Eigen::Matrix4d::Identity();
    M_A.block<3, 3>(0, 0) = R.matrix();
    M_A.block<3, 1>(0, 3) = t;

    // ── Round-trip : B = log(A), A' = exp(B)
    const Vec6d B    = A.computeLog();
    const SE3d  Ap   = SE3d::computeExp(B);

    Eigen::Matrix4d M_Ap = Eigen::Matrix4d::Identity();
    M_Ap.block<3, 3>(0, 0) = Ap.matrix().block<3, 3>(0, 0);
    M_Ap.block<3, 1>(0, 3) = Ap.computeAction(Vec3d::Zero());

    // ── A == A' at matrix level
    const double err = (M_A - M_Ap).norm();
    EXPECT_LT(err, 1e-13)
        << "‖A − exp(log(A))‖_F = " << err
        << "\nA   =\n" << M_A
        << "\nA'  =\n" << M_Ap
        << "\nlog(A) = " << B.transpose();
}
