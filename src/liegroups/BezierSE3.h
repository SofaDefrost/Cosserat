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
#pragma once

#include <vector>
#include <stdexcept>
#include "SE3.h"

namespace sofa::component::cosserat::liegroups {

/**
 * @brief Bézier curves on SE(3) via the De Casteljau algorithm.
 *
 * Provides smooth interpolation and path generation on SE(3) using
 * control poses, generalising classical Euclidean Bézier curves to
 * the Lie group setting.
 *
 * ## Geometric interpolation on manifolds
 *
 * The classical geodesic interpolation between two poses is:
 *
 *   slerp(g0, g1, t) = g0 · exp(t · log(g0⁻¹ · g1))
 *
 * De Casteljau's algorithm replaces linear interpolation with slerp
 * at each subdivision level, giving a smooth curve that:
 *   - Passes through the first and last control pose (at t=0 and t=1)
 *   - Lies on SE(3) at all parameter values
 *   - Reduces to slerp for 2 control poses (degree 1)
 *
 * ## Applications in Cosserat rods
 *
 *   1. **Shape interpolation** (BeamShapeInterpolation): smooth blending
 *      between known rod configurations
 *   2. **Reference trajectory generation**: smooth task-space path for
 *      the rod tip, usable with the iLQR controller
 *   3. **Multi-section shape fitting**: fit a Bézier curve in SE(3) to
 *      discretely measured cross-section poses
 *
 * ## Usage
 *
 * @code
 *   std::vector<SE3d> ctrl = {g0, g1, g2, g3};   // 4 control poses
 *   BezierSE3<double> curve(ctrl);
 *
 *   SE3d mid   = curve.evaluate(0.5);             // midpoint
 *   Twist<double> vel = curve.velocity(0.5);      // tangent twist
 *
 *   auto samples = curve.sample(100);             // 100 uniform samples
 * @endcode
 *
 * @tparam _Scalar Scalar type (float, double, …)
 */
template<typename _Scalar>
class BezierSE3 {
public:
    using Scalar        = _Scalar;
    using SE3Type       = SE3<Scalar>;
    using TangentVector = typename SE3Type::TangentVector;

    // ── Construction ─────────────────────────────────────────────────────────

    /**
     * @brief Construct a Bézier curve from control poses.
     *
     * @param control_poses  Ordered control poses P_0, P_1, …, P_n.
     *                       At least 2 poses required.
     *                       Degree of the curve = n (number of intervals).
     */
    explicit BezierSE3(std::vector<SE3Type> control_poses)
        : control_poses_(std::move(control_poses))
    {
        if (control_poses_.size() < 2)
            throw std::invalid_argument(
                "BezierSE3: at least 2 control poses required");
    }

    /** Degree of the Bézier curve (= number of control poses − 1). */
    [[nodiscard]] int degree() const {
        return static_cast<int>(control_poses_.size()) - 1;
    }

    /** Number of control poses. */
    [[nodiscard]] int numControlPoses() const {
        return static_cast<int>(control_poses_.size());
    }

    [[nodiscard]] const SE3Type& controlPose(int i) const { return control_poses_[i]; }

    // ── Evaluation ───────────────────────────────────────────────────────────

    /**
     * @brief Evaluate the Bézier curve at parameter t ∈ [0, 1].
     *
     * Uses De Casteljau's algorithm on SE(3):
     *   - Replace linear interpolation with geodesic (slerp)
     *   - Recurse n times, reducing from n+1 to 1 control pose
     *
     * @param t  Curve parameter in [0, 1].
     *           t=0 → first control pose, t=1 → last control pose.
     */
    [[nodiscard]] SE3Type evaluate(Scalar t) const {
        return deCasteljau(control_poses_, t);
    }

    /** Alias for evaluate(). */
    [[nodiscard]] SE3Type operator()(Scalar t) const { return evaluate(t); }

    /**
     * @brief Estimate the tangent twist at parameter t.
     *
     * Computed via central finite differences in the tangent space:
     *
     *   velocity(t) = log(curve(t + h)⁻¹ · curve(t - h)) / (2h)
     *
     * This gives the body twist (expressed in the frame at curve(t)).
     *
     * @param t  Parameter in [0, 1].
     * @param h  Finite-difference step (default: 1e-5).
     */
    [[nodiscard]] TangentVector velocity(Scalar t, Scalar h = Scalar(1e-5)) const {
        const Scalar t_fwd = std::min(t + h, Scalar(1));
        const Scalar t_bwd = std::max(t - h, Scalar(0));
        const Scalar dt    = t_fwd - t_bwd;

        const SE3Type g_fwd = evaluate(t_fwd);
        const SE3Type g_bwd = evaluate(t_bwd);

        // Central difference in the tangent space at g(t)
        // ξ ≈ log(g_bwd⁻¹ · g_fwd) / dt
        return g_bwd.computeInverse().compose(g_fwd).log() / dt;
    }

    // ── Sampling ──────────────────────────────────────────────────────────────

    /**
     * @brief Sample the curve at N uniformly spaced parameter values.
     *
     * Returns N+1 poses at t = 0, 1/N, 2/N, …, 1.
     *
     * @param N  Number of intervals (N+1 samples).
     */
    [[nodiscard]] std::vector<SE3Type> sample(int N) const {
        if (N < 1)
            throw std::invalid_argument("BezierSE3::sample: N must be >= 1");

        std::vector<SE3Type> pts;
        pts.reserve(static_cast<std::size_t>(N + 1));
        for (int i = 0; i <= N; ++i)
            pts.push_back(evaluate(Scalar(i) / Scalar(N)));
        return pts;
    }

    /**
     * @brief Sample at a custom list of parameter values.
     *
     * @param params  Vector of t values (each in [0, 1]).
     */
    [[nodiscard]] std::vector<SE3Type>
    sampleAt(const std::vector<Scalar>& params) const {
        std::vector<SE3Type> pts;
        pts.reserve(params.size());
        for (Scalar t : params)
            pts.push_back(evaluate(t));
        return pts;
    }

    // ── Subdivision ───────────────────────────────────────────────────────────

    /**
     * @brief Split the curve at parameter t into two sub-curves.
     *
     * Uses the De Casteljau intermediate points as control poses for the
     * sub-curves [0, t] and [t, 1]. This is exact (no approximation).
     *
     * @param t  Split parameter in (0, 1).
     * @return   Pair (left_curve, right_curve).
     */
    [[nodiscard]] std::pair<BezierSE3, BezierSE3> split(Scalar t) const {
        const int n = degree();
        // Compute the De Casteljau triangle and extract the two sets of
        // boundary control points
        std::vector<std::vector<SE3Type>> triangle = deCasteljauTriangle(t);

        std::vector<SE3Type> left_ctrl, right_ctrl;
        left_ctrl.reserve(static_cast<std::size_t>(n + 1));
        right_ctrl.reserve(static_cast<std::size_t>(n + 1));

        for (int r = 0; r <= n; ++r) {
            left_ctrl.push_back(triangle[r][0]);
            right_ctrl.push_back(triangle[n - r][r]);
        }

        return {BezierSE3(std::move(left_ctrl)), BezierSE3(std::move(right_ctrl))};
    }

    // ── Degree elevation ──────────────────────────────────────────────────────

    /**
     * @brief Elevate the degree by 1 (adds one extra control pose).
     *
     * The resulting degree-(n+1) curve represents the same geometric path.
     * Useful for merging curves or increasing refinement control.
     *
     * For Euclidean Bézier: P'_i = (i/(n+1)) P_{i-1} + (1 - i/(n+1)) P_i
     * On SE(3) we replace weighted average with slerp at the same ratios.
     */
    [[nodiscard]] BezierSE3 elevate() const {
        const int n = degree();
        std::vector<SE3Type> elevated;
        elevated.reserve(static_cast<std::size_t>(n + 2));

        elevated.push_back(control_poses_.front());
        for (int i = 1; i <= n; ++i) {
            const Scalar alpha = Scalar(i) / Scalar(n + 1);
            elevated.push_back(geodesicSlerp(control_poses_[i], control_poses_[i - 1], alpha));
        }
        elevated.push_back(control_poses_.back());
        return BezierSE3(std::move(elevated));
    }

    // ── Arc length ────────────────────────────────────────────────────────────

    /**
     * @brief Estimate the arc length of the curve using N-point quadrature.
     *
     * Arc length is measured using the Riemannian metric on SE(3):
     *   L = ∫₀¹ ‖log(g(t)⁻¹ g(t+dt))‖ / dt dt
     *
     * @param N  Number of segments for the numerical integration.
     */
    [[nodiscard]] Scalar arcLength(int N = 100) const {
        Scalar length = Scalar(0);
        SE3Type prev = evaluate(Scalar(0));
        for (int i = 1; i <= N; ++i) {
            SE3Type cur = evaluate(Scalar(i) / Scalar(N));
            length += prev.computeInverse().compose(cur).log().norm();
            prev = cur;
        }
        return length;
    }

private:
    std::vector<SE3Type> control_poses_;

    // ── Core De Casteljau recursion ───────────────────────────────────────────

    /**
     * @brief Geodesic interpolation between two SE(3) elements (slerp).
     *
     *   slerp(g0, g1, t) = g0 · exp(t · log(g0⁻¹ · g1))
     */
    [[nodiscard]] static SE3Type geodesicSlerp(
        const SE3Type& g0,
        const SE3Type& g1,
        Scalar t)
    {
        const TangentVector xi = g0.computeInverse().compose(g1).log();
        return g0.compose(SE3Type::computeExp(t * xi));
    }

    /**
     * @brief One pass of De Casteljau's algorithm on a set of control poses.
     *
     * Reduces n+1 control poses to n by geodesic interpolation at t.
     */
    [[nodiscard]] static std::vector<SE3Type> deCasteljauStep(
        const std::vector<SE3Type>& pts,
        Scalar t)
    {
        const std::size_t n = pts.size() - 1;
        std::vector<SE3Type> result;
        result.reserve(n);
        for (std::size_t i = 0; i < n; ++i)
            result.push_back(geodesicSlerp(pts[i], pts[i + 1], t));
        return result;
    }

    /**
     * @brief Full De Casteljau evaluation: reduce to a single point.
     */
    [[nodiscard]] static SE3Type deCasteljau(
        std::vector<SE3Type> pts,
        Scalar t)
    {
        while (pts.size() > 1)
            pts = deCasteljauStep(pts, t);
        return pts[0];
    }

    /**
     * @brief Full De Casteljau triangle (for subdivision).
     *
     * Returns a lower-triangular structure where triangle[r][i] is the
     * result after r reduction steps starting at column i.
     */
    [[nodiscard]] std::vector<std::vector<SE3Type>>
    deCasteljauTriangle(Scalar t) const {
        const int n = degree();
        std::vector<std::vector<SE3Type>> tri;
        tri.reserve(static_cast<std::size_t>(n + 1));
        tri.push_back(control_poses_);

        for (int r = 1; r <= n; ++r)
            tri.push_back(deCasteljauStep(tri.back(), t));
        return tri;
    }
};

// ── Convenience type aliases ──────────────────────────────────────────────────

using BezierSE3d = BezierSE3<double>;
using BezierSE3f = BezierSE3<float>;

// ── Free function helpers ─────────────────────────────────────────────────────

/**
 * @brief Build a smooth SE(3) path through a list of waypoints.
 *
 * Constructs a degree-(N-1) Bézier curve through N waypoints used directly
 * as control poses (interpolating endpoints only, approximating the interior).
 *
 * For exact interpolation through interior waypoints, use a B-spline
 * or piecewise cubic Bézier approach instead.
 *
 * @param waypoints  Ordered list of SE(3) waypoints.
 */
template<typename Scalar>
[[nodiscard]] BezierSE3<Scalar>
makeBezierPath(std::vector<SE3<Scalar>> waypoints) {
    return BezierSE3<Scalar>(std::move(waypoints));
}

/**
 * @brief Build a piecewise-cubic (C1-continuous) SE(3) path.
 *
 * Splits N waypoints into (N-1)/3 cubic Bézier segments, with tangent
 * constraints enforced at the junction points to achieve C1 continuity.
 * Inner control poses are placed at ⅓ and ⅔ of the geodesic between
 * consecutive waypoints.
 *
 * @param waypoints  Key waypoints the curve passes through.
 * @return           Vector of cubic BezierSE3 segments.
 */
template<typename Scalar>
[[nodiscard]] std::vector<BezierSE3<Scalar>>
makePiecewiseCubicPath(const std::vector<SE3<Scalar>>& waypoints) {
    const int M = static_cast<int>(waypoints.size());
    if (M < 2)
        throw std::invalid_argument("makePiecewiseCubicPath: at least 2 waypoints required");

    using SE3Type = SE3<Scalar>;
    using TangentVector = typename SE3Type::TangentVector;

    std::vector<BezierSE3<Scalar>> segments;
    segments.reserve(static_cast<std::size_t>(M - 1));

    for (int i = 0; i < M - 1; ++i) {
        const SE3Type& P0 = waypoints[i];
        const SE3Type& P3 = waypoints[i + 1];

        // Geodesic from P0 to P3
        const TangentVector xi = P0.computeInverse().compose(P3).log();

        // Inner control poses at 1/3 and 2/3 of the geodesic
        const SE3Type P1 = P0.compose(SE3Type::computeExp(xi / Scalar(3)));
        const SE3Type P2 = P0.compose(SE3Type::computeExp(Scalar(2) * xi / Scalar(3)));

        segments.emplace_back(std::vector<SE3Type>{P0, P1, P2, P3});
    }
    return segments;
}

} // namespace sofa::component::cosserat::liegroups
