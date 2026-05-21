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

#include <functional>
#include "SE3.h"
#include "Twist.h"

namespace sofa::component::cosserat::liegroups {

/**
 * @brief Geometric integrators for ODEs on SE(3).
 *
 * Solves the right-invariant ODE describing Cosserat rod kinematics:
 *
 *   g'(s) = g(s) · hat(strain(s))
 *
 * where g(s) ∈ SE(3) is the cross-section pose at arc-length s, and
 * strain(s) ∈ se(3) is the 6D strain vector (angular, linear).
 *
 * ## Why geometric integrators matter
 *
 * The standard expCosserat() is an Euler step: exact only when strain is
 * constant on the element. For varying strains (e.g. Legendre polynomial
 * parameterisation), higher-order methods give significantly better
 * accuracy without reducing the step size.
 *
 * ## Available methods
 *
 * | Method         | Order | Notes                                  |
 * |----------------|-------|----------------------------------------|
 * | Euler          | 1     | Equivalent to expCosserat()            |
 * | Midpoint       | 2     | One extra strain evaluation            |
 * | RKMK4          | 4     | Magnus expansion + commutator correction|
 *
 * ## Usage (piecewise-constant strain — existing use case)
 * @code
 *   // Constant strain over the element
 *   auto strain_field = [&strain](double) { return strain; };
 *   SE3d g1 = SE3Integrator<double>::integrateRKMK4(g0, strain_field, s0, ds);
 * @endcode
 *
 * ## Usage (Legendre polynomial strain)
 * @code
 *   auto strain_field = [&coeffs, order](double s) {
 *       return evalLegendreStrain(coeffs, order, s);
 *   };
 *   SE3d g1 = SE3Integrator<double>::integrateRKMK4(g0, strain_field, s0, ds);
 * @endcode
 *
 * @tparam _Scalar Scalar type (float, double, autodiff dual, …)
 */
template<typename _Scalar>
class SE3Integrator {
public:
    using Scalar        = _Scalar;
    using SE3Type       = SE3<Scalar>;
    using TangentVector = typename SE3Type::TangentVector;
    using Vector3       = typename SE3Type::Vector3;
    using Matrix6       = typename Types<Scalar>::Matrix6;

    /** Callable: arc-length → strain vector ∈ se(3). */
    using StrainField = std::function<TangentVector(Scalar)>;

    /**
     * @brief Integration method selector.
     */
    enum class Method {
        Euler,    ///< 1st order — equivalent to expCosserat()
        Midpoint, ///< 2nd order (RKMK2)
        RKMK4     ///< 4th order (RK4 + Magnus commutator correction)
    };

    // ── Single-step integrators ──────────────────────────────────────────────

    /**
     * @brief Euler step (order 1) — equivalent to expCosserat().
     *
     * g(s+ds) = g(s) · Exp(ds · strain(s))
     *
     * Exact for piecewise-constant strain. Retained as a reference and
     * for backward-compatibility comparisons.
     *
     * @param g0     Initial pose at s0
     * @param strain Strain field s → ξ(s)
     * @param s0     Starting arc-length
     * @param ds     Step size (element length)
     */
    [[nodiscard]] static SE3Type integrateEuler(
        const SE3Type&   g0,
        const StrainField& strain,
        Scalar s0,
        Scalar ds)
    {
        const TangentVector k1 = strain(s0);
        return g0.compose(SE3Type::expCosserat(k1, ds));
    }

    /**
     * @brief Midpoint rule (order 2, RKMK2).
     *
     * σ = ds · strain(s0 + ds/2)
     * g(s+ds) = g(s) · Exp(σ)
     *
     * Two strain evaluations. Good balance of cost vs accuracy for
     * moderate curvatures.
     *
     * @param g0     Initial pose at s0
     * @param strain Strain field s → ξ(s)
     * @param s0     Starting arc-length
     * @param ds     Step size (element length)
     */
    [[nodiscard]] static SE3Type integrateMidpoint(
        const SE3Type&    g0,
        const StrainField& strain,
        Scalar s0,
        Scalar ds)
    {
        const TangentVector sigma = ds * strain(s0 + ds * Scalar(0.5));
        return g0.compose(SE3Type::computeExp(sigma));
    }

    /**
     * @brief 4th-order Runge-Kutta-Munthe-Kaas (RKMK4) step.
     *
     * Uses classical RK4 quadrature weights combined with the leading
     * Magnus commutator correction:
     *
     *   k1 = strain(s0)
     *   k2 = strain(s0 + ds/2)
     *   k3 = strain(s0 + ds/2)
     *   k4 = strain(s0 + ds)
     *
     *   σ = ds·(k1 + 2·k2 + 2·k3 + k4)/6 + (ds²/12)·[k1, k4]
     *
     *   g(s+ds) = g(s) · Exp(σ)
     *
     * The commutator term [k1, k4] accounts for the non-commutativity of
     * SE(3) and lifts the method from 4th order (commutative) to 4th order
     * (non-commutative / Magnus sense).
     *
     * For constant strain, k1=k2=k3=k4 and [k1,k4]=0, so RKMK4 = Euler.
     *
     * Reference: Munthe-Kaas (1998), Iserles et al. (2000) §IV.
     *
     * @param g0     Initial pose at s0
     * @param strain Strain field s → ξ(s)
     * @param s0     Starting arc-length
     * @param ds     Step size (element length)
     */
    [[nodiscard]] static SE3Type integrateRKMK4(
        const SE3Type&    g0,
        const StrainField& strain,
        Scalar s0,
        Scalar ds)
    {
        const TangentVector k1 = strain(s0);
        const TangentVector k2 = strain(s0 + ds * Scalar(0.5));
        const TangentVector k3 = strain(s0 + ds * Scalar(0.5));
        const TangentVector k4 = strain(s0 + ds);

        // RK4 weighted average in the algebra
        const TangentVector rk4_avg =
            ds * (k1 + Scalar(2)*k2 + Scalar(2)*k3 + k4) / Scalar(6);

        // Leading Magnus commutator correction (makes method 4th-order for
        // non-commutative algebras)
        const TangentVector comm = (ds * ds / Scalar(12)) * lieBracket(k1, k4);

        const TangentVector sigma = rk4_avg + comm;
        return g0.compose(SE3Type::computeExp(sigma));
    }

    // ── Multi-step integration ───────────────────────────────────────────────

    /**
     * @brief Integrate the Cosserat ODE over [s0, s0+L] with N uniform steps.
     *
     * Chains N single steps of the chosen method to produce the final pose.
     * Useful for:
     *   - Sub-stepping a single element for accuracy
     *   - Integrating along the full rod from base to tip
     *
     * @param g0      Initial pose at s0
     * @param strain  Strain field s → ξ(s)
     * @param s0      Starting arc-length
     * @param L       Total length to integrate over
     * @param N       Number of sub-steps (N=1 ≡ single element step)
     * @param method  Integration method (default: RKMK4)
     * @return Final pose g(s0 + L)
     */
    [[nodiscard]] static SE3Type integrate(
        const SE3Type&    g0,
        const StrainField& strain,
        Scalar s0,
        Scalar L,
        int    N      = 1,
        Method method = Method::RKMK4)
    {
        const Scalar ds = L / Scalar(N);
        SE3Type g = g0;

        for (int i = 0; i < N; ++i) {
            const Scalar s = s0 + Scalar(i) * ds;
            switch (method) {
                case Method::Euler:
                    g = integrateEuler(g, strain, s, ds);
                    break;
                case Method::Midpoint:
                    g = integrateMidpoint(g, strain, s, ds);
                    break;
                case Method::RKMK4:
                default:
                    g = integrateRKMK4(g, strain, s, ds);
                    break;
            }
        }
        return g;
    }

    /**
     * @brief Integrate and record all intermediate poses.
     *
     * Same as integrate() but stores the pose at each sub-step boundary.
     * Useful when you need poses at intermediate arc-lengths (e.g. when
     * populating output frame arrays from a LegendrePolynomialsMapping).
     *
     * @param g0      Initial pose at s0
     * @param strain  Strain field s → ξ(s)
     * @param s0      Starting arc-length
     * @param L       Total length to integrate over
     * @param N       Number of sub-steps
     * @param method  Integration method (default: RKMK4)
     * @return Vector of N+1 poses: [g(s0), g(s0+ds), …, g(s0+L)]
     */
    [[nodiscard]] static std::vector<SE3Type> integratePath(
        const SE3Type&    g0,
        const StrainField& strain,
        Scalar s0,
        Scalar L,
        int    N      = 1,
        Method method = Method::RKMK4)
    {
        std::vector<SE3Type> path;
        path.reserve(static_cast<std::size_t>(N + 1));
        path.push_back(g0);

        const Scalar ds = L / Scalar(N);
        SE3Type g = g0;

        for (int i = 0; i < N; ++i) {
            const Scalar s = s0 + Scalar(i) * ds;
            switch (method) {
                case Method::Euler:
                    g = integrateEuler(g, strain, s, ds);
                    break;
                case Method::Midpoint:
                    g = integrateMidpoint(g, strain, s, ds);
                    break;
                case Method::RKMK4:
                default:
                    g = integrateRKMK4(g, strain, s, ds);
                    break;
            }
            path.push_back(g);
        }
        return path;
    }

private:
    /**
     * @brief Lie bracket [xi, eta] in se(3).
     *
     * For ξ = [φ, ρ] and η = [ψ, σ] (Cosserat convention):
     *   [ξ, η] = [φ×ψ,  φ×σ - ψ×ρ]
     *
     * This is the adjoint representation: ad_ξ · η.
     */
    [[nodiscard]] static TangentVector lieBracket(
        const TangentVector& xi,
        const TangentVector& eta)
    {
        const Vector3 phi = xi.template head<3>();
        const Vector3 rho = xi.template tail<3>();
        const Vector3 psi = eta.template head<3>();
        const Vector3 sig = eta.template tail<3>();

        TangentVector result;
        result.template head<3>() = phi.cross(psi);
        result.template tail<3>() = phi.cross(sig) - psi.cross(rho);
        return result;
    }
};

// ── Convenience type aliases ──────────────────────────────────────────────────

using SE3Integratord = SE3Integrator<double>;
using SE3Integratorf = SE3Integrator<float>;

} // namespace sofa::component::cosserat::liegroups
