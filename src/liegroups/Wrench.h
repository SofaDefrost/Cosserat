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

#include <Eigen/Dense>
#include <iostream>
#include "Types.h"
#include "Twist.h"

namespace sofa::component::cosserat::liegroups {

/**
 * @brief Wrench — element of the dual Lie algebra se(3)*.
 *
 * Represents a generalised force (resultant force + moment) acting on a
 * rigid body cross-section. In Cosserat rod mechanics these are the
 * internal stress resultants transmitted through a cross-section.
 *
 * Storage convention (dual to Cosserat strain / Twist convention):
 *   w = [τ, f]ᵀ ∈ ℝ⁶
 *   τ = torque / moment (head<3>)  — dual to angular strain φ
 *   f = force           (tail<3>)  — dual to linear  strain ρ
 *
 * The duality pairing (virtual power) is:
 *   ⟨w, ξ⟩ = w^T ξ = τ·φ + f·ρ    (scalar)
 *
 * Frame-change rule (co-adjoint action):
 *   w_a = Ad_{g_ab}^T · w_b
 *
 * This preserves virtual power: ⟨w_b, ξ_b⟩ = ⟨w_a, ξ_a⟩
 * where ξ_b = Ad_{g_ab} · ξ_a.
 *
 * Use transformBy(g.adjoint()) which applies Ad^T internally.
 *
 * @tparam _Scalar Scalar type (float, double, autodiff dual, …)
 */
template<typename _Scalar>
class Wrench {
public:
    using Scalar  = _Scalar;
    using LGTypes = Types<Scalar>;

    using Vector3 = typename LGTypes::Vector3;
    using Vector6 = typename LGTypes::Vector6;
    using Matrix6 = typename LGTypes::Matrix6;

    // ── Constructors ────────────────────────────────────────────────────────

    Wrench() : m_data(Vector6::Zero()) {}

    explicit Wrench(const Vector6& data) : m_data(data) {}

    Wrench(const Vector3& torque, const Vector3& force) {
        m_data.template head<3>() = torque;
        m_data.template tail<3>() = force;
    }

    static Wrench Zero() { return Wrench(Vector6::Zero()); }

    // ── Accessors ───────────────────────────────────────────────────────────

    /** Torque / moment τ — first 3 components, dual to angular strain φ. */
    [[nodiscard]] Vector3 torque() const { return m_data.template head<3>(); }

    /** Force f — last 3 components, dual to linear strain ρ. */
    [[nodiscard]] Vector3 force()  const { return m_data.template tail<3>(); }

    [[nodiscard]] const Vector6& data() const { return m_data; }
    [[nodiscard]]       Vector6& data()       { return m_data; }

    [[nodiscard]] Scalar  operator[](int i) const { return m_data[i]; }
    [[nodiscard]] Scalar& operator[](int i)       { return m_data[i]; }

    // ── Duality pairing ──────────────────────────────────────────────────────

    /**
     * @brief Virtual power: ⟨w, ξ⟩ = τ·φ + f·ρ.
     *
     * This is the fundamental duality pairing between se(3)* and se(3).
     * In mechanics it equals the power delivered by the wrench w
     * to a body moving with twist ξ.
     */
    [[nodiscard]] Scalar dot(const Twist<Scalar>& twist) const {
        return m_data.dot(twist.data());
    }

    // ── Frame change ─────────────────────────────────────────────────────────

    /**
     * @brief Express this wrench in a new frame (co-adjoint transport).
     *
     * Given g = SE3 element mapping frame a to frame b:
     *   w_a = Ad_{g_ab}^T · w_b
     *
     * Equivalently, to push a wrench from frame a to frame b:
     *   w_b = Ad_{g_ab}^{-T} · w_a  =  (Ad_{g_ab}^{-1})^T · w_a
     *
     * This method applies Ad^T, so:
     *   wrench_in_a = wrench_in_b.transformBy(g_ab.adjoint())
     *
     * Virtual power is preserved:
     *   ⟨w_b, ξ_b⟩ = ⟨w_a, ξ_a⟩  when ξ_b = Ad * ξ_a
     *
     * @param Ad 6×6 adjoint matrix from SE3::adjoint() or SE3::computeAdjoint()
     */
    [[nodiscard]] Wrench transformBy(const Matrix6& Ad) const {
        return Wrench(Ad.transpose() * m_data);
    }

    /**
     * @brief Inverse co-adjoint transport: push wrench from a to b.
     *
     *   w_b = Ad_{g_ab}^{-T} · w_a
     *
     * Usage: wrench_in_b = wrench_in_a.transformByInverse(g_ab.adjoint())
     *
     * @param Ad 6×6 adjoint matrix from SE3::adjoint()
     */
    [[nodiscard]] Wrench transformByInverse(const Matrix6& Ad) const {
        return Wrench(Ad.transpose().lu().solve(m_data));
    }

    // ── Norms ────────────────────────────────────────────────────────────────

    [[nodiscard]] Scalar norm()        const { return m_data.norm(); }
    [[nodiscard]] Scalar torqueNorm()  const { return torque().norm(); }
    [[nodiscard]] Scalar forceNorm()   const { return force().norm(); }
    [[nodiscard]] Scalar squaredNorm() const { return m_data.squaredNorm(); }

    // ── Arithmetic ───────────────────────────────────────────────────────────

    [[nodiscard]] Wrench operator+(const Wrench& other) const { return Wrench(m_data + other.m_data); }
    [[nodiscard]] Wrench operator-(const Wrench& other) const { return Wrench(m_data - other.m_data); }
    [[nodiscard]] Wrench operator-()                    const { return Wrench(-m_data); }
    [[nodiscard]] Wrench operator*(Scalar s)            const { return Wrench(m_data * s); }
    [[nodiscard]] Wrench operator/(Scalar s)            const { return Wrench(m_data / s); }

    Wrench& operator+=(const Wrench& other) { m_data += other.m_data; return *this; }
    Wrench& operator-=(const Wrench& other) { m_data -= other.m_data; return *this; }
    Wrench& operator*=(Scalar s)            { m_data *= s;            return *this; }

    [[nodiscard]] friend Wrench operator*(Scalar s, const Wrench& w) { return w * s; }

    // ── Comparison ───────────────────────────────────────────────────────────

    [[nodiscard]] bool isApprox(const Wrench& other,
                                Scalar eps = LGTypes::epsilon()) const {
        return m_data.isApprox(other.m_data, eps);
    }

    [[nodiscard]] bool isZero(Scalar eps = LGTypes::epsilon()) const {
        return m_data.isZero(eps);
    }

    // ── I/O ──────────────────────────────────────────────────────────────────

    friend std::ostream& operator<<(std::ostream& os, const Wrench& w) {
        os << "Wrench[τ=(" << w.torque().transpose()
           << "), f=(" << w.force().transpose() << ")]";
        return os;
    }

private:
    Vector6 m_data;  ///< [τ(torque), f(force)]
};

// ── Convenience type aliases ──────────────────────────────────────────────────

using Wrenchd = Wrench<double>;
using Wrenchf = Wrench<float>;

// ── Free function helpers ─────────────────────────────────────────────────────

/**
 * @brief Virtual power delivered by a wrench to a body with a given twist.
 *
 * P = ⟨w, ξ⟩ = τ·φ + f·ρ
 */
template<typename Scalar>
[[nodiscard]] Scalar virtualPower(const Wrench<Scalar>& w, const Twist<Scalar>& xi) {
    return w.dot(xi);
}

} // namespace sofa::component::cosserat::liegroups
