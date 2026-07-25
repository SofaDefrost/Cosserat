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

namespace sofa::component::cosserat::liegroups {

/**
 * @brief Body twist — element of the Lie algebra se(3).
 *
 * Represents the instantaneous velocity (or strain rate) of a rigid body,
 * expressed in the body frame. In Cosserat rod mechanics this is the
 * strain vector along the rod centreline.
 *
 * Storage convention (matches Cosserat strain convention):
 *   ξ = [φ, ρ]ᵀ ∈ ℝ⁶
 *   φ = angular part (head<3>) : torsion / bending strains
 *   ρ = linear part  (tail<3>) : elongation / shearing strains
 *
 * Frame-change rule (adjoint action):
 *   ξ_b = Ad_{g_ab} · ξ_a      (g transforms frame a → frame b)
 *
 * Use transformBy(g.adjoint()) to change the expression frame.
 *
 * @tparam _Scalar Scalar type (float, double, autodiff dual, …)
 */
template<typename _Scalar>
class Twist {
public:
    using Scalar  = _Scalar;
    using LGTypes = Types<Scalar>;

    using Vector3    = typename LGTypes::Vector3;
    using Vector6    = typename LGTypes::Vector6;
    using Matrix3    = typename LGTypes::Matrix3;
    using Matrix4    = typename LGTypes::template Matrix<4, 4>;
    using Matrix6    = typename LGTypes::Matrix6;

    // ── Constructors ────────────────────────────────────────────────────────

    Twist() : m_data(Vector6::Zero()) {}

    explicit Twist(const Vector6& data) : m_data(data) {}

    Twist(const Vector3& angular, const Vector3& linear) {
        m_data.template head<3>() = angular;
        m_data.template tail<3>() = linear;
    }

    static Twist Zero()   { return Twist(Vector6::Zero()); }

    // ── Accessors ───────────────────────────────────────────────────────────

    /** Angular part φ (torsion / bending) — first 3 components. */
    [[nodiscard]] Vector3 angular() const { return m_data.template head<3>(); }

    /** Linear part ρ (elongation / shearing) — last 3 components. */
    [[nodiscard]] Vector3 linear()  const { return m_data.template tail<3>(); }

    [[nodiscard]] const Vector6& data() const { return m_data; }
    [[nodiscard]]       Vector6& data()       { return m_data; }

    [[nodiscard]] Scalar operator[](int i) const { return m_data[i]; }
    [[nodiscard]] Scalar& operator[](int i)      { return m_data[i]; }

    // ── Lie algebra structure ────────────────────────────────────────────────

    /**
     * @brief Hat operator: vector → 4×4 se(3) matrix.
     *
     * For ξ = [φ, ρ]:
     *   ξ̂ = [ φ×  ρ ]  (top-left 3×3 is skew(φ), top-right col is ρ)
     *        [  0  0 ]
     */
    [[nodiscard]] Matrix4 hat() const {
        Matrix4 Xi = Matrix4::Zero();
        const Vector3 phi = angular();
        const Vector3 rho = linear();
        // skew-symmetric block for angular part
        Xi(0,1) = -phi[2];  Xi(0,2) =  phi[1];
        Xi(1,0) =  phi[2];  Xi(1,2) = -phi[0];
        Xi(2,0) = -phi[1];  Xi(2,1) =  phi[0];
        // linear part in last column
        Xi(0,3) = rho[0];
        Xi(1,3) = rho[1];
        Xi(2,3) = rho[2];
        return Xi;
    }

    /**
     * @brief Lie bracket [this, other] = ad_ξ · η.
     *
     * For ξ = [φ, ρ] and η = [ψ, σ]:
     *   [ξ, η] = [φ×ψ,  φ×σ - ψ×ρ]
     */
    [[nodiscard]] Twist bracket(const Twist& other) const {
        const Vector3 phi = angular();
        const Vector3 rho = linear();
        const Vector3 psi = other.angular();
        const Vector3 sig = other.linear();
        return Twist(phi.cross(psi),
                     phi.cross(sig) - psi.cross(rho));
    }

    /**
     * @brief Small adjoint matrix ad_ξ such that [ξ, η] = ad_ξ · η.
     *
     * ad_ξ = [ φ×   0  ]    (6×6, acts on twist vectors)
     *         [ ρ×   φ× ]
     */
    [[nodiscard]] Matrix6 smallAdjoint() const {
        Matrix6 ad = Matrix6::Zero();
        const Vector3 phi = angular();
        const Vector3 rho = linear();

        // skew(phi) blocks on diagonal
        const Matrix3 phi_x = skew(phi);
        const Matrix3 rho_x = skew(rho);

        ad.template block<3,3>(0,0) = phi_x;
        ad.template block<3,3>(3,3) = phi_x;
        ad.template block<3,3>(3,0) = rho_x;
        return ad;
    }

    // ── Frame change ─────────────────────────────────────────────────────────

    /**
     * @brief Express this twist in a new frame using the adjoint matrix.
     *
     * Given g = SE3 element mapping frame a to frame b:
     *   ξ_b = Ad_g · ξ_a
     *
     * Usage:
     *   Twist xi_b = xi_a.transformBy(g_ab.adjoint());
     *
     * @param Ad 6×6 adjoint matrix from SE3::adjoint() or SE3::computeAdjoint()
     */
    [[nodiscard]] Twist transformBy(const Matrix6& Ad) const {
        return Twist(Ad * m_data);
    }

    // ── Norms ────────────────────────────────────────────────────────────────

    [[nodiscard]] Scalar norm()        const { return m_data.norm(); }
    [[nodiscard]] Scalar angularNorm() const { return angular().norm(); }
    [[nodiscard]] Scalar linearNorm()  const { return linear().norm(); }
    [[nodiscard]] Scalar squaredNorm() const { return m_data.squaredNorm(); }

    // ── Arithmetic ───────────────────────────────────────────────────────────

    [[nodiscard]] Twist operator+(const Twist& other) const { return Twist(m_data + other.m_data); }
    [[nodiscard]] Twist operator-(const Twist& other) const { return Twist(m_data - other.m_data); }
    [[nodiscard]] Twist operator-()                   const { return Twist(-m_data); }
    [[nodiscard]] Twist operator*(Scalar s)           const { return Twist(m_data * s); }
    [[nodiscard]] Twist operator/(Scalar s)           const { return Twist(m_data / s); }

    Twist& operator+=(const Twist& other) { m_data += other.m_data; return *this; }
    Twist& operator-=(const Twist& other) { m_data -= other.m_data; return *this; }
    Twist& operator*=(Scalar s)           { m_data *= s;            return *this; }

    [[nodiscard]] friend Twist operator*(Scalar s, const Twist& t) { return t * s; }

    // ── Comparison ───────────────────────────────────────────────────────────

    [[nodiscard]] bool isApprox(const Twist& other,
                                Scalar eps = LGTypes::epsilon()) const {
        return m_data.isApprox(other.m_data, eps);
    }

    [[nodiscard]] bool isZero(Scalar eps = LGTypes::epsilon()) const {
        return m_data.isZero(eps);
    }

    // ── I/O ──────────────────────────────────────────────────────────────────

    friend std::ostream& operator<<(std::ostream& os, const Twist& t) {
        os << "Twist[φ=(" << t.angular().transpose()
           << "), ρ=(" << t.linear().transpose() << ")]";
        return os;
    }

private:
    Vector6 m_data;  ///< [φ(angular), ρ(linear)]

    /** Build 3×3 skew-symmetric matrix from a vector. */
    static Matrix3 skew(const Vector3& v) {
        Matrix3 S;
        S <<       Scalar(0), -v[2],  v[1],
              v[2],       Scalar(0), -v[0],
             -v[1],  v[0],       Scalar(0);
        return S;
    }
};

// ── Convenience type aliases ──────────────────────────────────────────────────

using Twistd = Twist<double>;
using Twistf = Twist<float>;

} // namespace sofa::component::cosserat::liegroups
