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
#include <cmath>
#include "SE3.h"
#include "GaussianOnManifold.h"

namespace sofa::component::cosserat::liegroups {

/**
 * @brief Compositional uncertainty propagation along a Cosserat rod.
 *
 * Extends GaussianOnManifold with the physically correct step that
 * propagates pose uncertainty through one rod element, accounting for
 * both:
 *   1. Transport of existing pose uncertainty  (via Ad_{g⁻¹})
 *   2. Injection of new strain noise           (via the tangent-exp map)
 *
 * ## Mathematical derivation
 *
 * Consider the right-perturbation model:
 *   G_k = μ_k · exp(ε_k)   where ε_k ~ N(0, Σ_k)  (pose uncertainty)
 *
 * After composing with a noisy element g_k = exp(ξ_k + δξ_k),
 * where δξ_k ~ N(0, Σ_ξk) (strain uncertainty):
 *
 *   G_{k+1} = G_k · g_k
 *
 * First-order linearisation gives:
 *
 *   ε_{k+1} ≈ Ad_{g_k⁻¹} · ε_k  +  T(L, ξ_k) · δξ_k
 *
 * where T(L, ξ) = ∫₀ᴸ exp(s · ad_ξ) ds is the tangent-exponential matrix
 * (computed by computeTangExpImplementation in CosseratGeometryMapping).
 *
 * The propagated covariance is therefore:
 *
 *   Σ_{k+1} = Ad_{g_k⁻¹} · Σ_k · Ad_{g_k⁻¹}ᵀ  +  T · Σ_ξk · Tᵀ
 *              ↑ transport of existing uncertainty   ↑ injected strain noise
 *
 * ## Usage: filtering along the rod
 *
 * @code
 *   using Propagator = CosseratUncertaintyPropagator<double>;
 *
 *   GaussianOnManifold<SE3d> state(base_pose, Sigma_0);
 *
 *   for (int k = 0; k < N; ++k) {
 *       Propagator::Section sec;
 *       sec.strain        = strain_k;
 *       sec.length        = L_k;
 *       sec.strain_cov    = Sigma_xi_k;  // 6×6 strain noise covariance
 *       state = Propagator::propagateStep(state, sec);
 *       rod_gaussians.push_back(state);
 *   }
 * @endcode
 *
 * @tparam _Scalar Scalar type (float, double, …)
 */
template<typename _Scalar>
class CosseratUncertaintyPropagator {
public:
    using Scalar    = _Scalar;
    using SE3Type   = SE3<Scalar>;
    using GaussianSE3 = GaussianOnManifold<SE3Type>;

    using Matrix6  = typename Types<Scalar>::Matrix6;
    using Vector6  = typename Types<Scalar>::Vector6;
    using TangentVector = typename SE3Type::TangentVector;
    using Vector3  = typename SE3Type::Vector3;

    /**
     * @brief Data for one rod section (one integration step).
     */
    struct Section {
        TangentVector strain    = TangentVector::Zero(); ///< Nominal strain ξ_k ∈ se(3)
        Scalar        length    = Scalar(0);             ///< Arc-length L_k of this section
        Matrix6       strain_cov = Matrix6::Zero();      ///< Strain noise covariance Σ_ξk
    };

    // ── Single-step propagation ───────────────────────────────────────────────

    /**
     * @brief Propagate pose uncertainty through one rod section.
     *
     * Σ_{k+1} = Ad_{g⁻¹} · Σ_k · Ad_{g⁻¹}ᵀ  +  T · Σ_ξ · Tᵀ
     *
     * where:
     *   g = SE3::expCosserat(section.strain, section.length)
     *   T = tangentExp(section.length, section.strain)  [6×6]
     *
     * @param prior    Gaussian pose distribution at the start of the section.
     * @param section  Section geometry and strain noise.
     * @return         Gaussian pose distribution at the end of the section.
     */
    [[nodiscard]] static GaussianSE3 propagateStep(
        const GaussianSE3& prior,
        const Section&     section)
    {
        // ── Mean propagation ─────────────────────────────────────────────────
        const SE3Type g = SE3Type::expCosserat(section.strain, section.length);
        const SE3Type new_mean = prior.getMean().compose(g);

        // ── Covariance propagation ───────────────────────────────────────────
        // 1. Transport of existing pose uncertainty
        const Matrix6 Ad_g_inv = g.computeAdjoint().inverse();
        const typename GaussianSE3::CovarianceMatrix& Sigma_k = prior.getCovariance();
        const typename GaussianSE3::CovarianceMatrix transport =
            Ad_g_inv * Sigma_k * Ad_g_inv.transpose();

        // 2. Injected strain noise via the tangent-exponential map T(L, ξ)
        const Matrix6 T = computeTangentExp(section.length, section.strain);
        const typename GaussianSE3::CovarianceMatrix process_noise =
            T * section.strain_cov * T.transpose();

        const typename GaussianSE3::CovarianceMatrix new_cov = transport + process_noise;

        return GaussianSE3(new_mean, new_cov);
    }

    // ── Full-rod propagation ──────────────────────────────────────────────────

    /**
     * @brief Propagate uncertainty along the full rod.
     *
     * Returns the Gaussian distribution at each section endpoint,
     * starting from a base distribution.
     *
     * @param base      Gaussian at the rod base (root pose + uncertainty).
     * @param sections  Ordered list of sections (root → tip).
     * @return          N+1 Gaussian distributions: [base, g1, g2, …, gN].
     */
    [[nodiscard]] static std::vector<GaussianSE3> propagateAlongRod(
        const GaussianSE3&         base,
        const std::vector<Section>& sections)
    {
        std::vector<GaussianSE3> result;
        result.reserve(sections.size() + 1);
        result.push_back(base);

        for (const auto& sec : sections)
            result.push_back(propagateStep(result.back(), sec));

        return result;
    }

    // ── Marginalisation utilities ─────────────────────────────────────────────

    /**
     * @brief Extract the translational (position) uncertainty at a section.
     *
     * Returns the 3×3 covariance block corresponding to the linear part
     * of the perturbation (last 3 components of the 6D tangent, following
     * the Cosserat convention [φ, ρ]).
     *
     * @param g  Gaussian SE(3) distribution at a section.
     */
    [[nodiscard]] static Eigen::Matrix<Scalar, 3, 3>
    translationCovariance(const GaussianSE3& g) {
        return g.getCovariance().template block<3, 3>(3, 3);
    }

    /**
     * @brief Extract the rotational (orientation) uncertainty at a section.
     *
     * Returns the 3×3 covariance block for the angular part
     * (first 3 components of the 6D tangent: [φ_x, φ_y, φ_z]).
     *
     * @param g  Gaussian SE(3) distribution at a section.
     */
    [[nodiscard]] static Eigen::Matrix<Scalar, 3, 3>
    rotationCovariance(const GaussianSE3& g) {
        return g.getCovariance().template block<3, 3>(0, 0);
    }

    /**
     * @brief 3σ confidence ellipsoid radii for the translational part.
     *
     * Returns the square roots of the 3 eigenvalues of the 3×3 translation
     * covariance block, scaled by 3 (i.e. the 3σ semi-axes of the ellipsoid).
     *
     * Useful for visualising the rod tip uncertainty in SOFA.
     *
     * @param g  Gaussian SE(3) distribution at a section.
     */
    [[nodiscard]] static Eigen::Matrix<Scalar, 3, 1>
    tipConfidenceRadii(const GaussianSE3& g) {
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<Scalar, 3, 3>> eig(
            translationCovariance(g));
        return Scalar(3) * eig.eigenvalues().cwiseSqrt();
    }

private:
    /**
     * @brief Tangent-exponential map T(L, ξ): se(3) → se(3).
     *
     * Computes T = ∫₀ᴸ exp(s · ad_ξ) ds via the closed-form series:
     *
     * For θ = ‖φ‖ (angular part of ξ):
     *
     *   Small θ: T = L·I + (L²/2)·ad_ξ
     *
     *   General: T = L·I + α₁·ad_ξ + α₂·ad_ξ² + α₃·ad_ξ³ + α₄·ad_ξ⁴
     *
     * These are the same scalar coefficients as computeTangExpImplementation()
     * in CosseratGeometryMapping, centralised here in the liegroups layer.
     *
     * @param L      Arc-length of the section.
     * @param strain Strain vector ξ = [φ, ρ] (angular, linear).
     */
    [[nodiscard]] static Matrix6 computeTangentExp(Scalar L, const TangentVector& strain) {
        const Vector3 phi = strain.template head<3>();
        const Scalar theta = phi.norm();

        // Build ad_ξ (small adjoint matrix of the strain)
        const Matrix6 ad_xi = computeSmallAdjoint(strain);
        const Matrix6 I6    = Matrix6::Identity();

        if (theta <= std::numeric_limits<Scalar>::epsilon() * Scalar(100)) {
            // Small-angle approximation: T ≈ L·I + (L²/2)·ad_ξ
            return L * I6 + (L * L / Scalar(2)) * ad_xi;
        }

        const Scalar Lt     = L * theta;
        const Scalar cos_Lt = std::cos(Lt);
        const Scalar sin_Lt = std::sin(Lt);
        const Scalar t2     = theta * theta;
        const Scalar t3     = t2 * theta;
        const Scalar t4     = t2 * t2;
        const Scalar t5     = t4 * theta;

        const Scalar a1 = (Scalar(4) - Scalar(4)*cos_Lt - Lt*sin_Lt) / (Scalar(2)*t2);
        const Scalar a2 = (Scalar(4)*Lt*theta + Lt*theta*cos_Lt - Scalar(5)*sin_Lt) / (Scalar(2)*t3);
        const Scalar a3 = (Scalar(2) - Scalar(2)*cos_Lt - Lt*sin_Lt) / (Scalar(2)*t4);
        const Scalar a4 = (Scalar(2)*Lt*theta + Lt*theta*cos_Lt - Scalar(3)*sin_Lt) / (Scalar(2)*t5);

        const Matrix6 ad2 = ad_xi * ad_xi;
        const Matrix6 ad3 = ad2   * ad_xi;
        const Matrix6 ad4 = ad3   * ad_xi;

        return L * I6 + a1 * ad_xi + a2 * ad2 + a3 * ad3 + a4 * ad4;
    }

    /**
     * @brief Small adjoint matrix ad_ξ for ξ = [φ, ρ] ∈ se(3).
     *
     * ad_ξ = [ φ×   0  ]  (6×6)
     *         [ ρ×   φ× ]
     */
    [[nodiscard]] static Matrix6 computeSmallAdjoint(const TangentVector& xi) {
        Matrix6 ad = Matrix6::Zero();
        const Vector3 phi = xi.template head<3>();
        const Vector3 rho = xi.template tail<3>();

        const Eigen::Matrix<Scalar, 3, 3> phi_x = skew(phi);
        const Eigen::Matrix<Scalar, 3, 3> rho_x = skew(rho);

        ad.template block<3,3>(0,0) = phi_x;
        ad.template block<3,3>(3,3) = phi_x;
        ad.template block<3,3>(3,0) = rho_x;
        return ad;
    }

    static Eigen::Matrix<Scalar,3,3> skew(const Vector3& v) {
        Eigen::Matrix<Scalar,3,3> S;
        S <<  Scalar(0), -v[2],  v[1],
               v[2],  Scalar(0), -v[0],
              -v[1],   v[0],  Scalar(0);
        return S;
    }
};

// ── Convenience type aliases ──────────────────────────────────────────────────

using CosseratUncertaintyPropagator_d = CosseratUncertaintyPropagator<double>;
using CosseratUncertaintyPropagator_f = CosseratUncertaintyPropagator<float>;

} // namespace sofa::component::cosserat::liegroups
