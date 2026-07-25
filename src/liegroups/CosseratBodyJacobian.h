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
#include <vector>
#include <stdexcept>
#include "SE3.h"
#include "Twist.h"
#include "Wrench.h"

namespace sofa::component::cosserat::liegroups {

/**
 * @brief Body Jacobian of a discretised Cosserat rod.
 *
 * Centralises the kinematic computation that is currently scattered across
 * BaseCosseratMapping::applyJ / applyJT.
 *
 * ## Mathematical background
 *
 * For a rod with N elements, the body twist at section k follows the
 * incremental recurrence:
 *
 *   η_0 = base_twist
 *   η_k = Ad_{g_k⁻¹} · (η_{k-1} + J_local_k · ξ̇_k)
 *
 * where:
 *   - g_k  = expCosserat(ξ_k, L_k)  is the local SE(3) transform for element k
 *   - J_local_k ∈ ℝ^{6×6} is the tangent-exponential matrix for element k,
 *     i.e. what computeTangExpImplementation() returns
 *
 * Unrolling the recurrence, the full body twist at section k is:
 *
 *   η_k = Σ_{j=1}^{k} B_{j,k} · ξ̇_j  +  Ad_{G_k⁻¹} · η_0
 *
 * where:
 *   B_{j,k} = Ad_{G_k⁻¹ · G_{j-1}} · J_local_j     (G_k = g_1·…·g_k)
 *
 * The full Body Jacobian at section k is therefore:
 *
 *   J_body(k) = [B_{1,k} | B_{2,k} | ... | B_{k,k} | 0 | ... | 0]  ∈ ℝ^{6 × 6N}
 *
 * ## Usage
 *
 * @code
 *   CosseratBodyJacobian<double> J(N);
 *
 *   // 1. Feed section data (in order, root → tip)
 *   for (int k = 0; k < N; ++k)
 *       J.pushSection(g_k, J_local_k);
 *
 *   // 2. Forward pass: body twists from strain rates
 *   auto twists = J.applyForward(strain_rates, base_twist);
 *
 *   // 3. Backward pass: strain forces from external wrenches (applyJT)
 *   auto strain_forces = J.applyTranspose(external_wrenches, base_wrench_out);
 *
 *   // 4. Explicit matrix at section k (for control / IK)
 *   Eigen::MatrixXd Jk = J.getJacobianAtSection(k);
 * @endcode
 *
 * @tparam _Scalar Scalar type (float, double, …)
 */
template<typename _Scalar>
class CosseratBodyJacobian {
public:
    using Scalar      = _Scalar;
    using SE3Type     = SE3<Scalar>;
    using TwistType   = Twist<Scalar>;
    using WrenchType  = Wrench<Scalar>;

    using Matrix6  = typename Types<Scalar>::Matrix6;
    using Vector6  = typename Types<Scalar>::Vector6;
    using MatrixXd = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using VectorXd = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    // ── Construction ─────────────────────────────────────────────────────────

    explicit CosseratBodyJacobian(int expected_sections = 0) {
        if (expected_sections > 0) reserve(expected_sections);
    }

    void reserve(int n) {
        g_.reserve(static_cast<std::size_t>(n));
        J_local_.reserve(static_cast<std::size_t>(n));
        Ad_inv_.reserve(static_cast<std::size_t>(n));
    }

    /**
     * @brief Clear all stored section data.
     */
    void clear() {
        g_.clear();
        J_local_.clear();
        Ad_inv_.clear();
    }

    // ── Building the Jacobian ─────────────────────────────────────────────────

    /**
     * @brief Append one section to the rod, from root to tip.
     *
     * @param g_k        Local SE(3) transform for this element:
     *                   SE3::expCosserat(strain, length)
     * @param J_local_k  6×6 tangent-exponential matrix for this element:
     *                   result of computeTangExpImplementation(length, strain)
     */
    void pushSection(const SE3Type& g_k, const Matrix6& J_local_k) {
        Ad_inv_.push_back(g_k.computeAdjoint().inverse());
        g_.push_back(g_k);
        J_local_.push_back(J_local_k);
    }

    /** Number of sections pushed so far. */
    [[nodiscard]] int size() const { return static_cast<int>(g_.size()); }

    // ── Forward pass: applyJ ──────────────────────────────────────────────────

    /**
     * @brief Propagate body twists from strain rates (applyJ).
     *
     * Implements:
     *   η_0 = base_twist
     *   η_k = Ad_{g_k⁻¹} · (η_{k-1} + J_local_k · ξ̇_k)
     *
     * @param strain_rates  N body strain rates ξ̇_k ∈ se(3), one per section.
     * @param base_twist    Body twist at the root (η_0).
     * @return              N+1 body twists [η_0, η_1, …, η_N].
     */
    [[nodiscard]] std::vector<TwistType> applyForward(
        const std::vector<TwistType>& strain_rates,
        const TwistType& base_twist = TwistType::Zero()) const
    {
        const int N = size();
        if (static_cast<int>(strain_rates.size()) != N)
            throw std::invalid_argument(
                "CosseratBodyJacobian::applyForward: "
                "strain_rates.size() must equal number of sections");

        std::vector<TwistType> twists;
        twists.reserve(static_cast<std::size_t>(N + 1));
        twists.push_back(base_twist);

        for (int k = 0; k < N; ++k) {
            // η_k = Ad_{g_k⁻¹} · (η_{k-1} + J_local_k · ξ̇_k)
            const Vector6 eta_prev = twists.back().data();
            const Vector6 contribution = J_local_[k] * strain_rates[k].data();
            const Vector6 eta_k = Ad_inv_[k] * (eta_prev + contribution);
            twists.emplace_back(eta_k);
        }
        return twists;
    }

    // ── Backward pass: applyJT ────────────────────────────────────────────────

    /**
     * @brief Propagate wrenches back to strain forces (applyJT).
     *
     * The transpose of applyForward. Implements virtual-work duality:
     *   ⟨f_out, δq⟩ = ⟨w_ext, η⟩   for all δq
     *
     * Backward recurrence (dual of applyForward):
     *
     *   f_N = w_N                                      (init at tip)
     *   q_k = J_local_k^T · Ad_{g_k⁻¹}ᵀ · f_{k+1}   (strain force)
     *   f_k = w_k + Ad_{g_k⁻¹}ᵀ · f_{k+1}            (co-twist transport)
     *
     * where sections are 0-indexed (section k maps position k → k+1).
     * The accumulated co-twist f_k carries only DOWNSTREAM wrenches
     * (positions k+1 … N); w_k itself never enters q_k.
     *
     * @param external_wrenches  N+1 external wrenches at sections [w_0, …, w_N].
     * @param base_wrench_out    Accumulated wrench at the root (output).
     * @return                   N strain forces (generalised forces in strain space).
     */
    [[nodiscard]] std::vector<WrenchType> applyTranspose(
        const std::vector<WrenchType>& external_wrenches,
        WrenchType& base_wrench_out) const
    {
        const int N = size();
        if (static_cast<int>(external_wrenches.size()) != N + 1)
            throw std::invalid_argument(
                "CosseratBodyJacobian::applyTranspose: "
                "external_wrenches.size() must equal N+1");

        std::vector<WrenchType> strain_forces(static_cast<std::size_t>(N), WrenchType::Zero());

        // Backward accumulation — dual of the forward pass
        //
        // Forward:  η_k = Ad_{g_k⁻¹} · (η_{k-1} + J_k · ξ̇_k)
        //
        // Dual (virtual power):
        //   f_{k-1} = w_{k-1} + Ad_{g_k⁻¹}ᵀ · f_k        (co-twist transport)
        //   q_k     = J_kᵀ · Ad_{g_k⁻¹}ᵀ · f_k            (strain force)
        //
        // Key: q_k depends only on f_k (wrenches DOWNSTREAM of section k).
        // w_{k-1} is upstream of section k and must NOT appear in q_k.
        //
        // Here sections are 0-indexed: section k maps position k → position k+1.
        //   Ad_inv_[k] = Ad_{g_k⁻¹}
        //   f_{k+1} ≡ acc  before the update step

        Vector6 acc = external_wrenches[N].data();  // f_N = w_N  (tip)

        for (int k = N - 1; k >= 0; --k) {
            // Transport f_{k+1} back through section k:
            //   transported = Ad_{g_k⁻¹}ᵀ · f_{k+1}
            const Vector6 transported = Ad_inv_[k].transpose() * acc;

            // Strain force: q_k = J_kᵀ · Ad_{g_k⁻¹}ᵀ · f_{k+1}
            strain_forces[static_cast<std::size_t>(k)] =
                WrenchType(J_local_[k].transpose() * transported);

            // Update accumulator: f_k = w_k + Ad_{g_k⁻¹}ᵀ · f_{k+1}
            acc = external_wrenches[k].data() + transported;
        }

        // Base wrench: f_0 = w_0 + Ad_{g_0⁻¹}ᵀ · f_1
        base_wrench_out = WrenchType(acc);
        return strain_forces;
    }

    // ── Explicit matrix access ────────────────────────────────────────────────

    /**
     * @brief Assemble the explicit 6×6N Body Jacobian at section k.
     *
     * Column block j (for j ≤ k):
     *   B_{j,k} = Ad_{G_k⁻¹ · G_{j-1}} · J_local_j
     *           = Ad_{g_k⁻¹} · Ad_{g_{k-1}⁻¹} · … · Ad_{g_{j+1}⁻¹} · J_local_j
     *
     * Computed via backward accumulation from k to j.
     *
     * @param k  Section index (0-based, 0 = first element after root)
     * @return   6×6N matrix; columns j*6…j*6+5 are the contribution of section j.
     */
    [[nodiscard]] MatrixXd getJacobianAtSection(int k) const {
        const int N = size();
        if (k < 0 || k >= N)
            throw std::out_of_range("CosseratBodyJacobian::getJacobianAtSection: index out of range");

        MatrixXd J = MatrixXd::Zero(6, 6 * N);

        // Start from section k and transport backward toward the root
        // cumulative_Ad = Ad_{g_k⁻¹} · Ad_{g_{k-1}⁻¹} · … · Ad_{g_{j+1}⁻¹}
        Matrix6 cumulative_Ad = Matrix6::Identity();

        for (int j = k; j >= 0; --j) {
            // B_{j,k} = cumulative_Ad · J_local_j
            J.template block<6, 6>(0, j * 6) = cumulative_Ad * J_local_[j];
            // Extend transport one more step toward the root
            if (j > 0)
                cumulative_Ad = cumulative_Ad * Ad_inv_[j];
        }
        return J;
    }

    /**
     * @brief Assemble the full stacked Body Jacobian (6N × 6N).
     *
     * Row block k (rows k*6…k*6+5) = getJacobianAtSection(k).
     * The matrix is lower block-triangular.
     *
     * Useful for:
     *   - Computing the pseudoinverse for task-space IK
     *   - Analysing rod manipulability
     *   - Linearised dynamics matrices
     *
     * Complexity: O(N³) in the worst case (dominated by the accumulations).
     * For real-time control, prefer applyForward / applyTranspose.
     */
    [[nodiscard]] MatrixXd getFullJacobian() const {
        const int N = size();
        MatrixXd J_full = MatrixXd::Zero(6 * N, 6 * N);
        for (int k = 0; k < N; ++k)
            J_full.block(k * 6, 0, 6, 6 * N) = getJacobianAtSection(k);
        return J_full;
    }

    // ── Manipulability ────────────────────────────────────────────────────────

    /**
     * @brief Yoshikawa manipulability measure at section k.
     *
     *   w(k) = sqrt(det(J_k · J_k^T))
     *
     * Measures how isotropically the rod tip can move in task space
     * from section k. Zero indicates a singular configuration.
     *
     * @param k  Section index
     */
    [[nodiscard]] Scalar manipulability(int k) const {
        const MatrixXd Jk = getJacobianAtSection(k);
        const MatrixXd JJT = Jk * Jk.transpose();
        return std::sqrt(std::max(Scalar(0), JJT.determinant()));
    }

    /**
     * @brief Minimum singular value of J_k (condition number proxy).
     *
     * Small minimum singular value → near-singular configuration.
     *
     * @param k  Section index
     */
    [[nodiscard]] Scalar minSingularValue(int k) const {
        const MatrixXd Jk = getJacobianAtSection(k);
        Eigen::JacobiSVD<MatrixXd> svd(Jk);
        return svd.singularValues().minCoeff();
    }

    // ── Accessors ─────────────────────────────────────────────────────────────

    [[nodiscard]] const SE3Type&  getSE3(int k)      const { return g_[k]; }
    [[nodiscard]] const Matrix6&  getLocalJacobian(int k) const { return J_local_[k]; }
    [[nodiscard]] const Matrix6&  getAdjointInverse(int k) const { return Ad_inv_[k]; }

private:
    std::vector<SE3Type>  g_;        ///< Local SE3 transforms g_k = expCosserat(ξ_k, L_k)
    std::vector<Matrix6>  J_local_;  ///< Local tangent-exp matrices J_local_k
    std::vector<Matrix6>  Ad_inv_;   ///< Precomputed Ad_{g_k⁻¹} for each section
};

// ── Convenience type aliases ──────────────────────────────────────────────────

using CosseratBodyJacobiand = CosseratBodyJacobian<double>;
using CosseratBodyJacobianf = CosseratBodyJacobian<float>;

} // namespace sofa::component::cosserat::liegroups
