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

#include <Cosserat/config.h>
#include <Cosserat/mapping/Strain2FramesCosseratMapping.h>

#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/objectmodel/Event.h>
#include <sofa/core/objectmodel/Link.h>
#include <sofa/simulation/AnimateBeginEvent.h>

namespace Cosserat::controller {

	using namespace sofa::component::cosserat::liegroups;

	/**
	 * @brief Quasi-static iLQR tip-tracking controller for Cosserat rods.
	 *
	 * At each simulation step, computes optimal strain corrections Δξ that
	 * steer the rod tip toward a desired target pose, minimising:
	 *
	 *   J(ξ) = ½ ‖log(g_ref⁻¹ · g_tip)‖²_Q  +  ½ ‖ξ‖²_R
	 *
	 * where:
	 *   - g_tip ∈ SE(3)  is the current tip pose (from output frames)
	 *   - g_ref ∈ SE(3)  is the desired target pose (d_targetPose)
	 *   - Q  (d_Q_tip)   is the tip error weight (scalar × I₆)
	 *   - R  (d_R_strain) is the strain regularisation weight (scalar × I_{6N})
	 *
	 * ## Control modes
	 *
	 * | d_mode | Name           | Update rule                                   |
	 * |--------|----------------|-----------------------------------------------|
	 * | 0      | Gradient       | Δξ = -α · J_tip^T · Q · e                    |
	 * | 1      | Gauss-Newton   | (J^T·Q·J + R·I)·Δξ = -J^T·Q·e  (one Newton) |
	 *
	 * Both modes use the body Jacobian J_tip ∈ ℝ^{6×(N·DOF)} provided by
	 * `Strain2FramesCosseratMapping::getBodyJacobian().getJacobianAtSection(N-1)`.
	 *
	 * For piecewise-constant strain (the standard Cosserat discretisation),
	 * Gauss-Newton converges in one step for small errors.  Multiple iLQR
	 * iterations (d_maxIterations > 1) amortise linearisation errors for
	 * large deformations.
	 *
	 * ## SOFA integration
	 *
	 * The controller listens to AnimateBeginEvent and writes corrected strain
	 * values back to the strain mechanical state at the start of each step,
	 * so that the next apply() sees the updated configuration.
	 *
	 * @tparam TIn1 Strain DOF type (typically Vec3Types or Vec6Types)
	 * @tparam TIn2 Rigid base type (typically Rigid3Types)
	 * @tparam TOut Output frame type (typically Rigid3Types)
	 */
	template<class TIn1, class TIn2, class TOut>
	class CosseratILQRController : public sofa::core::objectmodel::BaseObject {
	public:
		SOFA_CLASS(SOFA_TEMPLATE3(CosseratILQRController, TIn1, TIn2, TOut),
				   sofa::core::objectmodel::BaseObject);

		using MappingType   = Cosserat::mapping::Strain2FramesCosseratMapping<TIn1, TIn2, TOut>;
		using SE3Types      = SE3<double>;
		using SO3Type       = typename SE3Types::SO3Type;
		using TangentVector = typename SE3Types::TangentVector;
		using BodyJacobian  = CosseratBodyJacobian<double>;
		using MatrixXd      = Eigen::MatrixXd;
		using VectorXd      = Eigen::VectorXd;

		using Coord1   = sofa::Coord_t<TIn1>;
		using OutCoord = sofa::Coord_t<TOut>;

		static constexpr int NStrainDOF = static_cast<int>(Coord1::total_size);

	public:
		//////////////////////////////////////////////////////////////////////
		/// @name Links
		/// @{

		/// The Cosserat mapping whose body Jacobian drives the controller.
		sofa::core::objectmodel::SingleLink<
			CosseratILQRController,
			MappingType,
			sofa::core::objectmodel::BaseLink::FLAG_STOREPATH |
			sofa::core::objectmodel::BaseLink::FLAG_STRONGLINK>
		l_mapping;

		/// @}
		//////////////////////////////////////////////////////////////////////

		//////////////////////////////////////////////////////////////////////
		/// @name Data fields
		/// @{

		/// Desired tip pose in world frame (SE(3) = position + quaternion).
		sofa::Data<OutCoord> d_targetPose;

		/// Tip error weight Q (scalar, applied isotropically on se(3)).
		sofa::Data<double> d_Q_tip;

		/// Strain regularisation weight R (scalar, Tikhonov on Δξ).
		sofa::Data<double> d_R_strain;

		/// Gradient step size α (only used in mode 0 — gradient descent).
		sofa::Data<double> d_stepSize;

		/// Maximum iLQR iterations per simulation step.
		sofa::Data<int> d_maxIterations;

		/// Convergence threshold on ‖Δξ‖ (stops early if update is small).
		sofa::Data<double> d_tolerance;

		/// Control mode:
		///   0 = gradient descent (Δξ = -α · J^T · Q · e)
		///   1 = Gauss-Newton     (solve (J^T·Q·J + R·I)·Δξ = -J^T·Q·e)
		sofa::Data<int> d_mode;

		/// Index of the output frame to track as the tip (default: last frame).
		/// -1 means "last frame".
		sofa::Data<int> d_tipIndex;

		/// Whether to apply the strain update every step (true = closed-loop).
		sofa::Data<bool> d_active;

		// ── Diagnostics (read-only) ──

		/// Current tip error ‖log(g_ref⁻¹ · g_tip)‖ (updated each step).
		sofa::Data<double> d_tipError;

		/// Current manipulability at the tip section.
		sofa::Data<double> d_manipulability;

		/// @}
		//////////////////////////////////////////////////////////////////////

	public:
		//////////////////////////////////////////////////////////////////////
		/// @name Inherited from BaseObject
		/// @{

		void init() override;
		void handleEvent(sofa::core::objectmodel::Event *event) override;

		/// @}
		//////////////////////////////////////////////////////////////////////

		//////////////////////////////////////////////////////////////////////
		/// @name Controller API
		/// @{

		/**
		 * @brief Run one iLQR optimisation sweep (up to d_maxIterations steps).
		 *
		 * @return Per-section strain corrections Δξ_k (size = N sections).
		 *         Each entry is a NStrainDOF-dimensional vector (first NStrainDOF
		 *         components of the 6D body tangent).
		 */
		[[nodiscard]] std::vector<Coord1> computeControl();

		/// @}
		//////////////////////////////////////////////////////////////////////

	protected:
		CosseratILQRController();
		~CosseratILQRController() override = default;

	private:
		/**
		 * @brief Extract the current tip SE(3) pose from the mapping output.
		 * @param tip_idx  Frame index to use as tip (-1 = last frame).
		 */
		[[nodiscard]] SE3Types getCurrentTipPose(int tip_idx) const;

		/**
		 * @brief Convert an OutCoord (Rigid3 position+quaternion) to SE3.
		 */
		[[nodiscard]] static SE3Types rigidToSE3(const OutCoord &coord);

		/**
		 * @brief One gradient-descent update step.
		 *
		 * Δξ = -α · J_tip^T · Q · e   (unpacked to NStrainDOF per section)
		 *
		 * @param J_tip   6×6N tip Jacobian
		 * @param e       6D tip error log(g_ref⁻¹·g_tip)
		 * @param alpha   step size
		 * @param Q       tip error weight
		 * @return        Δξ as a single 6N vector
		 */
		[[nodiscard]] static VectorXd
		gradientStep(const MatrixXd &J_tip, const TangentVector &e,
		             double alpha, double Q);

		/**
		 * @brief One Gauss-Newton update step.
		 *
		 * Solves (J^T·Q·J + R·I) · Δξ = -J^T·Q·e via Cholesky decomposition.
		 *
		 * @param J_tip   6×6N tip Jacobian
		 * @param e       6D tip error
		 * @param Q       tip weight
		 * @param R       regularisation weight
		 * @return        Δξ as a single 6N vector
		 */
		[[nodiscard]] static VectorXd
		gaussNewtonStep(const MatrixXd &J_tip, const TangentVector &e,
		                double Q, double R);

		/**
		 * @brief Unpack a 6N delta vector into N per-section Coord1 corrections.
		 *
		 * Only the first NStrainDOF components of each 6D block are retained
		 * (discards the linear part for Vec3 strain inputs).
		 *
		 * @param delta_xi  6N flattened Δξ vector
		 * @param N         number of sections
		 * @return          N Coord1 corrections
		 */
		[[nodiscard]] static std::vector<Coord1>
		unpackDelta(const VectorXd &delta_xi, int N);
	};

#if !defined(SOFA_COSSERAT_CPP_CosseratILQRController)
	extern template class SOFA_COSSERAT_API CosseratILQRController<
		sofa::defaulttype::Vec3Types,
		sofa::defaulttype::Rigid3Types,
		sofa::defaulttype::Rigid3Types>;
#endif

} // namespace Cosserat::controller
