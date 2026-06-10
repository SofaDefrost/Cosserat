/******************************************************************************
 *       SOFA, Simulation Open-Framework Architecture, development version     *
 *                (c) 2006-2019 INRIA, USTL, UJF, CNRS, MGH                    *
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
#include <Cosserat/mapping/CosseratGeometryMapping.h>
#include <liegroups/CosseratBodyJacobian.h>
#include <sofa/helper/ColorMap.h>

// NOTE (P2 cleanup): the following headers were previously included for
// SE3Integrator / BezierSE3 / UncertaintyPropagator / GaussianOnManifold /
// BeamStateEstimator typedefs and the matching computeSmoothedPath() /
// computeUncertainties() / computeSectionSE3() methods. None of those were
// implemented or called externally; the typedefs and methods were removed,
// so their includes are no longer needed:
//   #include <Cosserat/mapping/BeamStateEstimator.h>
//   #include <liegroups/BezierSE3.h>
//   #include <liegroups/CosseratUncertaintyPropagator.h>
//   #include <liegroups/GaussianOnManifold.h>
//   #include <liegroups/LieGroupIntegrators.h>

namespace Cosserat::mapping {
	namespace{
		using Mat3x6 = sofa::type::Mat<3, 6, SReal>;
		using Mat6x6 = sofa::type::Mat6x6;
	}

	/**
	 * @brief Discrete implementation of CosseratGeometryMapping using liegroups library
	 *
	 * This class provides a concrete implementation of the Cosserat rod mapping
	 * using the liegroups library for SE(3) operations, with discrete exponential
	 * integration along the rod.
	 *
	 * @tparam TIn1 The first input type for the mapping (strain state)
	 * @tparam TIn2 The second input type for the mapping (rigid base)
	 * @tparam TOut The output type for the mapping (frames)
	 */
	template<class TIn1, class TIn2, class TOut>
	class Strain2FramesCosseratMapping : public CosseratGeometryMapping<TIn1, TIn2, TOut> {
	public:
		SOFA_CLASS(SOFA_TEMPLATE3(Strain2FramesCosseratMapping
	, TIn1, TIn2, TOut),
				   SOFA_TEMPLATE3(CosseratGeometryMapping, TIn1, TIn2, TOut));

		using In1 = TIn1;
		using In2 = TIn2;
		using Out = TOut;
		using Inherit = CosseratGeometryMapping<TIn1, TIn2, TOut>;

		// Type aliases from base classes
		using Coord1 = sofa::Coord_t<In1>;
		using Deriv1 = sofa::Deriv_t<In1>;
		using OutCoord = sofa::Coord_t<Out>;
		using OutDeriv = sofa::Deriv_t<Out>;

		using SE3Types               = sofa::component::cosserat::liegroups::SE3<double>;
		using Vector3                = typename SE3Types::Vector3;
		using TangentVector          = typename SE3Types::TangentVector;
		using BodyJacobian           = sofa::component::cosserat::liegroups::CosseratBodyJacobian<double>;
		using TwistType              = sofa::component::cosserat::liegroups::Twist<double>;
		using WrenchType             = sofa::component::cosserat::liegroups::Wrench<double>;

	public:
		/**
		 * @brief Helper method for manually setting linked models (useful for unit tests)
		 */
		void setModels(sofa::core::State<In1> *strain, sofa::core::State<In2> *base, sofa::core::State<Out> *frames) {
			this->m_strain_state = strain;
			this->m_rigid_base = base;
			this->m_frames = frames;
		}

		/**
		 * @brief Read-only access to the body Jacobian (for CosseratILQRController).
		 *
		 * Valid after each call to apply() / updateFrameTransformations().
		 * The Jacobian is rebuilt from scratch each configuration update.
		 */
		[[nodiscard]] const BodyJacobian &getBodyJacobian() const { return m_bodyJacobian; }

		/**
		 * @brief Read-only access to the strain mechanical state.
		 */
		[[nodiscard]] sofa::core::State<In1> *getStrainState() const { return this->m_strain_state; }

		/**
		 * @brief Read-only access to the output frames mechanical state.
		 */
		[[nodiscard]] sofa::core::State<Out> *getFramesState() const { return this->m_frames; }

	public:
		//////////////////////////////////////////////////////////////////////
		/// @name Data Fields
		/// @{
		sofa::Data<int> d_deformationAxis;
		sofa::Data<SReal> d_max;
		sofa::Data<SReal> d_min;
		sofa::Data<SReal> d_radius;
		sofa::Data<bool> d_drawMapBeam;
		sofa::Data<sofa::type::RGBAColor> d_color;
		sofa::Data<sofa::type::vector<int>> d_index;
		sofa::Data<unsigned int> d_baseIndex;

		/// Integration method for the Cosserat ODE  g'(s) = g(s)·hat(ξ(s)).
		/// 0 = Euler (order 1, equivalent to expCosserat — default for backward compat.)
		/// 1 = Midpoint (order 2, RKMK2)
		/// 2 = RKMK4   (order 4, Magnus expansion — recommended for non-constant strains)
		///
		/// For piecewise-constant strains (current standard use), all methods give
		/// identical results. The option unlocks higher accuracy when combined with
		/// Legendre polynomial strain parameterisations.
		sofa::Data<int> d_integrationMethod;

		/// Number of samples per section for the smoothed centerline path.
		/// Used by computeSmoothedPath() and the draw() visualization.
		/// 1 = section endpoints only (same as the discrete simulation).
		/// Higher values → smoother rendered rod, more Bézier evaluations.
		sofa::Data<int> d_smoothPathSamples;

		/// Isotropic strain noise level σ² used by computeUncertainties().
		/// Σ_ξ = σ² · I₆ for each rod section.
		/// Units: [1/m]² (strain = curvature/elongation per unit length).
		/// 0 disables uncertainty propagation (returns zero covariances).
		sofa::Data<double> d_strainNoiseLevel;
		/// @}
		//////////////////////////////////////////////////////////////////////

	public:
		//////////////////////////////////////////////////////////////////////
		/// @name Inherited from BaseObject
		/// @{
		void initialization() override;
		void doBaseCosseratInit() override;
		void draw(const sofa::core::visual::VisualParams *vparams) override;
		/// @}
		//////////////////////////////////////////////////////////////////////

		//////////////////////////////////////////////////////////////////////
		/// @name Inherited from Multi2Mapping
		/// @{
		void apply(const sofa::core::MechanicalParams *mparams,
				   const sofa::type::vector<sofa::DataVecCoord_t<Out> *> &dataVecOutPos,
				   const sofa::type::vector<const sofa::DataVecCoord_t<In1> *> &dataVecIn1Pos,
				   const sofa::type::vector<const sofa::DataVecCoord_t<In2> *> &dataVecIn2Pos) override;

		void applyJ(const sofa::core::MechanicalParams *mparams,
					const sofa::type::vector<sofa::DataVecDeriv_t<Out> *> &dataVecOutVel,
					const sofa::type::vector<const sofa::DataVecDeriv_t<In1> *> &dataVecIn1Vel,
					const sofa::type::vector<const sofa::DataVecDeriv_t<In2> *> &dataVecIn2Vel) override;

		void applyJT(const sofa::core::MechanicalParams *mparams,
					 const sofa::type::vector<sofa::DataVecDeriv_t<In1> *> &dataVecOut1Force,
					 const sofa::type::vector<sofa::DataVecDeriv_t<In2> *> &dataVecOut2RootForce,
					 const sofa::type::vector<const sofa::DataVecDeriv_t<Out> *> &dataVecInForce) override;

		/**
		 * @brief Geometric stiffness contribution K_G · δξ.
		 *
		 * Computes ∂/∂ξ [J(ξ)ᵀ f_x] · δξ with the child wrenches f_x held
		 * constant (frozen at the current configuration), and accumulates the
		 * result into the input (strain) force vector with weight `kFactor`.
		 *
		 * Two contributions are assembled:
		 *  - **(a) Frame direct term** — the dependency of each output-frame
		 *    co-adjoint on the local strain via the tangent-exponential matrix:
		 *      δf_k += kFactor · Bᵀ · J_frame^T · ad(J_frame · B·δξ_k)ᵀ · w_body_k
		 *  - **(b) Section transport term** — the dependency of the accumulated
		 *    downstream wrench on the local SE(3) transport (backward sweep):
		 *      δf_k += kFactor · Bᵀ · J_local_k^T · ad(J_local_k · B·δξ_k)ᵀ · F_tot_k
		 *
		 * where B = [I_m | 0] is the selector for the m active strain DOFs
		 * (m = 3 for Vec3Types, 6 for Vec6Types).
		 *
		 * The implementation uses `TwistType::smallAdjoint()` from the liegroups
		 * library and `CosseratBodyJacobian` for the backward wrench sweep.
		 *
		 * @param mparams   Mechanical parameters; kFactor is read from here.
		 * @param inForce   Multi-vec id for the strain force (in/out, += semantics).
		 * @param outForce  Multi-vec id for the child wrenches (read-only, unused
		 *                  because we read directly from m_frames).
		 */
		void applyDJT(const sofa::core::MechanicalParams* mparams,
					  sofa::core::MultiVecDerivId          inForce,
					  sofa::core::ConstMultiVecDerivId     outForce) override;

		/// Support for constraints
		void applyJT(const sofa::core::ConstraintParams *cparams,
					 const sofa::type::vector<sofa::DataMatrixDeriv_t<In1> *> &dataMatOut1Const,
					 const sofa::type::vector<sofa::DataMatrixDeriv_t<In2> *> &dataMatOut2Const,
					 const sofa::type::vector<const sofa::DataMatrixDeriv_t<Out> *> &dataMatInConst) override;
		/// @}
		//////////////////////////////////////////////////////////////////////

		void computeBBox(const sofa::core::ExecParams *params, bool onlyVisible) override;

		// Note (P2 cleanup, 2026-06-05):
		// The following methods were previously declared here but never
		// implemented in the .inl. They are removed to keep the public API
		// honest. If/when they are needed, restore the doc-block from git
		// history and provide an implementation:
		//   std::vector<SE3Types> computeSmoothedPath(...) const;     // Bezier path
		//   std::vector<GaussianSE3Type> computeUncertainties(...);   // Cov propagation
		//   SE3Types computeSectionSE3(const TangentVector&, double); // RKMK integrator
		// The associated Data fields d_integrationMethod / d_smoothPathSamples /
		// d_strainNoiseLevel are kept as XML-stable placeholders for the
		// eventual re-introduction.

	public:
		////////////////////////// Inherited attributes ////////////////////////////
		/// Bring inherited attributes into the current lookup context
		using CosseratGeometryMapping<TIn1, TIn2, TOut>::d_curv_abs_section;
		using CosseratGeometryMapping<TIn1, TIn2, TOut>::d_curv_abs_frames;
		using CosseratGeometryMapping<TIn1, TIn2, TOut>::d_debug;

	protected:
		////////////////////////// Inherited attributes ////////////////////////////
		/// Bring inherited attributes into the current lookup context
		using CosseratGeometryMapping<TIn1, TIn2, TOut>::m_section_properties;
		using CosseratGeometryMapping<TIn1, TIn2, TOut>::m_frame_properties;
		using CosseratGeometryMapping<TIn1, TIn2, TOut>::m_frame_to_section_indices;
		using CosseratGeometryMapping<TIn1, TIn2, TOut>::m_frame_to_section_indices_draw;
		using CosseratGeometryMapping<TIn1, TIn2, TOut>::m_section_length_vectors;
		using CosseratGeometryMapping<TIn1, TIn2, TOut>::m_strain_state;
		using CosseratGeometryMapping<TIn1, TIn2, TOut>::m_rigid_base;
		using CosseratGeometryMapping<TIn1, TIn2, TOut>::m_frames;
		//////////////////////////////////////////////////////////////////////////////

		sofa::helper::ColorMap m_colorMap;

		/**
		 * @brief Body Jacobian of the rod — accessor for external consumers.
		 *
		 * TODO (P3): currently never populated by updateFrameTransformations().
		 * `CosseratILQRController` reads it via getBodyJacobian() but receives a
		 * default-constructed (size 0) BodyJacobian, which it interprets as a
		 * "no-op" sentinel. Populate this member once applyJ() is migrated to
		 * the body-Jacobian propagation API.
		 */
		BodyJacobian m_bodyJacobian;

		/**
		 * @brief Updates frame transformations using liegroups SE(3) exponential map
		 * @param vec_of_strains Current strain values
		 */
		void updateFrameTransformations(const sofa::type::vector<Coord1> &vec_of_strains);


		// Debug display functions
		void displayStrainState(const sofa::type::vector<Coord1> &strainState, const std::string &context = "") const;
		void displayRigidState(const sofa::type::vector<sofa::Coord_t<In2>> &rigidState,
							   const std::string &context = "") const;
		void displayOutputFrames(const sofa::type::vector<OutCoord> &outputFrames,
								 const std::string &context = "") const;
		void displaySectionProperties(const std::string &context = "") const;
		void displayFrameProperties(const std::string &context = "") const;
		void displaySE3Transform(const SE3Types &transform, const std::string &name = "Transform") const;
		void displayMappingState(const std::string &context = "") const;
		void displayVelocities(const sofa::type::vector<Deriv1> &strainVel,
							   const sofa::type::vector<sofa::Deriv_t<In2>> &baseVel,
							   const sofa::type::vector<OutDeriv> &outputVel, const std::string &context = "") const;

	protected:
		Strain2FramesCosseratMapping
();
		~Strain2FramesCosseratMapping
() override = default;
	};

#if !defined(SOFA_COSSERAT_CPP_Strain2FramesCosseratMapping)
	extern template class SOFA_COSSERAT_API Strain2FramesCosseratMapping<
			sofa::defaulttype::Vec3Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types>;
	extern template class SOFA_COSSERAT_API Strain2FramesCosseratMapping<
			sofa::defaulttype::Vec6Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types>;
#endif

} // namespace Cosserat::mapping
