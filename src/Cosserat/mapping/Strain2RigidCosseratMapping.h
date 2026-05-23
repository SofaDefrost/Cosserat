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
#include <Cosserat/mapping/BeamStateEstimator.h>
#include <Cosserat/mapping/CosseratGeometryMapping.h>
#include <liegroups/BezierSE3.h>
#include <liegroups/CosseratBodyJacobian.h>
#include <liegroups/CosseratUncertaintyPropagator.h>
#include <liegroups/GaussianOnManifold.h>
#include <liegroups/LieGroupIntegrators.h>
#include <sofa/helper/ColorMap.h>

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
	class Strain2RigidCosseratMapping : public CosseratGeometryMapping<TIn1, TIn2, TOut> {
	public:
		SOFA_CLASS(SOFA_TEMPLATE3(Strain2RigidCosseratMapping
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

		// using SectionProperties = typename CosseratGeometryMapping<TIn1,TIn2,TOut>::SectionProperties;
		// using FrameInfo = typename FrameInfo;
		using SE3Types               = sofa::component::cosserat::liegroups::SE3<double>;
		using Vector3                = typename SE3Types::Vector3;
		using TangentVector          = typename SE3Types::TangentVector;
		using BodyJacobian           = sofa::component::cosserat::liegroups::CosseratBodyJacobian<double>;
		using TwistType              = sofa::component::cosserat::liegroups::Twist<double>;
		using WrenchType             = sofa::component::cosserat::liegroups::Wrench<double>;
		using SE3Integrator          = sofa::component::cosserat::liegroups::SE3Integrator<double>;
		using BezierSE3Type          = sofa::component::cosserat::liegroups::BezierSE3<double>;
		using UncertaintyPropagator  = sofa::component::cosserat::liegroups::CosseratUncertaintyPropagator<double>;
		using GaussianSE3Type        = sofa::component::cosserat::liegroups::GaussianOnManifold<SE3Types>;

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

		void applyDJT(const sofa::core::MechanicalParams * /*mparams*/, sofa::core::MultiVecDerivId /*inForce*/,
					  sofa::core::ConstMultiVecDerivId /*outForce*/) override {}

		/// Support for constraints
		void applyJT(const sofa::core::ConstraintParams *cparams,
					 const sofa::type::vector<sofa::DataMatrixDeriv_t<In1> *> &dataMatOut1Const,
					 const sofa::type::vector<sofa::DataMatrixDeriv_t<In2> *> &dataMatOut2Const,
					 const sofa::type::vector<const sofa::DataMatrixDeriv_t<Out> *> &dataMatInConst) override;
		/// @}
		//////////////////////////////////////////////////////////////////////

		void computeBBox(const sofa::core::ExecParams *params, bool onlyVisible) override;

		//////////////////////////////////////////////////////////////////////
		/// @name BezierSE3 utilities
		/// @{

		/**
		 * @brief Build a smooth centerline path through the computed section poses.
		 *
		 * Constructs a piecewise-cubic Bézier curve (C1-continuous) in SE(3)
		 * through the current section node poses, then samples it uniformly.
		 *
		 * For each consecutive pair of node poses (g_k, g_{k+1}):
		 *   - Two inner control poses are placed at ⅓ and ⅔ of the geodesic,
		 *     giving a degree-3 segment with C1 continuity at junctions.
		 *   - The curve is sampled at `samplesPerSection` arc-length steps.
		 *
		 * This path is geometrically smoother than the piecewise-exponential
		 * result from apply(), and is useful for:
		 *   - High-quality visualization (draw())
		 *   - Reference trajectory for CosseratILQRController
		 *   - Shape comparison / error analysis
		 *
		 * @note apply() must have been called first (node poses come from
		 *       m_section_properties built by updateFrameTransformations()).
		 *
		 * @param rigidBase        Current rigid base transformation (In2 input).
		 * @param samplesPerSection Number of uniform samples per section
		 *                         (default = d_smoothPathSamples.getValue()).
		 * @return Vector of SE3 poses along the centerline, from base to tip.
		 */
		[[nodiscard]] std::vector<SE3Types>
		computeSmoothedPath(const sofa::type::vector<sofa::Coord_t<In2>> &rigidBase,
		                    int samplesPerSection = -1) const;

		/// @}
		//////////////////////////////////////////////////////////////////////

		//////////////////////////////////////////////////////////////////////
		/// @name Uncertainty propagation
		/// @{

		/**
		 * @brief Propagate pose uncertainty along the rod using CosseratUncertaintyPropagator.
		 *
		 * Starting from the base with zero initial covariance, propagates a
		 * Gaussian distribution over SE(3) through each rod section.  The
		 * strain noise covariance is taken as σ²·I₆ where σ² = d_strainNoiseLevel.
		 *
		 * Returns N+1 Gaussian distributions: [base, after-sec-0, …, after-sec-N-1].
		 * Each Gaussian carries a 6×6 covariance in the body-frame tangent space.
		 *
		 * @note apply() must have been called first so that m_section_properties
		 *       is populated with current strains and lengths.
		 *
		 * @param rigidBase  Current rigid base transformation (In2 input).
		 * @return           N+1 Gaussian SE(3) distributions, root to tip.
		 */
		[[nodiscard]] std::vector<GaussianSE3Type>
		computeUncertainties(const sofa::type::vector<sofa::Coord_t<In2>> &rigidBase) const;

		/// @}
		//////////////////////////////////////////////////////////////////////

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
		using CosseratGeometryMapping<TIn1, TIn2, TOut>::m_frameProperties;
		using CosseratGeometryMapping<TIn1, TIn2, TOut>::m_indices_vectors;
		using CosseratGeometryMapping<TIn1, TIn2, TOut>::m_indices_vectors_draw;
		using CosseratGeometryMapping<TIn1, TIn2, TOut>::m_beam_length_vectors;
		using CosseratGeometryMapping<TIn1, TIn2, TOut>::m_strain_state;
		using CosseratGeometryMapping<TIn1, TIn2, TOut>::m_rigid_base;
		using CosseratGeometryMapping<TIn1, TIn2, TOut>::m_frames;
		//////////////////////////////////////////////////////////////////////////////

		sofa::helper::ColorMap m_colorMap;

		/**
		 * @brief Body Jacobian of the rod — centralises applyJ / applyJT kinematics.
		 *
		 * Built by updateFrameTransformations() after each configuration update.
		 * Sections are stored 0-indexed (section k ↔ m_section_properties[k+1]).
		 *
		 * - applyJ  uses m_bodyJacobian.applyForward()  for the node-level propagation.
		 * - applyJT uses m_bodyJacobian.applyTranspose() for the node-level backward pass.
		 */
		BodyJacobian m_bodyJacobian;

		/**
		 * @brief Left-invariant EKF estimator for the rod base pose.
		 *
		 * Updated by the host SOFA component or an external observer via:
		 *   m_stateEstimator.predict(strain, length, strain_cov)
		 *   m_stateEstimator.update(measurement, R)
		 *
		 * computeUncertainties() uses its current mean as the base pose for
		 * full-rod uncertainty propagation via UncertaintyPropagator.
		 */
		BeamStateEstimator m_stateEstimator;

		/**
		 * @brief Updates frame transformations using liegroups SE(3) exponential map
		 * @param vec_of_strains Current strain values
		 */
		void updateFrameTransformations(const sofa::type::vector<Coord1> &vec_of_strains);

		/**
		 * @brief Single-element SE(3) integration, method chosen by d_integrationMethod.
		 *
		 * @param strain  6D body-frame strain vector
		 * @param length  Element arc-length
		 * @return        Local SE(3) step (not composed with previous elements)
		 */
		[[nodiscard]] SE3Types computeSectionSE3(const TangentVector &strain, double length) const;


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
		Strain2RigidCosseratMapping
();
		~Strain2RigidCosseratMapping
() override = default;
	};

#if !defined(SOFA_COSSERAT_CPP_Strain2RigidCosseratMappin)
	extern template class SOFA_COSSERAT_API Strain2RigidCosseratMapping<
			sofa::defaulttype::Vec3Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types>;
	// Vec6 instantiation is currently disabled
	// extern template class SOFA_COSSERAT_API Strain2RigidCosseratMapping<
	// 		sofa::defaulttype::Vec6Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types>;
#endif

} // namespace Cosserat::mapping
