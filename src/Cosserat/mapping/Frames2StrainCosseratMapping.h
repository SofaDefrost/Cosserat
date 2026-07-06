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
 *******************************************************************************
 * Authors: The SOFA Team and external contributors (see Authors.txt)          *
 *                                                                             *
 * Contact information: contact@sofa-framework.org                             *
 ******************************************************************************/
#pragma once

#include <Cosserat/config.h>
#include <Cosserat/mapping/CosseratBeamGeometry.h>
#include <sofa/core/Mapping.h>
#include <sofa/helper/ColorMap.h>

namespace Cosserat::mapping {

	namespace {
		using Mat3x6 = sofa::type::Mat<3, 6, SReal>;
		using Mat6x6 = sofa::type::Mat6x6;
	}

	/**
	 * @brief SOFA Mapping computing per-section Cosserat strains from absolute
	 *        Rigid3 frames along a rod.
	 *
	 * Mathematically the mapping is 1-to-1:
	 *
	 *   ξ_k = log( g_a⁻¹ · g_b ) / h_k          for each section k
	 *
	 * where g_a, g_b are the two boundary frames of the k-th section. The
	 * rigid base that used to be a second input is geometrically *invariant*
	 * (strain is local), so this class inherits from a single-input
	 * `sofa::core::Mapping<TIn, TOut>` rather than the historical
	 * `Multi2Mapping`. See FRAMES2STRAIN_ANALYSIS.md §7 for the rationale.
	 *
	 * Inheritance model (Z1):
	 *   - sofa::core::Mapping<TIn, TOut> provides the SOFA Mapping interface
	 *     and the BaseObject lineage.
	 *   - CosseratBeamGeometry<TOut> is the non-BaseObject mixin holding the
	 *     section / frame topology and tangent-exp matrices. Templated on the
	 *     strain (output) type. Its Data members are registered with this
	 *     BaseObject via addData() in the constructor.
	 *
	 * @tparam TIn  Frames input type — Rigid3Types.
	 * @tparam TOut Strain output type — Vec3Types or Vec6Types.
	 */
	template<class TIn, class TOut>
	class Frames2StrainCosseratMapping
		: public sofa::core::Mapping<TIn, TOut>
		, public CosseratBeamGeometry<TOut>
	{
	public:
		SOFA_CLASS(SOFA_TEMPLATE2(Frames2StrainCosseratMapping, TIn, TOut),
				   SOFA_TEMPLATE2(sofa::core::Mapping, TIn, TOut));

		using In       = TIn;
		using Out      = TOut;
		using Inherit  = sofa::core::Mapping<TIn, TOut>;
		using Geometry = CosseratBeamGeometry<TOut>;

		using Coord    = sofa::Coord_t<In>;
		using Deriv    = sofa::Deriv_t<In>;
		using OutCoord = sofa::Coord_t<Out>;
		using OutDeriv = sofa::Deriv_t<Out>;

		using SE3Types      = sofa::component::cosserat::liegroups::SE3<double>;
		using Vector3       = typename SE3Types::Vector3;
		using TangentVector = typename SE3Types::TangentVector;

		// Re-export the geometry mixin's Data members
		using Geometry::d_curv_abs_section;
		using Geometry::d_curv_abs_frames;
		using Geometry::d_debug;

	public:
		// ── Data fields (visualization) ──────────────────────────────────
		sofa::Data<int>                     d_deformationAxis;
		sofa::Data<SReal>                   d_max;
		sofa::Data<SReal>                   d_min;
		sofa::Data<SReal>                   d_radius;
		sofa::Data<bool>                    d_drawMapBeam;
		sofa::Data<sofa::type::RGBAColor>   d_color;
		sofa::Data<sofa::type::vector<int>> d_index;

	public:
		// ── Inherited from BaseObject ─────────────────────────────────────
		void init() override;
		void draw(const sofa::core::visual::VisualParams *vparams) override;

		// ── Inherited from Mapping<TIn, TOut> ─────────────────────────────
		void apply(const sofa::core::MechanicalParams *mparams,
				   sofa::DataVecCoord_t<Out> &out,
				   const sofa::DataVecCoord_t<In> &in) override;

		void applyJ(const sofa::core::MechanicalParams *mparams,
					sofa::DataVecDeriv_t<Out> &out,
					const sofa::DataVecDeriv_t<In> &in) override;

		void applyJT(const sofa::core::MechanicalParams *mparams,
					 sofa::DataVecDeriv_t<In> &out,
					 const sofa::DataVecDeriv_t<Out> &in) override;

		/// Geometric stiffness K_G · δξ — NOT YET IMPLEMENTED (one-time warning).
		void applyDJT(const sofa::core::MechanicalParams *mparams,
					  sofa::core::MultiVecDerivId          inForce,
					  sofa::core::ConstMultiVecDerivId     outForce) override;

		/// Constraint version
		void applyJT(const sofa::core::ConstraintParams *cparams,
					 sofa::DataMatrixDeriv_t<In> &out,
					 const sofa::DataMatrixDeriv_t<Out> &in) override;

		void computeBBox(const sofa::core::ExecParams *params, bool onlyVisible) override;

		// Note: the exact SE(3) exp-map Jacobians (Jr⁻¹ etc.) now live in
		// liegroups/SE3.h as SE3::computeRightJacobianInverse(...); the mapping
		// calls them directly (no local dexp helpers).

	public:
		// Re-export geometry state for convenience inside .inl
		using Geometry::m_section_properties;
		using Geometry::m_frame_properties;
		using Geometry::m_frame_to_section_indices;
		using Geometry::m_frame_to_section_indices_draw;
		using Geometry::m_section_length_vectors;

	protected:
		sofa::helper::ColorMap m_colorMap;

		// Mechanical-state pointers cached at init().
		sofa::core::State<In>  *m_input  = nullptr; ///< Rigid3 frames
		sofa::core::State<Out> *m_output = nullptr; ///< Vec3 / Vec6 strains

	protected:
		Frames2StrainCosseratMapping();
		~Frames2StrainCosseratMapping() override = default;
	};

#if !defined(SOFA_COSSERAT_CPP_Frames2StrainCosseratMapping)
	extern template class SOFA_COSSERAT_API Frames2StrainCosseratMapping<
			sofa::defaulttype::Rigid3Types, sofa::defaulttype::Vec3Types>;
	extern template class SOFA_COSSERAT_API Frames2StrainCosseratMapping<
			sofa::defaulttype::Rigid3Types, sofa::defaulttype::Vec6Types>;
#endif

} // namespace Cosserat::mapping
