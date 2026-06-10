#pragma once

#include <Cosserat/config.h>
#include <Cosserat/mapping/CosseratBeamGeometry.h>
#include <sofa/core/Multi2Mapping.h>
#include <sofa/core/objectmodel/BaseContext.h>

namespace Cosserat::mapping {

	/**
	 * @brief SOFA Multi2Mapping carrying the mechanical states wired up by the
	 *        scene graph (strain, rigid base, output frames) plus the geometric
	 *        machinery inherited from CosseratBeamGeometry<TIn1>.
	 *
	 * Inheritance model (Z1):
	 *   - Multi2Mapping<TIn1, TIn2, TOut> provides the SOFA Mapping interface
	 *     and the BaseObject lineage.
	 *   - CosseratBeamGeometry<TIn1> is a non-BaseObject mixin that brings the
	 *     section / frame / topology state. It exposes Data MEMBERS but does
	 *     not register them — CosseratGeometryMapping's constructor calls
	 *     initData(...) on them so they appear as Data of this BaseObject.
	 *
	 * This design avoids the diamond that would arise if CosseratBeamGeometry
	 * also inherited from BaseObject.
	 */
	template<class TIn1, class TIn2, class TOut>
	class CosseratGeometryMapping
		: public sofa::core::Multi2Mapping<TIn1, TIn2, TOut>
		, public CosseratBeamGeometry<TIn1>
	{
	public:
		SOFA_ABSTRACT_CLASS(SOFA_TEMPLATE3(CosseratGeometryMapping, TIn1, TIn2, TOut),
							SOFA_TEMPLATE3(sofa::core::Multi2Mapping, TIn1, TIn2, TOut));

		using In1     = TIn1;
		using In2     = TIn2;
		using Out     = TOut;
		using Inherit = sofa::core::Multi2Mapping<TIn1, TIn2, TOut>;
		using Geometry = CosseratBeamGeometry<TIn1>;

		// Re-export the Data members owned by the geometry mixin so that
		// scene code can find them via the usual `d_curv_abs_section.…` form.
		using Geometry::d_curv_abs_section;
		using Geometry::d_curv_abs_frames;
		using Geometry::d_debug;

	public:
		void init() override;

		/// Hook called by init() AFTER the mechanical states are wired up but
		/// BEFORE updateGeometryInfo() / initializeSectionProperties() run.
		virtual void doBaseCosseratInit() = 0;

		/// Hook for derived classes to seed per-frame initial transforms after
		/// the mechanical states are bound but before geometry is built.
		virtual void initialization() = 0;

	protected:
		// Mechanical-state pointers. Names retain the mapping-direction-neutral
		// labels for now (legacy compatibility with Strain2Frames and
		// Frames2Strain consumers); renaming to a fully neutral convention is
		// deferred to a dedicated rename pass.
		sofa::core::State<In1> *m_strain_state = nullptr; ///< In1 input (strains for Strain2Frames)
		sofa::core::State<In2> *m_rigid_base   = nullptr; ///< In2 input (rigid base)
		sofa::core::State<Out> *m_frames       = nullptr; ///< Out (output frames for Strain2Frames)

	protected:
		CosseratGeometryMapping();
		~CosseratGeometryMapping() override = default;

		CosseratGeometryMapping &operator=(const CosseratGeometryMapping &) = delete;
	};

#if !defined(SOFA_COSSERAT_CPP_CosseratGeometryMapping)
	extern template class SOFA_COSSERAT_API CosseratGeometryMapping<
			sofa::defaulttype::Vec3Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types>;
#endif

} // namespace Cosserat::mapping
