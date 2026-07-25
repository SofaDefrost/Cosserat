/*
 * CosseratGeometryMapping.inl
 * Implementation details for the CosseratGeometryMapping class.
 *
 * The bulk of the geometry helpers (SectionInfo, FrameInfo, BeamTopology,
 * updateGeometryInfo, initializeSectionProperties, initializeFrameProperties,
 * updateTangExpSE3, validateJacobianAccuracy, …) was moved to
 * CosseratBeamGeometry.{h,inl} in commit "extract CosseratBeamGeometry mixin".
 *
 * What remains here:
 *   - Constructor: registers the Data members owned by the geometry mixin
 *     with this BaseObject via initData().
 *   - init(): wires the 3 mechanical states, calls the derived-class hooks,
 *     then delegates the geometry build to the mixin.
 */
#pragma once

#include <Cosserat/config.h>
#include <Cosserat/mapping/CosseratGeometryMapping.h>
#include <iostream>
#include <sofa/core/Multi2Mapping.inl>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/helper/logging/Message.h>

namespace Cosserat::mapping {

	using namespace sofa::component::cosserat::liegroups;
	using sofa::helper::getReadAccessor;
	using sofa::type::Quat;
	using sofa::type::Vec3;
	using sofa::type::Vec6;
	using sofa::type::vector;

	template<class TIn1, class TIn2, class TOut>
	CosseratGeometryMapping<TIn1, TIn2, TOut>::CosseratGeometryMapping()
		: sofa::core::Multi2Mapping<TIn1, TIn2, TOut>()
		, CosseratBeamGeometry<TIn1>()
	{
		// Register the Data members owned by the CosseratBeamGeometry mixin.
		// The standard SOFA `d_foo(initData(...))` pattern requires being in
		// the member-initializer-list, which we cannot use here because the
		// Data fields live in the mixin base (already default-constructed by
		// the time this constructor body runs). Instead we use addData() with
		// explicit name/help setup — the SOFA Base.h says this is the
		// supported alternative when the InitData pattern is not applicable.
		this->d_curv_abs_section.setName("curv_abs_input");
		this->d_curv_abs_section.setHelp("Curvilinear abscissa of the input sections along the rod");
		this->addData(&this->d_curv_abs_section);

		this->d_curv_abs_frames.setName("curv_abs_output");
		this->d_curv_abs_frames.setHelp("Curvilinear abscissa of the output frames along the rod");
		this->addData(&this->d_curv_abs_frames);

		this->d_debug.setName("debug");
		this->d_debug.setHelp("Enable debug output");
		this->d_debug.setValue(false);
		this->addData(&this->d_debug);

		msg_info("CosseratGeometryMapping") << "CosseratGeometryMapping constructor called";
	}

	template<class TIn1, class TIn2, class TOut>
	void CosseratGeometryMapping<TIn1, TIn2, TOut>::init() {
		msg_info("CosseratGeometryMapping") << "Initializing CosseratGeometryMapping...";

		m_strain_state = nullptr;
		m_rigid_base   = nullptr;
		m_frames       = nullptr;

		if (this->fromModels1.empty()) {
			msg_error() << "Input1 (strain state) not found";
			return;
		}
		if (this->fromModels2.empty()) {
			msg_error() << "Input2 (rigid base) not found";
			return;
		}
		if (this->toModels.empty()) {
			msg_error() << "Output (frames) missing";
			return;
		}

		m_strain_state = this->fromModels1[0];
		m_rigid_base   = this->fromModels2[0];
		m_frames       = this->toModels[0];

		// Derived-class hook (seed per-frame transforms etc.)
		initialization();

		// Build geometry from the curvilinear abscissas + current strain.
		this->updateGeometryInfo();

		const auto &strain = m_strain_state->read(sofa::core::vec_id::read_access::position)->getValue();
		this->initializeSectionProperties(strain);

		this->initializeFrameProperties();

		if (!this->validateSectionProperties()) {
			msg_error() << "Invalid section properties detected";
			return;
		}

		// Pre-compute adjoint matrices for performance (caches them in SectionInfo)
		for (const auto &section : this->m_section_properties) {
			section.getAdjoint();
		}

		// Parent (Multi2Mapping) init
		Inherit::init();
	}

} // namespace Cosserat::mapping
