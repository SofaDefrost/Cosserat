/*
 * HookeSeratBaseMapping.inl
 * Implementation details for the HookeSeratBaseMapping class.
 * This file contains implementations for functions inline in the HookeSeratBaseMapping class.
 */
#pragma once

#include <Cosserat/config.h>
#include <Cosserat/mapping/HookeSeratBaseMapping.h>

namespace Cosserat::mapping {


	template <class TIn1, class TIn2, class TOut>
 HookeSeratBaseMapping<TIn1, TIn2, TOut>::HookeSeratBaseMapping()
	 // TODO(dmarchal: 2024/06/12): please add the help comments !
	 : d_curv_abs_section(initData(&d_curv_abs_section, "curv_abs_input",
								   " need to be com....")),
	   d_curv_abs_frames(initData(&d_curv_abs_frames, "curv_abs_output",
								  " need to be com....")),
	   d_debug(initData(&d_debug, false, "debug", "printf for the debug")),
	   m_indexInput(0) {}


using namespace sofa::component::cosserat::liegroups;

// Implementation of the constructor for the HookeSeratBaseMapping class
// Initialize the HookeSeratBaseMapping
template <class TIn1, class TIn2, class TOut>
void HookeSeratBaseMapping<TIn1, TIn2, TOut>::init(){
	try {
		// Validation des données
		if (!validateSectionProperties()) {
			throw std::runtime_error("Invalid section properties");
		}

		// Vérification de la continuité
		if (!checkContinuity()) {
			msg_warning() << "Rod sections are not continuous";
		}

		// Pré-calcul des matrices adjointes (se fait automatiquement avec le cache)
		for (const auto& section : m_sectionProperties) {
			section.getAdjoint(); // Force le calcul
		}

		// Initialisation des états
		m_strain_state = dynamic_cast<sofa::core::State<In1>*>(this->fromModels1[0].get());
		m_rigid_base = dynamic_cast<sofa::core::State<In2>*>(this->fromModels2[0].get());
		m_frames = dynamic_cast<sofa::core::State<Out>*>(this->toModels[0].get());

		if (!m_strain_state || !m_rigid_base || !m_frames) {
			throw std::runtime_error("Failed to initialize mechanical states");
		}

		Inherit::init();

	} catch (const std::exception& e) {
		msg_error() << "Initialization failed: " << e.what();
		throw;
	}
}

} // namespace Cosserat::mapping

