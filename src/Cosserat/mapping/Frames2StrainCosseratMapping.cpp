/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, development version     *
*                (c) 2006-2019 INRIA, USTL, UJF, CNRS, MGH                    *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
******************************************************************************/
#define SOFA_COSSERAT_CPP_Frames2StrainCosseratMapping
#include <Cosserat/mapping/Frames2StrainCosseratMapping.inl>
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/VecTypes.h>

namespace Cosserat::mapping {
	// Two-argument template instantiation: TIn = Rigid3, TOut = Vec3 / Vec6.
	template class SOFA_COSSERAT_API Frames2StrainCosseratMapping<
			sofa::defaulttype::Rigid3Types, sofa::defaulttype::Vec3Types>;
	template class SOFA_COSSERAT_API Frames2StrainCosseratMapping<
			sofa::defaulttype::Rigid3Types, sofa::defaulttype::Vec6Types>;
} // namespace Cosserat::mapping

namespace Cosserat {
	// Register in the Factory
	void registerFrames2StrainCosseratMapping(sofa::core::ObjectFactory* factory) {
		factory->registerObjects(
				sofa::core::ObjectRegistrationData(
						"Cosserat strain mapping: computes per-section strains (Vec3/Vec6) from absolute "
						"Rigid3 frames along a beam. Inherits from sofa::core::Mapping (single input) — the "
						"rigid base used to be a second input but is geometrically invariant and was removed "
						"in the option-B refactor; see FRAMES2STRAIN_ANALYSIS.md §7.")
						.add<mapping::Frames2StrainCosseratMapping<sofa::defaulttype::Rigid3Types,
																   sofa::defaulttype::Vec3Types>>(true)
						.add<mapping::Frames2StrainCosseratMapping<sofa::defaulttype::Rigid3Types,
																   sofa::defaulttype::Vec6Types>>());
	}
} // namespace Cosserat
