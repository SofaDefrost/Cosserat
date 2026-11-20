/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, development version     *
*                (c) 2006-2019 INRIA, USTL, UJF, CNRS, MGH                    *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
******************************************************************************/
#define SOFA_COSSERAT_CPP_LinearStrainElementMapping
#include <Cosserat/mapping/LinearStrainElementMapping.inl>
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/VecTypes.h>

namespace Cosserat::mapping {
	template class SOFA_COSSERAT_API LinearStrainElementMapping<
			sofa::defaulttype::Rigid3Types, sofa::defaulttype::Vec6Types, sofa::defaulttype::Vec6Types>;
} // namespace Cosserat::mapping

namespace Cosserat {
	void registerLinearStrainElementMapping(sofa::core::ObjectFactory *factory) {
		factory->registerObjects(
				sofa::core::ObjectRegistrationData(
						"Linear Strain Element mapping (Xun & Tamadazte, Comput. Struct. 325 (2026) 108181): "
						"reconstructs the mean strain of each rod element from the two nodal Rigid3 poses and a "
						"strain slope, xi_bar = A^-1 Log(g_a^-1 g_b) with A = h I - h^3/12 ad_beta "
						"(fourth-order Magnus). Elements are node pairs, so rod networks and closed loops are "
						"supported. beta = 0 gives the constant-strain element of Frames2StrainCosseratMapping.")
						.add<mapping::LinearStrainElementMapping<sofa::defaulttype::Rigid3Types,
																 sofa::defaulttype::Vec6Types,
																 sofa::defaulttype::Vec6Types>>(true));
	}
} // namespace Cosserat
