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
#define SOFA_COSSERAT_CPP_Frames2StrainCosseratMapping
#include <Cosserat/mapping/Frames2StrainCosseratMapping.inl>
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/VecTypes.h>

namespace Cosserat::mapping {
	template class SOFA_COSSERAT_API Frames2StrainCosseratMapping<
			sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Vec3Types>;
	// template class SOFA_COSSERAT_API Frames2StrainCosseratMapping<
	// 		sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Vec6Types>;
} // namespace Cosserat::mapping

namespace Cosserat{
	// Register in the Factory
	void registerFrames2StrainCosseratMapping(sofa::core::ObjectFactory* factory)
	{
		factory->registerObjects(
				sofa::core::ObjectRegistrationData(
						"This component computes Cosserat strains from input frames and rigid base. "
						"It takes two mechanical" 
						"objects as inputs; the rigid base of the beam (with 6 degree of freedom) and the frames (Rigid3) "
						"in global coordinates." 
						"Using these inputs, the component computes and outputs the strains (Vec3) of the beam.")
						.add<mapping::Frames2StrainCosseratMapping<sofa::defaulttype::Rigid3Types, 
																sofa::defaulttype::Rigid3Types, 
																sofa::defaulttype::Vec3Types>>(true));
			//.add< mapping::Frames2StrainCosseratMapping< sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Vec6Types> >());
	}
}// namespace Cosserat