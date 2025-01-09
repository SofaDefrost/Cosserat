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
#define SOFA_COSSERAT_CPP_DifferenceMultiMapping
#include <Cosserat/mapping/DifferenceMultiMapping.inl>

#include <sofa/defaulttype/VecTypes.h>
#include <sofa/core/ObjectFactory.h>

namespace Cosserat
{

using namespace sofa::defaulttype;

template class SOFA_COSSERAT_API mapping::DifferenceMultiMapping< sofa::defaulttype::Vec3Types, sofa::defaulttype::Vec3Types, sofa::defaulttype::Vec3Types >;

// Register in the Factory
void registerDifferenceMultiMapping(sofa::core::ObjectFactory* factory)
{
  factory->registerObjects(sofa::core::ObjectRegistrationData("Set the positions and velocities of points attached to a rigid parent")
        .add< mapping::DifferenceMultiMapping< sofa::defaulttype::Vec3Types, sofa::defaulttype::Vec3Types, sofa::defaulttype::Vec3Types > >()) ;
}

} // namespace Cosserat
