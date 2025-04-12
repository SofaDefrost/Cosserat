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
using namespace sofa::defaulttype;
namespace Cosserat::mapping
{

template class SOFA_COSSERAT_API mapping::DifferenceMultiMapping< sofa::defaulttype::Vec3Types, sofa::defaulttype::Vec3Types, sofa::defaulttype::Vec3Types >;

} // namespace Cosserat::mapping

namespace Cosserat
{
// Register in the Factory
void registerDifferenceMultiMapping(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData(
      "he DifferenceMultiMapping class is a component used to enforce constraints between two sets of points. "
      "Takes two sets of coordinates (x1 and x2) as input, then finds the closest points between these sets of coordinates. "
      "Then compute the force to keep the points in a desired relationship to each other. Particularly useful for "
      "scenarios like cable-structure interactions. ")
        .add< mapping::DifferenceMultiMapping< Vec3Types, Vec3Types, Vec3Types > >()) ;
}

} // namespace sofa.
