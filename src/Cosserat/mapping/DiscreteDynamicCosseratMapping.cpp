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
#define SOFA_COSSERAT_CPP_DiscreteDynamicCosseratMapping
#include <Cosserat/mapping/DiscreteDynamicCosseratMapping.inl>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/ObjectFactory.h>
using namespace sofa::defaulttype;
namespace Cosserat::mapping
{
template class SOFA_COSSERAT_API DiscreteDynamicCosseratMapping< Vec3Types, Rigid3Types, Rigid3Types >;
} // namespace sofa::component::mapping

namespace Cosserat
{
// Register in the Factory
void registerDiscretDynamicCosseratMapping(sofa::core::ObjectFactory* factory)
{
  factory->registerObjects( sofa::core::ObjectRegistrationData(
  "This component is designed for dynamic simulations, in contrast to the DiscreteCosseratMapping "
  "component, which is primarily used for quasi-static scenarios. ")
  .add<mapping::DiscreteDynamicCosseratMapping<Vec3Types, Rigid3Types, Rigid3Types>>());
}
}
