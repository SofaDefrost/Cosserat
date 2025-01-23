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
#define SOFA_COSSERAT_CPP_RigidDistanceMapping
#include <Cosserat/mapping/RigidDistanceMapping.inl>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/ObjectFactory.h>

using namespace sofa::defaulttype;
namespace Cosserat::mapping
{
template class SOFA_COSSERAT_API RigidDistanceMapping< Rigid3Types, Rigid3Types, Rigid3Types >;
} // namespace Cosserat::mapping

namespace Cosserat
{

// Register in the Factory
void registerRigidDistanceMapping(sofa::core::ObjectFactory* factory)
{
  factory->registerObjects(sofa::core::ObjectRegistrationData(
      "This is a mapping class that handles relationships between two input models (In1, In2) and an "
      "output model (Out). It specifically computes and maintains the spatial relationships (both position and "
      "orientation) between pairs of rigid bodies. It handles the mapping between the two input models and an the "
      "output model, calculating relative transformations using quaternions and providing mechanisms for force "
      "feedback and constraint handling. ")
      .add< mapping::RigidDistanceMapping< Rigid3Types, Rigid3Types, Rigid3Types > >()) ;
}
}