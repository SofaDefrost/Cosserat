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
#define SOFA_COMPONENT_MAPPING_POE_MAPING_CPP
#include "RigidRigidCosseratMapping.inl"
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/ObjectFactory.h>

namespace sofa
{

namespace component
{

namespace mapping
{
using namespace sofa::defaulttype;

// Register in the Factory
int RigidRigidCosseratMappingClass = core::RegisterObject("Set the positions and velocities of points attached to a rigid parent")
        .add< RigidRigidCosseratMapping< sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types > >() ;


template class SOFA_COSSERAT_MAPPING_API RigidRigidCosseratMapping< sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types >;


} // namespace mapping

} // namespace component

} // namespace sofa.
