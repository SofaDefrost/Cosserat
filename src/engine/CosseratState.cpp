/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 beta 4      *
*                (c) 2006-2009 MGH, INRIA, USTL, UJF, CNRS                    *
*                                                                             *
* This library is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This library is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this library; if not, write to the Free Software Foundation,     *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
*******************************************************************************
*                               SOFA :: Modules                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/

#define SOFA_COMPONENT_CONTAINER_TEMPERATURESTATE_CPP
#include "TemperatureState.inl"
#include <sofa/core/ObjectFactory.h>
#include <SofaBaseMechanics/MechanicalObject.h>

namespace sofa
{

namespace component
{

namespace container
{

using namespace core::behavior;
using namespace defaulttype;

SOFA_DECL_CLASS(TemperatureState)

int TemperatureStateClass = core::RegisterObject("electrical state vectors")
.add< TemperatureState<Vec3Types> >(true) // default template
.add< TemperatureState<Vec2Types> >()
.add< TemperatureState<Vec1Types> >()
.add< TemperatureState<Vec6Types> >()
.add< TemperatureState<Rigid3Types> >()
.add< TemperatureState<Rigid2Types> >()

;

// template specialization must be in the same namespace as original namespace for GCC 4.1
// g++ 4.1 requires template instantiations to be declared on a parent namespace from the template class.
template class TemperatureState<Vec3Types>;
template class TemperatureState<Vec2Types>;
template class TemperatureState<Vec1Types>;
template class TemperatureState<Vec6Types>;
template class TemperatureState<Rigid3Types>;
template class TemperatureState<Rigid2Types>;



}

} // namespace component

} // namespace sofa
