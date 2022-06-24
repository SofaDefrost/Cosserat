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

#pragma once
#include "PhysicalState.inl"
#include <sofa/core/ObjectFactory.h>
#include <SofaBaseMechanics/MechanicalObject.h>

namespace sofa::component::container
{

using namespace core::behavior;
using namespace defaulttype;

SOFA_DECL_CLASS(PhysicalState)

int PhysicalStateClass = core::RegisterObject("Cosserat state vectors: torsion, bending on Y, bensing on z !")
.add< PhysicalState<Vec3Types> >(true) // default template
.add< PhysicalState<Vec2Types> >()
.add< PhysicalState<Vec1Types> >()
.add< PhysicalState<Vec6Types> >()
.add< PhysicalState<Rigid3Types> >()
.add< PhysicalState<Rigid2Types> >()
;

// template specialization must be in the same namespace as original namespace for GCC 4.1
// g++ 4.1 requires template instantiations to be declared on a parent namespace from the template class.
template class PhysicalState<Vec3Types>;
template class PhysicalState<Vec2Types>;
template class PhysicalState<Vec1Types>;
template class PhysicalState<Vec6Types>;
template class PhysicalState<Rigid3Types>;
template class PhysicalState<Rigid2Types>;

} // namespace sofa::component::container
