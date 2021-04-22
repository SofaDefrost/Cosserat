/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
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
#define SOFA_COMPONENT_FORCEFIELD_UNIFORMVELOCITYDAMPINGFORCEFIELD_CPP

#include "MyUniformVelocityDampingForceField.inl"
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>

namespace sofa::component::forcefield
{

using namespace sofa::defaulttype;

int MyUniformVelocityDampingForceFieldClass = core::RegisterObject("Uniform velocity damping")
.add< MyUniformVelocityDampingForceField<Vec3Types> >()
.add< MyUniformVelocityDampingForceField<Vec2Types> >()
.add< MyUniformVelocityDampingForceField<Vec1Types> >()
.add< MyUniformVelocityDampingForceField<Vec6Types> >()
.add< MyUniformVelocityDampingForceField<Rigid3Types> >()
.add< MyUniformVelocityDampingForceField<Rigid2Types> >()

;


    template class SOFA_SOFABOUNDARYCONDITION_API MyUniformVelocityDampingForceField<Vec3Types>;
    template class SOFA_SOFABOUNDARYCONDITION_API MyUniformVelocityDampingForceField<Vec2Types>;
    template class SOFA_SOFABOUNDARYCONDITION_API MyUniformVelocityDampingForceField<Vec1Types>;
    template class SOFA_SOFABOUNDARYCONDITION_API MyUniformVelocityDampingForceField<Vec6Types>;
    template class SOFA_SOFABOUNDARYCONDITION_API MyUniformVelocityDampingForceField<Rigid3Types>;
    template class SOFA_SOFABOUNDARYCONDITION_API MyUniformVelocityDampingForceField<Rigid2Types>;


} // namespace sofa::component::forcefield
