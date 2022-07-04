
/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 MGH, INRIA, USTL, UJF, CNRS                    *
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
* This component is not open-source                                           *
*                                                                             *
* Authors: Yinoussa Adagolodjo/adagolodjo@protonamil.com                      *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#include<sofa/defaulttype/VecTypes.h>
#include <sofa/core/ObjectFactory.h>
#include "CosseratNeedleSlidingConstraint.inl"



namespace sofa::component::constraintset
{

using sofa::defaulttype::Rigid3Types;
using namespace sofa::helper;
using namespace sofa::core;

////////////////////////////////////////////    FACTORY    //////////////////////////////////////////////
// Registering the component
// see: http://wiki.sofa-framework.org/wiki/ObjectFactory
// 1-RegisterObject("description") + .add<> : Register the component
// 2-.add<>(true) : Set default template

int CosseratNeedleSlidingConstraintClass = RegisterObject("Simulate sliding contraints for needle insertion.")
.add< CosseratNeedleSlidingConstraint<sofa::defaulttype::Vec3Types> >(true)
;

} // namespace sofa


