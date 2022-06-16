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
#ifndef SOFA_COMPONENT_CONTAINER_TEMPERATURESTATE_INL
#define SOFA_COMPONENT_CONTAINER_TEMPERATURESTATE_INL

#include "TemperatureState.h"


namespace sofa
{

namespace component
{

namespace container
{
using namespace std;
using namespace sofa::core;
using namespace sofa::core::topology;
using namespace sofa::defaulttype;

template <class DataTypes>
TemperatureState<DataTypes>::TemperatureState()
    : computeBoundingBox ( initData(&computeBoundingBox, (bool) true, "computeBoundingBox","Boolean to activate the computation of BoundingBox") )
{
}


template <class DataTypes>
void TemperatureState<DataTypes>::computeBBox(const core::ExecParams* params)
{
    if ( computeBoundingBox.getValue())
    {
         std::cout << "Compute BBox" << std::endl;
         Inherited::computeBBox(params);
    }
}



} // namespace container

} // namespace component

} // namespace sofa

#endif //SOFA_COMPONENT_CONTAINER_TEMPERATURESTATE_INL
