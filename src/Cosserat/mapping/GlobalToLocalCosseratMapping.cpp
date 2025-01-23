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
#define SOFA_COSSERAT_CPP_GlobalToLocalCosseratMapping
#include <Cosserat/mapping/GlobalToLocalCosseratMapping.inl>

#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/ObjectFactory.h>

namespace Cosserat::mapping
{
    using namespace sofa::defaulttype;
template class SOFA_COSSERAT_API GlobalToLocalCosseratMapping< sofa::defaulttype::Rigid3Types,
    sofa::defaulttype::Vec3Types >;
} // namespace sofa::component::mapping

namespace Cosserat
{
// Register in the Factory
void registerGlobalToLocalCosseratMapping(sofa::core::ObjectFactory* factory)
{
  factory->registerObjects(sofa::core::ObjectRegistrationData(
  "This component facilitates the creation of Cosserat Cables in SOFA simulations. It takes two mechanical"
  "objects as inputs: the rigid base of the beam (with 6 degrees of freedom) and the local coordinates of the beam. Using "
    "these inputs, the component computes and outputs the mechanical positions of the beam in global coordinates. "
    "Like any mapping, it updates the positions and velocities of the outputs based on the inputs. "
    "Additionally, forces applied to the outputs are propagated back to the inputs, ensuring bidirectional coupling.")
    .add< mapping::GlobalToLocalCosseratMapping< sofa::defaulttype::Rigid3Types, sofa::defaulttype::Vec3Types>>(true));
}
}
