/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture                          *
*                (c) 2006-2018 INRIA, USTL, UJF, CNRS, MGH                    *
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
*                           Plugin Cosserat    v1.0                           *
*				                                              *
* This plugin is also distributed under the GNU LGPL (Lesser General          *
* Public License) license with the same conditions than SOFA.                 *
*                                                                             *
* Contributors: Defrost team  (INRIA, University of Lille, CNRS,              *
*               Ecole Centrale de Lille)                                      *
*                                                                             *
* Contact information: https://project.inria.fr/softrobot/contact/            *
*                                                                             *
******************************************************************************/
#define SOFA_COSSERAT_CPP_CosseratInternalActuation
#include <Cosserat/forcefield/CosseratInternalActuation.inl>
#include <sofa/core/ObjectFactory.h>

namespace Cosserat
{

    using namespace sofa::defaulttype;
void registerCosseratInternalActuation(sofa::core::ObjectFactory *factory) {
  factory->registerObjects(
     sofa::core::ObjectRegistrationData(
          "This component is used to compute internal stress for torsion (along x) and bending (y and z)")
          .add<sofa::component::forcefield::CosseratInternalActuation<Vec3Types>>(true));
}

}
namespace sofa::component::forcefield
{
    using namespace sofa::defaulttype;
    template class SOFA_COSSERAT_API CosseratInternalActuation<defaulttype::Vec3Types>;
}
