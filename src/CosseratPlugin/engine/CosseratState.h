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
#include <SofaBaseMechanics/MechanicalObject.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/helper/vector.h>

namespace sofa::component::container
{

    using namespace core::behavior;
    using namespace core::objectmodel;
    using namespace sofa::defaulttype;
    using namespace sofa::core::topology;
    using sofa::defaulttype::Vector3;

    template <class DataTypes>
    class TemperatureState : public sofa::component::container::MechanicalObject<DataTypes>
    {
    public:
        SOFA_CLASS(SOFA_TEMPLATE(TemperatureState, DataTypes), SOFA_TEMPLATE(sofa::core::behavior::MechanicalState, DataTypes));

        typedef sofa::core::behavior::MechanicalState<DataTypes> Inherited;
        typedef typename DataTypes::VecCoord VecCoord;
        typedef typename DataTypes::VecDeriv VecDeriv;
        typedef typename DataTypes::Coord Coord;
        typedef typename DataTypes::Deriv Deriv;
        typedef typename Coord::value_type Real;

        /// assumes the mechanical object type (3D)
        typedef Vec<3, Real> Vec3;
        typedef StdVectorTypes<Vec3, Vec3, Real> MechanicalTypes;

        TemperatureState();

        /// compute the bounding box of the current state
        Data<bool> computeBoundingBox;
        void computeBBox(const core::ExecParams *params);
    };

} // namespace container
