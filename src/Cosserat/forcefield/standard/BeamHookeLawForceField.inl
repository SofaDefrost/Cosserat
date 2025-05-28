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
#pragma once
#include <Cosserat/forcefield/BeamHookeLawForceField.h>

namespace sofa::component::forcefield
{

template<typename DataTypes>
BeamHookeLawForceField<DataTypes>::BeamHookeLawForceField()
    : Inherit1()
{
    // Constructor initializes base class
}

template<typename DataTypes>
typename BeamHookeLawForceField<DataTypes>::Vector 
BeamHookeLawForceField<DataTypes>::getPosition(const Coord& coord) const
{
    if constexpr (DataTypes::spatial_dimensions > 3) {
        // For Vec6Types, position is the first three components
        return Vector(coord[0], coord[1], coord[2]);
    } else {
        // For Vec3Types, use the entire vector as position
        return coord;
    }
}

template<typename DataTypes>
typename BeamHookeLawForceField<DataTypes>::SO3Type 
BeamHookeLawForceField<DataTypes>::getRotation(const Coord& coord) const
{
    if constexpr (DataTypes::spatial_dimensions > 3) {
        // For Vec6Types, extract rotation from last three components
        // Using exponential map to convert rotation vector to SO3
        Vector3 rotVec(coord[3], coord[4], coord[5]);
        return SO3Type::exp(rotVec);
    } else {
        // For Vec3Types, return identity rotation as there's no rotation component
        return SO3Type::identity();
    }
}

template<typename DataTypes>
typename BeamHookeLawForceField<DataTypes>::Vector 
BeamHookeLawForceField<DataTypes>::getForce(const Deriv& deriv) const
{
    if constexpr (DataTypes::spatial_dimensions > 3) {
        // For Vec6Types, force is the first three components
        return Vector(deriv[0], deriv[1], deriv[2]);
    } else {
        // For Vec3Types, use the entire vector as force
        return deriv;
    }
}

template<typename DataTypes>
typename BeamHookeLawForceField<DataTypes>::Vector 
BeamHookeLawForceField<DataTypes>::getMoment(const Deriv& deriv) const
{
    if constexpr (DataTypes::spatial_dimensions > 3) {
        // For Vec6Types, moment is the last three components
        return Vector(deriv[3], deriv[4], deriv[5]);
    } else {
        // For Vec3Types, return zero moment as there's no moment component
        return Vector(0, 0, 0);
    }
}

template<typename DataTypes>
typename BeamHookeLawForceField<DataTypes>::Deriv 
BeamHookeLawForceField<DataTypes>::createDeriv(const Vector& force, const Vector& moment) const
{
    if constexpr (DataTypes::spatial_dimensions > 3) {
        // For Vec6Types, combine force and moment into a 6D vector
        return Deriv(force[0], force[1], force[2], moment[0], moment[1], moment[2]);
    } else {
        // For Vec3Types, only use force component
        return Deriv(force[0], force[1], force[2]);
    }
}

} // namespace sofa::component::forcefield
