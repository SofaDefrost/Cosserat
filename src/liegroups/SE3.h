// This file contains a forward declaration for the SE3 (Special Euclidean group
// in 3D) class.

/******************************************************************************
 *                 SOFA, Simulation Open-Framework Architecture * (c) 2006
 *INRIA, USTL, UJF, CNRS, MGH                     *
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
 * along with this program. If not, see <http://www.gnu.org/licenses/\>. *
 ******************************************************************************/

// #ifndef SOFA_COMPONENT_COSSERAT_LIEGROUPS_SE3_H
// #define SOFA_COMPONENT_COSSERAT_LIEGROUPS_SE3_H
#pragma once

#include "LieGroupBase.h" // Then the base class interface
#include "SO3.h"          // Then other dependencies
#include "Types.h"        // Then our type system
#include <Eigen/Geometry> // Include Eigen first

// Forward declaration outside the namespace
namespace sofa::component::cosserat::liegroups {
template <typename Scalar> class SE3;
}

namespace sofa::component::cosserat::liegroups {