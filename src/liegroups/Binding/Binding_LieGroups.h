/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2021 INRIA, USTL, UJF, CNRS, MGH                     *
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
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/

#pragma once

#include <memory>
#include <type_traits>
#include <array>
#include <tuple>
#include <random>
#include <iostream>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

namespace sofapython3 {

// Add SO2 class bindings to the module
void moduleAddSO2(pybind11::module &m);

// Add SO3 class bindings to the module
void moduleAddSO3(pybind11::module &m);

// Add SE2 class bindings to the module
void moduleAddSE2(pybind11::module &m);

// Add SE3 class bindings to the module
void moduleAddSE3(pybind11::module &m);

// Add SGal3 class bindings to the module
void moduleAddSGal3(pybind11::module &m);

// Add SE23 class bindings to the module
void moduleAddSE23(pybind11::module &m);

// Add Bundle class bindings to the module
void moduleAddBundle(pybind11::module &m);

// Add Utility functions for Lie groups
void moduleAddLieGroupUtils(pybind11::module &m);

// Add all Lie group bindings to the module
void moduleAddLieGroups(pybind11::module &m);

} // namespace sofapython3

