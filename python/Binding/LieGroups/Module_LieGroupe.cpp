// This file defines the main pybind11 module for the Cosserat plugin, including bindings for Lie groups and PointsManager.

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

#include <memory>
#include <functional>
#include <type_traits>
#include <iterator>

// Fix namespace issues with a wrapper to provide std namespace prefixes
namespace std {
    using ::std::allocator;
    using ::std::forward_iterator_tag;
    using ::std::pointer_traits;
    using ::std::__remove_cv_t;
    using ::std::__rebind_pointer_t;
    using ::std::__conditional_t;
    using ::std::is_same;
    using ::std::is_pointer;
    using ::std::__forward_list_node;
    using ::std::__begin_node_of;
    using ::std::__forward_begin_node;
}

#include <pybind11/pybind11.h>
#include "Binding_LieGroups.h"


namespace py { using namespace pybind11; }

namespace sofapython3
{

PYBIND11_MODULE(Cosserat, m)
{
    // Only add Lie groups related functionality
    moduleAddLieGroups(m);
}

} // namespace sofapython3