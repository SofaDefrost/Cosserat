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
#include "Binding_HookeSeratMapping.h"
#include <Cosserat/mapping/HookeSeratBaseMapping.h>
#include <Cosserat/mapping/HookeSeratDiscretMapping.h>
#include <SofaPython3/PythonFactory.h>
#include <SofaPython3/Sofa/Core/Binding_Base.h>
#include <SofaPython3/Sofa/Core/Binding_BaseContext.h>
#include <pybind11/stl.h>

namespace py {
using namespace pybind11;
}

using namespace sofa::core::objectmodel;
using namespace Cosserat::mapping;

namespace sofapython3 {

void moduleAddHookeSeratMapping(py::module &m) {

    using namespace sofa::defaulttype;

    py::class_<SectionInfo>(m, "SectionInfo")
        .def(py::init<>()) // Add other constructors as needed
        .def_property("length", &SectionInfo::getLength, &SectionInfo::setLength)
        .def("getStrainsVec", &SectionInfo::getStrainsVec, py::return_value_policy::reference_internal);

    py::class_<FrameInfo>(m, "FrameInfo")
        .def(py::init<>()) // Add other constructors as needed
        .def_property("length", &FrameInfo::getLength, &FrameInfo::setLength)
        .def_property("kappa", &FrameInfo::getKappa, &FrameInfo::setKappa)
        .def_property("transformation", &FrameInfo::getTransformation, &FrameInfo::setTransformation);


    // Explicit instantiation for Vec3Types
    using HookeSeratDiscretMapping3 = HookeSeratDiscretMapping<Vec3Types, Rigid3Types, Rigid3Types>;
    py::class_<HookeSeratDiscretMapping3, HookeSeratBaseMapping<Vec3Types, Rigid3Types, Rigid3Types>, py_shared_ptr<HookeSeratDiscretMapping3>> c3(m, "HookeSeratDiscretMapping3");

    PythonFactory::registerType<HookeSeratDiscretMapping3>(
        [](sofa::core::objectmodel::Base *object) {
            return py::cast(dynamic_cast<HookeSeratDiscretMapping3 *>(object));
        });


    // Explicit instantiation for Vec6Types
    using HookeSeratDiscretMapping6 = HookeSeratDiscretMapping<Vec6Types, Rigid3Types, Rigid3Types>;
    py::class_<HookeSeratDiscretMapping6, HookeSeratBaseMapping<Vec6Types, Rigid3Types, Rigid3Types>, py_shared_ptr<HookeSeratDiscretMapping6>> c6(m, "HookeSeratDiscretMapping6");

    PythonFactory::registerType<HookeSeratDiscretMapping6>(
        [](sofa::core::objectmodel::Base *object) {
            return py::cast(dynamic_cast<HookeSeratDiscretMapping6 *>(object));
        });
}

} // namespace sofapython3
