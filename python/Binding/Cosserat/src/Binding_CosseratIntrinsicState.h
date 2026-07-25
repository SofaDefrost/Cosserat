#pragma once
// pybind11 binding declaration for CosseratIntrinsicState and PainlessBeamForceField.
// Activated in Module_Cosserat.cpp.

#include <pybind11/pybind11.h>

namespace sofapython3 {

/// Registers CosseratIntrinsicState in the Python module and in SOFA's PythonFactory.
void moduleAddCosseratIntrinsicState(pybind11::module& m);

/// Registers PainlessBeamForceField in the Python module and in SOFA's PythonFactory.
void moduleAddPainlessBeamForceField(pybind11::module& m);

}  // namespace sofapython3
