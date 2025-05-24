// pybind11 bindings for the Cosserat Lie group library.
// Exposes SO2 and SO3 (double precision) to SofaPython3.
//
// Author  : Younes Adagolodjo (DEFROST / INRIA / Polytech Lille)
// Branch  : painless/base-geometry

/******************************************************************************
 *                 SOFA, Simulation Open-Framework Architecture                *
 *                    (c) 2021 INRIA, USTL, UJF, CNRS, MGH                     *
 *                                                                             *
 * This library is free software; you can redistribute it and/or modify it    *
 * under the terms of the GNU Lesser General Public License as published by   *
 * the Free Software Foundation; either version 2.1 of the License, or (at    *
 * your option) any later version.                                             *
 ******************************************************************************/

#include <array>
#include <random>
#include <sstream>

#include "Binding_LieGroups.h"
#include <liegroups/LieGroupBase.h>
#include <liegroups/LieGroupBase.inl>
#include <liegroups/SO2.h>
#include <liegroups/SO3.h>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>    // Eigen <-> numpy automatic conversions
#include <pybind11/stl.h>      // std::vector <-> Python list

namespace py = pybind11;
namespace lg = sofa::component::cosserat::liegroups;

namespace sofapython3 {

// ── SO2 ───────────────────────────────────────────────────────────────────────

void moduleAddSO2(py::module& m) {
    py::class_<lg::SO2<double>>(m, "SO2",
        "2D rotation group SO(2).  Internally stored as a unit complex number.")
        .def(py::init<>(), "Identity rotation.")
        .def(py::init<double>(), "Construct from angle [rad].", py::arg("angle"))
        .def("__mul__", &lg::SO2<double>::operator*, "Group composition.")
        .def("inverse", &lg::SO2<double>::inverse,  "Inverse rotation.")
        .def("matrix",  &lg::SO2<double>::matrix,   "2×2 rotation matrix (numpy).")
        .def("angle",   &lg::SO2<double>::angle,    "Rotation angle [rad].")
        .def("setAngle",&lg::SO2<double>::setAngle, py::arg("angle"))
        .def("complex", &lg::SO2<double>::complex,  "Complex (cos θ, sin θ) as 2-vector.")
        .def("direction",    &lg::SO2<double>::direction)
        .def("perpendicular",&lg::SO2<double>::perpendicular)
        .def_static("exp", &lg::SO2<double>::exp,  py::arg("omega"))
        .def("log",        &lg::SO2<double>::log)
        .def("adjoint",    &lg::SO2<double>::adjoint)
        .def("isApprox",   &lg::SO2<double>::isApprox,
             py::arg("other"), py::arg("eps") = 1e-12)
        .def_static("identity", []() { return lg::SO2<double>::Identity(); })
        .def_static("hat", &lg::SO2<double>::hat, py::arg("omega"))
        .def_static("vee", &lg::SO2<double>::vee, py::arg("matrix"))
        .def("act", [](const lg::SO2<double>& self, const Eigen::Vector2d& p) {
                return self.act(p);
             }, py::arg("point"))
        .def("__repr__", [](const lg::SO2<double>& self) {
                std::ostringstream ss;
                ss << "SO2(angle=" << self.angle() << " rad)";
                return ss.str();
             })
        .def_static("random", [](py::object seed) {
                static std::random_device rd;
                static std::mt19937 gen(rd());
                if (!seed.is_none()) gen.seed(seed.cast<unsigned>());
                return lg::SO2<double>::computeRandom(gen);
             }, py::arg("seed") = py::none());
}

// ── SO3 ───────────────────────────────────────────────────────────────────────

void moduleAddSO3(py::module& m) {
    py::class_<lg::SO3<double>>(m, "SO3",
        R"doc(
3D rotation group SO(3).  Internally stored as a unit quaternion.

Constructors
------------
SO3()                              # identity
SO3(angle, axis)                   # angle [rad] + axis (auto-normalised)
SO3(Quaterniond)                   # Eigen quaternion (also accepts numpy [x,y,z,w])
SO3.from_quat_xyzw(x, y, z, w)    # explicit SOFA-convention components
SO3.from_quat_array([x,y,z,w])    # from numpy array

Key methods
-----------
R * R2                → SO3      composition
R.inverse()           → SO3
R.act(v)              → R³       rotate 3-vector (numpy in, numpy out)
R.log()               → R³       logarithm ∈ so(3)  [rad]
SO3.exp(omega)        → SO3      exp map (omega ∈ R³)
R.toRotationMatrix()  → 3×3 numpy
R.quaternion()        → 4-numpy  [qx, qy, qz, qw]  (SOFA convention)
)doc")

        // ── Constructors ────────────────────────────────────────────────────────
        .def(py::init<>(), "Identity rotation.")
        .def(py::init<double, const Eigen::Vector3d&>(),
             py::arg("angle"), py::arg("axis"),
             "Construct from angle [rad] and axis (normalised internally).")
        .def(py::init<const Eigen::Quaterniond&>(),
             py::arg("quat"),
             "Construct from Eigen Quaternion (normalised).")
        .def(py::init<const Eigen::Matrix3d&>(),
             py::arg("matrix"),
             "Construct from 3×3 rotation matrix.")

        // ── From Python-friendly quaternion representations ─────────────────────
        .def_static("from_quat_xyzw",
             [](double qx, double qy, double qz, double qw) {
                 return lg::SO3<double>(
                     Eigen::Quaterniond(qw, qx, qy, qz).normalized());
             },
             py::arg("qx"), py::arg("qy"), py::arg("qz"), py::arg("qw"),
             "Construct from quaternion components in SOFA order (x, y, z, w).")
        .def_static("from_quat_array",
             [](const Eigen::Vector4d& q) {
                 return lg::SO3<double>(
                     Eigen::Quaterniond(q[3], q[0], q[1], q[2]).normalized());
             },
             py::arg("q_xyzw"),
             "Construct from numpy array [qx, qy, qz, qw].")

        // ── Group operations ────────────────────────────────────────────────────
        .def("__mul__", &lg::SO3<double>::operator*, "Group composition R * R2.")
        .def("inverse", &lg::SO3<double>::inverse,   "Inverse rotation R^{-1}.")

        // ── Action on vectors ───────────────────────────────────────────────────
        .def("act",
             [](const lg::SO3<double>& self, const Eigen::Vector3d& v) {
                 return self.act(v);
             },
             py::arg("v"), "Rotate a 3-vector: R·v  (numpy in, numpy out).")

        // ── Lie group maps ──────────────────────────────────────────────────────
        .def_static("exp",
             [](const Eigen::Vector3d& omega) {
                 return lg::SO3<double>::exp(omega);
             },
             py::arg("omega"),
             "Exponential map so(3) → SO(3).  omega is a 3-vector [rad].")
        .def("log",
             [](const lg::SO3<double>& self) -> Eigen::Vector3d {
                 return self.log();
             },
             "Logarithm map SO(3) → so(3).  Returns a 3-vector [rad].")
        .def("adjoint", &lg::SO3<double>::adjoint,
             "Adjoint = rotation matrix for SO(3).")

        // ── Conversions ─────────────────────────────────────────────────────────
        .def("toRotationMatrix",
             [](const lg::SO3<double>& self) -> Eigen::Matrix3d {
                 return self.toRotationMatrix();
             },
             "3×3 rotation matrix as numpy array (row-major).")
        .def("quaternion",
             [](const lg::SO3<double>& self) -> Eigen::Vector4d {
                 const auto& q = self.quaternion();
                 return Eigen::Vector4d(q.x(), q.y(), q.z(), q.w());
             },
             "Unit quaternion as numpy [qx, qy, qz, qw] (SOFA convention, w last).")

        // ── Utilities ───────────────────────────────────────────────────────────
        .def("isApprox", &lg::SO3<double>::isApprox,
             py::arg("other"), py::arg("eps") = 1e-12)
        .def_static("identity",
             []() { return lg::SO3<double>(); },
             "Identity rotation.")
        .def_static("hat", &lg::SO3<double>::hat,
             py::arg("omega"),
             "Hat operator: ω ∈ R³  →  skew-symmetric 3×3.")
        .def_static("vee", &lg::SO3<double>::vee,
             py::arg("matrix"),
             "Vee operator: skew-symmetric 3×3  →  ω ∈ R³.")
        .def("__repr__",
             [](const lg::SO3<double>& self) {
                 const auto& q = self.quaternion();
                 std::ostringstream ss;
                 ss << "SO3(quat=[" << q.x() << ", " << q.y() << ", "
                    << q.z() << ", " << q.w() << "])";
                 return ss.str();
             })
        .def_static("random",
             [](py::object seed) {
                 static std::random_device rd;
                 static std::mt19937 gen(rd());
                 if (!seed.is_none()) gen.seed(seed.cast<unsigned>());
                 return lg::SO3<double>::computeRandom(gen);
             },
             py::arg("seed") = py::none(),
             "Uniformly random SO(3) element (Shoemake method).");
}

// ── moduleAddLieGroups ────────────────────────────────────────────────────────

void moduleAddLieGroups(py::module& m) {
    moduleAddSO2(m);
    moduleAddSO3(m);
}

}  // namespace sofapython3
