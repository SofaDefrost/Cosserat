// pybind11 bindings for CosseratIntrinsicState and PainlessBeamForceField.
//
// Exposed API
// -----------
//
// CosseratIntrinsicState
//   .getPositions()              → numpy (N+1, 3)  — node positions
//   .setPositions(arr)           ← numpy (N+1, 3)
//   .getOrientations()           → list[SO3]       — segment orientations
//   .setOrientations(list[SO3])  ← list of SO3 objects
//   .getOrientationsQuat()       → numpy (N, 4)    — quaternions [x,y,z,w]
//   .setOrientationsQuat(arr)    ← numpy (N, 4)    — quaternions [x,y,z,w]
//   .getRestLengths()            → list[float]
//   .getLinearStrain(i)          → numpy 3
//   .getAngularStrain(i)         → numpy 3
//   .getNbSegments()             → int
//
// PainlessBeamForceField
//   .computeDForcesFromData(kFactor=1.0)   — fills d_df_positions / d_df_angles
//   .getDfPositions()            → numpy (N+1, 3)
//   .getDfAngles()               → numpy (N, 3)

/******************************************************************************
 *                 SOFA, Simulation Open-Framework Architecture                *
 * (c) 2021 INRIA, USTL, UJF, CNRS, MGH                                       *
 * LGPL v2.1 or later                                                          *
 ******************************************************************************/

#include "Binding_CosseratIntrinsicState.h"

#include <Cosserat/engine/CosseratIntrinsicState.h>
#include <Cosserat/engine/PainlessBeamForceField.h>

#include <SofaPython3/PythonFactory.h>
#include <SofaPython3/Sofa/Core/Binding_Base.h>

#include <pybind11/eigen.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using CIS  = sofa::component::cosserat::engine::CosseratIntrinsicState;
using PBFF = sofa::component::cosserat::engine::PainlessBeamForceField;
using SO3  = sofa::component::cosserat::liegroups::SO3<double>;

namespace sofapython3 {

// ── helpers ───────────────────────────────────────────────────────────────────

static Eigen::MatrixXd vec3_to_matrix(
        const sofa::type::vector<sofa::type::Vec3d>& v) {
    Eigen::MatrixXd m(static_cast<Eigen::Index>(v.size()), 3);
    for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(v.size()); ++i) {
        m(i, 0) = v[i][0];
        m(i, 1) = v[i][1];
        m(i, 2) = v[i][2];
    }
    return m;
}

static sofa::type::vector<sofa::type::Vec3d> matrix_to_vec3(
        const Eigen::MatrixXd& m) {
    sofa::type::vector<sofa::type::Vec3d> v(static_cast<size_t>(m.rows()));
    for (Eigen::Index i = 0; i < m.rows(); ++i)
        v[static_cast<size_t>(i)] = sofa::type::Vec3d(m(i,0), m(i,1), m(i,2));
    return v;
}

// ── CosseratIntrinsicState binding ───────────────────────────────────────────

void moduleAddCosseratIntrinsicState(py::module& m) {

    py::class_<CIS, sofa::core::objectmodel::BaseObject,
               sofapython3::py_shared_ptr<CIS>>(m, "CosseratIntrinsicState",
        R"doc(
Staggered Cosserat beam mechanical state.

Stores N+1 node positions (Vec3d) and N segment orientations (SO3).

Python API
----------
getPositions()                  → numpy (N+1, 3)
setPositions(arr)               ← numpy (N+1, 3)
getOrientations()               → list[SO3]
setOrientations(list[SO3])
getOrientationsQuat()           → numpy (N, 4)  [qx,qy,qz,qw]
setOrientationsQuat(arr)        ← numpy (N, 4)  [qx,qy,qz,qw]
getRestLengths()                → list[float]
getLinearStrain(i)              → numpy 3
getAngularStrain(i)             → numpy 3
getNbSegments()                 → int
)doc")

        // ── Positions ──────────────────────────────────────────────────────────
        .def("getPositions",
             [](CIS& self) -> Eigen::MatrixXd {
                 return vec3_to_matrix(self.getPositions());
             },
             "Node positions as numpy array of shape (N+1, 3).")

        .def("setPositions",
             [](CIS& self, const Eigen::MatrixXd& arr) {
                 if (arr.cols() != 3) throw py::value_error("Expected shape (N+1, 3).");
                 auto& data = *self.d_positions.beginEdit();
                 data.resize(static_cast<size_t>(arr.rows()));
                 for (Eigen::Index i = 0; i < arr.rows(); ++i)
                     data[static_cast<size_t>(i)] =
                         sofa::type::Vec3d(arr(i,0), arr(i,1), arr(i,2));
                 self.d_positions.endEdit();
             },
             py::arg("positions"),
             "Set node positions from numpy array of shape (N+1, 3).")

        // ── Orientations (as SO3 list) ─────────────────────────────────────────
        .def("getOrientations",
             [](CIS& self) -> std::vector<SO3> {
                 const auto& R = self.getOrientations();
                 return std::vector<SO3>(R.begin(), R.end());
             },
             "Segment orientations as a list of SO3 objects (length N).")

        .def("setOrientations",
             [](CIS& self, const std::vector<SO3>& R_list) {
                 auto& data = *self.d_orientations.beginEdit();
                 data.resize(R_list.size());
                 for (size_t i = 0; i < R_list.size(); ++i)
                     data[i] = R_list[i];
                 self.d_orientations.endEdit();
             },
             py::arg("orientations"),
             "Set segment orientations from a list of SO3 objects.")

        // ── Orientations (as quaternion array) ────────────────────────────────
        .def("getOrientationsQuat",
             [](CIS& self) -> Eigen::MatrixXd {
                 const auto& R = self.getOrientations();
                 Eigen::MatrixXd q(static_cast<Eigen::Index>(R.size()), 4);
                 for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(R.size()); ++i) {
                     const auto& ei_q = R[static_cast<size_t>(i)].quaternion();
                     q(i, 0) = ei_q.x();
                     q(i, 1) = ei_q.y();
                     q(i, 2) = ei_q.z();
                     q(i, 3) = ei_q.w();
                 }
                 return q;
             },
             "Segment orientations as numpy (N, 4) array [qx, qy, qz, qw].")

        .def("setOrientationsQuat",
             [](CIS& self, const Eigen::MatrixXd& q) {
                 if (q.cols() != 4) throw py::value_error("Expected shape (N, 4).");
                 auto& data = *self.d_orientations.beginEdit();
                 data.resize(static_cast<size_t>(q.rows()));
                 for (Eigen::Index i = 0; i < q.rows(); ++i) {
                     Eigen::Quaterniond eq(q(i,3), q(i,0), q(i,1), q(i,2));
                     data[static_cast<size_t>(i)] = SO3(eq.normalized());
                 }
                 self.d_orientations.endEdit();
             },
             py::arg("quats"),
             "Set segment orientations from numpy (N, 4) array [qx, qy, qz, qw].")

        // ── Rest lengths ───────────────────────────────────────────────────────
        .def("getRestLengths",
             [](CIS& self) -> std::vector<double> {
                 return self.getRestLengths();
             },
             "Rest lengths h_i of the N segments (computed during init).")

        // ── Strain accessors ───────────────────────────────────────────────────
        .def("getLinearStrain",
             [](CIS& self, size_t i) -> Eigen::Vector3d {
                 const auto v = self.getLinearStrain(i);
                 return Eigen::Vector3d(v[0], v[1], v[2]);
             },
             py::arg("i"),
             "Linear strain Γ_i = R_i^{-1}(x_{i+1}-x_i)/h_i - e1 for segment i.")

        .def("getAngularStrain",
             [](CIS& self, size_t i) -> Eigen::Vector3d {
                 const auto v = self.getAngularStrain(i);
                 return Eigen::Vector3d(v[0], v[1], v[2]);
             },
             py::arg("i"),
             "Angular strain Ω_i = log(R_{i-1}^T R_i)/h_dual at node i.")

        // ── Misc ───────────────────────────────────────────────────────────────
        .def("getNbSegments", &CIS::getNbSegments,
             "Number of segments N (= number of SO3 DOFs).");

    // Register for automatic downcasting in SofaPython3
    sofapython3::PythonFactory::registerType<CIS>(
        [](sofa::core::objectmodel::Base* obj) {
            return py::cast(dynamic_cast<CIS*>(obj));
        });
}

// ── PainlessBeamForceField binding ────────────────────────────────────────────

void moduleAddPainlessBeamForceField(py::module& m) {

    py::class_<PBFF, sofa::core::objectmodel::BaseObject,
               sofapython3::py_shared_ptr<PBFF>>(m, "PainlessBeamForceField",
        R"doc(
Staggered Cosserat beam elastic force field (Painless discretisation).

Python integration API
-----------------------
After writing d_dx_positions / d_dx_angles Data fields (or using the numpy
setters below), call computeDForcesFromData() to fill the df outputs.

computeDForcesFromData(kFactor=1.0)
getDfPositions()   → numpy (N+1, 3)
getDfAngles()      → numpy (N, 3)
getNodalForces()   → numpy (N+1, 3)   (filled by addForce / the animation loop)
getSegmentTorques()→ numpy (N, 3)
)doc")

        // ── Direct differential force computation (Python Euler integrator) ────
        .def("computeDForcesFromData",
             &PBFF::computeDForcesFromData,
             py::arg("kFactor") = 1.0,
             "Compute df from d_dx_positions/d_dx_angles and write to d_df_* Data fields.")

        // ── Set displacement inputs ────────────────────────────────────────────
        .def("setDxPositions",
             [](PBFF& self, const Eigen::MatrixXd& arr) {
                 if (arr.cols() != 3) throw py::value_error("Expected shape (N+1, 3).");
                 auto& data = *self.d_dx_positions.beginEdit();
                 data = matrix_to_vec3(arr);
                 self.d_dx_positions.endEdit();
             },
             py::arg("dx"),
             "Set position displacement input (N+1, 3) for computeDForcesFromData.")

        .def("setDxAngles",
             [](PBFF& self, const Eigen::MatrixXd& arr) {
                 if (arr.cols() != 3) throw py::value_error("Expected shape (N, 3).");
                 auto& data = *self.d_dx_angles.beginEdit();
                 data = matrix_to_vec3(arr);
                 self.d_dx_angles.endEdit();
             },
             py::arg("dw"),
             "Set angular displacement input (N, 3) for computeDForcesFromData.")

        // ── Read force outputs ─────────────────────────────────────────────────
        .def("getNodalForces",
             [](PBFF& self) -> Eigen::MatrixXd {
                 return vec3_to_matrix(self.d_nodalForces.getValue());
             },
             "Elastic forces on N+1 nodes (world frame), filled by addForce.")

        .def("getSegmentTorques",
             [](PBFF& self) -> Eigen::MatrixXd {
                 return vec3_to_matrix(self.d_segmentTorques.getValue());
             },
             "Elastic torques on N segments (body frame), filled by addForce.")

        .def("getDfPositions",
             [](PBFF& self) -> Eigen::MatrixXd {
                 return vec3_to_matrix(self.d_df_positions.getValue());
             },
             "Differential forces on N+1 nodes, filled by computeDForcesFromData.")

        .def("getDfAngles",
             [](PBFF& self) -> Eigen::MatrixXd {
                 return vec3_to_matrix(self.d_df_angles.getValue());
             },
             "Differential torques on N segments, filled by computeDForcesFromData.");

    sofapython3::PythonFactory::registerType<PBFF>(
        [](sofa::core::objectmodel::Base* obj) {
            return py::cast(dynamic_cast<PBFF*>(obj));
        });
}

}  // namespace sofapython3
