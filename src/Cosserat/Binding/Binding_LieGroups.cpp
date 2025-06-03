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

#include "Binding_LieGroups.h"
#include <Cosserat/liegroups/SO2.h>
#include <Cosserat/liegroups/SO3.h>
#include <Cosserat/liegroups/SE2.h>
#include <Cosserat/liegroups/SE3.h>
#include <Cosserat/liegroups/SGal3.h>
#include <Cosserat/liegroups/Bundle.h>
#include <Cosserat/liegroups/Utils.h>
#include <pybind11/eigen.h>
#include <pybind11/operators.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>

using namespace sofa::component::cosserat::liegroups;

namespace py = pybind11;

namespace sofapython3 {

template <typename T>
void addCommonMethods(py::class_<T>& cls) {
    cls.def("inverse", &T::computeInverse, "Returns the inverse element");
    cls.def("log", &T::computeLog, "Returns the logarithm mapping to the Lie algebra");
    cls.def("adjoint", &T::computeAdjoint, "Returns the adjoint representation");
    cls.def("act", py::overload_cast<const typename T::Vector&>(&T::computeAction, py::const_), "Apply group action on a vector");
    cls.def("matrix", &T::matrix, "Returns the matrix representation");
    cls.def("__mul__", &T::compose, "Compose two group elements");
    cls.def("isApprox", &T::computeIsApprox, "Check if approximately equal to another element",
            py::arg("other"), py::arg("eps") = 1e-10);
}

void moduleAddSO2(py::module& m) {
    using Scalar = double; // Use double as the default scalar type
    using SO2d = SO2<Scalar>;
    using Vector = typename SO2d::Vector;
    using TangentVector = typename SO2d::TangentVector;

    py::class_<SO2d> so2(m, "SO2", "Special Orthogonal group in 2D - rotations in the plane");
    
    // Constructors
    so2.def(py::init<>(), "Default constructor - identity rotation");
    so2.def(py::init<Scalar>(), "Construct from angle (in radians)", py::arg("angle"));
    
    // Static methods
    so2.def_static("Identity", &SO2d::computeIdentity, "Returns the identity element");
    so2.def_static("exp", &SO2d::computeExp, "Exponential map from Lie algebra", py::arg("algebra_element"));
    
    // Properties
    so2.def_property_readonly("angle", &SO2d::angle, "The rotation angle in radians");
    
    // Common methods
    addCommonMethods<SO2d>(so2);
}

void moduleAddSO3(py::module& m) {
    using Scalar = double; // Use double as the default scalar type
    using SO3d = SO3<Scalar>;
    using Vector = typename SO3d::Vector;
    using Matrix = typename SO3d::Matrix;
    using TangentVector = typename SO3d::TangentVector;
    using Quaternion = typename SO3d::Quaternion;
    
    py::class_<SO3d> so3(m, "SO3", "Special Orthogonal group in 3D - rotations in 3D space");
    
    // Constructors
    so3.def(py::init<>(), "Default constructor - identity rotation");
    so3.def(py::init<const Scalar&, const Vector&>(), "Construct from angle-axis representation",
            py::arg("angle"), py::arg("axis"));
    so3.def(py::init<const Quaternion&>(), "Construct from quaternion", py::arg("quat"));
    so3.def(py::init<const Matrix&>(), "Construct from rotation matrix", py::arg("mat"));
    
    // Static methods
    so3.def_static("Identity", &SO3d::computeIdentity, "Returns the identity element");
    so3.def_static("exp", &SO3d::computeExp, "Exponential map from Lie algebra", py::arg("omega"));
    so3.def_static("hat", &SO3d::hat, "Convert vector to skew-symmetric matrix", py::arg("v"));
    so3.def_static("vee", &SO3d::vee, "Convert skew-symmetric matrix to vector", py::arg("Omega"));
    
    // Properties
    so3.def_property_readonly("quaternion", &SO3d::quaternion, "The quaternion representation");
    
    // Common methods
    addCommonMethods<SO3d>(so3);
    
    // Additional methods
    so3.def("angleAxis", &SO3d::angleAxis, "Get the angle-axis representation");
}

void moduleAddSE2(py::module& m) {
    using Scalar = double; // Use double as the default scalar type
    using SE2d = SE2<Scalar>;
    using SO2d = SO2<Scalar>;
    using Vector = typename SE2d::Vector;
    using Vector2 = typename SE2d::Vector2;
    using TangentVector = typename SE2d::TangentVector;
    
    py::class_<SE2d> se2(m, "SE2", "Special Euclidean group in 2D - rigid transformations in the plane");
    
    // Constructors
    se2.def(py::init<>(), "Default constructor - identity transformation");
    se2.def(py::init<const SO2d&, const Vector2&>(), "Construct from rotation and translation",
            py::arg("rotation"), py::arg("translation"));
    se2.def(py::init<const Scalar&, const Vector2&>(), "Construct from angle and translation",
            py::arg("angle"), py::arg("translation"));
    
    // Static methods
    se2.def_static("Identity", &SE2d::identity, "Returns the identity element");
    se2.def_static("exp", &SE2d::exp, "Exponential map from Lie algebra", py::arg("algebra_element"));
    
    // Properties
    se2.def_property_readonly("rotation", &SE2d::rotation, "The rotation component");
    se2.def_property_readonly("translation", &SE2d::translation, "The translation component");
    
    // Common methods
    addCommonMethods<SE2d>(se2);
}

void moduleAddSE3(py::module& m) {
    using Scalar = double; // Use double as the default scalar type
    using SE3d = SE3<Scalar>;
    using SO3d = SO3<Scalar>;
    using Vector = typename SE3d::Vector;
    using Vector3 = typename SE3d::Vector3;
    using Matrix4 = typename SE3d::Matrix4;
    using TangentVector = typename SE3d::TangentVector;
    
    py::class_<SE3d> se3(m, "SE3", "Special Euclidean group in 3D - rigid transformations in 3D space");
    
    // Constructors
    se3.def(py::init<>(), "Default constructor - identity transformation");
    se3.def(py::init<const SO3d&, const Vector3&>(), "Construct from rotation and translation",
            py::arg("rotation"), py::arg("translation"));
    se3.def(py::init<const Matrix4&>(), "Construct from homogeneous transformation matrix",
            py::arg("T"));
    
    // Static methods
    se3.def_static("Identity", &SE3d::computeIdentity, "Returns the identity element");
    se3.def_static("exp", &SE3d::computeExp, "Exponential map from Lie algebra", py::arg("twist"));
    se3.def_static("BCH", &SE3d::BCH, "Baker-Campbell-Hausdorff formula for SE(3)",
                   py::arg("X"), py::arg("Y"));
    
    // Properties
    se3.def_property_readonly("rotation", py::overload_cast<>(&SE3d::rotation, py::const_), "The rotation component");
    se3.def_property_readonly("translation", py::overload_cast<>(&SE3d::translation, py::const_), "The translation component");
    
    // Common methods
    addCommonMethods<SE3d>(se3);
}

void moduleAddSGal3(py::module& m) {
    using Scalar = double; // Use double as the default scalar type
    using SGal3d = SGal3<Scalar>;
    using SE3d = SE3<Scalar>;
    using SO3d = SO3<Scalar>;
    using Vector = typename SGal3d::Vector;
    using Vector3 = typename SGal3d::Vector3;
    using Matrix = typename SGal3d::Matrix;
    using TangentVector = typename SGal3d::TangentVector;
    
    py::class_<SGal3d> sgal3(m, "SGal3", "Special Galilean group in 3D - spacetime transformations including velocity");
    
    // Constructors
    sgal3.def(py::init<>(), "Default constructor - identity transformation");
    sgal3.def(py::init<const SE3d&, const Vector3&, Scalar>(), "Construct from SE3, velocity, and time",
             py::arg("pose"), py::arg("velocity"), py::arg("time") = 0.0);
    
    // Factory functions
    m.def("fromComponents", &fromComponents<Scalar>, 
          "Create SGal3 transformation from position, rotation, velocity, and time",
          py::arg("position"), py::arg("rotation"), py::arg("velocity"), py::arg("time") = 0.0);
    
    m.def("fromPositionEulerVelocityTime", &fromPositionEulerVelocityTime<Scalar>,
          "Create SGal3 transformation from position, Euler angles (roll, pitch, yaw), velocity, and time",
          py::arg("position"), py::arg("roll"), py::arg("pitch"), py::arg("yaw"), 
          py::arg("velocity"), py::arg("time") = 0.0);
    
    m.def("toPositionEulerVelocityTime", &toPositionEulerVelocityTime<Scalar>,
          "Convert SGal3 transformation to position, Euler angles, velocity, and time vector",
          py::arg("transform"));
    
    m.def("interpolate", &interpolate<Scalar>,
          "Interpolate between two Galilean transformations",
          py::arg("from"), py::arg("to"), py::arg("t"));
    
    // Properties
    sgal3.def_property_readonly("pose", py::overload_cast<>(&SGal3d::pose, py::const_), 
                              "The SE3 component (pose)");
    sgal3.def_property_readonly("velocity", py::overload_cast<>(&SGal3d::velocity, py::const_), 
                              "The linear velocity component");
    sgal3.def_property_readonly("time", &SGal3d::time, 
                              "The time component");
    
    // Common methods
    addCommonMethods<SGal3d>(sgal3);
    
    // Additional methods
    sgal3.def("transform", &SGal3d::transform, 
             "Transform a point-velocity-time tuple",
             py::arg("point_vel_time"));
    
    // Static methods
    sgal3.def_static("Identity", &SGal3d::Identity, "Returns the identity element");
    sgal3.def_static("exp", &SGal3d::computeExp, "Exponential map from Lie algebra", 
                    py::arg("algebra_element"));
    sgal3.def_static("BCH", &SGal3d::BCH, "Baker-Campbell-Hausdorff formula for SGal3",
                    py::arg("X"), py::arg("Y"));
    
    // Operators
    sgal3.def("__truediv__", [](const SGal3d& g, const Vector& point_vel_time) {
        return point_vel_time / g;
    }, "Apply inverse transformation to a point-velocity-time tuple", 
       py::arg("point_vel_time"));
}

void moduleAddSE23(py::module& m) {
    using Scalar = double; // Use double as the default scalar type
    using SE23d = SE23<Scalar>;
    using SE3d = SE3<Scalar>;
    using SO3d = SO3<Scalar>;
    using Vector = typename SE23d::Vector;
    using Vector3 = typename SE23d::Vector3;
    using TangentVector = typename SE23d::TangentVector;
    
    py::class_<SE23d> se23(m, "SE23", "Special Euclidean group in (2+1)D - rigid transformations with velocity");
    
    // Constructors
    se23.def(py::init<>(), "Default constructor - identity transformation with zero velocity");
    se23.def(py::init<const SE3d&, const Vector3&>(), "Construct from SE3 pose and velocity",
             py::arg("pose"), py::arg("velocity"));
    se23.def(py::init<const SO3d&, const Vector3&, const Vector3&>(), "Construct from rotation, position, and velocity",
             py::arg("rotation"), py::arg("position"), py::arg("velocity"));
    
    // Factory functions
    m.def("fromComponents", &fromComponents<Scalar>, 
          "Create SE23 transformation from position, rotation, and velocity",
          py::arg("position"), py::arg("rotation"), py::arg("velocity"));
    
    m.def("fromPositionEulerVelocity", &fromPositionEulerVelocity<Scalar>,
          "Create SE23 transformation from position, Euler angles (roll, pitch, yaw), and velocity",
          py::arg("position"), py::arg("roll"), py::arg("pitch"), py::arg("yaw"), 
          py::arg("velocity"));
    
    m.def("toPositionEulerVelocity", &toPositionEulerVelocity<Scalar>,
          "Convert SE23 transformation to position, Euler angles, and velocity vector",
          py::arg("transform"));
    
    m.def("interpolate", &interpolate<Scalar>,
          "Interpolate between two SE23 transformations",
          py::arg("from"), py::arg("to"), py::arg("t"));
    
    // Properties
    se23.def_property_readonly("pose", py::overload_cast<>(&SE23d::pose, py::const_), 
                              "The SE3 component (pose)");
    se23.def_property_readonly("velocity", py::overload_cast<>(&SE23d::velocity, py::const_), 
                              "The linear velocity component");
    
    // Common methods
    addCommonMethods<SE23d>(se23);
    
    // Transform a point-velocity pair
    se23.def("transform", [](const SE23d& self, const Vector& point_vel) {
        return self.act(point_vel);
    }, "Transform a point-velocity pair", py::arg("point_vel"));
    
    // Static methods
    se23.def_static("Identity", &SE23d::identity, "Returns the identity element");
    se23.def_static("exp", &SE23d::exp, "Exponential map from Lie algebra", 
                    py::arg("algebra_element"));
    se23.def_static("BCH", &SE23d::BCH, "Baker-Campbell-Hausdorff formula for SE23",
                    py::arg("X"), py::arg("Y"));
    
    // Operators
    se23.def("__truediv__", [](const SE23d& g, const Vector& point_vel) {
        return point_vel / g;
    }, "Apply inverse transformation to a point-velocity pair", 
       py::arg("point_vel"));
}

void moduleAddBundle(py::module& m) {
    using Scalar = double; // Use double as the default scalar type
    
    // Specialization for SE3 + R3 (pose + velocity)
    using SE3d = SE3<Scalar>;
    using Vector3 = typename SE3d::Vector3;
    using PoseVel = Bundle<SE3<Scalar>, Vector3>;
    
    py::class_<PoseVel> pose_vel(m, "PoseVel", "Bundle of SE3 pose and R3 velocity");
    
    // Constructors
    pose_vel.def(py::init<>(), "Default constructor - identity transformation with zero velocity");
    pose_vel.def(py::init<const SE3d&, const Vector3&>(), 
                "Construct from SE3 pose and R3 velocity",
                py::arg("pose"), py::arg("velocity"));
    
    // Access components
    pose_vel.def("get_pose", &PoseVel::template get<0>, "Get the SE3 pose component");
    pose_vel.def("get_velocity", &PoseVel::template get<1>, "Get the R3 velocity component");
    
    // Common methods
    addCommonMethods<PoseVel>(pose_vel);
    
    // Specialization for SO3 + R3 (rotation + translation)
    using SO3d = SO3<Scalar>;
    using RotTrans = Bundle<SO3<Scalar>, Vector3>;
    
    py::class_<RotTrans> rot_trans(m, "RotTrans", "Bundle of SO3 rotation and R3 translation");
    
    // Constructors
    rot_trans.def(py::init<>(), "Default constructor - identity rotation with zero translation");
    rot_trans.def(py::init<const SO3d&, const Vector3&>(), 
                "Construct from SO3 rotation and R3 translation",
                py::arg("rotation"), py::arg("translation"));
    
    // Access components
    rot_trans.def("get_rotation", &RotTrans::template get<0>, "Get the SO3 rotation component");
    rot_trans.def("get_translation", &RotTrans::template get<1>, "Get the R3 translation component");
    
    // Common methods
    addCommonMethods<RotTrans>(rot_trans);
}

void moduleAddLieGroupUtils(py::module& m) {
    using LieUtils = Cosserat::LieGroupUtils<double>;
    
    // Angle utilities
    m.def("normalize_angle", &LieUtils::normalizeAngle, 
          "Normalize angle to [-π, π]", 
          py::arg("angle"));
    
    m.def("angle_difference", &LieUtils::angleDifference, 
          "Compute the difference between two angles with proper wrapping", 
          py::arg("a"), py::arg("b"));
    
    m.def("angle_distance", &LieUtils::angleDistance, 
          "Compute bi-invariant distance between two angles (as SO(2) elements)", 
          py::arg("a"), py::arg("b"));
    
    m.def("is_angle_near_zero", &LieUtils::isAngleNearZero, 
          "Check if an angle is near zero (within epsilon)", 
          py::arg("angle"));
    
    m.def("are_angles_nearly_equal", &LieUtils::areAnglesNearlyEqual, 
          "Check if two angles are nearly equal (within epsilon, considering wrapping)", 
          py::arg("a"), py::arg("b"));
    
    // Interpolation utilities
    m.def("lerp", &LieUtils::lerp, 
          "Linear interpolation between two scalars", 
          py::arg("a"), py::arg("b"), py::arg("t"));
    
    m.def("slerp_angle", &LieUtils::slerpAngle, 
          "Spherical linear interpolation (SLERP) between two angles", 
          py::arg("a"), py::arg("b"), py::arg("t"));
    
    // Numerical utilities
    m.def("sinc", &LieUtils::sinc, 
          "Compute sinc(x) = sin(x)/x with numerical stability for small x", 
          py::arg("x"));
    
    m.def("one_minus_cos", &LieUtils::oneMinusCos, 
          "Numerically stable computation of 1 - cos(x) for small x", 
          py::arg("x"));
    
    // Vector utilities
    m.def("safe_normalize", [](const Eigen::VectorXd& v) {
        return LieUtils::safeNormalize(v);
    }, "Safe vector normalization that handles near-zero vectors", 
       py::arg("v"));
    
    m.def("project_vector", [](const Eigen::VectorXd& v, const Eigen::VectorXd& onto) {
        return LieUtils::projectVector(v, onto);
    }, "Project a vector onto another vector", 
       py::arg("v"), py::arg("onto"));
    
    // SE(2) utilities
    m.def("interpolate_se2_path", [](const Eigen::Vector3d& start, 
                                   const Eigen::Vector3d& end, 
                                   double t) {
        return LieUtils::interpolateSE2Path(start, end, t);
    }, "Path interpolation between two SE(2) elements represented as [angle, x, y]", 
       py::arg("start"), py::arg("end"), py::arg("t"));
    
    // Constants
    m.attr("epsilon") = LieUtils::epsilon;
    m.attr("pi") = LieUtils::pi;
    m.attr("two_pi") = LieUtils::two_pi;
}

void moduleAddLieGroups(py::module& m) {
    // Create a submodule for Lie groups
    py::module liegroups = m.def_submodule("LieGroups", "Lie group implementations");
    
    // Add all Lie group classes to the submodule
    moduleAddSO2(liegroups);
    moduleAddSO3(liegroups);
    moduleAddSE2(liegroups);
    moduleAddSE3(liegroups);
    moduleAddSGal3(liegroups);
    moduleAddSE23(liegroups);
    moduleAddBundle(liegroups);
    moduleAddLieGroupUtils(liegroups);
}

} // namespace sofapython3

