// This file contains the pybind11 bindings for the Lie groups, exposing C++ Lie group functionalities to Python.

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
 ******************************************************************************/

// Standard library includes
#include <array>
#include <random>
#include <sstream>

// Lie group headers (via CMake include paths)
#include "Binding_LieGroups.h"
#include <liegroups/LieGroupBase.h>
#include <liegroups/LieGroupBase.inl>
#include <liegroups/RealSpace.h>
#include <liegroups/SE2.h>
#include <liegroups/SE23.h>
#include <liegroups/SE3.h>
#include <liegroups/SE3.inl>
#include <liegroups/SGal3.h>
#include <liegroups/SO2.h>
#include <liegroups/SO3.h>
#include <liegroups/Bundle.h>



// pybind11 includes
#include <pybind11/pybind11.h>
namespace py = pybind11;

namespace lg = sofa::component::cosserat::liegroups;


namespace sofapython3 {

    void moduleAddSO2(py::module &m) {
		py::class_<sofa::component::cosserat::liegroups::SO2<double>>(m, "SO2", "2D rotation group SO(2)")
				.def(py::init<>(), "Default constructor (identity rotation)")
				.def(py::init<double>(), "Construct from angle in radians", py::arg("angle"))
				.def("__mul__", &sofa::component::cosserat::liegroups::SO2<double>::operator*, "Group composition")
				.def("inverse", &sofa::component::cosserat::liegroups::SO2<double>::inverse, "Compute inverse rotation")
				.def("matrix", &sofa::component::cosserat::liegroups::SO2<double>::matrix, "Get 2x2 rotation matrix")
				.def("angle", &sofa::component::cosserat::liegroups::SO2<double>::angle, "Get rotation angle in radians")
				.def("setAngle", &sofa::component::cosserat::liegroups::SO2<double>::setAngle, "Set rotation angle in radians", py::arg("angle"))
				.def("complex", &sofa::component::cosserat::liegroups::SO2<double>::complex, "Get complex representation (cos θ, sin θ)")
				.def("direction", &sofa::component::cosserat::liegroups::SO2<double>::direction, "Get unit direction vector")
				.def("perpendicular", &sofa::component::cosserat::liegroups::SO2<double>::perpendicular, "Get perpendicular unit vector")
				.def_static("exp", &sofa::component::cosserat::liegroups::SO2<double>::exp, "Exponential map from so(2) to SO(2)", py::arg("omega"))
				.def("log", &sofa::component::cosserat::liegroups::SO2<double>::log, "Logarithmic map from SO(2) to so(2)")
				.def("adjoint", &lg::SO2<double>::adjoint, "Adjoint representation (identity for SO(2))")
				.def("isApprox", &lg::SO2<double>::isApprox, "Check approximate equality", py::arg("other"),
					 py::arg("eps") = 1e-12)
				.def_static(
						"identity", []() { return lg::SO2<double>::Identity(); }, "Get identity element")
				.def_static("hat", &lg::SO2<double>::hat, "Hat operator: R -> so(2) matrix", py::arg("omega"))
				.def_static("vee", &lg::SO2<double>::vee, "Vee operator: so(2) matrix -> R", py::arg("matrix"))
				.def(
						"act", [](const lg::SO2<double> &self, const Eigen::Vector2d &point) { return self.act(point); },
						"Apply rotation to a 2D point", py::arg("point"))
				.def("__repr__",
					 [](const lg::SO2<double> &self) {
						 std::ostringstream oss;
						 oss << "SO2(angle=" << self.angle() << " rad, " << (self.angle() * 180.0 / M_PI) << " deg)";
						 return oss.str();
					 })
				.def_static(
						"random",
						[](py::object seed = py::none()) {
							static std::random_device rd;
							static std::mt19937 gen(rd());
							if (!seed.is_none()) {
								gen.seed(seed.cast<unsigned int>());
							}
							return lg::SO2<double>::computeRandom(gen);
						},
						"Generate random SO(2) element", py::arg("seed") = py::none());
	}

    void moduleAddLieGroups(py::module &m) {
        // Add all Lie group bindings
        //moduleAddSO2(m);
        /*moduleAddSO3(m);
        moduleAddSE2(m);
        moduleAddSE3(m);
        moduleAddSGal3(m);
        moduleAddSE23(m);
        moduleAddRealSpace(m);
        moduleAddBundle(m);
        moduleAddLieGroupUtils(m); */
    }

}// namespace sofapython3

	/**/

	/*void moduleAddSO3(py::module &m) {
		py::class_<lg::SO3<double>>(m, "SO3", "3D rotation group SO(3)")
				.def(py::init<>(), "Default constructor (identity rotation)")
				.def(py::init<double, const Eigen::Vector3d &>(), "Construct from angle-axis representation",
					 py::arg("angle"), py::arg("axis"))
				.def(py::init<const Eigen::Quaterniond &>(), "Construct from quaternion", py::arg("quaternion"))
				.def(py::init<const Eigen::Matrix3d &>(), "Construct from rotation matrix", py::arg("matrix"))
				.def("__mul__", &lg::SO3<double>::operator*, "Group composition")
				.def("inverse", &lg::SO3<double>::inverse, "Compute inverse rotation")
				.def("matrix", &lg::SO3<double>::matrix, "Get 3x3 rotation matrix")
				.def("quaternion", &lg::SO3<double>::quaternion, "Get quaternion representation")
				.def_static("exp", &lg::SO3<double>::exp, "Exponential map from so(3) to SO(3)", py::arg("omega"))
				.def("log", &lg::SO3<double>::log, "Logarithmic map from SO(3) to so(3)")
				.def("adjoint", &lg::SO3<double>::adjoint, "Adjoint representation (rotation matrix for SO(3))")
				.def("isApprox", &lg::SO3<double>::isApprox, "Check approximate equality", py::arg("other"),
					 py::arg("eps") = 1e-12)
				.def_static("identity", &lg::SO3<double>::identity, "Get identity element")
				.def_static("hat", &lg::SO3<double>::hat, "Hat operator: R^3 -> so(3) matrix", py::arg("omega"))
				.def_static("vee", &lg::SO3<double>::vee, "Vee operator: so(3) matrix -> R^3", py::arg("matrix"))
				.def(
						"act", [](const lg::SO3<double> &self, const Eigen::Vector3d &point) { return self.act(point); },
						"Apply rotation to a 3D point", py::arg("point"))
				.def("__repr__",
					 [](const lg::SO3<double> &self) {
						 std::ostringstream oss;
						 const auto quat = self.quaternion();
						 oss << "SO3(quat=[" << quat.w() << ", " << quat.x() << ", " << quat.y() << ", " << quat.z()
							 << "])";
						 return oss.str();
					 })
				.def_static(
						"random",
						[](py::object seed = py::none()) {
							static std::random_device rd;
							static std::mt19937 gen(rd());
							if (!seed.is_none()) {
								gen.seed(seed.cast<unsigned int>());
							}
							return lg::SO3<double>::computeRandom(gen);
						},
						"Generate random SO(3) element", py::arg("seed") = py::none());
	}*/

	/*void moduleAddSE2(py::module &m) {
		py::class_<lg::SE2<double>>(m, "SE2", "2D Euclidean group SE(2)")
				.def(py::init<>(), "Default constructor (identity transformation)")
				.def(py::init<const lg::SO2<double> &, const Eigen::Vector2d &>(),
					 "Construct from rotation and translation", py::arg("rotation"), py::arg("translation"))
				.def(py::init<const Eigen::Matrix3d &>(), "Construct from 3x3 transformation matrix", py::arg("matrix"))
				.def("__mul__", &lg::SE2<double>::operator*, "Group composition")
				.def("inverse", &lg::SE2<double>::inverse, "Compute inverse transformation")
				.def("matrix", &lg::SE2<double>::matrix, "Get 3x3 transformation matrix")
				.def("rotation", static_cast<const lg::SO2<double> &(lg::SE2<double>::*) () const>(&lg::SE2<double>::rotation),
					 "Get rotation part")
				.def("translation",
					 static_cast<const typename lg::SE2<double>::Vector2 &(lg::SE2<double>::*) () const>(
							 &lg::SE2<double>::translation),
					 "Get translation part")
				.def_static("exp", &lg::SE2<double>::exp, "Exponential map from se(2) to SE(2)", py::arg("xi"))
				.def("log", &lg::SE2<double>::log, "Logarithmic map from SE(2) to se(2)")
				.def("adjoint", &lg::SE2<double>::adjoint, "Adjoint representation")
				.def("isApprox", &lg::SE2<double>::isApprox, "Check approximate equality", py::arg("other"),
					 py::arg("eps") = 1e-12)
				.def_static(
						"identity", []() { return lg::SE2<double>::Identity(); }, "Get identity element")
				.def(
						"act", [](const lg::SE2<double> &self, const Eigen::Vector2d &point) { return self.act(point); },
						"Apply transformation to a 2D point", py::arg("point"))
				.def("__repr__",
					 [](const lg::SE2<double> &self) {
						 std::ostringstream oss;
						 const auto &rot = self.rotation();
						 const auto &trans = self.translation();
						 oss << "SE2(angle=" << rot.angle() << " rad, translation=[" << trans.x() << ", " << trans.y()
							 << "])";
						 return oss.str();
					 })
				.def_static(
						"random",
						[](py::object seed = py::none()) {
							static std::random_device rd;
							static std::mt19937 gen(rd());
							if (!seed.is_none()) {
								gen.seed(seed.cast<unsigned int>());
							}
							return lg::SE2<double>::computeRandom(gen);
						},
						"Generate random SE(2) element", py::arg("seed") = py::none());
	}*/

	/*void moduleAddSE3(py::module &m) {
		py::class_<lg::SE3<double>>(m, "SE3", "3D Euclidean group SE(3)")
				.def(py::init<>(), "Default constructor (identity transformation)")
				.def(py::init<const lg::SO3<double> &, const Eigen::Vector3d &>(),
					 "Construct from rotation and translation", py::arg("rotation"), py::arg("translation"))
				.def(py::init<const Eigen::Matrix4d &>(), "Construct from 4x4 transformation matrix", py::arg("matrix"))
				.def("__mul__", &lg::SE3<double>::operator*, "Group composition")
				.def("inverse", &lg::SE3<double>::inverse, "Compute inverse transformation")
				.def("matrix", &lg::SE3<double>::matrix, "Get 4x4 transformation matrix")
				.def("rotation", static_cast<const lg::SO3<double> &(lg::SE3<double>::*) () const>(&lg::SE3<double>::rotation),
					 "Get rotation part")
				.def("translation",
					 static_cast<const typename lg::SE3<double>::Vector3 &(lg::SE3<double>::*) () const>(
							 &lg::SE3<double>::translation),
					 "Get translation part")
				.def_static("exp", &lg::SE3<double>::exp, "Exponential map from se(3) to SE(3)", py::arg("xi"))
				.def("log", &lg::SE3<double>::log, "Logarithmic map from SE(3) to se(3)")
				.def("adjoint", &lg::SE3<double>::adjoint, "Adjoint representation")
				.def("isApprox", &lg::SE3<double>::isApprox, "Check approximate equality", py::arg("other"),
					 py::arg("eps") = 1e-12)
				.def_static(
						"identity", []() { return lg::SE3<double>::Identity(); }, "Get identity element")
				.def(
						"act", [](const lg::SE3<double> &self, const Eigen::Vector3d &point) { return self.act(point); },
						"Apply transformation to a 3D point", py::arg("point"))
				.def("__repr__",
					 [](const lg::SE3<double> &self) {
						 std::ostringstream oss;
						 const auto &rot = self.rotation().quaternion();
						 const auto &trans = self.translation();
						 oss << "SE3(quat=[" << rot.w() << ", " << rot.x() << ", " << rot.y() << ", " << rot.z()
							 << "], translation=[" << trans.x() << ", " << trans.y() << ", " << trans.z() << "])";
						 return oss.str();
					 })
				.def_static(
						"random",
						[](py::object seed = py::none()) {
							static std::random_device rd;
							static std::mt19937 gen(rd());
							if (!seed.is_none()) {
								gen.seed(seed.cast<unsigned int>());
							}
							return lg::SE3<double>::computeRandom(gen);
						},
						"Generate random SE(3) element", py::arg("seed") = py::none());
	}*/

	/*void moduleAddSGal3(py::module &m) {
		py::class_<lg::SGal3<double>>(m, "SGal3", "Special Galilean group SGal(3)")
				.def(py::init<>(), "Default constructor (identity transformation)")
				.def(py::init<const lg::SE3<double> &, const Eigen::Vector3d &, double>(),
					 "Construct from pose, velocity, and time", py::arg("pose"), py::arg("velocity"), py::arg("time"))
				.def(py::init<const lg::SO3<double> &, const Eigen::Vector3d &, const Eigen::Vector3d &, double>(),
					 "Construct from rotation, position, velocity, and time", py::arg("rotation"), py::arg("position"),
					 py::arg("velocity"), py::arg("time"))
				.def("inverse", &lg::SGal3<double>::inverse, "Compute inverse transformation")
				.def("matrix", &lg::SGal3<double>::extendedMatrix, "Get 6x6 extended transformation matrix")
				.def("pose", static_cast<const lg::SE3<double> &(lg::SGal3<double>::*) () const>(&lg::SGal3<double>::pose),
					 "Get pose part")
				.def("velocity",
					 static_cast<const Eigen::Vector3d &(lg::SGal3<double>::*) () const>(&lg::SGal3<double>::velocity),
					 "Get velocity part")
				.def("time", static_cast<const double &(lg::SGal3<double>::*) () const>(&lg::SGal3<double>::time),
					 "Get time part")
				.def_static("exp", &lg::SGal3<double>::exp, "Exponential map from sgal(3) to SGal(3)", py::arg("xi"))
				.def("log", &lg::SGal3<double>::log, "Logarithmic map from SGal(3) to sgal(3)")
				.def("adjoint", &lg::SGal3<double>::adjoint, "Adjoint representation")
				.def("isApprox", &lg::SGal3<double>::isApprox, "Check approximate equality", py::arg("other"),
					 py::arg("eps") = 1e-12)
				.def_static(
						"identity", []() { return lg::SGal3<double>::Identity(); }, "Get identity element")
				.def(
						"act",
						[](const lg::SGal3<double> &self, const Eigen::Matrix<double, 10, 1> &point_vel_time) {
							return self.act(point_vel_time);
						},
						"Apply transformation to a 10D point-velocity-time vector", py::arg("point_vel_time"))
				.def("__repr__",
					 [](const lg::SGal3<double> &self) {
						 std::ostringstream oss;
						 const auto &pose = self.pose();
						 const auto &vel = self.velocity();
						 oss << "SGal3(pose=" << pose << ", velocity=[" << vel.x() << ", " << vel.y() << ", " << vel.z()
							 << "], time=" << self.time() << ")";
						 return oss.str();
					 })
				.def_static(
						"random",
						[](py::object seed = py::none()) {
							static std::random_device rd;
							static std::mt19937 gen(rd());
							if (!seed.is_none()) {
								gen.seed(seed.cast<unsigned int>());
							}
							return lg::SGal3<double>::computeRandom(gen);
						},
						"Generate random SGal3 element", py::arg("seed") = py::none());
	}*/

	/*void moduleAddSE23(py::module &m) {
		py::class_<lg::SE23<double>>(m, "SE23", "Extended 3D Euclidean group SE_2(3)")
				.def(py::init<>(), "Default constructor (identity transformation)")
				.def(py::init<const lg::SE3<double> &, const Eigen::Vector3d &>(), "Construct from pose and velocity",
					 py::arg("pose"), py::arg("velocity"))
				.def(py::init<const lg::SO3<double> &, const Eigen::Vector3d &, const Eigen::Vector3d &>(),
					 "Construct from rotation, position, and velocity", py::arg("rotation"), py::arg("position"),
					 py::arg("velocity"))
				.def("inverse", &lg::SE23<double>::inverse, "Compute inverse transformation")
				.def("matrix", &lg::SE23<double>::extendedMatrix, "Get 5x5 extended transformation matrix")
				.def("pose", static_cast<const lg::SE3<double> &(lg::SE23<double>::*) () const>(&lg::SE23<double>::pose),
					 "Get pose part")
				.def("velocity",
					 static_cast<const Eigen::Vector3d &(lg::SE23<double>::*) () const>(&lg::SE23<double>::velocity),
					 "Get velocity part")
				.def_static("exp", &lg::SE23<double>::exp, "Exponential map from se_2(3) to SE_2(3)", py::arg("xi"))
				.def("log", &lg::SE23<double>::log, "Logarithmic map from SE_2(3) to se_2(3)")
				.def("adjoint", &lg::SE23<double>::adjoint, "Adjoint representation")
				.def("isApprox", &lg::SE23<double>::isApprox, "Check approximate equality", py::arg("other"),
					 py::arg("eps") = 1e-12)
				.def_static(
						"identity", []() { return lg::SE23<double>::Identity(); }, "Get identity element")
				.def(
						"act",
						[](const lg::SE23<double> &self, const Eigen::Matrix<double, 6, 1> &point_vel) {
							return self.act(point_vel);
						},
						"Apply transformation to a 6D point-velocity vector", py::arg("point_vel"))
				.def("__repr__",
					 [](const lg::SE23<double> &self) {
						 std::ostringstream oss;
						 const auto &pose = self.pose();
						 const auto &vel = self.velocity();
						 oss << "SE23(pose=" << pose << ", velocity=[" << vel.x() << ", " << vel.y() << ", " << vel.z()
							 << "])";
						 return oss.str();
					 })
				.def_static(
						"random",
						[](py::object seed = py::none()) {
							static std::random_device rd;
							static std::mt19937 gen(rd());
							if (!seed.is_none()) {
								gen.seed(seed.cast<unsigned int>());
							}
							return lg::SE23<double>::computeRandom(gen);
						},
						"Generate random SE23 element", py::arg("seed") = py::none());
	}*/

	/*void moduleAddBundle(py::module &m) {
		// SE3_Velocity Bundle (SE3 + R^6 velocity)
		using SE3_Velocity_double = lg::Bundle<lg::SE3<double>, lg::RealSpace<double, 6>>;
		py::class_<SE3_Velocity_double>(m, "SE3_Velocity", "SE(3) x R^6 bundle for pose with 6D velocity")
				.def(py::init<>(), "Default constructor (identity)")
				.def(py::init<const lg::SE3<double> &, const lg::RealSpace<double, 6> &>(), "Construct from SE(3) and R^6",
					 py::arg("pose"), py::arg("velocity"))
				.def("__mul__", &SE3_Velocity_double::operator*, "Group composition")
				.def("inverse", &SE3_Velocity_double::inverse, "Compute inverse")
				.def("log", &SE3_Velocity_double::log, "Logarithmic map to algebra")
				.def_static(
						"exp",
						[](const SE3_Velocity_double::TangentVector &xi) {
							return SE3_Velocity_double::identity().exp(xi);
						},
						"Exponential map from algebra", py::arg("xi"))
				.def("adjoint", &SE3_Velocity_double::adjoint, "Adjoint representation")
				.def("isApprox", &SE3_Velocity_double::isApprox, "Check approximate equality", py::arg("other"),
					 py::arg("eps") = 1e-12)
				.def_static(
						"identity", []() { return SE3_Velocity_double::identity(); }, "Get identity element")
				.def(
						"get_pose", [](const SE3_Velocity_double &self) { return self.get<0>(); },
						"Get the SE(3) component")
				.def(
						"get_velocity", [](const SE3_Velocity_double &self) { return self.get<1>(); },
						"Get the R^6 component")
				.def("__repr__", [](const SE3_Velocity_double &self) {
					std::ostringstream oss;
					oss << "SE3_Velocity(pose=" << self.get<0>() << ", velocity=" << self.get<1>() << ")";
					return oss.str();
				});*/

		// SE2_Velocity Bundle (SE2 + R^3 velocity)
		/*using SE2_Velocity_double = lg::Bundle<lg::SE2<double>, lg::RealSpace<double, 3>>;
		py::class_<SE2_Velocity_double>(m, "SE2_Velocity", "SE(2) x R^3 bundle for 2D pose with 3D velocity")
				.def(py::init<>(), "Default constructor (identity)")
				.def(py::init<const lg::SE2<double> &, const lg::RealSpace<double, 3> &>(), "Construct from SE(2) and R^3",
					 py::arg("pose"), py::arg("velocity"))
				.def("__mul__", &SE2_Velocity_double::operator*, "Group composition")
				.def("inverse", &SE2_Velocity_double::inverse, "Compute inverse")
				.def("log", &SE2_Velocity_double::log, "Logarithmic map to algebra")
				.def_static(
						"exp",
						[](const SE2_Velocity_double::TangentVector &xi) {
							return SE2_Velocity_double::identity().exp(xi);
						},
						"Exponential map from algebra", py::arg("xi"))
				.def("adjoint", &SE2_Velocity_double::adjoint, "Adjoint representation")
				.def("isApprox", &SE2_Velocity_double::isApprox, "Check approximate equality", py::arg("other"),
					 py::arg("eps") = 1e-12)
				.def_static(
						"identity", []() { return SE2_Velocity_double::identity(); }, "Get identity element")
				.def(
						"get_pose", [](const SE2_Velocity_double &self) { return self.get<0>(); },
						"Get the SE(2) component")
				.def(
						"get_velocity", [](const SE2_Velocity_double &self) { return self.get<1>(); },
						"Get the R^3 component")
				.def("__repr__", [](const SE2_Velocity_double &self) {
					std::ostringstream oss;
					oss << "SE2_Velocity(pose=" << self.get<0>() << ", velocity=" << self.get<1>() << ")";
					return oss.str();
				});
	}*/

	/*void moduleAddRealSpace(py::module &m) {
		// lg::RealSpace<3> bindings
		py::class_<lg::RealSpace<double, 3>>(m, "R3", "3D real space R^3")
				.def(py::init<>(), "Default constructor (zero vector)")
				.def(py::init<const Eigen::Vector3d &>(), "Construct from 3D vector", py::arg("vector"))
				.def(
						"__mul__",
						[](const lg::RealSpace<double, 3> &self, const lg::RealSpace<double, 3> &other) {
							return self.compose(other);
						},
						"Group composition (vector addition)")
				.def(
						"inverse", [](const lg::RealSpace<double, 3> &self) { return self.inverse(); },
						"Inverse (vector negation)")
				.def(
						"vector", [](const lg::RealSpace<double, 3> &self) { return self.computeLog(); },
						"Get the underlying vector");

		// lg::RealSpace<6> bindings
		py::class_<lg::RealSpace<double, 6>>(m, "R6", "6D real space R^6")
				.def(py::init<>(), "Default constructor (zero vector)")
				.def(py::init<const Eigen::Matrix<double, 6, 1> &>(), "Construct from 6D vector", py::arg("vector"))
				.def(
						"__mul__",
						[](const lg::RealSpace<double, 6> &self, const lg::RealSpace<double, 6> &other) {
							return self.compose(other);
						},
						"Group composition (vector addition)")
				.def(
						"inverse", [](const lg::RealSpace<double, 6> &self) { return self.inverse(); },
						"Inverse (vector negation)")
				.def(
						"vector", [](const lg::RealSpace<double, 6> &self) { return self.computeLog(); },
						"Get the underlying vector");
	}*/

	/*void moduleAddLieGroupUtils(py::module &m) {
		// Utility functions for interpolation, etc.
		// Note: slerp function would need to be implemented in the Lie group classes
		// m.def("slerp", [](const lg::SO3<double>& a, const lg::SO3<double>& b, double t) {
		//     return a.interpolate(b, t);
		// }, "Spherical linear interpolation for SO3");
		(void)m; // suppress unused parameter warning
	}*/