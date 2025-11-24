// This file contains the pybind11 bindings for the Lie groups, exposing C++ Lie group functionalities to Python.

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

// Standard library includes first
#include <array>
#include <functional>
#include <iostream>
#include <iterator>
#include <memory>
#include <random>
#include <tuple>
#include <type_traits>

// Eigen/SOFA includes before pybind11
#include "../../../../src/liegroups/LieGroupBase.h"
#include "../../../../src/liegroups/LieGroupBase.inl"
#include "../../../../src/liegroups/RealSpace.h"
#include "../../../../src/liegroups/SE2.h"
#include "../../../../src/liegroups/SE3.h"
#include "../../../../src/liegroups/SGal3.h"
#include "../../../../src/liegroups/SO2.h"
#include "../../../../src/liegroups/SO3.h"
#include "../../../../src/liegroups/Sim3.h"
#include "../../../../src/liegroups/Uncertainty.h"
#include "../../../../src/liegroups/Types.h"
#include "../../../../src/liegroups/Bundle.h"

// pybind11 includes last to avoid template conflicts
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "Binding_LieGroups.h"

namespace sofa::component::cosserat::liegroups {

	namespace detail {

		// Helper to compute total dimension of all groups
		template<typename... Groups>
		struct TotalDimension;

		template<>
		struct TotalDimension<> {
			static constexpr int value = 0;
		};

		template<typename Group, typename... Rest>
		struct TotalDimension<Group, Rest...> {
			static constexpr int value = Group::Dim + TotalDimension<Rest...>::value;
		};

		// Helper to compute total action dimension
		template<typename... Groups>
		struct TotalActionDimension;

		template<>
		struct TotalActionDimension<> {
			static constexpr int value = 0;
		};

		template<typename Group, typename... Rest>
		struct TotalActionDimension<Group, Rest...> {
			// This assumes we know actionDimension at compile time
			// For runtime computation, we'll use a different approach
			static constexpr int value = -1; // Indicates runtime computation needed
		};

		// Helper to check if all types are LieGroupBase derivatives with same scalar
		// type
		template<typename FirstScalar, typename... Groups>
		struct AllAreLieGroups;

		template<typename FirstScalar>
		struct AllAreLieGroups<FirstScalar> : std::true_type {};

		template<typename FirstScalar, typename Group, typename... Rest>
		struct AllAreLieGroups<FirstScalar, Group, Rest...> {
			static constexpr bool value =
					std::is_base_of_v<LieGroupBase<typename Group::Scalar, std::integral_constant<int, Group::Dim>,
												   Group::Dim, Group::Dim>,
									  Group> &&
					std::is_same_v<FirstScalar, typename Group::Scalar> && AllAreLieGroups<FirstScalar, Rest...>::value;
		};

		// Compile-time offset computation
		template<int Index, typename... Groups>
		struct OffsetAt;

		template<int Index>
		struct OffsetAt<Index> {
			static constexpr int value = 0;
		};

		template<int Index, typename Group, typename... Rest>
		struct OffsetAt<Index, Group, Rest...> {
			static constexpr int value = (Index == 0) ? 0 : Group::Dim + OffsetAt<Index - 1, Rest...>::value;
		};

		// Runtime offset computation for action dimensions
		template<typename... Groups>
		class ActionOffsets {
		public:
			explicit ActionOffsets(const std::tuple<Groups...> &groups) {
				computeOffsets(groups, std::index_sequence_for<Groups...>());
			}

			int operator[](std::size_t i) const { return offsets_[i]; }
			int total() const { return total_; }

		private:
			std::array<int, sizeof...(Groups)> offsets_;
			int total_ = 0;

			template<std::size_t... Is>
			void computeOffsets(const std::tuple<Groups...> &groups, std::index_sequence<Is...>) {
				int offset = 0;
				((offsets_[Is] = offset, offset += std::get<Is>(groups).actionDimension()), ...);
				total_ = offset;
			}
		};

	} // namespace detail



// Python bindings implementation
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "../../../../src/liegroups/SE2.h"
#include "../../../../src/liegroups/SE23.h"
#include "../../../../src/liegroups/SE3.h"
#include "../../../../src/liegroups/SE3.inl"
#include "../../../../src/liegroups/SGal3.h"
#include "../../../../src/liegroups/SO2.h"
#include "../../../../src/liegroups/SO3.h"
#include "Binding_LieGroups.h"

namespace py = pybind11;
using namespace sofa::component::cosserat::liegroups;

namespace sofapython3 {

	void moduleAddSO2(py::module &m) {
		// SO2 bindings
		py::class_<SO2<double>>(m, "SO2", "2D rotation group SO(2)")
				.def(py::init<>(), "Default constructor (identity rotation)")
				.def(py::init<double>(), "Construct from angle in radians", py::arg("angle"))
				.def("__mul__", &SO2<double>::operator*, "Group composition")
				.def("inverse", &SO2<double>::inverse, "Compute inverse rotation")
				.def("matrix", &SO2<double>::matrix, "Get 2x2 rotation matrix")
				.def("angle", &SO2<double>::angle, "Get rotation angle in radians")
				.def("setAngle", &SO2<double>::setAngle, "Set rotation angle in radians", py::arg("angle"))
				.def("complex", &SO2<double>::complex, "Get complex representation (cos θ, sin θ)")
				.def("direction", &SO2<double>::direction, "Get unit direction vector")
				.def("perpendicular", &SO2<double>::perpendicular, "Get perpendicular unit vector")
				.def_static("exp", &SO2<double>::exp, "Exponential map from so(2) to SO(2)", py::arg("omega"))
				.def("log", &SO2<double>::log, "Logarithmic map from SO(2) to so(2)")
				.def("adjoint", &SO2<double>::adjoint, "Adjoint representation (identity for SO(2))")
				.def("isApprox", &SO2<double>::isApprox, "Check approximate equality", py::arg("other"), py::arg("eps") = 1e-12)
				.def_static("identity", &SO2<double>::Identity, "Get identity element")
				.def("__add__", &SO2<double>::operator+, "Right-plus operator: X + τ = X ⊕ τ", py::arg("tau"))
				.def("__sub__", [](const SO2<double>& self, const SO2<double>& other) {
					return self - other;
				}, "Right-minus operator: Y - X = Log(X⁻¹ ◦ Y)", py::arg("other"))
				.def("distance", &SO2<double>::distance, "Compute geodesic distance to another element", py::arg("other"))
				.def("lerp", &SO2<double>::lerp, "Linear interpolation between elements", py::arg("other"), py::arg("t"))
				.def("slerp", &SO2<double>::slerp, "Spherical linear interpolation between elements", py::arg("other"), py::arg("t"))
				.def_static("BCH", &SO2<double>::BCH, "Baker-Campbell-Hausdorff formula", py::arg("v"), py::arg("w"), py::arg("order") = 3)
				.def("actionJacobian", &SO2<double>::actionJacobian, "Jacobian of group action", py::arg("point"))
				.def("dlog", &SO2<double>::dlog, "Differential of logarithm map")
				.def_static("hat", &SO2<double>::hat, "Hat operator: R -> so(2) matrix", py::arg("omega"))
				.def_static("vee", &SO2<double>::vee, "Vee operator: so(2) matrix -> R", py::arg("matrix"))
				.def("act", [](const SO2<double>& self, const Eigen::Vector2d& point) {
					return self.act(point);
				}, "Apply rotation to a 2D point", py::arg("point"))
				.def("__repr__", [](const SO2<double>& self) {
					std::ostringstream oss;
					oss << "SO2(angle=" << self.angle() << " rad, " << (self.angle() * 180.0 / M_PI) << " deg)";
					return oss.str();
				})
				.def_static("random", [](py::object seed = py::none()) {
					static std::random_device rd;
					static std::mt19937 gen(rd());
					if (!seed.is_none()) {
						gen.seed(seed.cast<unsigned int>());
					}
					return SO2<double>::computeRandom(gen);
				}, "Generate random SO(2) element", py::arg("seed") = py::none());
	}

	void moduleAddSO3(py::module &m) {
		// SO3 bindings
		py::class_<SO3<double>>(m, "SO3", "3D rotation group SO(3)")
				.def(py::init<>(), "Default constructor (identity rotation)")
				.def(py::init<double, const Eigen::Vector3d &>(), "Construct from angle-axis representation", py::arg("angle"), py::arg("axis"))
				.def(py::init<const Eigen::Quaterniond &>(), "Construct from quaternion", py::arg("quaternion"))
				.def(py::init<const Eigen::Matrix3d &>(), "Construct from rotation matrix", py::arg("matrix"))
				.def("__mul__", &SO3<double>::operator*, "Group composition")
				.def("inverse", &SO3<double>::inverse, "Compute inverse rotation")
				.def("matrix", &SO3<double>::matrix, "Get 3x3 rotation matrix")
				.def("quaternion", &SO3<double>::quaternion, "Get quaternion representation")
				.def_static("exp", &SO3<double>::exp, "Exponential map from so(3) to SO(3)", py::arg("omega"))
				.def("log", &SO3<double>::log, "Logarithmic map from SO(3) to so(3)")
				.def("adjoint", &SO3<double>::adjoint, "Adjoint representation (rotation matrix for SO(3))")
				.def("isApprox", &SO3<double>::isApprox, "Check approximate equality", py::arg("other"), py::arg("eps") = 1e-12)
				.def_static("identity", &SO3<double>::Identity, "Get identity element")
				.def("__add__", &SO3<double>::operator+, "Right-plus operator: X + τ = X ⊕ τ", py::arg("tau"))
				.def("__sub__", [](const SO3<double>& self, const SO3<double>& other) {
					return self - other;
				}, "Right-minus operator: Y - X = Log(X⁻¹ ◦ Y)", py::arg("other"))
				.def("distance", &SO3<double>::distance, "Compute geodesic distance to another element", py::arg("other"))
				.def("lerp", &SO3<double>::lerp, "Linear interpolation between elements", py::arg("other"), py::arg("t"))
				.def("slerp", &SO3<double>::slerp, "Spherical linear interpolation between elements", py::arg("other"), py::arg("t"))
				.def_static("BCH", &SO3<double>::BCH, "Baker-Campbell-Hausdorff formula", py::arg("v"), py::arg("w"), py::arg("order") = 3)
				.def("actionJacobian", &SO3<double>::actionJacobian, "Jacobian of group action", py::arg("point"))
				.def("dlog", &SO3<double>::dlog, "Differential of logarithm map")
				.def_static("hat", &SO3<double>::hat, "Hat operator: R^3 -> so(3) matrix", py::arg("omega"))
				.def_static("vee", &SO3<double>::vee, "Vee operator: so(3) matrix -> R^3", py::arg("matrix"))
				.def("act", [](const SO3<double>& self, const Eigen::Vector3d& point) {
					return self.act(point);
				}, "Apply rotation to a 3D point", py::arg("point"))
				.def("__repr__", [](const SO3<double>& self) {
					std::ostringstream oss;
					const auto quat = self.quaternion();
					oss << "SO3(quat=[" << quat.w() << ", " << quat.x() << ", " << quat.y() << ", " << quat.z() << "])";
					return oss.str();
				})
				.def_static("random", [](py::object seed = py::none()) {
					static std::random_device rd;
					static std::mt19937 gen(rd());
					if (!seed.is_none()) {
						gen.seed(seed.cast<unsigned int>());
					}
					return SO3<double>::computeRandom(gen);
				}, "Generate random SO(3) element", py::arg("seed") = py::none());
	}

	void moduleAddSE2(py::module &m) {
		// SE2 bindings
		py::class_<SE2<double>>(m, "SE2", "2D Euclidean group SE(2)")
				.def(py::init<>(), "Default constructor (identity transformation)")
				.def(py::init<const SO2<double> &, const Eigen::Vector2d &>(), "Construct from rotation and translation", py::arg("rotation"), py::arg("translation"))
				.def(py::init<const Eigen::Matrix3d &>(), "Construct from 3x3 transformation matrix", py::arg("matrix"))
				.def("__mul__", &SE2<double>::operator*, "Group composition")
				.def("inverse", &SE2<double>::inverse, "Compute inverse transformation")
				.def("matrix", &SE2<double>::matrix, "Get 3x3 transformation matrix")
				.def("rotation", static_cast<const SO2<double> &(SE2<double>::*) () const>(&SE2<double>::rotation), "Get rotation part")
				.def("translation", static_cast<const typename SE2<double>::Vector2 &(SE2<double>::*) () const>(
											&SE2<double>::translation), "Get translation part")
				.def_static("exp", &SE2<double>::exp, "Exponential map from se(2) to SE(2)", py::arg("xi"))
				.def("log", &SE2<double>::log, "Logarithmic map from SE(2) to se(2)")
				.def("adjoint", &SE2<double>::adjoint, "Adjoint representation")
				.def("isApprox", &SE2<double>::isApprox, "Check approximate equality", py::arg("other"), py::arg("eps") = 1e-12)
				.def_static("identity", &SE2<double>::Identity, "Get identity element")
				.def("__add__", &SE2<double>::operator+, "Right-plus operator: X + τ = X ⊕ τ", py::arg("tau"))
				.def("__sub__", [](const SE2<double>& self, const SE2<double>& other) {
					return self - other;
				}, "Right-minus operator: Y - X = Log(X⁻¹ ◦ Y)", py::arg("other"))
				.def("distance", &SE2<double>::distance, "Compute geodesic distance to another element", py::arg("other"))
				.def("lerp", &SE2<double>::lerp, "Linear interpolation between elements", py::arg("other"), py::arg("t"))
				.def("slerp", &SE2<double>::slerp, "Spherical linear interpolation between elements", py::arg("other"), py::arg("t"))
				.def_static("BCH", &SE2<double>::BCH, "Baker-Campbell-Hausdorff formula", py::arg("v"), py::arg("w"), py::arg("order") = 3)
				.def("actionJacobian", &SE2<double>::actionJacobian, "Jacobian of group action", py::arg("point"))
				.def("dlog", &SE2<double>::dlog, "Differential of logarithm map")
				.def("act", [](const SE2<double>& self, const Eigen::Vector2d& point) {
					return self.act(point);
				}, "Apply transformation to a 2D point", py::arg("point"))
				.def("__repr__", [](const SE2<double>& self) {
					std::ostringstream oss;
					const auto& rot = self.rotation();
					const auto& trans = self.translation();
					oss << "SE2(angle=" << rot.angle() << " rad, translation=[" << trans.x() << ", " << trans.y() << "])";
					return oss.str();
				})
				.def_static("random", [](py::object seed = py::none()) {
					static std::random_device rd;
					static std::mt19937 gen(rd());
					if (!seed.is_none()) {
						gen.seed(seed.cast<unsigned int>());
					}
					return SE2<double>::computeRandom(gen);
				}, "Generate random SE(2) element", py::arg("seed") = py::none());
	}

	void moduleAddSE3(py::module &m) {
		// SE3 bindings with enhanced functionality including Jacobians
		py::class_<SE3<double>>(m, "SE3", "3D Euclidean group SE(3)")
				.def(py::init<>(), "Default constructor (identity transformation)")
				.def(py::init<const SO3<double> &, const Eigen::Vector3d &>(), "Construct from rotation and translation", py::arg("rotation"), py::arg("translation"))
				.def(py::init<const Eigen::Matrix4d &>(), "Construct from 4x4 transformation matrix", py::arg("matrix"))
				.def("__mul__", &SE3<double>::operator*, "Group composition")
				.def("inverse", &SE3<double>::inverse, "Compute inverse transformation")
				.def("matrix", &SE3<double>::matrix, "Get 4x4 transformation matrix")
				.def("rotation", static_cast<const SO3<double> &(SE3<double>::*) () const>(&SE3<double>::rotation), "Get rotation part")
				.def("translation", static_cast<const typename SE3<double>::Vector3 &(SE3<double>::*) () const>(
											&SE3<double>::translation), "Get translation part")
				.def_static("exp", &SE3<double>::exp, "Exponential map from se(3) to SE(3)", py::arg("xi"))
				.def_static("expCosserat", &SE3<double>::expCosserat, "Cosserat-style exponential map", py::arg("strain"), py::arg("length"))
				.def("log", &SE3<double>::log, "Logarithmic map from SE(3) to se(3)")
				.def("adjoint", &SE3<double>::adjoint, "Adjoint representation")
				.def("isApprox", &SE3<double>::isApprox, "Check approximate equality", py::arg("other"), py::arg("eps") = 1e-12)
				.def_static("identity", &SE3<double>::Identity, "Get identity element")
				.def("__add__", &SE3<double>::operator+, "Right-plus operator: X + τ = X ⊕ τ", py::arg("tau"))
				.def("__sub__", [](const SE3<double>& self, const SE3<double>& other) {
					return self - other;
				}, "Right-minus operator: Y - X = Log(X⁻¹ ◦ Y)", py::arg("other"))
				.def("distance", static_cast<double (SE3<double>::*)(const SE3<double>&, double, double) const>(&SE3<double>::distance), "Compute weighted distance to another element", py::arg("other"), py::arg("w_rot") = 1.0, py::arg("w_trans") = 1.0)
				.def("lerp", &SE3<double>::lerp, "Linear interpolation between elements", py::arg("other"), py::arg("t"))
				.def("slerp", &SE3<double>::slerp, "Spherical linear interpolation between elements", py::arg("other"), py::arg("t"))
				.def_static("BCH", &SE3<double>::BCH, "Baker-Campbell-Hausdorff formula", py::arg("v"), py::arg("w"), py::arg("order") = 3)
				.def("actionJacobian", &SE3<double>::actionJacobian, "Jacobian of group action", py::arg("point"))
				.def("dlog", &SE3<double>::dlog, "Differential of logarithm map")
				.def("act", [](const SE3<double>& self, const Eigen::Vector3d& point) {
					return self.act(point);
				}, "Apply transformation to a 3D point", py::arg("point"))
				// Jacobian methods for uncertainty propagation
				.def_static("rightJacobian", &SE3<double>::rightJacobian, "Right Jacobian for uncertainty propagation", py::arg("xi"))
				.def_static("leftJacobian", &SE3<double>::leftJacobian, "Left Jacobian for uncertainty propagation", py::arg("xi"))
				.def_static("rightJacobianInverse", &SE3<double>::rightJacobianInverse, "Inverse right Jacobian", py::arg("xi"))
				.def_static("leftJacobianInverse", &SE3<double>::leftJacobianInverse, "Inverse left Jacobian", py::arg("xi"))
				// Advanced features
				.def("interpolate", &SE3<double>::interpolate, "Interpolate between transformations", py::arg("other"), py::arg("t"))
				.def_static("generateTrajectory", &SE3<double>::generateTrajectory, "Generate trajectory between waypoints", py::arg("waypoints"), py::arg("num_points") = 10)
				.def("__repr__", [](const SE3<double>& self) {
					std::ostringstream oss;
					const auto& rot = self.rotation().quaternion();
					const auto& trans = self.translation();
					oss << "SE3(quat=[" << rot.w() << ", " << rot.x() << ", " << rot.y() << ", " << rot.z() << "], translation=[" << trans.x() << ", " << trans.y() << ", " << trans.z() << "])";
					return oss.str();
				})
				.def_static("random", [](py::object seed = py::none()) {
					static std::random_device rd;
					static std::mt19937 gen(rd());
					if (!seed.is_none()) {
						gen.seed(seed.cast<unsigned int>());
					}
					return SE3<double>::computeRandom(gen);
				}, "Generate random SE(3) element", py::arg("seed") = py::none());
	}

	void moduleAddSGal3(py::module &m) {
		// SGal3 bindings
		py::class_<SGal3<double>>(m, "SGal3", "Special Galilean group SGal(3)")
				.def(py::init<>(), "Default constructor (identity transformation)")
				.def(py::init<const SE3<double> &, const Eigen::Vector3d &, double>(), "Construct from pose, velocity, and time", py::arg("pose"), py::arg("velocity"), py::arg("time"))
				.def(py::init<const SO3<double> &, const Eigen::Vector3d &, const Eigen::Vector3d &, double>(), "Construct from rotation, position, velocity, and time", py::arg("rotation"), py::arg("position"), py::arg("velocity"), py::arg("time"))
				.def("inverse", &SGal3<double>::inverse, "Compute inverse transformation")
				.def("matrix", &SGal3<double>::extendedMatrix, "Get 6x6 extended transformation matrix")
				.def("pose", static_cast<const SE3<double> &(SGal3<double>::*)() const>(&SGal3<double>::pose), "Get pose part")
				.def("velocity", static_cast<const Eigen::Vector3d &(SGal3<double>::*)() const>(&SGal3<double>::velocity), "Get velocity part")
				.def("time", static_cast<const double &(SGal3<double>::*)() const>(&SGal3<double>::time), "Get time part")
				.def_static("exp", &SGal3<double>::exp, "Exponential map from sgal(3) to SGal(3)", py::arg("xi"))
				.def("log", &SGal3<double>::log, "Logarithmic map from SGal(3) to sgal(3)")
				.def("adjoint", &SGal3<double>::adjoint, "Adjoint representation")
				.def("isApprox", &SGal3<double>::isApprox, "Check approximate equality", py::arg("other"), py::arg("eps") = 1e-12)
				.def_static("identity", &SGal3<double>::Identity, "Get identity element")
				.def("act", [](const SGal3<double>& self, const Eigen::Matrix<double, 10, 1>& point_vel_time) {
					return self.act(point_vel_time);
				}, "Apply transformation to a 10D point-velocity-time vector", py::arg("point_vel_time"))
				.def("__repr__", [](const SGal3<double>& self) {
					std::ostringstream oss;
					const auto& pose = self.pose();
					const auto& vel = self.velocity();
					oss << "SGal3(pose=" << pose << ", velocity=[" << vel.x() << ", " << vel.y() << ", " << vel.z() << "], time=" << self.time() << ")";
					return oss.str();
				})
				.def_static("random", [](py::object seed = py::none()) {
					static std::random_device rd;
					static std::mt19937 gen(rd());
					if (!seed.is_none()) {
						gen.seed(seed.cast<unsigned int>());
					}
					return SGal3<double>::computeRandom(gen);
				}, "Generate random SGal3 element", py::arg("seed") = py::none());
	}

	void moduleAddSim3(py::module &m) {
		// Sim3 bindings for similarity transformations
		py::class_<Sim3<double>>(m, "Sim3", "3D Similarity group Sim(3) - rotations, translations, and scaling")
				.def(py::init<>(), "Default constructor (identity transformation)")
				.def(py::init<const SO3<double> &, const Eigen::Vector3d &, double>(), "Construct from rotation, translation, and scale", py::arg("rotation"), py::arg("translation"), py::arg("scale"))
				.def(py::init<const Eigen::Quaterniond &, const Eigen::Vector3d &, double>(), "Construct from quaternion, translation, and scale", py::arg("quaternion"), py::arg("translation"), py::arg("scale"))
				.def(py::init<const Eigen::Matrix4d &>(), "Construct from 4x4 similarity matrix", py::arg("matrix"))
				.def("__mul__", &Sim3<double>::operator*, "Group composition")
				.def("inverse", &Sim3<double>::inverse, "Compute inverse transformation")
				.def("matrix", &Sim3<double>::matrix, "Get 4x4 transformation matrix")
				.def("rotation", static_cast<const SO3<double> &(Sim3<double>::*) () const>(&Sim3<double>::rotation), "Get rotation part")
				.def("translation", static_cast<const typename Sim3<double>::Vector3 &(Sim3<double>::*) () const>(
											&Sim3<double>::translation), "Get translation part")
				.def("scale", static_cast<const double &(Sim3<double>::*) () const>(&Sim3<double>::scale), "Get scale factor")
				.def_static("exp", &Sim3<double>::exp, "Exponential map from sim(3) to Sim(3)", py::arg("xi"))
				.def("log", &Sim3<double>::log, "Logarithmic map from Sim(3) to sim(3)")
				.def("adjoint", &Sim3<double>::adjoint, "Adjoint representation")
				.def("isApprox", &Sim3<double>::isApprox, "Check approximate equality", py::arg("other"), py::arg("eps") = 1e-12)
				.def_static("identity", &Sim3<double>::Identity, "Get identity element")
				.def("act", [](const Sim3<double>& self, const Eigen::Vector3d& point) {
					return self.act(point);
				}, "Apply similarity transformation to a 3D point", py::arg("point"))
				.def("__repr__", [](const Sim3<double>& self) {
					std::ostringstream oss;
					const auto& rot = self.rotation().quaternion();
					const auto& trans = self.translation();
					oss << "Sim3(quat=[" << rot.w() << ", " << rot.x() << ", " << rot.y() << ", " << rot.z() << "], translation=[" << trans.x() << ", " << trans.y() << ", " << trans.z() << "], scale=" << self.scale() << ")";
					return oss.str();
				})
				.def_static("random", [](py::object seed = py::none()) {
					static std::random_device rd;
					static std::mt19937 gen(rd());
					if (!seed.is_none()) {
						gen.seed(seed.cast<unsigned int>());
					}
					return Sim3<double>::computeRandom(gen);
				}, "Generate random Sim(3) element", py::arg("seed") = py::none());
	}

	void moduleAddUncertainty(py::module &m) {
		// GaussianOnManifold bindings for uncertainty propagation
		py::class_<GaussianOnManifold<SE3<double>>>(m, "GaussianSE3", "Gaussian distribution on SE(3) manifold for uncertainty propagation")
				.def(py::init<const SE3<double> &, const Eigen::Matrix<double, 6, 6> &>(),
					 "Construct from mean and covariance", py::arg("mean"), py::arg("covariance"))
				.def("mean", &GaussianOnManifold<SE3<double>>::mean, "Get the mean element")
				.def("covariance", &GaussianOnManifold<SE3<double>>::covariance, "Get the covariance matrix")
				.def("toGlobalFrame", &GaussianOnManifold<SE3<double>>::toGlobalFrame, "Transform covariance to global frame")
				.def("composeWith", &GaussianOnManifold<SE3<double>>::composeWith,
					 "Compose with tangent space increment", py::arg("delta"), py::arg("delta_cov"))
				.def("mahalanobisDistance", &GaussianOnManifold<SE3<double>>::mahalanobisDistance,
					 "Compute Mahalanobis distance to another distribution", py::arg("other"))
				.def("__repr__", [](const GaussianOnManifold<SE3<double>>& self) {
					std::ostringstream oss;
					oss << "GaussianSE3(mean=" << self.mean() << ")";
					return oss.str();
				});

		py::class_<GaussianOnManifold<SO3<double>>>(m, "GaussianSO3", "Gaussian distribution on SO(3) manifold for uncertainty propagation")
				.def(py::init<const SO3<double> &, const Eigen::Matrix3d &>(),
					 "Construct from mean and covariance", py::arg("mean"), py::arg("covariance"))
				.def("mean", &GaussianOnManifold<SO3<double>>::mean, "Get the mean element")
				.def("covariance", &GaussianOnManifold<SO3<double>>::covariance, "Get the covariance matrix")
				.def("toGlobalFrame", &GaussianOnManifold<SO3<double>>::toGlobalFrame, "Transform covariance to global frame")
				.def("mahalanobisDistance", &GaussianOnManifold<SO3<double>>::mahalanobisDistance,
					 "Compute Mahalanobis distance to another distribution", py::arg("other"))
				.def("__repr__", [](const GaussianOnManifold<SO3<double>>& self) {
					std::ostringstream oss;
					oss << "GaussianSO3(mean=" << self.mean() << ")";
					return oss.str();
				});

		// Uncertainty propagation utilities
		py::class_<UncertaintyPropagation>(m, "UncertaintyPropagation", "Utilities for uncertainty propagation on Lie groups")
				.def_static("isotropic", &UncertaintyPropagation::isotropic<SE3<double>>,
						"Create isotropic Gaussian distribution", py::arg("mean"), py::arg("std_dev"))
				.def_static("diagonal", &UncertaintyPropagation::diagonal<SE3<double>>,
						"Create diagonal covariance Gaussian distribution", py::arg("mean"), py::arg("variances"));
	}

	void moduleAddSE23(py::module &m) {
		// SE23 bindings
		py::class_<SE23<double>>(m, "SE23", "Extended 3D Euclidean group SE_2(3)")
				.def(py::init<>(), "Default constructor (identity transformation)")
				.def(py::init<const SE3<double> &, const Eigen::Vector3d &>(), "Construct from pose and velocity", py::arg("pose"), py::arg("velocity"))
				.def(py::init<const SO3<double> &, const Eigen::Vector3d &, const Eigen::Vector3d &>(), "Construct from rotation, position, and velocity", py::arg("rotation"), py::arg("position"), py::arg("velocity"))
				.def("inverse", &SE23<double>::inverse, "Compute inverse transformation")
				.def("matrix", &SE23<double>::extendedMatrix, "Get 5x5 extended transformation matrix")
				.def("pose", static_cast<const SE3<double> &(SE23<double>::*)() const>(&SE23<double>::pose), "Get pose part")
				.def("velocity", static_cast<const Eigen::Vector3d &(SE23<double>::*)() const>(&SE23<double>::velocity), "Get velocity part")
				.def_static("exp", &SE23<double>::exp, "Exponential map from se_2(3) to SE_2(3)", py::arg("xi"))
				.def("log", &SE23<double>::log, "Logarithmic map from SE_2(3) to se_2(3)")
				.def("adjoint", &SE23<double>::adjoint, "Adjoint representation")
				.def("isApprox", &SE23<double>::isApprox, "Check approximate equality", py::arg("other"), py::arg("eps") = 1e-12)
				.def_static("identity", &SE23<double>::Identity, "Get identity element")
				.def("act", [](const SE23<double>& self, const Eigen::Matrix<double, 6, 1>& point_vel) {
					return self.act(point_vel);
				}, "Apply transformation to a 6D point-velocity vector", py::arg("point_vel"))
				.def("__repr__", [](const SE23<double>& self) {
					std::ostringstream oss;
					const auto& pose = self.pose();
					const auto& vel = self.velocity();
					oss << "SE23(pose=" << pose << ", velocity=[" << vel.x() << ", " << vel.y() << ", " << vel.z() << "])";
					return oss.str();
				})
				.def_static("random", [](py::object seed = py::none()) {
					static std::random_device rd;
					static std::mt19937 gen(rd());
					if (!seed.is_none()) {
						gen.seed(seed.cast<unsigned int>());
					}
					return SE23<double>::computeRandom(gen);
				}, "Generate random SE23 element", py::arg("seed") = py::none());
	}

	void moduleAddBundle(py::module &m) {
		// SE3_Velocity Bundle (SE3 + R^6 velocity)
		using SE3_Velocity_double = Bundle<SE3<double>, RealSpace<double, 6>>;
		py::class_<SE3_Velocity_double>(m, "SE3_Velocity", "SE(3) x R^6 bundle for pose with 6D velocity")
				.def(py::init<>(), "Default constructor (identity)")
				.def(py::init<const SE3<double> &, const RealSpace<double, 6> &>(), "Construct from SE(3) and R^6", py::arg("pose"), py::arg("velocity"))
				.def("__mul__", &SE3_Velocity_double::operator*, "Group composition")
				.def("inverse", &SE3_Velocity_double::inverse, "Compute inverse")
				.def("log", &SE3_Velocity_double::log, "Logarithmic map to algebra")
				.def_static("exp", [](const SE3_Velocity_double::TangentVector& xi) {
					return SE3_Velocity_double::identity().exp(xi);
				}, "Exponential map from algebra", py::arg("xi"))
				.def("adjoint", &SE3_Velocity_double::adjoint, "Adjoint representation")
				.def("isApprox", &SE3_Velocity_double::isApprox, "Check approximate equality", py::arg("other"), py::arg("eps") = 1e-12)
				.def_static("identity", []() { return SE3_Velocity_double::Identity(); }, "Get identity element")
				.def("get_pose", [](const SE3_Velocity_double& self) { return self.get<0>(); }, "Get the SE(3) component")
				.def("get_velocity", [](const SE3_Velocity_double& self) { return self.get<1>(); }, "Get the R^6 component")
				.def("__repr__", [](const SE3_Velocity_double& self) {
					std::ostringstream oss;
					oss << "SE3_Velocity(pose=" << self.get<0>() << ", velocity=" << self.get<1>() << ")";
					return oss.str();
				});

		// SE2_Velocity Bundle (SE2 + R^3 velocity)
		using SE2_Velocity_double = Bundle<SE2<double>, RealSpace<double, 3>>;
		py::class_<SE2_Velocity_double>(m, "SE2_Velocity", "SE(2) x R^3 bundle for 2D pose with 3D velocity")
				.def(py::init<>(), "Default constructor (identity)")
				.def(py::init<const SE2<double> &, const RealSpace<double, 3> &>(), "Construct from SE(2) and R^3", py::arg("pose"), py::arg("velocity"))
				.def("__mul__", &SE2_Velocity_double::operator*, "Group composition")
				.def("inverse", &SE2_Velocity_double::inverse, "Compute inverse")
				.def("log", &SE2_Velocity_double::log, "Logarithmic map to algebra")
				.def_static("exp", [](const SE2_Velocity_double::TangentVector& xi) {
					return SE2_Velocity_double::identity().exp(xi);
				}, "Exponential map from algebra", py::arg("xi"))
				.def("adjoint", &SE2_Velocity_double::adjoint, "Adjoint representation")
				.def("isApprox", &SE2_Velocity_double::isApprox, "Check approximate equality", py::arg("other"), py::arg("eps") = 1e-12)
				.def_static("identity", []() { return SE2_Velocity_double::Identity(); }, "Get identity element")
				.def("get_pose", [](const SE2_Velocity_double& self) { return self.get<0>(); }, "Get the SE(2) component")
				.def("get_velocity", [](const SE2_Velocity_double& self) { return self.get<1>(); }, "Get the R^3 component")
				.def("__repr__", [](const SE2_Velocity_double& self) {
					std::ostringstream oss;
					oss << "SE2_Velocity(pose=" << self.get<0>() << ", velocity=" << self.get<1>() << ")";
					return oss.str();
				});
	}

	void moduleAddRealSpace(py::module &m) {
		// RealSpace bindings
		py::class_<RealSpace<double, 3>>(m, "R3", "3D real space R^3")
				.def(py::init<>(), "Default constructor (zero vector)")
				.def(py::init<const Eigen::Vector3d &>(), "Construct from 3D vector", py::arg("vector"))
				.def("__add__", [](const RealSpace<double, 3>& self, const RealSpace<double, 3>& other) {
					return self + other;
				}, "Vector addition")
				.def("__neg__", [](const RealSpace<double, 3>& self) {
					return -self;
				}, "Vector negation")
				.def("vector", [](const RealSpace<double, 3>& self) {
					return self.computeLog();
				}, "Get the underlying vector")
				.def_static("hat", &RealSpace<double, 3>::hat, "Hat operator: R^3 -> 3x3 matrix", py::arg("v"))
				.def_static("vee", &RealSpace<double, 3>::vee, "Vee operator: 3x3 matrix -> R^3", py::arg("X"));

		py::class_<RealSpace<double, 6>>(m, "R6", "6D real space R^6")
				.def(py::init<>(), "Default constructor (zero vector)")
				.def(py::init<const Eigen::Matrix<double, 6, 1> &>(), "Construct from 6D vector", py::arg("vector"))
				.def("__add__", [](const RealSpace<double, 6>& self, const RealSpace<double, 6>& other) {
					return self + other;
				}, "Vector addition")
				.def("__neg__", [](const RealSpace<double, 6>& self) {
					return -self;
				}, "Vector negation")
				.def("vector", [](const RealSpace<double, 6>& self) {
					return self.computeLog();
				}, "Get the underlying vector")
				.def_static("hat", &RealSpace<double, 6>::hat, "Hat operator: R^6 -> 6x6 matrix", py::arg("v"))
				.def_static("vee", &RealSpace<double, 6>::vee, "Vee operator: 6x6 matrix -> R^6", py::arg("X"));
	}

	void moduleAddLieGroupUtils(py::module &m) {
		// Utility functions for interpolation, etc.
		// Note: slerp function would need to be implemented in the Lie group classes
		// m.def("slerp", [](const SO3<double>& a, const SO3<double>& b, double t) {
		//     return a.interpolate(b, t);
		// }, "Spherical linear interpolation for SO3");
	}

	void moduleAddLieGroups(py::module &m) {
		// Add all Lie group bindings
		moduleAddSO2(m);
		moduleAddSO3(m);
		moduleAddSE2(m);
		moduleAddSE3(m);
		moduleAddSGal3(m);
		moduleAddSE23(m);
		moduleAddRealSpace(m);
		moduleAddBundle(m);
		moduleAddLieGroupUtils(m);
	}

} // namespace sofapython3
