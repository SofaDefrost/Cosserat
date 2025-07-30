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
#include "../../../../src/liegroups/Types.h"
#include "../../../../src/liegroups/LieGroupBase.h"
#include "../../../../src/liegroups/LieGroupBase.inl"
#include "../../../../src/liegroups/RealSpace.h"
#include "../../../../src/liegroups/SO2.h"
#include "../../../../src/liegroups/SO3.h"
#include "../../../../src/liegroups/SE2.h"
#include "../../../../src/liegroups/SE3.h"
#include "../../../../src/liegroups/SGal3.h"

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

	/**
	 * @brief Implementation of product manifold bundle of Lie groups
	 *
	 * This class implements a Cartesian product of multiple Lie groups, allowing
	 * them to be treated as a single composite Lie group. The bundle maintains the
	 * product structure while providing all necessary group operations.
	 *
	 * Mathematical Background:
	 * For Lie groups G₁, G₂, ..., Gₙ, the product manifold G = G₁ × G₂ × ... × Gₙ
	 * is also a Lie group with:
	 * - Group operation: (g₁, g₂, ..., gₙ) * (h₁, h₂, ..., hₙ) = (g₁h₁, g₂h₂, ...,
	 * gₙhₙ)
	 * - Identity: (e₁, e₂, ..., eₙ) where eᵢ is identity of Gᵢ
	 * - Inverse: (g₁, g₂, ..., gₙ)⁻¹ = (g₁⁻¹, g₂⁻¹, ..., gₙ⁻¹)
	 * - Lie algebra: 𝔤 = 𝔤₁ ⊕ 𝔤₂ ⊕ ... ⊕ 𝔤ₙ (direct sum)
	 * - Adjoint: block diagonal with Adᵢ on diagonal blocks
	 *
	 * Usage Examples:
	 * ```cpp
	 * // Rigid body pose with velocity
	 * using RigidBodyState = Bundle<SE3<double>, RealSpace<double, 6>>;
	 *
	 * // 2D robot with multiple joints
	 * using Robot2D = Bundle<SE2<double>, SO2<double>, SO2<double>>;
	 *
	 * // Create and manipulate
	 * auto state1 = RigidBodyState(SE3<double>::identity(),
	 *                             RealSpace<double, 6>::zero());
	 * auto state2 = RigidBodyState(pose, velocity);
	 * auto combined = state1 * state2;
	 * ```
	 *
	 * @tparam Groups The Lie groups to bundle together (must have same scalar type)
	 */
	template<typename... Groups>
	class Bundle
		: public LieGroupBase<typename std::tuple_element_t<0, std::tuple<Groups...>>::Scalar,
							  std::integral_constant<int, detail::TotalDimension<Groups...>::value>,
							  detail::TotalDimension<Groups...>::value, detail::TotalDimension<Groups...>::value> {

		// Compile-time validation
		static_assert(sizeof...(Groups) > 0, "Bundle must contain at least one group");

		using FirstGroup = std::tuple_element_t<0, std::tuple<Groups...>>;
		using FirstScalar = typename FirstGroup::Scalar;

		static_assert(detail::AllAreLieGroups<FirstScalar, Groups...>::value,
					  "All template parameters must be Lie groups with the same scalar type");

	public:
		using Base = LieGroupBase<FirstScalar, std::integral_constant<int, detail::TotalDimension<Groups...>::value>,
								  detail::TotalDimension<Groups...>::value, detail::TotalDimension<Groups...>::value>;
		using Scalar = typename Base::Scalar;
		using Vector = typename Base::Vector;
		using Matrix = typename Base::Matrix;
		using TangentVector = typename Base::TangentVector;
		using AdjointMatrix = typename Base::AdjointMatrix;

		static constexpr int Dim = Base::Dim;
		static constexpr std::size_t NumGroups = sizeof...(Groups);

		using GroupTuple = std::tuple<Groups...>;

		// Compile-time offset table for algebra elements
		template<std::size_t I>
		static constexpr int AlgebraOffset = detail::OffsetAt<I, Groups...>::value;

		template<std::size_t I>
		using GroupType = std::tuple_element_t<I, GroupTuple>;

		// ========== Constructors ==========

		/**
		 * @brief Default constructor creates identity bundle
		 */
		Bundle() : m_groups(), m_action_offsets(m_groups) {}

		/**
		 * @brief Construct from individual group elements
		 */
		explicit Bundle(const Groups &...groups) : m_groups(groups...), m_action_offsets(m_groups) {}

		/**
		 * @brief Construct from tuple of groups
		 */
		explicit Bundle(const GroupTuple &groups) : m_groups(groups), m_action_offsets(m_groups) {}

		/**
		 * @brief Construct from Lie algebra vector
		 */
		explicit Bundle(const TangentVector &algebra_element) : Bundle() { *this = exp(algebra_element); }

		/**
		 * @brief Copy constructor
		 */
		Bundle(const Bundle &other) = default;

		/**
		 * @brief Move constructor
		 */
		Bundle(Bundle &&other) noexcept = default;

		/**
		 * @brief Copy assignment
		 */
		Bundle &operator=(const Bundle &other) = default;

		/**
		 * @brief Move assignment
		 */
		Bundle &operator=(Bundle &&other) noexcept = default;

		// ========== Group Operations ==========

		/**
		 * @brief Group composition (component-wise)
		 * Implements: (g₁, ..., gₙ) * (h₁, ..., hₙ) = (g₁h₁, ..., gₙhₙ)
		 */
		Bundle operator*(const Bundle &other) const {
			return multiply_impl(other, std::index_sequence_for<Groups...>());
		}

		/**
		 * @brief In-place group composition
		 */
		Bundle &operator*=(const Bundle &other) {
			multiply_assign_impl(other, std::index_sequence_for<Groups...>());
			return *this;
		}

		/**
		 * @brief Inverse element (component-wise)
		 * Implements: (g₁, ..., gₙ)⁻¹ = (g₁⁻¹, ..., gₙ⁻¹)
		 */
		Bundle inverse() const { return inverse_impl(std::index_sequence_for<Groups...>()); }

		// ========== Lie Algebra Operations ==========

		/**
		 * @brief Exponential map from Lie algebra to bundle
		 * The Lie algebra of the product is the direct sum of individual algebras
		 */
		Bundle exp(const TangentVector &algebra_element) const {
			validateAlgebraElement(algebra_element);
			return exp_impl(algebra_element, std::index_sequence_for<Groups...>());
		}

		/**
		 * @brief Logarithmic map from bundle to Lie algebra
		 * Maps to the direct sum of individual Lie algebras
		 */
		TangentVector log() const { return log_impl(std::index_sequence_for<Groups...>()); }

		/**
		 * @brief Adjoint representation (block diagonal structure)
		 * Ad_{(g₁,...,gₙ)} = diag(Ad_{g₁}, ..., Ad_{gₙ})
		 */
		AdjointMatrix adjoint() const { return adjoint_impl(std::index_sequence_for<Groups...>()); }

		// ========== Group Actions ==========

		/**
		 * @brief Group action on a point (component-wise on appropriate subspaces)
		 * Each group acts on its corresponding portion of the input vector
		 */
		Vector act(const Vector &point) const {
			validateActionInput(point);
			return act_impl(point, std::index_sequence_for<Groups...>());
		}

		/**
		 * @brief Batch group action on multiple points
		 */
		template<int N>
		Eigen::Matrix<Scalar, Eigen::Dynamic, N> act(const Eigen::Matrix<Scalar, Eigen::Dynamic, N> &points) const {

			if (points.rows() != actionDimension()) {
				throw std::invalid_argument("Point matrix has wrong dimension");
			}

			Eigen::Matrix<Scalar, Eigen::Dynamic, N> result(actionDimension(), N);

			for (int i = 0; i < N; ++i) {
				result.col(i) = act(points.col(i));
			}

			return result;
		}

		// ========== Utility Functions ==========

		/**
		 * @brief Check if approximately equal to another bundle
		 */
		bool isApprox(const Bundle &other, const Scalar &eps = Types<Scalar>::epsilon()) const {
			return isApprox_impl(other, eps, std::index_sequence_for<Groups...>());
		}

		/**
		 * @brief Equality operator
		 */
		bool operator==(const Bundle &other) const { return isApprox(other); }

		/**
		 * @brief Inequality operator
		 */
		bool operator!=(const Bundle &other) const { return !(*this == other); }

		/**
		 * @brief Linear interpolation between two bundles (geodesic in product space)
		 * @param other Target bundle
		 * @param t Interpolation parameter [0,1]
		 */
		Bundle interpolate(const Bundle &other, const Scalar &t) const {
			if (t < 0 || t > 1) {
				throw std::invalid_argument("Interpolation parameter must be in [0,1]");
			}

			TangentVector delta = (inverse() * other).log();
			return *this * exp(t * delta);
		}

		/**
		 * @brief Generate random bundle element
		 */
		template<typename Generator>
		static Bundle random(Generator &gen, const Scalar &scale = Scalar(1)) {
			return random_impl(gen, scale, std::index_sequence_for<Groups...>());
		}

		// ========== Accessors ==========

		/**
		 * @brief Get the identity element
		 */
		static const Bundle &identity() {
			static const Bundle id;
			return id;
		}

		/**
		 * @brief Get the dimension of the Lie algebra
		 */
		int algebraDimension() const { return Dim; }

		/**
		 * @brief Get the dimension of the space the group acts on (computed at
		 * runtime)
		 */
		int actionDimension() const { return m_action_offsets.total(); }

		/**
		 * @brief Access individual group elements (const)
		 */
		template<std::size_t I>
		const auto &get() const {
			static_assert(I < NumGroups, "Index out of bounds");
			return std::get<I>(m_groups);
		}

		/**
		 * @brief Access individual group elements (mutable)
		 */
		template<std::size_t I>
		auto &get() {
			static_assert(I < NumGroups, "Index out of bounds");
			// Need to recompute action offsets if groups are modified
			auto &result = std::get<I>(m_groups);
			// In practice, you might want to make this const and provide setters
			return result;
		}

		/**
		 * @brief Set individual group element
		 */
		template<std::size_t I>
		void set(const GroupType<I> &group) {
			std::get<I>(m_groups) = group;
			m_action_offsets = detail::ActionOffsets<Groups...>(m_groups);
		}

		/**
		 * @brief Get the underlying tuple
		 */
		const GroupTuple &groups() const { return m_groups; }

		/**
		 * @brief Get algebra element for specific group
		 */
		template<std::size_t I>
		typename GroupType<I>::TangentVector getAlgebraElement() const {
			return std::get<I>(m_groups).log();
		}

		/**
		 * @brief Set from algebra element for specific group
		 */
		template<std::size_t I>
		void setFromAlgebra(const typename GroupType<I>::TangentVector &algebra) {
			std::get<I>(m_groups) = GroupType<I>().exp(algebra);
			m_action_offsets = detail::ActionOffsets<Groups...>(m_groups);
		}

		// ========== Stream Output ==========

		/**
		 * @brief Output stream operator
		 */
		friend std::ostream &operator<<(std::ostream &os, const Bundle &bundle) {
			os << "Bundle<" << NumGroups << ">(";
			bundle.print_impl(os, std::index_sequence_for<Groups...>());
			os << ")";
			return os;
		}

	private:
		GroupTuple m_groups; ///< Tuple of group elements
		detail::ActionOffsets<Groups...> m_action_offsets; ///< Cached action dimension offsets

		// ========== Implementation Helpers ==========

		// Validation helpers
		void validateAlgebraElement(const TangentVector &element) const {
			if (element.size() != Dim) {
				throw std::invalid_argument("Algebra element has wrong dimension");
			}
		}

		void validateActionInput(const Vector &point) const {
			if (point.size() != actionDimension()) {
				throw std::invalid_argument("Action input has wrong dimension");
			}
		}

		// Multiplication implementation
		template<std::size_t... Is>
		Bundle multiply_impl(const Bundle &other, std::index_sequence<Is...>) const {
			return Bundle((std::get<Is>(m_groups) * std::get<Is>(other.m_groups))...);
		}

		template<std::size_t... Is>
		void multiply_assign_impl(const Bundle &other, std::index_sequence<Is...>) {
			((std::get<Is>(m_groups) *= std::get<Is>(other.m_groups)), ...);
		}

		// Inverse implementation
		template<std::size_t... Is>
		Bundle inverse_impl(std::index_sequence<Is...>) const {
			return Bundle((std::get<Is>(m_groups).inverse())...);
		}

		// Exponential map implementation
		template<std::size_t... Is>
		Bundle exp_impl(const TangentVector &algebra_element, std::index_sequence<Is...>) const {
			return Bundle(
					(GroupType<Is>().exp(algebra_element.template segment<GroupType<Is>::Dim>(AlgebraOffset<Is>)))...);
		}

		// Logarithmic map implementation
		template<std::size_t... Is>
		TangentVector log_impl(std::index_sequence<Is...>) const {
			TangentVector result;
			((result.template segment<GroupType<Is>::Dim>(AlgebraOffset<Is>) = std::get<Is>(m_groups).log()), ...);
			return result;
		}

		// Adjoint implementation (block diagonal)
		template<std::size_t... Is>
		AdjointMatrix adjoint_impl(std::index_sequence<Is...>) const {
			AdjointMatrix result = AdjointMatrix::Zero();
			((result.template block<GroupType<Is>::Dim, GroupType<Is>::Dim>(AlgebraOffset<Is>, AlgebraOffset<Is>) =
					  std::get<Is>(m_groups).adjoint()),
			 ...);
			return result;
		}

		// Group action implementation
		template<std::size_t... Is>
		Vector act_impl(const Vector &point, std::index_sequence<Is...>) const {
			Vector result(actionDimension());

			// Apply each group's action to its corresponding subspace
			((applyGroupAction<Is>(result, point)), ...);

			return result;
		}

		template<std::size_t I>
		void applyGroupAction(Vector &result, const Vector &point) const {
			const auto &group = std::get<I>(m_groups);
			int in_offset = m_action_offsets[I];
			int out_offset = in_offset; // Same offset for output
			int dim = group.actionDimension();

			auto input_segment = point.segment(in_offset, dim);
			auto output_segment = result.segment(out_offset, dim);
			output_segment = group.act(input_segment);
		}

		// Approximate equality implementation
		template<std::size_t... Is>
		bool isApprox_impl(const Bundle &other, const Scalar &eps, std::index_sequence<Is...>) const {
			return (std::get<Is>(m_groups).isApprox(std::get<Is>(other.m_groups), eps) && ...);
		}

		// Random generation implementation
		template<typename Generator, std::size_t... Is>
		static Bundle random_impl(Generator &gen, const Scalar &scale, std::index_sequence<Is...>) {
			return Bundle((GroupType<Is>::random(gen, scale))...);
		}

		// Stream output implementation
		template<std::size_t... Is>
		void print_impl(std::ostream &os, std::index_sequence<Is...>) const {
			bool first = true;
			((os << (first ? (first = false, "") : ", ") << std::get<Is>(m_groups)), ...);
		}
	};

	// ========== Type Aliases ==========

	// Common bundles for robotics applications
	template<typename Scalar>
	using SE3_Velocity = Bundle<SE3<Scalar>, RealSpace<Scalar, 6>>;

	template<typename Scalar>
	using SE2_Velocity = Bundle<SE2<Scalar>, RealSpace<Scalar, 3>>;

	template<typename Scalar, int N>
	using SE3_Joints = Bundle<SE3<Scalar>, RealSpace<Scalar, N>>;

	// Convenience aliases
	template<typename... Groups>
	using Bundlef = Bundle<Groups...>;

	template<typename... Groups>
	using Bundled = Bundle<Groups...>;

} // namespace sofa::component::cosserat::liegroups

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
		py::class_<SO2<double>>(m, "SO2")
				.def(py::init<>())
				.def(py::init<double>())
				.def("__mul__", &SO2<double>::operator*)
				.def("inverse", &SO2<double>::inverse)
				.def("matrix", &SO2<double>::matrix)
				.def("angle", &SO2<double>::angle)
				.def("exp", &SO2<double>::exp)
				.def("log", &SO2<double>::log)
				.def("adjoint", &SO2<double>::adjoint)
				.def("isApprox", &SO2<double>::isApprox)
				.def_static("identity", &SO2<double>::identity)
				.def_static("hat", &SO2<double>::hat)
				.def("act", static_cast<typename SO2<double>::Vector (SO2<double>::*)(
									const typename SO2<double>::Vector &) const>(&SO2<double>::act));
	}

	void moduleAddSO3(py::module &m) {
		// SO3 bindings
		py::class_<SO3<double>>(m, "SO3")
				.def(py::init<>())
				.def(py::init<double, const Eigen::Vector3d &>())
				.def(py::init<const Eigen::Quaterniond &>())
				.def(py::init<const Eigen::Matrix3d &>())
				.def("__mul__", &SO3<double>::operator*)
				.def("inverse", &SO3<double>::inverse)
				.def("matrix", &SO3<double>::matrix)
				.def("quaternion", &SO3<double>::quaternion)
				.def("exp", &SO3<double>::exp)
				.def("log", &SO3<double>::log)
				.def("adjoint", &SO3<double>::adjoint)
				.def("isApprox", &SO3<double>::isApprox)
				.def_static("identity", &SO3<double>::identity)
				.def_static("hat", &SO3<double>::hat)
				.def_static("vee", &SO3<double>::vee)
				.def("act", static_cast<typename SO3<double>::Vector (SO3<double>::*)(
									const typename SO3<double>::Vector &) const>(&SO3<double>::computeAction));
	}

	void moduleAddSE2(py::module &m) {
		// SE2 bindings
		py::class_<SE2<double>>(m, "SE2")
				.def(py::init<>())
				.def(py::init<const SO2<double> &, const Eigen::Vector2d &>())
				.def("__mul__", &SE2<double>::operator*)
				.def("inverse", &SE2<double>::inverse)
				.def("matrix", &SE2<double>::matrix)
				.def("rotation", static_cast<const SO2<double> &(SE2<double>::*) () const>(&SE2<double>::rotation))
				.def("translation", static_cast<const typename SE2<double>::Vector2 &(SE2<double>::*) () const>(
											&SE2<double>::translation))
				.def("exp", &SE2<double>::exp)
				.def("log", &SE2<double>::log)
				.def("adjoint", &SE2<double>::adjoint)
				.def("isApprox", &SE2<double>::isApprox)
				.def_static("identity", &SE2<double>::identity)
				.def("act", static_cast<typename SE2<double>::Vector2 (SE2<double>::*)(
									const typename SE2<double>::Vector2 &) const>(&SE2<double>::act));
	}

	void moduleAddSE3(py::module &m) {
		// SE3 bindings with enhanced functionality
		py::class_<SE3<double>>(m, "SE3")
				.def(py::init<>())
				.def(py::init<const SO3<double> &, const Eigen::Vector3d &>())
				.def(py::init<const Eigen::Matrix4d &>())
				.def("__mul__", &SE3<double>::operator*)
				.def("inverse", &SE3<double>::inverse)
				.def("matrix", &SE3<double>::matrix)
				.def("rotation", static_cast<const SO3<double> &(SE3<double>::*) () const>(&SE3<double>::rotation))
				.def("translation", static_cast<const typename SE3<double>::Vector3 &(SE3<double>::*) () const>(
											&SE3<double>::translation))
				.def("exp", &SE3<double>::exp)
				.def("log", &SE3<double>::log)
				.def("adjoint", &SE3<double>::adjoint)
				.def("isApprox", &SE3<double>::isApprox)
				.def_static("identity", &SE3<double>::Identity)
				// Add the hat operator (maps 6D vector to 4x4 matrix)
				.def_static(
						"hat", [](const Eigen::Matrix<double, 6, 1> &xi) { return dualMatrix<double>(xi); },
						"Map 6D vector to 4x4 matrix representation (hat operator)")
				// Add co-adjoint (transpose of adjoint)
				.def(
						"co_adjoint", [](const SE3<double> &self) { return self.adjoint().transpose(); },
						"Co-adjoint representation (transpose of adjoint)")
				.def(
						"coadjoint", [](const SE3<double> &self) { return self.adjoint().transpose(); },
						"Co-adjoint representation (alias for co_adjoint)")
				.def("act", static_cast<typename SE3<double>::Vector3 (SE3<double>::*)(
									const typename SE3<double>::Vector3 &) const>(&SE3<double>::act))
				// Baker-Campbell-Hausdorff formula
				.def_static("BCH", &SE3<double>::BCH, "Baker-Campbell-Hausdorff formula");
	}

	void moduleAddSGal3(py::module &m) {
		// SGal3 bindings (placeholder for now)
		// Implementation depends on the actual SGal3 class structure
	}

	void moduleAddSE23(py::module &m) {
		// SE23 bindings (placeholder for now)
		// Implementation depends on the actual SE23 class structure
	}

	void moduleAddBundle(py::module &m) {
		// Bundle bindings (placeholder for now)
		// This would require template instantiation for specific Bundle types
		// For example: Bundle<SE3<double>, RealSpace<double, 6>>
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
		moduleAddBundle(m);
		moduleAddLieGroupUtils(m);
	}

} // namespace sofapython3

