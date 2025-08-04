// This file defines the base class template for all Lie group implementations,
// using the Curiously Recurring Template Pattern (CRTP).

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
 * along with this program. If not, see <http://www.gnu.org/licenses/>.        *
 ******************************************************************************/

#pragma once

#include <array>
#include <concepts>
#include <iostream>
#include <memory>
#include <random>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include "Types.h"

namespace sofa::component::cosserat::liegroups {

	// Forward declaration for matrix representation of Lie algebra
	template<typename _Scalar, int _Dim>
	class LieAlgebra;

	// Modern C++20 concepts for better error messages and type safety
	template<typename T>
	concept FloatingPoint = std::is_floating_point_v<T>;

	template<typename T, typename Scalar>
	concept HasScalarType = std::same_as<typename T::Scalar, Scalar>;

	/**
	 * @brief Exception hierarchy for Lie group operations
	 */
	struct LieGroupException : public std::runtime_error {
		using std::runtime_error::runtime_error;
	};

	struct InverseFailureException : public LieGroupException {
		InverseFailureException() : LieGroupException("Failed to compute inverse: singular element") {}
	};

	struct LogarithmException : public LieGroupException {
		LogarithmException() : LieGroupException("Logarithm undefined at this point") {}
	};

	struct NumericalInstabilityException : public LieGroupException {
		NumericalInstabilityException(const std::string &msg) : LieGroupException("Numerical instability: " + msg) {}
	};

	/**
	 * @brief Base class template for all Lie group implementations using CRTP
	 *
	 * This class defines the interface that all Lie group implementations must
	 * satisfy. It uses the Curiously Recurring Template Pattern (CRTP) to allow
	 * static polymorphism, providing better performance than virtual functions.
	 *
	 * @tparam Derived The derived class implementing the specific Lie group
	 * @tparam _Scalar The scalar type used for computations (must be a
	 * floating-point type)
	 * @tparam _Dim The dimension of the group representation
	 * @tparam _AlgebraDim The dimension of the Lie algebra (default: same as _Dim)
	 * @tparam _ActionDim The dimension of vectors the group acts on (default: same
	 * as _Dim)
	 */
	template<typename Derived, FloatingPoint _Scalar, int _Dim, int _AlgebraDim = _Dim, int _ActionDim = _Dim>
		requires(_Dim > 0 && _AlgebraDim > 0 && _ActionDim > 0)
	class LieGroupBase {
	public:
		// Type aliases for scalar and types
		using Scalar = _Scalar;
		using Types = Types<Scalar>;

		// Dimensions as static constants
		static constexpr int Dim = _Dim;
		static constexpr int AlgebraDim = _AlgebraDim;
		static constexpr int ActionDim = _ActionDim;

		// Define commonly used types
		using Vector = typename Types::template Vector<Dim>;
		using Matrix = typename Types::template Matrix<Dim, Dim>;

		using TangentVector = typename Types::template Vector<AlgebraDim>;
		using AdjointMatrix = typename Types::template Matrix<AlgebraDim, AlgebraDim>;
		using AlgebraMatrix = typename Types::template Matrix<Dim, Dim>;

		using ActionVector = typename Types::template Vector<ActionDim>;
		using JacobianMatrix = typename Types::template Matrix<ActionDim, AlgebraDim>;

		using DerivedType = Derived;

		// Type aliases for common patterns
		using value_type = Scalar;
		using tangent_vector_type = TangentVector;
		using adjoint_matrix_type = AdjointMatrix;

	protected:
		/**
		 * @brief Protected constructor to prevent direct instantiation
		 * Only derived classes can construct this base
		 */
		LieGroupBase() = default;

	public:
		// Rule of 5 - modern C++ best practices
		LieGroupBase(const LieGroupBase &) = default;
		LieGroupBase(LieGroupBase &&) noexcept = default;
		LieGroupBase &operator=(const LieGroupBase &) = default;
		LieGroupBase &operator=(LieGroupBase &&) noexcept = default;

		/**
		 * @brief Access to the derived object (const version)
		 * @return Reference to the derived object
		 */
		[[nodiscard]] constexpr const Derived &derived() const noexcept { return static_cast<const Derived &>(*this); }

		/**
		 * @brief Access to the derived object (non-const version)
		 * @return Reference to the derived object
		 */
		[[nodiscard]] constexpr Derived &derived() noexcept { return static_cast<Derived &>(*this); }

		/**
		 * @brief In-place group composition with perfect forwarding
		 * @param other Another element of the same Lie group
		 * @return Reference to this after composition
		 */
		template<typename T>
			requires std::same_as<std::decay_t<T>, Derived>
		constexpr Derived &operator*=(T &&other) noexcept {
			derived() = derived().compose(std::forward<T>(other));
			return derived();
		}

		/**
		 * @brief Group composition operation
		 * @param other Another element of the same Lie group
		 * @return The composition this * other
		 */
		[[nodiscard]] constexpr Derived operator*(const Derived &other) const noexcept {
			return derived().compose(other);
		}

		/**
		 * @brief Compute the inverse element
		 * @return The inverse element such that this * inverse() = identity()
		 * @throws InverseFailureException if the element is singular (not invertible)
		 */
		[[nodiscard]] Derived inverse() const { return derived().computeInverse(); }

		/**
		 * @brief Exponential map from Lie algebra to Lie group
		 * @param algebra_element Element of the Lie algebra (tangent space at
		 * identity)
		 * @return The corresponding element in the Lie group
		 */
		[[nodiscard]] static constexpr Derived exp(const TangentVector &algebra_element) noexcept {
			return Derived::computeExp(algebra_element);
		}

		/**
		 * @brief Logarithmic map from Lie group to Lie algebra
		 * @return The corresponding element in the Lie algebra
		 * @throws LogarithmException if the element is outside the domain of the log
		 * map
		 */
		[[nodiscard]] TangentVector log() const { return derived().computeLog(); }

		/**
		 * @brief Adjoint representation of the group element
		 * @return The adjoint matrix representing the action on the Lie algebra
		 */
		[[nodiscard]] AdjointMatrix adjoint() const noexcept { return derived().computeAdjoint(); }

		/**
		 * @brief Group action on a point
		 * @param point The point to transform
		 * @return The transformed point
		 */
		[[nodiscard]] ActionVector act(const ActionVector &point) const noexcept {
			return derived().computeAction(point);
		}

		/**
		 * @brief Overloaded function call operator for group action
		 * @param point The point to transform
		 * @return The transformed point
		 */
		[[nodiscard]] ActionVector operator()(const ActionVector &point) const noexcept { return act(point); }

		/**
		 * @brief Check if this element is approximately equal to another
		 * @param other Another element of the same Lie group
		 * @param eps Tolerance for comparison (optional)
		 * @return true if elements are approximately equal
		 */
		[[nodiscard]] bool isApprox(const Derived &other, const Scalar &eps = Types::epsilon()) const noexcept {
			return derived().computeIsApprox(other, eps);
		}

		/**
		 * @brief Get the identity element of the group
		 * @return The identity element
		 */
		[[nodiscard]] static constexpr Derived Identity() noexcept { return Derived::computeIdentity(); }

		/**
		 * @brief Compute distance between two group elements
		 * @param other Another element of the same Lie group
		 * @return A scalar representing the distance
		 */
		[[nodiscard]] Scalar distance(const Derived &other) const noexcept;

		/**
		 * @brief Interpolate between two group elements
		 * @param other Target group element
		 * @param t Interpolation parameter between 0 and 1
		 * @return Interpolated group element
		 */
		[[nodiscard]] Derived interpolate(const Derived &other, const Scalar &t) const noexcept;

		/**
		 * @brief Linear interpolation (alias for interpolate)
		 * @param other Target group element
		 * @param t Interpolation parameter between 0 and 1
		 * @return Interpolated group element
		 */
		[[nodiscard]] Derived lerp(const Derived &other, const Scalar &t) const noexcept {
			return interpolate(other, t);
		}

		/**
		 * @brief Spherical linear interpolation with better numerical properties
		 * @param other Target group element
		 * @param t Interpolation parameter between 0 and 1
		 * @return Interpolated group element
		 */
		[[nodiscard]] Derived slerp(const Derived &other, const Scalar &t) const noexcept;

		/**
		 * @brief Get the Jacobian of the group action
		 * @param point The point at which to compute the Jacobian
		 * @return The Jacobian matrix
		 */
		[[nodiscard]] JacobianMatrix actionJacobian(const ActionVector &point) const noexcept;

		/**
		 * @brief Convert a tangent vector to a Lie algebra matrix representation (hat
		 * operator)
		 * @param v Vector representation of Lie algebra element
		 * @return Matrix representation in the Lie algebra
		 */
		[[nodiscard]] static AlgebraMatrix hat(const TangentVector &v) noexcept;

		/**
		 * @brief Convert a Lie algebra matrix to its vector representation (vee
		 * operator)
		 * @param X Matrix representation in the Lie algebra
		 * @return Vector representation of the Lie algebra element
		 */
		[[nodiscard]] static TangentVector vee(const AlgebraMatrix &X) noexcept;

		/**
		 * @brief Baker-Campbell-Hausdorff formula for Lie algebra elements
		 * @param v First tangent vector
		 * @param w Second tangent vector
		 * @param order Order of approximation (1-5, default: 3)
		 * @return Tangent vector approximating log(exp(v)*exp(w))
		 */
		[[nodiscard]] static TangentVector BCH(const TangentVector &v, const TangentVector &w, int order = 3);

		/**
		 * @brief Get the differential of the exponential map
		 * @param v Tangent vector
		 * @return Matrix representing the differential of exp at v
		 */
		[[nodiscard]] static AdjointMatrix dexp(const TangentVector &v);

		/**
		 * @brief Get the inverse of the differential of the exponential map
		 * @param v Tangent vector
		 * @return Matrix representing the inverse differential of exp at v
		 */
		[[nodiscard]] static AdjointMatrix dexpInv(const TangentVector &v);

		/**
		 * @brief Get the differential of the logarithm map
		 * @return Matrix representing the differential of log at the current point
		 */
		[[nodiscard]] AdjointMatrix dlog() const;

		/**
		 * @brief Create a group element from parameters with perfect forwarding
		 * @tparam Args Parameter types
		 * @param args Constructor arguments
		 * @return New group element
		 */
		template<typename... Args>
		[[nodiscard]] static constexpr Derived create(Args &&...args) {
			return Derived(std::forward<Args>(args)...);
		}

		/**
		 * @brief Compute the adjoint representation of a Lie algebra element
		 *
		 * For matrix Lie algebras, this computes the matrix representation of
		 * the adjoint action ad_X(Y) = [X,Y] = XY - YX
		 *
		 * @param v Element of the Lie algebra in vector form
		 * @return Matrix representing the adjoint action
		 */
		[[nodiscard]] static AdjointMatrix ad(const TangentVector &v);

		/**
		 * @brief Compute the group manifold dimension
		 * @return The dimension of the manifold
		 */
		[[nodiscard]] static constexpr int manifoldDim() noexcept { return AlgebraDim; }

		/**
		 * @brief Check if the current element is close to the identity
		 * @param eps Tolerance for comparison
		 * @return true if close to identity
		 */
		[[nodiscard]] bool isNearIdentity(const Scalar &eps = Types::epsilon()) const noexcept {
			return isApprox(Identity(), eps);
		}

		/**
		 * @brief Compute the matrix logarithm with improved numerical stability
		 * @return Matrix logarithm of the current element
		 */
		[[nodiscard]] AlgebraMatrix matrixLog() const;

		/**
		 * @brief Get a random element of the group (for testing/sampling)
		 * @param generator Random number generator
		 * @return Random group element
		 */
		template<typename Generator>
		[[nodiscard]] static Derived Random(Generator &gen) {
			return Derived::computeRandom(gen);
		}

		/**
		 * @brief Stream output operator
		 */
		friend std::ostream &operator<<(std::ostream &os, const LieGroupBase &g) { return g.derived().print(os); }

		/**
		 * @brief Get type information string for debugging
		 */
		[[nodiscard]] static constexpr std::string_view typeName() noexcept { return Derived::getTypeName(); }

		/**
		 * @brief Validate that the current element is a valid group element
		 * @return true if valid, false otherwise
		 */
		[[nodiscard]] bool isValid() const noexcept { return derived().computeIsValid(); }

		/**
		 * @brief Normalize the group element to ensure it remains on the manifold
		 * This is useful for numerical stability after many operations
		 */
		void normalize() noexcept { derived().computeNormalize(); }

		/**
		 * @brief Get normalized copy of the group element
		 * @return Normalized copy
		 */
		[[nodiscard]] Derived normalized() const noexcept {
			Derived result = derived();
			result.normalize();
			return result;
		}

		/**
		 * @brief Compute the squared distance (more efficient than distance when only
		 * comparison is needed)
		 * @param other Another element of the same Lie group
		 * @return Squared distance
		 */
		[[nodiscard]] Scalar squaredDistance(const Derived &other) const noexcept;

	protected:
		/**
		 * @brief Helper for derived classes to validate construction parameters
		 */
		template<typename... Args>
		static constexpr bool validateConstructorArgs(Args &&...args) noexcept {
			return Derived::isValidConstruction(std::forward<Args>(args)...);
		}
	};

	/**
	 * @brief Utility traits for Lie group types
	 */
	template<typename T>
	struct is_lie_group : std::false_type {};

	template<typename Derived, FloatingPoint Scalar, int Dim, int AlgDim, int ActDim>
	struct is_lie_group<LieGroupBase<Derived, Scalar, Dim, AlgDim, ActDim>> : std::true_type {};

	template<typename T>
	inline constexpr bool is_lie_group_v = is_lie_group<T>::value;

	/**
	 * @brief Concept to check if a type is a Lie group
	 */
	template<typename T>
	concept LieGroup = is_lie_group_v<T>;

} // namespace sofa::component::cosserat::liegroups

// #endif // SOFA_COMPONENT_COSSERAT_LIEGROUPS_LIEGROUPBASE_H

