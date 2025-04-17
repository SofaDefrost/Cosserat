/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                  *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
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

#ifndef SOFA_COMPONENT_COSSERAT_LIEGROUPS_LIEGROUPBASE_H
#define SOFA_COMPONENT_COSSERAT_LIEGROUPS_LIEGROUPBASE_H

#include "Types.h"
#include <type_traits>
#include <stdexcept>
#include <memory>

namespace sofa::component::cosserat::liegroups
{

// Forward declaration for matrix representation of Lie algebra
template<typename _Scalar, int _Dim>
class LieAlgebra;

/**
 * @brief Base class template for all Lie group implementations using CRTP
 *
 * This class defines the interface that all Lie group implementations must satisfy.
 * It uses the Curiously Recurring Template Pattern (CRTP) to allow static polymorphism,
 * providing better performance than virtual functions.
 *
 * @tparam Derived The derived class implementing the specific Lie group
 * @tparam _Scalar The scalar type used for computations (must be a floating-point type)
 * @tparam _Dim The dimension of the group representation
 * @tparam _AlgebraDim The dimension of the Lie algebra (default: same as _Dim)
 * @tparam _ActionDim The dimension of vectors the group acts on (default: same as _Dim)
 */
template<typename Derived, typename _Scalar,
         int _Dim, int _AlgebraDim = _Dim, int _ActionDim = _Dim>
class LieGroupBase {
public:
    static_assert(std::is_floating_point<_Scalar>::value,
                 "Scalar type must be a floating-point type");
    static_assert(_Dim > 0, "Dimension must be positive");
    static_assert(_AlgebraDim > 0, "Algebra dimension must be positive");
    static_assert(_ActionDim > 0, "Action dimension must be positive");

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

  // Tag for exception safety
  struct InverseFailureException : public std::runtime_error {
    InverseFailureException() : std::runtime_error("Failed to compute inverse: singular element") {}
  };

  /**
     * @brief Access to the derived object (const version)
     * @return Reference to the derived object
     */
  inline const Derived& derived() const noexcept {
    return static_cast<const Derived&>(*this);
  }

  /**
     * @brief Access to the derived object (non-const version)
     * @return Reference to the derived object
     */
  inline Derived& derived() noexcept {
    return static_cast<Derived&>(*this);
  }

  /**
     * @brief In-place group composition
     * @param other Another element of the same Lie group
     * @return Reference to this after composition
     */
  inline Derived& operator*=(const Derived& other) noexcept {
    // In-place group composition (g = g * h)
    // Uses the derived class's compose implementation
    derived() = derived().compose(other);
    return derived();
  }


  /**
     * @brief Group composition operation
     * @param other Another element of the same Lie group
     * @return The composition this * other
     */
  inline Derived operator*(const Derived& other) const noexcept {
    // Binary group composition (g * h)
    // Returns a new group element without modifying the operands
    return derived().compose(other);
  }

  /**
     * @brief Compute the inverse element
     * @return The inverse element such that this * inverse() = identity()
     * @throws InverseFailureException if the element is singular (not invertible)
     */
  inline Derived inverse() const
  {
    // Computes the inverse element g^(-1)
    // Such that g * g^(-1) = g^(-1) * g = identity
    return derived().computeInverse();
  }
  //
  /**
     * @brief Exponential map from Lie algebra to Lie group
     * @param algebra_element Element of the Lie algebra (tangent space at identity)
     * @return The corresponding element in the Lie group
     */
  static inline Derived exp(const TangentVector& algebra_element) noexcept {
    return Derived::computeExp(algebra_element);
  }

  /**
     * @brief Logarithmic map from Lie group to Lie algebra
     * @return The corresponding element in the Lie algebra
     * @throws std::domain_error if the element is outside the domain of the log map
     */
  inline TangentVector log() const {
    return derived().computeLog();
  }

  /**
     * @brief Adjoint representation of the group element
     * @return The adjoint matrix representing the action on the Lie algebra
     */
  inline AdjointMatrix adjoint() const noexcept {
    return derived().computeAdjoint();
  }

  /**
     * @brief Group action on a point
     * @param point The point to transform
     * @return The transformed point
     */
  inline ActionVector act(const ActionVector& point) const noexcept {
    return derived().computeAction(point);
  }

  /**
     * @brief Check if this element is approximately equal to another
     * @param other Another element of the same Lie group
     * @param eps Tolerance for comparison (optional)
     * @return true if elements are approximately equal
     */
  inline bool isApprox(const Derived& other,
                       const Scalar& eps = Types::epsilon()) const noexcept {
    return derived().computeIsApprox(other, eps);
  }

  /**
     * @brief Get the identity element of the group
     * @return The identity element
     */
  static inline Derived Identity() noexcept {
    return Derived::computeIdentity();
  }

  /**
     * @brief Compute distance between two group elements
     * @param other Another element of the same Lie group
     * @return A scalar representing the distance
     */
  Scalar distance(const Derived& other) const noexcept;

  /**
     * @brief Interpolate between two group elements
     * @param other Target group element
     * @param t Interpolation parameter between 0 and 1
     * @return Interpolated group element
     */
  Derived interpolate(const Derived& other, const Scalar& t) const noexcept;

  /**
     * @brief Get the Jacobian of the group action
     * @param point The point at which to compute the Jacobian
     * @return The Jacobian matrix
     */
  JacobianMatrix actionJacobian(const ActionVector& point) const noexcept;

  /**
     * @brief Convert a tangent vector to a Lie algebra matrix representation (hat operator)
     * @param v Vector representation of Lie algebra element
     * @return Matrix representation in the Lie algebra
     */
  static AlgebraMatrix hat(const TangentVector& v) noexcept;

  /**
     * @brief Convert a Lie algebra matrix to its vector representation (vee operator)
     * @param X Matrix representation in the Lie algebra
     * @return Vector representation of the Lie algebra element
     */
  static TangentVector vee(const AlgebraMatrix& X) noexcept;

  /**
     * @brief Baker-Campbell-Hausdorff formula for Lie algebra elements
     * @param v First tangent vector
     * @param w Second tangent vector
     * @param order Order of approximation (1-3)
     * @return Tangent vector approximating log(exp(v)*exp(w))
     */
  static TangentVector BCH(const TangentVector& v,
                          const TangentVector& w, int order = 2);

  /**
     * @brief Get the differential of the exponential map
     * @param v Tangent vector
     * @return Matrix representing the differential of exp at v
     */
  static AdjointMatrix dexp(const TangentVector& v);

  /**
     * @brief Get the differential of the logarithm map
     * @return Matrix representing the differential of log at the current point
     */
  AdjointMatrix dlog() const;

  /**
     * @brief Create a group element from parameters
     * @tparam Args Parameter types
     * @param args Constructor arguments
     * @return New group element
     */
  template<typename... Args>
    static inline Derived create(Args&&... args) {
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
    static AdjointMatrix ad(const TangentVector& v);
};

} // namespace sofa::component::cosserat::liegroups

#endif // SOFA_COMPONENT_COSSERAT_LIEGROUPS_LIEGROUPBASE_H
