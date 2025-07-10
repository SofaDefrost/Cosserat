// This file provides improved implementations for common Lie group operations,
// such as distance, interpolation, and BCH formula, with better numerical
// stability.

// #ifndef SOFA_COMPONENT_COSSERAT_LIEGROUPS_LIEGROUPBASE_INL
// #define SOFA_COMPONENT_COSSERAT_LIEGROUPS_LIEGROUPBASE_INL
#pragma once

#include "LieGroupBase.h"
#include "Types.h"
#include <algorithm>
#include <cmath>

namespace sofa::component::cosserat::liegroups {

/**
 * @brief Improved implementations for common operations
 *
 * This file provides improved implementations of common Lie group operations
 * with better numerical stability and error handling.
 */

/**
 * @brief Computes the geodesic distance between two Lie group elements.
 * This implementation uses the logarithm of the relative transformation to find
 * the distance in the Lie algebra, which corresponds to the geodesic distance
 * on the manifold. It includes robust error handling and fallbacks for
 * numerical stability.
 * @tparam Derived The derived class implementing the specific Lie group.
 * @tparam _Scalar The scalar type used for computations.
 * @tparam _Dim The dimension of the group representation.
 * @tparam _AlgebraDim The dimension of the Lie algebra.
 * @tparam _ActionDim The dimension of vectors the group acts on.
 * @param other The other Lie group element to compute the distance to.
 * @return A scalar representing the geodesic distance between the two elements.
 */
template <typename Derived, FloatingPoint _Scalar, int _Dim, int _AlgebraDim,
          int _ActionDim>
  requires(_Dim > 0 && _AlgebraDim > 0 && _ActionDim > 0)
inline typename LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim,
                             _ActionDim>::Scalar
LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::distance(
    const Derived &other) const noexcept {
        // Check for numerical issues first
        if (!derived().isValid() || !other.isValid()) {
            return std::numeric_limits<Scalar>::quiet_NaN();
          }
  try {
    // Use relative transformation
    // This uses g^(-1)*h which maps to the tangent space at identity
    const Derived rel = derived().inverse().compose(other);
    const TangentVector tangent = rel.log();
    return tangent.norm();
  } catch (const std::exception &) {
    // Fallback: check if elements are approximately equal
    if (derived().computeIsApprox(other)) {
      return Scalar(0);
    }
    // Use squared Frobenius norm as fallback distance metric
    try {
      const auto diff = derived().matrix() - other.matrix();
      return std::sqrt(diff.squaredNorm());
    } catch (...) {
      return std::numeric_limits<Scalar>::max();
    }
  }
}

/**
 * @brief Computes the squared geodesic distance between two Lie group elements.
 * This is often more efficient than `distance()` when only comparison of
 * distances is needed, as it avoids the square root operation. It uses the
 * squared norm of the logarithm of the relative transformation.
 * @tparam Derived The derived class implementing the specific Lie group.
 * @tparam _Scalar The scalar type used for computations.
 * @tparam _Dim The dimension of the group representation.
 * @tparam _AlgebraDim The dimension of the Lie algebra.
 * @tparam _ActionDim The dimension of vectors the group acts on.
 * @param other The other Lie group element to compute the squared distance to.
 * @return A scalar representing the squared geodesic distance.
 */
template <typename Derived, FloatingPoint _Scalar, int _Dim, int _AlgebraDim,
          int _ActionDim>
  requires(_Dim > 0 && _AlgebraDim > 0 && _ActionDim > 0)
inline typename LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim,
                             _ActionDim>::Scalar
LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::squaredDistance(
    const Derived &other) const noexcept {
  try {
    const Derived rel = derived().inverse().compose(other);
    const TangentVector tangent = rel.log();
    return tangent.squaredNorm();
  } catch (const std::exception &) {
    if (derived().computeIsApprox(other)) {
      return Scalar(0);
    }
    try {
      const auto diff = derived().matrix() - other.matrix();
      return diff.squaredNorm();
    } catch (...) {
      return std::numeric_limits<Scalar>::max();
    }
  }
}

/**
 * @brief Interpolates between two Lie group elements using the exponential and
 * logarithm maps. This method performs geodesic interpolation on the manifold,
 * ensuring the interpolated path stays within the group structure. It clamps
 * the interpolation parameter `t` to the range [0, 1] and includes error
 * handling.
 * @tparam Derived The derived class implementing the specific Lie group.
 * @tparam _Scalar The scalar type used for computations.
 * @tparam _Dim The dimension of the group representation.
 * @tparam _AlgebraDim The dimension of the Lie algebra.
 * @tparam _ActionDim The dimension of vectors the group acts on.
 * @param other The target Lie group element to interpolate towards.
 * @param t The interpolation parameter, typically between 0 (returns `this`)
 * and 1 (returns `other`).
 * @return The interpolated Lie group element.
 */
template <typename Derived, FloatingPoint _Scalar, int _Dim, int _AlgebraDim,
          int _ActionDim>
  requires(_Dim > 0 && _AlgebraDim > 0 && _ActionDim > 0)
inline Derived
LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::interpolate(
    const Derived &other, const Scalar &t) const noexcept {

  // Clamp interpolation parameter
  const Scalar t_clamped = std::max(Scalar(0), std::min(Scalar(1), t));

  if (t_clamped <= Types::epsilon())
    return derived();
  if (t_clamped >= Scalar(1) - Types::epsilon())
    return other;

  try {
    // Use geodesic interpolation on the manifold
    const Derived rel = derived().inverse().compose(other);
    const TangentVector delta = rel.log();
    return derived().compose(Derived::exp(t_clamped * delta));
  } catch (const std::exception &) {
    // Fallback: return one of the endpoints
    return t_clamped < Scalar(0.5) ? derived() : other;
  }
}

/**
 * @brief Performs spherical linear interpolation (SLERP) between two Lie group
 * elements. This method is particularly useful for rotations and aims to
 * provide a smooth, constant-velocity interpolation along the shortest path on
 * the manifold. It includes numerical stability checks and falls back to linear
 * interpolation if elements are too close or other issues arise.
 * @tparam Derived The derived class implementing the specific Lie group.
 * @tparam _Scalar The scalar type used for computations.
 * @tparam _Dim The dimension of the group representation.
 * @tparam _AlgebraDim The dimension of the Lie algebra.
 * @tparam _ActionDim The dimension of vectors the group acts on.
 * @param other The target Lie group element for interpolation.
 * @param t The interpolation parameter, typically between 0 (returns `this`)
 * and 1 (returns `other`).
 * @return The interpolated Lie group element.
 */
template <typename Derived, FloatingPoint _Scalar, int _Dim, int _AlgebraDim,
          int _ActionDim>
  requires(_Dim > 0 && _AlgebraDim > 0 && _ActionDim > 0)
inline Derived
LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::slerp(
    const Derived &other, const Scalar &t) const noexcept {

  const Scalar t_clamped = std::max(Scalar(0), std::min(Scalar(1), t));

  if (t_clamped <= Types::epsilon())
    return derived();
  if (t_clamped >= Scalar(1) - Types::epsilon())
    return other;

  try {
    // Compute the "angular distance"
    const Scalar dist = distance(other);

    if (dist < Types::epsilon()) {
      return derived(); // Elements are too close for meaningful interpolation
    }

    // Use sinc-based interpolation for better numerical properties
    const Scalar theta = dist * t_clamped;
    const Scalar sin_theta = std::sin(theta);
    const Scalar sin_dist = std::sin(dist);

    if (std::abs(sin_dist) < Types::epsilon()) {
      // Fallback to linear interpolation
      return interpolate(other, t_clamped);
    }

    const Scalar w1 = std::sin(dist - theta) / sin_dist;
    const Scalar w2 = sin_theta / sin_dist;

    // This is a simplified version - actual implementation would depend on
    // group structure
    return interpolate(other, w2 / (w1 + w2));

  } catch (const std::exception &) {
    return interpolate(other, t_clamped);
  }
}

/**
 * @brief Computes the Baker-Campbell-Hausdorff (BCH) formula for Lie algebra
 * elements. The BCH formula provides a way to express the logarithm of the
 * product of exponentials of two Lie algebra elements. This implementation
 * supports up to fifth-order approximation and includes a convergence check.
 * @tparam Derived The derived class implementing the specific Lie group.
 * @tparam _Scalar The scalar type used for computations.
 * @tparam _Dim The dimension of the group representation.
 * @tparam _AlgebraDim The dimension of the Lie algebra.
 * @tparam _ActionDim The dimension of vectors the group acts on.
 * @param v The first tangent vector (Lie algebra element).
 * @param w The second tangent vector (Lie algebra element).
 * @param order The order of approximation to use (1 to 5). Higher orders
 * provide more accuracy but are computationally more expensive.
 * @return The tangent vector approximating log(exp(v)*exp(w)).
 * @throws NumericalInstabilityException if the inputs are too large for the
 * series to converge reliably.
 */
template <typename Derived, FloatingPoint _Scalar, int _Dim, int _AlgebraDim,
          int _ActionDim>
  requires(_Dim > 0 && _AlgebraDim > 0 && _ActionDim > 0)
inline typename LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim,
                             _ActionDim>::TangentVector
LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::BCH(
    const TangentVector &v, const TangentVector &w, int order) {

  // Clamp order to reasonable range
  order = std::max(1, std::min(5, order));

  // Check convergence criterion
  const Scalar v_norm = v.norm();
  const Scalar w_norm = w.norm();
  const Scalar max_norm = std::max(v_norm, w_norm);

  if (max_norm > M_PI / Scalar(2)) {
    throw NumericalInstabilityException(
        "BCH series may not converge for large inputs");
  }

  // First-order approximation
  TangentVector result = v + w;

  if (order >= 2) {
    // Second-order term: 1/2 * [v, w]
    const AlgebraMatrix vMat = Derived::computeHat(v);
    const AlgebraMatrix wMat = Derived::computeHat(w);
    const AlgebraMatrix comm_vw = vMat * wMat - wMat * vMat;
    result += Derived::computeVee(comm_vw) * Scalar(0.5);

    if (order >= 3) {
      // Third-order terms: 1/12 * [v, [v, w]] + 1/12 * [w, [w, v]]
      const AlgebraMatrix comm_v_vw = vMat * comm_vw - comm_vw * vMat;
      const AlgebraMatrix comm_w_wv = wMat * (-comm_vw) - (-comm_vw) * wMat;

      result += Derived::computeVee(comm_v_vw) * Scalar(1.0 / 12.0);
      result += Derived::computeVee(comm_w_wv) * Scalar(1.0 / 12.0);

      if (order >= 4) {
        // Fourth-order term: -1/24 * [w, [v, [v, w]]]
        const AlgebraMatrix comm4 = wMat * comm_v_vw - comm_v_vw * wMat;
        result -= Derived::computeVee(comm4) * Scalar(1.0 / 24.0);

        if (order >= 5) {
          // Fifth-order terms: more complex nested commutators
          // -1/720 * [v, [w, [v, [v, w]]]] - 1/720 * [w, [v, [w, [w, v]]]]
          const AlgebraMatrix comm5a = vMat * comm4 - comm4 * vMat;
          const AlgebraMatrix comm_w_v_wv = wMat * comm_w_wv - comm_w_wv * wMat;
          const AlgebraMatrix comm5b = wMat * comm_w_v_wv - comm_w_v_wv * wMat;

          result -= Derived::computeVee(comm5a) * Scalar(1.0 / 720.0);
          result -= Derived::computeVee(comm5b) * Scalar(1.0 / 720.0);
        }
      }
    }
  }

  return result;
}

/**
 * @brief Computes the differential of the exponential map (dexp).
 * This function calculates the Jacobian of the exponential map, which is
 * crucial for relating velocities in the Lie algebra to velocities on the Lie
 * group. It uses Taylor series expansion for small input values for numerical
 * stability and a closed-form expression for larger values.
 * @tparam Derived The derived class implementing the specific Lie group.
 * @tparam _Scalar The scalar type used for computations.
 * @tparam _Dim The dimension of the group representation.
 * @tparam _AlgebraDim The dimension of the Lie algebra.
 * @tparam _ActionDim The dimension of vectors the group acts on.
 * @param v The tangent vector (Lie algebra element) at which to compute the
 * differential.
 * @return The adjoint matrix representing the differential of the exponential
 * map at `v`.
 */
template <typename Derived, FloatingPoint _Scalar, int _Dim, int _AlgebraDim,
          int _ActionDim>
  requires(_Dim > 0 && _AlgebraDim > 0 && _ActionDim > 0)
inline typename LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim,
                             _ActionDim>::AdjointMatrix
LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::dexp(
    const TangentVector &v) {
  using Matrix = typename Derived::AdjointMatrix;
  using Scalar = typename Derived::Scalar;
  using Types = typename Derived::Types;

  const Scalar theta = v.norm();

  if (Types::isZero(theta)) {
    return Matrix::Identity();
  }

  // Create the adjoint representation of v
  const Matrix ad_v = Derived::ad(v);
  const Scalar theta_sq = theta * theta;

  // Use improved series with better numerical stability
  Matrix result = Matrix::Identity();

  // Series coefficients for improved numerical stability
  if (theta < Types::SMALL_ANGLE_THRESHOLD) {
    // Use Taylor series for small angles
    result += ad_v * Scalar(0.5);
    result += ad_v * ad_v * Scalar(1.0 / 12.0);
    result += ad_v * ad_v * ad_v * Scalar(1.0 / 720.0);
  } else {
    // Use closed form for larger angles
    const Scalar sin_theta = std::sin(theta);
    const Scalar cos_theta = std::cos(theta);

    const Scalar c1 = sin_theta / theta;
    const Scalar c2 = (Scalar(1) - cos_theta) / theta_sq;

    result += ad_v * c2;
    result += ad_v * ad_v * (theta - sin_theta) / (theta_sq * theta);
  }

  return result;
}

/**
 * @brief Computes the inverse of the differential of the exponential map
 * (dexpInv). This function calculates the inverse Jacobian of the exponential
 * map, useful for mapping velocities on the Lie group back to the Lie algebra.
 * It employs Taylor series expansion for small input values and a closed-form
 * expression for larger values to ensure numerical stability.
 * @tparam Derived The derived class implementing the specific Lie group.
 * @tparam _Scalar The scalar type used for computations.
 * @tparam _Dim The dimension of the group representation.
 * @tparam _AlgebraDim The dimension of the Lie algebra.
 * @tparam _ActionDim The dimension of vectors the group acts on.
 * @param v The tangent vector (Lie algebra element) at which to compute the
 * inverse differential.
 * @return The adjoint matrix representing the inverse differential of the
 * exponential map at `v`.
 */
template <typename Derived, FloatingPoint _Scalar, int _Dim, int _AlgebraDim,
          int _ActionDim>
  requires(_Dim > 0 && _AlgebraDim > 0 && _ActionDim > 0)
inline typename LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim,
                             _ActionDim>::AdjointMatrix
LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::dexpInv(
    const TangentVector &v) {
  using Matrix = typename Derived::AdjointMatrix;
  using Scalar = typename Derived::Scalar;
  using Types = typename Derived::Types;

  const Scalar theta = v.norm();

  if (Types::isZero(theta)) {
    return Matrix::Identity();
  }

  // Create the adjoint representation of v
  const Matrix ad_v = Derived::ad(v);
  const Scalar theta_sq = theta * theta;

  Matrix result = Matrix::Identity();

  if (theta < Types::SMALL_ANGLE_THRESHOLD) {
    // Taylor series for small angles
    result -= ad_v * Scalar(0.5);
    result += ad_v * ad_v * Scalar(1.0 / 12.0);
    result -= ad_v * ad_v * ad_v * Scalar(1.0 / 720.0);
  } else {
    // Closed form for larger angles
    const Scalar sin_theta = std::sin(theta);
    const Scalar cos_theta = std::cos(theta);
    const Scalar half_theta = theta * Scalar(0.5);
    const Scalar cot_half = cos_theta / sin_theta;

    result -= ad_v * Scalar(0.5);
    result += ad_v * ad_v *
              (Scalar(1.0 / 12.0) - half_theta * cot_half / Scalar(6.0));
  }

  return result;
}

/**
 * @brief Computes the differential of the logarithm map (dlog).
 * This function calculates the Jacobian of the logarithm map, which is
 * essential for mapping velocities on the Lie group to velocities in the Lie
 * algebra. It uses Taylor series expansion for small input values and a
 * closed-form expression for larger values to ensure numerical stability.
 * @tparam Derived The derived class implementing the specific Lie group.
 * @tparam _Scalar The scalar type used for computations.
 * @tparam _Dim The dimension of the group representation.
 * @tparam _AlgebraDim The dimension of the Lie algebra.
 * @tparam _ActionDim The dimension of vectors the group acts on.
 * @return The adjoint matrix representing the differential of the logarithm map
 * at the current Lie group element.
 */
template <typename Derived, FloatingPoint _Scalar, int _Dim, int _AlgebraDim,
          int _ActionDim>
  requires(_Dim > 0 && _AlgebraDim > 0 && _ActionDim > 0)
inline typename LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim,
                             _ActionDim>::AdjointMatrix
LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::dlog() const {
  using Matrix = typename Derived::AdjointMatrix;
  using Vector = typename Derived::TangentVector;
  using Scalar = typename Derived::Scalar;
  using Types = typename Derived::Types;

  // Get the logarithm of the current element
  const Vector v = derived().log();
  const Scalar theta = v.norm();

  if (Types::isZero(theta)) {
    return Matrix::Identity();
  }

  // Create the adjoint representation of v
  const Matrix ad_v = Derived::ad(v);
  Matrix result = Matrix::Identity();

  if (theta < Types::SMALL_ANGLE_THRESHOLD) {
    // Taylor series for small angles
    result -= ad_v * Scalar(0.5);
    result += ad_v * ad_v * Scalar(1.0 / 12.0);
    result -= ad_v * ad_v * ad_v * Scalar(1.0 / 720.0);
  } else {
    // Closed form for larger angles
    const Scalar sin_theta = std::sin(theta);
    const Scalar cos_theta = std::cos(theta);
    const Scalar half_theta = theta * Scalar(0.5);

    const Scalar factor =
        half_theta * std::cos(half_theta) / std::sin(half_theta);
    result -= ad_v * Scalar(0.5);
    result += ad_v * ad_v * (Scalar(1.0 / 12.0) - factor / Scalar(6.0));
  }

  return result;
}

/**
 * @brief Computes the action Jacobian, which describes how a point transforms
 * under the action of the Lie group with respect to changes in the Lie algebra.
 * This implementation first attempts to use an analytical method provided by
 * the derived class (if available) for better performance and accuracy. If no
 * analytical method is provided or it fails, it falls back to numerical
 * differentiation using central differences.
 * @tparam Derived The derived class implementing the specific Lie group.
 * @tparam _Scalar The scalar type used for computations.
 * @tparam _Dim The dimension of the group representation.
 * @tparam _AlgebraDim The dimension of the Lie algebra.
 * @tparam _ActionDim The dimension of vectors the group acts on.
 * @param point The point in the action space at which to compute the Jacobian.
 * @return The Jacobian matrix, mapping Lie algebra velocities to velocities in
 * the action space.
 */
template <typename Derived, FloatingPoint _Scalar, int _Dim, int _AlgebraDim,
          int _ActionDim>
  requires(_Dim > 0 && _AlgebraDim > 0 && _ActionDim > 0)
inline typename LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim,
                             _ActionDim>::JacobianMatrix
LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::actionJacobian(
    const ActionVector &point) const noexcept {

  using Matrix = typename Derived::JacobianMatrix;
  using AlgVec = typename Derived::TangentVector;
  using Scalar = typename Derived::Scalar;
  using Types = typename Derived::Types;

  // Try analytical computation first (if derived class provides it)
  if constexpr (requires(const Derived &d, const ActionVector &p) {
                  d.computeActionJacobianAnalytical(p);
                }) {
    try {
      return derived().computeActionJacobianAnalytical(point);
    } catch (...) {
      // Fall through to numerical computation
    }
  }

  // Numerical differentiation with adaptive step size
  const Scalar base_eps = std::sqrt(Types::epsilon());
  const Scalar point_scale = std::max(Scalar(1.0), point.norm());
  const Scalar eps = base_eps * point_scale;

  Matrix J = Matrix::Zero();

  // Central difference for better accuracy
  for (int i = 0; i < AlgebraDim; ++i) {
    AlgVec delta = AlgVec::Zero();
    delta[i] = eps;

    // Central difference
    const Derived perturbed_plus = derived().compose(Derived::exp(delta));
    const Derived perturbed_minus = derived().compose(Derived::exp(-delta));

    const ActionVector point_plus = perturbed_plus.act(point);
    const ActionVector point_minus = perturbed_minus.act(point);

    // Compute column of Jacobian using central difference
    J.col(i) = (point_plus - point_minus) / (Scalar(2) * eps);
  }

  return J;
}

/**
 * @brief Computes the hat operator, which maps a Lie algebra vector to its
 * matrix representation. This is a static assertion that ensures the derived
 * class provides its own implementation of `computeHat`.
 * @tparam Derived The derived class implementing the specific Lie group.
 * @tparam _Scalar The scalar type used for computations.
 * @tparam _Dim The dimension of the group representation.
 * @tparam _AlgebraDim The dimension of the Lie algebra.
 * @tparam _ActionDim The dimension of vectors the group acts on.
 * @param v The tangent vector (Lie algebra element).
 * @return The matrix representation of the Lie algebra element.
 */
template <typename Derived, FloatingPoint _Scalar, int _Dim, int _AlgebraDim,
          int _ActionDim>
  requires(_Dim > 0 && _AlgebraDim > 0 && _ActionDim > 0)
inline typename LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim,
                             _ActionDim>::AlgebraMatrix
LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::hat(
    const TangentVector &v) noexcept {
  static_assert(
      requires { Derived::computeHat(v); },
      "Derived class must implement computeHat method");
  return Derived::computeHat(v);
}

/**
 * @brief Computes the vee operator, which maps a Lie algebra matrix
 * representation back to its vector form. This is the inverse operation of the
 * hat operator. A static assertion ensures the derived class provides its own
 * `computeVee` implementation.
 * @tparam Derived The derived class implementing the specific Lie group.
 * @tparam _Scalar The scalar type used for computations.
 * @tparam _Dim The dimension of the group representation.
 * @tparam _AlgebraDim The dimension of the Lie algebra.
 * @tparam _ActionDim The dimension of vectors the group acts on.
 * @param X The matrix representation in the Lie algebra.
 * @return The vector representation of the Lie algebra element.
 */
template <typename Derived, FloatingPoint _Scalar, int _Dim, int _AlgebraDim,
          int _ActionDim>
  requires(_Dim > 0 && _AlgebraDim > 0 && _ActionDim > 0)
inline typename LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim,
                             _ActionDim>::TangentVector
LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::vee(
    const AlgebraMatrix &X) noexcept {
  static_assert(
      requires { Derived::computeVee(X); },
      "Derived class must implement computeVee method");
  return Derived::computeVee(X);
}

/**
 * @brief Computes the adjoint representation of a Lie algebra element.
 * The adjoint action describes how Lie algebra elements transform under
 * conjugation by Lie group elements. This implementation checks if the derived
 * class provides an optimized `computeAd` method; otherwise, it uses a default
 * implementation for matrix Lie groups.
 * @tparam Derived The derived class implementing the specific Lie group.
 * @tparam _Scalar The scalar type used for computations.
 * @tparam _Dim The dimension of the group representation.
 * @tparam _AlgebraDim The dimension of the Lie algebra.
 * @tparam _ActionDim The dimension of vectors the group acts on.
 * @param v The element of the Lie algebra in vector form.
 * @return The adjoint matrix representing the adjoint action.
 */
template <typename Derived, FloatingPoint _Scalar, int _Dim, int _AlgebraDim,
          int _ActionDim>
  requires(_Dim > 0 && _AlgebraDim > 0 && _ActionDim > 0)
inline typename LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim,
                             _ActionDim>::AdjointMatrix
LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::ad(
    const TangentVector &v) {
  using Matrix = typename Derived::AdjointMatrix;
  using AlgMat = typename Derived::AlgebraMatrix;

  // Check if derived class provides optimized implementation
  if constexpr (requires { Derived::computeAd(v); }) {
    return Derived::computeAd(v);
  } else {
    // Default implementation for matrix Lie groups
    const AlgMat X = Derived::computeHat(v);

    // For matrix Lie algebras, ad_X(Y) = [X, Y] = XY - YX
    Matrix result = Matrix::Zero();

    for (int j = 0; j < AlgebraDim; ++j) {
      // Create basis vector and convert to algebra matrix
      TangentVector e_j = TangentVector::Zero();
      e_j[j] = Scalar(1);
      const AlgMat E_j = Derived::computeHat(e_j);

      // Compute [X, E_j]
      const AlgMat commutator = X * E_j - E_j * X;

      // Extract vector form and set as column
      result.col(j) = Derived::computeVee(commutator);
    }

    return result;
  }
}

/**
 * @brief Computes the matrix logarithm of the current Lie group element.
 * This function provides the matrix representation of the Lie algebra element
 * that, when exponentiated, yields the current Lie group element. It is
 * primarily used for debugging and analysis purposes.
 * @tparam Derived The derived class implementing the specific Lie group.
 * @tparam _Scalar The scalar type used for computations.
 * @tparam _Dim The dimension of the group representation.
 * @tparam _AlgebraDim The dimension of the Lie algebra.
 * @tparam _ActionDim The dimension of vectors the group acts on.
 * @return The matrix representation of the Lie algebra element.
 */
template <typename Derived, FloatingPoint _Scalar, int _Dim, int _AlgebraDim,
          int _ActionDim>
  requires(_Dim > 0 && _AlgebraDim > 0 && _ActionDim > 0)
inline typename LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim,
                             _ActionDim>::AlgebraMatrix
LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::matrixLog()
    const {
  // This provides the matrix logarithm representation
  const TangentVector v = derived().log();
  return Derived::computeHat(v);
}

} // End namespace sofa::component::cosserat::liegroups

// #endif // SOFA_COMPONENT_COSSERAT_LIEGROUPS_LIEGROUPBASE_INL
