// This file contains the implementation details for the SO3 (Special Orthogonal
// group in 3D) class.

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

namespace sofa::component::cosserat::liegroups {

/**
 * @brief Computes the geodesic distance between two rotations.
 * Uses the fact that the geodesic distance between two rotations is
 * twice the magnitude of the rotation angle of their difference.
 * @tparam _Scalar The scalar type.
 * @param other The other SO3 rotation.
 * @return The geodesic distance.
 */
template <typename _Scalar>
typename SO3<_Scalar>::Scalar
SO3<_Scalar>::distance(const SO3 &other) const noexcept {
  // Directly compute the log of the relative rotation
  const Eigen::Quaternion<Scalar> diff = m_quat.inverse() * other.m_quat;
  return Scalar(2.0) * std::atan2(diff.vec().norm(), std::abs(diff.w()));
}

/**
 * @brief Performs spherical linear interpolation between two rotations.
 * Uses quaternion SLERP which gives the geodesic path between
 * two rotations. The parameter t is clamped to [0,1].
 * @tparam _Scalar The scalar type.
 * @param other The target SO3 rotation.
 * @param t The interpolation parameter (between 0 and 1).
 * @return The interpolated SO3 rotation.
 */
template <typename _Scalar>
SO3<_Scalar> SO3<_Scalar>::computeInterpolate(const SO3 &other,
                                              const Scalar &t) const noexcept {
  // Clamp t to [0,1] for safety
  const Scalar t_clamped = std::max(Scalar(0), std::min(Scalar(1), t));
  // Use quaternion SLERP for optimal interpolation
  return SO3(m_quat.slerp(t_clamped, other.m_quat));
}

/**
 * @brief Implements the Baker-Campbell-Hausdorff formula for SO(3).
 * Computes log(exp(v)exp(w)) up to the specified order.
 * @tparam _Scalar The scalar type.
 * @param v The first tangent vector.
 * @param w The second tangent vector.
 * @param order The order of approximation (e.g., 1, 2, 3).
 * @return The resulting tangent vector.
 */
template <typename _Scalar>
typename SO3<_Scalar>::TangentVector
SO3<_Scalar>::BCH(const TangentVector &v, const TangentVector &w, int order) {
  // First-order approximation: v + w
  TangentVector result = v + w;

  if (order >= 2) {
    // Compute [v,w] once and store it
    const Matrix Vhat = hat(v);
    const Matrix What = hat(w);
    const Matrix VW = Vhat * What - What * Vhat;
    const TangentVector vw = vee(VW);

    // Second-order term: 1/2[v,w]
    result += vw * Scalar(0.5);

    if (order >= 3) {
      // Third-order term using stored [v,w]
      result += (Derived::computeVee(Vhat * Derived::computeHat(vw)) -
                 Derived::computeVee(What * Derived::computeHat(vw))) *
                Scalar(1.0 / 12.0);
    }
  }

  return result;
}

/**
 * @brief Computes the differential of the exponential map.
 * For small angles, uses a Taylor expansion. For larger angles, uses the
 * closed-form expression.
 * @tparam _Scalar The scalar type.
 * @param v The tangent vector.
 * @return The matrix representing the differential of exp at v.
 */
template <typename _Scalar>
typename SO3<_Scalar>::AdjointMatrix
SO3<_Scalar>::dexp(const TangentVector &v) {
  const Scalar theta = v.norm();

  if (theta < Types<Scalar>::SMALL_ANGLE_THRESHOLD) {
    return Matrix::Identity() + hat(v) * Scalar(0.5);
  }

  const Matrix V = hat(v);
  const Scalar theta2 = theta * theta;

  return Matrix::Identity() + (Scalar(1) - std::cos(theta)) / theta2 * V +
         (theta - std::sin(theta)) / (theta2 * theta) * (V * V);
}

/**
 * @brief Computes the differential of the logarithm map.
 * For small angles, uses a Taylor expansion. For larger angles, uses the
 * closed-form expression.
 * @tparam _Scalar The scalar type.
 * @return The matrix representing the differential of log at the current point.
 */
template <typename _Scalar>
typename SO3<_Scalar>::AdjointMatrix SO3<_Scalar>::dlog() const {
  const TangentVector omega = computeLog();
  const Scalar theta = omega.norm();

  if (theta < Types<Scalar>::SMALL_ANGLE_THRESHOLD) {
    return Matrix::Identity() - hat(omega) * Scalar(0.5);
  }

  const Matrix V = hat(omega);
  const Scalar theta2 = theta * theta;

  return Matrix::Identity() - Scalar(0.5) * V +
         (Scalar(1) / theta2 - (Scalar(1) + std::cos(theta)) /
                                   (Scalar(2) * theta * std::sin(theta))) *
             (V * V);
}

/**
 * @brief Computes the adjoint representation of the Lie algebra element.
 * For SO(3), this is equivalent to the hat map: ad(v)w = [v,w] = hat(v)w.
 * @tparam _Scalar The scalar type.
 * @param v The tangent vector.
 * @return The adjoint matrix.
 */
template <typename _Scalar>
typename SO3<_Scalar>::AdjointMatrix
SO3<_Scalar>::computeAd(const TangentVector &v) {
  // For SO(3), ad(v) is just the hat map
  return hat(v);
}

} // namespace sofa::component::cosserat::liegroups
