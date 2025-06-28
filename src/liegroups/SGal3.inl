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

// #ifndef SOFA_COMPONENT_COSSERAT_LIEGROUPS_SGAL3_INL
// #define SOFA_COMPONENT_COSSERAT_LIEGROUPS_SGAL3_INL
#pragma once

#include "SGal3.h"

namespace sofa::component::cosserat::liegroups {

template <typename _Scalar> class SGal3;

/**
 * @brief Additional operators and utility functions for SGal(3)
 */

/**
 * @brief Transform a point-velocity-time tuple by inverse transformation
 */
template <typename _Scalar>
Eigen::Matrix<_Scalar, 7, 1>
operator/(const Eigen::Matrix<_Scalar, 7, 1> &point_vel_time,
          const SGal3<_Scalar> &g) {
  return g.inverse().act(point_vel_time).template head<7>();
}

/**
 * @brief Create SGal(3) transformation from components
 */
template <typename _Scalar>
SGal3<_Scalar> fromComponents(const typename SGal3<_Scalar>::Vector3 &position,
                              const SO3<_Scalar> &rotation,
                              const typename SGal3<_Scalar>::Vector3 &velocity,
                              const _Scalar &time) {
  return SGal3<_Scalar>(SE3<_Scalar>(rotation, position), velocity, time);
}

/**
 * @brief Create SGal(3) transformation from position, Euler angles, velocity,
 * and time
 */
template <typename _Scalar>
SGal3<_Scalar> fromPositionEulerVelocityTime(
    const typename SGal3<_Scalar>::Vector3 &position, const _Scalar &roll,
    const _Scalar &pitch, const _Scalar &yaw,
    const typename SGal3<_Scalar>::Vector3 &velocity, const _Scalar &time) {
  return SGal3<_Scalar>(fromPositionEulerZYX(position, roll, pitch, yaw),
                        velocity, time);
}

/**
 * @brief Convert transformation to position, Euler angles, velocity, and time
 */
template <typename _Scalar>
typename SGal3<_Scalar>::Vector
toPositionEulerVelocityTime(const SGal3<_Scalar> &transform) {
  typename SGal3<_Scalar>::Vector result;
  result.resize(10);
  result.template segment<6>(0) = toPositionEulerZYX(transform.pose());
  result.template segment<3>(6) = transform.velocity();
  result[9] = transform.time();
  return result;
}

/**
 * @brief Interpolate between two Galilean transformations
 *
 * This implementation uses the exponential map to perform proper interpolation
 * in the Lie algebra space, including velocity and time interpolation.
 *
 * @param from Starting Galilean transformation
 * @param to Ending Galilean transformation
 * @param t Interpolation parameter in [0,1]
 * @return Interpolated Galilean transformation
 */
template <typename _Scalar>
SGal3<_Scalar> interpolate(const SGal3<_Scalar> &from, const SGal3<_Scalar> &to,
                           const _Scalar &t) {
  // Convert 'to' relative to 'from'
  SGal3<_Scalar> rel = from.inverse() * to;
  // Get the relative motion in the Lie algebra
  typename SGal3<_Scalar>::TangentVector delta = rel.log();
  // Scale it by t and apply it to 'from'
  return from * SGal3<_Scalar>().exp(t * delta);
}

/**
 * @brief Dual vector operator for sgal(3)
 * Maps a 10D vector to its dual 6x6 matrix representation
 */
template <typename _Scalar>
Eigen::Matrix<_Scalar, 6, 6>
dualMatrix(const typename SGal3<_Scalar>::TangentVector &xi) {
  Eigen::Matrix<_Scalar, 6, 6> xi_hat = Eigen::Matrix<_Scalar, 6, 6>::Zero();

  // Extract components
  const auto &v = xi.template segment<3>(0);    // Linear velocity
  const auto &w = xi.template segment<3>(3);    // Angular velocity
  const auto &beta = xi.template segment<3>(6); // Boost
  const _Scalar &tau = xi[9];                   // Time rate

  // Fill the matrix blocks
  xi_hat.template block<3, 3>(0, 0) = SO3<_Scalar>::hat(w);
  xi_hat.template block<3, 1>(0, 3) = v;
  xi_hat.template block<3, 1>(0, 4) = beta;
  xi_hat(4, 5) = tau;

  return xi_hat;
}

/**
 * @brief Specialization of the Baker-Campbell-Hausdorff formula for SGal(3)
 *
 * For SGal(3), the BCH formula has a closed form up to second order:
 * BCH(X,Y) = X + Y + 1/2[X,Y] + higher order terms
 * where [X,Y] is the Lie bracket for sgal(3).
 */
template <typename _Scalar>
typename SGal3<_Scalar>::TangentVector
SGal3<_Scalar>::BCH(const TangentVector &X, const TangentVector &Y) {
  // Extract components
  const auto &v1 = X.template segment<3>(0); // First linear velocity
  const auto &w1 = X.template segment<3>(3); // First angular velocity
  const auto &b1 = X.template segment<3>(6); // First boost
  const _Scalar &t1 = X[9];                  // First time rate

  const auto &v2 = Y.template segment<3>(0); // Second linear velocity
  const auto &w2 = Y.template segment<3>(3); // Second angular velocity
  const auto &b2 = Y.template segment<3>(6); // Second boost
  const _Scalar &t2 = Y[9];                  // Second time rate

  // Compute Lie bracket components
  const auto w1_cross_w2 = w1.cross(w2); // Angular x Angular
  const auto w1_cross_v2 = w1.cross(v2); // Angular x Linear
  const auto v1_cross_w2 = v1.cross(w2); // Linear x Angular
  const auto w1_cross_b2 = w1.cross(b2); // Angular x Boost
  const auto b1_cross_w2 = b1.cross(w2); // Boost x Angular

  // Combine terms for the BCH formula up to second order
  TangentVector result;
  result.template segment<3>(0) =
      v1 + v2 + _Scalar(0.5) * (w1_cross_v2 - v1_cross_w2);
  result.template segment<3>(3) = w1 + w2 + _Scalar(0.5) * w1_cross_w2;
  result.template segment<3>(6) =
      b1 + b2 + _Scalar(0.5) * (w1_cross_b2 - b1_cross_w2);
  result[9] = t1 + t2; // Time component adds linearly

  return result;
}

} // namespace sofa::component::cosserat::liegroups

// #endif // SOFA_COMPONENT_COSSERAT_LIEGROUPS_SGAL3_INL
