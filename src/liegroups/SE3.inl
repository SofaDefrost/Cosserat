// This file contains additional operators and utility functions for the SE3
// (Special Euclidean group in 3D) class.

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

// #ifndef SOFA_COMPONENT_COSSERAT_LIEGROUPS_SE3_INL
// #define SOFA_COMPONENT_COSSERAT_LIEGROUPS_SE3_INL
#pragma once

#include "SE3.h"

namespace sofa::component::cosserat::liegroups {

/**
 * @brief Additional operators and utility functions for SE3
 */

/**
 * @brief Transform a point by inverse transformation
 */
template <typename _Scalar>
typename SE3<_Scalar>::Vector3
operator/(const typename SE3<_Scalar>::Vector3 &v, const SE3<_Scalar> &g) {
  return g.inverse().act(v);
}

/**
 * @brief Create SE3 transformation from position and Euler angles (ZYX
 * convention)
 * @param position Translation vector
 * @param roll Rotation around X axis (in radians)
 * @param pitch Rotation around Y axis (in radians)
 * @param yaw Rotation around Z axis (in radians)
 */
template <typename _Scalar>
SE3<_Scalar>
fromPositionEulerZYX(const typename SE3<_Scalar>::Vector3 &position,
                     const _Scalar &roll, const _Scalar &pitch,
                     const _Scalar &yaw) {
  return SE3<_Scalar>(fromEulerZYX(roll, pitch, yaw), position);
}

/**
 * @brief Convert transformation to position and Euler angles (ZYX convention)
 * @return Vector containing (x, y, z, roll, pitch, yaw)
 */
template <typename _Scalar>
typename SE3<_Scalar>::Vector
toPositionEulerZYX(const SE3<_Scalar> &transform) {
  typename SE3<_Scalar>::Vector result;
  result.template head<3>() = transform.translation();
  result.template tail<3>() = toEulerZYX(transform.rotation());
  return result;
}

/**
 * @brief Interpolate between two transformations
 *
 * This implementation uses the exponential map to perform proper interpolation
 * in the Lie algebra space.
 *
 * @param from Starting transformation
 * @param to Ending transformation
 * @param t Interpolation parameter in [0,1]
 * @return Interpolated transformation
 */
template <typename _Scalar>
SE3<_Scalar> interpolate(const SE3<_Scalar> &from, const SE3<_Scalar> &to,
                         const _Scalar &t) {
  // Convert 'to' relative to 'from'
  SE3<_Scalar> rel = from.inverse() * to;
  // Get the relative motion in the Lie algebra
  typename SE3<_Scalar>::TangentVector delta = rel.log();
  // Scale it by t and apply it to 'from'
  return from * SE3<_Scalar>::exp(t * delta);
}

/**
 * @brief Dual vector operator for se(3)
 * Maps a 6D vector to its dual 4x4 matrix representation
 */
template <typename _Scalar>
Eigen::Matrix<_Scalar, 4, 4>
dualMatrix(const typename SE3<_Scalar>::TangentVector &xi) {
  Eigen::Matrix<_Scalar, 4, 4> xi_hat = Eigen::Matrix<_Scalar, 4, 4>::Zero();
  xi_hat.template block<3, 3>(0, 0) = SO3<_Scalar>::hat(xi.template tail<3>());
  xi_hat.template block<3, 1>(0, 3) = xi.template head<3>();
  return xi_hat;
}

// BCH implementation moved to header file

} // namespace sofa::component::cosserat::liegroups

// #endif // SOFA_COMPONENT_COSSERAT_LIEGROUPS_SE3_INL