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
	template<typename _Scalar>
	typename SE3<_Scalar>::Vector3 operator/(const typename SE3<_Scalar>::Vector3 &v, const SE3<_Scalar> &g) {
		return g.inverse().act(v);
	}

	/**
	 * @brief Rotation from Euler angles, ZYX convention: R = Rz(yaw) Ry(pitch) Rx(roll).
	 * (Used by fromPositionEulerZYX below; it was called but never defined anywhere.)
	 */
	template<typename _Scalar>
	SO3<_Scalar> fromEulerZYX(const _Scalar &roll, const _Scalar &pitch, const _Scalar &yaw) {
		using Q = Eigen::Quaternion<_Scalar>;
		const Q q = Q(Eigen::AngleAxis<_Scalar>(yaw, Eigen::Matrix<_Scalar, 3, 1>::UnitZ())) *
					Q(Eigen::AngleAxis<_Scalar>(pitch, Eigen::Matrix<_Scalar, 3, 1>::UnitY())) *
					Q(Eigen::AngleAxis<_Scalar>(roll, Eigen::Matrix<_Scalar, 3, 1>::UnitX()));
		return SO3<_Scalar>(q);
	}

	/**
	 * @brief Euler angles (roll, pitch, yaw), ZYX convention, of a rotation. At gimbal lock
	 * (|pitch| = pi/2) the roll is set to 0 and the yaw absorbs the remaining rotation.
	 */
	template<typename _Scalar>
	Eigen::Matrix<_Scalar, 3, 1> toEulerZYX(const SO3<_Scalar> &rotation) {
		const Eigen::Matrix<_Scalar, 3, 3> R = rotation.matrix();
		const _Scalar sp = -R(2, 0);
		_Scalar roll, pitch, yaw;
		if (std::abs(sp) < _Scalar(1) - Types<_Scalar>::epsilon() * _Scalar(16)) {
			pitch = std::asin(sp);
			roll = std::atan2(R(2, 1), R(2, 2));
			yaw = std::atan2(R(1, 0), R(0, 0));
		} else {
			pitch = sp > 0 ? _Scalar(M_PI) / _Scalar(2) : -_Scalar(M_PI) / _Scalar(2);
			roll = _Scalar(0);
			yaw = std::atan2(-R(0, 1), R(1, 1));
		}
		return Eigen::Matrix<_Scalar, 3, 1>(roll, pitch, yaw);
	}

	/**
	 * @brief Create SE3 transformation from position and Euler angles (ZYX
	 * convention)
	 * @param position Translation vector
	 * @param roll Rotation around X axis (in radians)
	 * @param pitch Rotation around Y axis (in radians)
	 * @param yaw Rotation around Z axis (in radians)
	 */
	template<typename _Scalar>
	SE3<_Scalar> fromPositionEulerZYX(const typename SE3<_Scalar>::Vector3 &position, const _Scalar &roll,
									  const _Scalar &pitch, const _Scalar &yaw) {
		return SE3<_Scalar>(fromEulerZYX(roll, pitch, yaw), position);
	}

	/**
	 * @brief Convert transformation to position and Euler angles (ZYX convention)
	 * @return Vector containing (x, y, z, roll, pitch, yaw)
	 */
	template<typename _Scalar>
	typename SE3<_Scalar>::TangentVector toPositionEulerZYX(const SE3<_Scalar> &transform) {
		// SE3::Vector is the 4-vector of the matrix representation; the result has 6 entries.
		typename SE3<_Scalar>::TangentVector result;
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
	template<typename _Scalar>
	SE3<_Scalar> interpolate(const SE3<_Scalar> &from, const SE3<_Scalar> &to, const _Scalar &t) {
		// Convert 'to' relative to 'from'
		SE3<_Scalar> rel = from.inverse() * to;
		// Get the relative motion in the Lie algebra
		typename SE3<_Scalar>::TangentVector delta = rel.log();
		// Scale it by t and apply it to 'from'
		return from * SE3<_Scalar>::exp(t * delta);
	}


	// BCH implementation moved to header file

} // namespace sofa::component::cosserat::liegroups

// #endif // SOFA_COMPONENT_COSSERAT_LIEGROUPS_SE3_INL
