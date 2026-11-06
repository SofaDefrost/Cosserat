// This file contains the implementation details for the SE23 (extended Special
// Euclidean group in 3D) class.

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

#ifndef SOFA_COMPONENT_COSSERAT_LIEGROUPS_SE23_INL
#define SOFA_COMPONENT_COSSERAT_LIEGROUPS_SE23_INL

#include "SE23.h"

namespace sofa::component::cosserat::liegroups {

	/**
	 * @brief Additional operators and utility functions for SE_2(3)
	 */

	/**
	 * @brief Transform a point-velocity pair by inverse transformation
	 */
	template<typename _Scalar>
	Eigen::Matrix<_Scalar, 6, 1> operator/(const Eigen::Matrix<_Scalar, 6, 1> &point_vel, const SE23<_Scalar> &g) {
		return g.inverse().act(point_vel).template head<6>();
	}

	/**
	 * @brief Create SE_2(3) transformation from pose and velocity components
	 */
	template<typename _Scalar>
	SE23<_Scalar> fromComponents(const typename SE23<_Scalar>::Vector3 &position, const SO3<_Scalar> &rotation,
								 const typename SE23<_Scalar>::Vector3 &velocity) {
		return SE23<_Scalar>(SE3<_Scalar>(rotation, position), velocity);
	}

	/**
	 * @brief Create SE_2(3) transformation from position, Euler angles, and
	 * velocity
	 */
	template<typename _Scalar>
	SE23<_Scalar> fromPositionEulerVelocity(const typename SE23<_Scalar>::Vector3 &position, const _Scalar &roll,
											const _Scalar &pitch, const _Scalar &yaw,
											const typename SE23<_Scalar>::Vector3 &velocity) {
		return SE23<_Scalar>(fromPositionEulerZYX(position, roll, pitch, yaw), velocity);
	}

	/**
	 * @brief Convert transformation to position, Euler angles, and velocity
	 */
	template<typename _Scalar>
	typename SE23<_Scalar>::Vector toPositionEulerVelocity(const SE23<_Scalar> &transform) {
		typename SE23<_Scalar>::Vector result;
		result.resize(9);
		result.template segment<6>(0) = toPositionEulerZYX(transform.pose());
		result.template segment<3>(6) = transform.velocity();
		return result;
	}

	/**
	 * @brief Interpolate between two extended poses
	 *
	 * This implementation uses the exponential map to perform proper interpolation
	 * in the Lie algebra space, including velocity interpolation.
	 *
	 * @param from Starting extended pose
	 * @param to Ending extended pose
	 * @param t Interpolation parameter in [0,1]
	 * @return Interpolated extended pose
	 */
	template<typename _Scalar>
	SE23<_Scalar> interpolate(const SE23<_Scalar> &from, const SE23<_Scalar> &to, const _Scalar &t) {
		// Convert 'to' relative to 'from'
		SE23<_Scalar> rel = from.inverse() * to;
		// Get the relative motion in the Lie algebra
		typename SE23<_Scalar>::TangentVector delta = rel.log();
		// Scale it by t and apply it to 'from'
		return from * SE23<_Scalar>().exp(t * delta);
	}

	/**
	 * @brief Dual vector operator for se_2(3): maps xi = [v ; w ; a] to the 5x5 matrix
	 * [[w^, v, a], [0, 0, 0], [0, 0, 0]] (the top-left block of SE23::computeHat).
	 */
	template<typename _Scalar>
	Eigen::Matrix<_Scalar, 5, 5> dualMatrix(const typename SE23<_Scalar>::TangentVector &xi) {
		return SE23<_Scalar>::computeHat(xi).template block<5, 5>(0, 0);
	}

	// BCH: the generic LieGroupBase::BCH(v, w, order) now works for SE23 (hat/vee are implemented).
	// The two-argument SE23::BCH that used to be defined here was never declared in the class.

} // namespace sofa::component::cosserat::liegroups

#endif // SOFA_COMPONENT_COSSERAT_LIEGROUPS_SE23_INL
