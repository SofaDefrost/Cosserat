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
	template<typename Scalar_>
	Eigen::Matrix<Scalar_, 6, 1> operator/(const Eigen::Matrix<Scalar_, 6, 1> &point_vel, const SE23<Scalar_> &g) {
		return g.inverse().act(point_vel).template head<6>();
	}

	/**
	 * @brief Create SE_2(3) transformation from pose and velocity components
	 */
	template<typename Scalar_>
	SE23<Scalar_> fromComponents(const typename SE23<Scalar_>::Vector3 &position, const SO3<Scalar_> &rotation,
								 const typename SE23<Scalar_>::Vector3 &velocity) {
		return SE23<Scalar_>(SE3<Scalar_>(rotation, position), velocity);
	}

	/**
	 * @brief Create SE_2(3) transformation from position, Euler angles, and
	 * velocity
	 */
	template<typename Scalar_>
	SE23<Scalar_> fromPositionEulerVelocity(const typename SE23<Scalar_>::Vector3 &position, const Scalar_ &roll,
											const Scalar_ &pitch, const Scalar_ &yaw,
											const typename SE23<Scalar_>::Vector3 &velocity) {
		return SE23<Scalar_>(fromPositionEulerZYX(position, roll, pitch, yaw), velocity);
	}

	/**
	 * @brief Convert transformation to position, Euler angles, and velocity
	 */
	template<typename Scalar_>
	typename SE23<Scalar_>::Vector toPositionEulerVelocity(const SE23<Scalar_> &transform) {
		typename SE23<Scalar_>::Vector result;
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
	template<typename Scalar_>
	SE23<Scalar_> interpolate(const SE23<Scalar_> &from, const SE23<Scalar_> &to, const Scalar_ &t) {
		// Convert 'to' relative to 'from'
		SE23<Scalar_> rel = from.inverse() * to;
		// Get the relative motion in the Lie algebra
		typename SE23<Scalar_>::TangentVector delta = rel.log();
		// Scale it by t and apply it to 'from'
		return from * SE23<Scalar_>().exp(t * delta);
	}


	// BCH: the generic LieGroup::BCH(v, w, order) now works for SE23 (hat/vee are implemented).
	// The two-argument SE23::BCH that used to be defined here was never declared in the class.

} // namespace sofa::component::cosserat::liegroups

#endif // SOFA_COMPONENT_COSSERAT_LIEGROUPS_SE23_INL
