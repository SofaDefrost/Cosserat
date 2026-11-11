// This file contains the implementation details for the SGal3 (Special Galilean
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
 * along with this program. If not, see <http://www.gnu.org/licenses/\>. *
 ******************************************************************************/

// #ifndef SOFA_COMPONENT_COSSERAT_LIEGROUPS_SGAL3_INL
// #define SOFA_COMPONENT_COSSERAT_LIEGROUPS_SGAL3_INL
#pragma once

#include "SGal3.h"

namespace sofa::component::cosserat::liegroups {

	/**
	 * @brief Additional operators and utility functions for SGal(3)
	 */

	/**
	 * @brief Transform a point-velocity-time tuple by inverse transformation
	 */
	template<typename Scalar_>
	Eigen::Matrix<Scalar_, 7, 1> operator/(const Eigen::Matrix<Scalar_, 7, 1> &point_vel_time,
										   const SGal3<Scalar_> &g) {
		return g.inverse().act(point_vel_time).template head<7>();
	}

	/**
	 * @brief Create SGal(3) transformation from components
	 */
	template<typename Scalar_>
	SGal3<Scalar_> fromComponents(const typename SGal3<Scalar_>::Vector3 &position, const SO3<Scalar_> &rotation,
								  const typename SGal3<Scalar_>::Vector3 &velocity, const Scalar_ &time) {
		return SGal3<Scalar_>(SE3<Scalar_>(rotation, position), velocity, time);
	}

	/**
	 * @brief Create SGal(3) transformation from position, Euler angles, velocity,
	 * and time
	 */
	template<typename Scalar_>
	SGal3<Scalar_> fromPositionEulerVelocityTime(const typename SGal3<Scalar_>::Vector3 &position, const Scalar_ &roll,
												 const Scalar_ &pitch, const Scalar_ &yaw,
												 const typename SGal3<Scalar_>::Vector3 &velocity,
												 const Scalar_ &time) {
		return SGal3<Scalar_>(fromPositionEulerZYX(position, roll, pitch, yaw), velocity, time);
	}

	/**
	 * @brief Convert transformation to position, Euler angles, velocity, and time
	 */
	template<typename Scalar_>
	typename SGal3<Scalar_>::Vector toPositionEulerVelocityTime(const SGal3<Scalar_> &transform) {
		typename SGal3<Scalar_>::Vector result;
		result.resize(10);
		result.template segment<6>(0) = toPositionEulerZYX(transform.pose());
		result.template segment<3>(6) = transform.velocity();
		result[9] = transform.time();
		return result;
	}


	// BCH: the generic LieGroup::BCH(v, w, order) is used; SGal3::BCH(X, Y) was declared and only
	// defined here, in a file nobody included, and with brackets that ignored the time terms.

} // namespace sofa::component::cosserat::liegroups

// #endif // SOFA_COMPONENT_COSSERAT_LIEGROUPS_SGAL3_INL
