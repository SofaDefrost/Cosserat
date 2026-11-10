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

	template<typename Scalar_>
	class SE23;

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

	/**
	 * @brief Dual vector operator for se_2(3)
	 * Maps a 9D vector to its dual 5x5 matrix representation
	 */
	template<typename Scalar_>
	Eigen::Matrix<Scalar_, 5, 5> dualMatrix(const typename SE23<Scalar_>::TangentVector &xi) {
		Eigen::Matrix<Scalar_, 5, 5> xi_hat = Eigen::Matrix<Scalar_, 5, 5>::Zero();

		// Extract components
		const auto &v = xi.template segment<3>(0); // Linear velocity
		const auto &w = xi.template segment<3>(3); // Angular velocity
		const auto &a = xi.template segment<3>(6); // Linear acceleration

		// Fill the matrix blocks
		xi_hat.template block<3, 3>(0, 0) = SO3<Scalar_>::hat(w);
		xi_hat.template block<3, 1>(0, 3) = v;
		xi_hat.template block<3, 1>(0, 4) = a;

		return xi_hat;
	}

	/**
	 * @brief Specialization of the Baker-Campbell-Hausdorff formula for SE_2(3)
	 *
	 * For SE_2(3), the BCH formula has a closed form up to second order:
	 * BCH(X,Y) = X + Y + 1/2[X,Y] + higher order terms
	 * where [X,Y] is the Lie bracket for se_2(3).
	 */
	template<typename Scalar_>
	typename SE23<Scalar_>::TangentVector SE23<Scalar_>::BCH(const TangentVector &X, const TangentVector &Y) {
		// Extract components
		const auto &v1 = X.template segment<3>(0); // First linear velocity
		const auto &w1 = X.template segment<3>(3); // First angular velocity
		const auto &a1 = X.template segment<3>(6); // First linear acceleration

		const auto &v2 = Y.template segment<3>(0); // Second linear velocity
		const auto &w2 = Y.template segment<3>(3); // Second angular velocity
		const auto &a2 = Y.template segment<3>(6); // Second linear acceleration

		// Compute Lie bracket components
		const auto w1_cross_w2 = w1.cross(w2); // Angular x Angular
		const auto w1_cross_v2 = w1.cross(v2); // Angular x Linear
		const auto v1_cross_w2 = v1.cross(w2); // Linear x Angular
		const auto w1_cross_a2 = w1.cross(a2); // Angular x Acceleration
		const auto a1_cross_w2 = a1.cross(w2); // Acceleration x Angular

		// Combine terms for the BCH formula up to second order
		TangentVector result;
		result.template segment<3>(0) = v1 + v2 + Scalar_(0.5) * (w1_cross_v2 - v1_cross_w2);
		result.template segment<3>(3) = w1 + w2 + Scalar_(0.5) * w1_cross_w2;
		result.template segment<3>(6) = a1 + a2 + Scalar_(0.5) * (w1_cross_a2 - a1_cross_w2);

		return result;
	}

} // namespace sofa::component::cosserat::liegroups

#endif // SOFA_COMPONENT_COSSERAT_LIEGROUPS_SE23_INL
