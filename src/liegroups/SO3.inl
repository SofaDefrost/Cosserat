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
	template<typename _Scalar>
	typename SO3<_Scalar>::Scalar SO3<_Scalar>::distance(const SO3 &other) const noexcept {
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
	template<typename _Scalar>
	SO3<_Scalar> SO3<_Scalar>::interpolate(const SO3 &other, const Scalar &t) const noexcept {
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
	template<typename _Scalar>
	typename SO3<_Scalar>::TangentVector SO3<_Scalar>::BCH(const TangentVector &v, const TangentVector &w, int order) {
		// First-order approximation: v + w
		TangentVector result = v + w;

		if (order >= 2) {
			// [v,w] = v x w for so(3) (vee of the matrix commutator)
			const TangentVector vw = v.cross(w);

			// Second-order term: 1/2[v,w]
			result += vw * Scalar(0.5);

			if (order >= 3) {
				// Third-order term: 1/12 ([v,[v,w]] + [w,[w,v]]) = 1/12 (v x vw - w x vw)
				result += (v.cross(vw) - w.cross(vw)) * Scalar(1.0 / 12.0);
			}
		}

		return result;
	}

	namespace detail {
	/**
	 * @brief Coefficients shared by the SO(3) and SE(3) Jacobians and exponentials.
	 *   a(t) = (1 - cos t) / t^2         = 1/2  - t^2/24  + t^4/720 - ...
	 *   b(t) = (t - sin t) / t^3         = 1/6  - t^2/120 + t^4/5040 - ...
	 *   c(t) = 1/t^2 - 1/(2 t tan(t/2))  = 1/12 + t^2/720 + t^4/30240 + ...   (inverse Jacobian)
	 *
	 * Below kSeriesAngle the Taylor series (6 terms, truncation < 1e-16) is used: the closed
	 * forms lose digits there, (1-cos t)/t^2 and (t-sin t)/t^3 cancel, and the previous switch at
	 * 1e-4 still left 2e-13 errors in the Jacobian just above it. Above, a uses 2 sin^2(t/2)/t^2
	 * (no cancellation) and c the half-angle form, which stays finite and accurate up to pi
	 * (the (1+cos t)/(2 t sin t) form loses everything near pi).
	 */
	template<typename Scalar>
	inline constexpr Scalar so3SeriesAngle() noexcept {
		return Types<Scalar>::seriesAngle();
	}

	template<typename Scalar>
	inline void so3JacobianCoefficients(const Scalar theta, Scalar &a, Scalar &b, Scalar &c) noexcept {
		// a and b are the scalar toolbox functions of Types.h (Taylor series / cancellation-free forms)
		a = Types<Scalar>::sinc2(theta);
		b = Types<Scalar>::sinc3(theta);
		const Scalar t2 = theta * theta;
		if (theta < so3SeriesAngle<Scalar>()) {
			const Scalar t4 = t2 * t2;
			c = Scalar(1) / Scalar(12) + t2 / Scalar(720) + t4 / Scalar(30240) + t4 * t2 / Scalar(1209600) +
				t4 * t4 / Scalar(47900160);
		} else {
			c = Scalar(1) / t2 - Scalar(1) / (Scalar(2) * theta * std::tan(Scalar(0.5) * theta));
		}
	}
	} // namespace detail

	/**
	 * @brief Computes the differential of the exponential map (the left Jacobian Jl).
	 *   Jl(w) = I + a(|w|) [w]x + b(|w|) [w]x^2
	 * @tparam _Scalar The scalar type.
	 * @param v The tangent vector.
	 * @return The matrix representing the differential of exp at v.
	 */
	template<typename _Scalar>
	typename SO3<_Scalar>::AdjointMatrix SO3<_Scalar>::dexp(const TangentVector &v) {
		Scalar a, b, c;
		detail::so3JacobianCoefficients(v.norm(), a, b, c);
		const Matrix V = SO3::computeHat(v);
		return Matrix::Identity() + a * V + b * (V * V);
	}

	/**
	 * @brief Computes the differential of the logarithm map (inverse of the left Jacobian,
	 * evaluated at log of the current element).
	 *   Jl^-1(w) = I - 1/2 [w]x + c(|w|) [w]x^2
	 * @tparam _Scalar The scalar type.
	 * @return The matrix representing the differential of log at the current point.
	 */
	template<typename _Scalar>
	typename SO3<_Scalar>::AdjointMatrix SO3<_Scalar>::dlog() const {
		const TangentVector omega = computeLog();
		Scalar a, b, c;
		detail::so3JacobianCoefficients(omega.norm(), a, b, c);
		const Matrix V = SO3::computeHat(omega);
		return Matrix::Identity() - Scalar(0.5) * V + c * (V * V);
	}

	/**
	 * @brief Computes the adjoint representation of the Lie algebra element.
	 * For SO(3), this is equivalent to the hat map: ad(v)w = [v,w] = hat(v)w.
	 * @tparam _Scalar The scalar type.
	 * @param v The tangent vector.
	 * @return The adjoint matrix.
	 */
	template<typename _Scalar>
	typename SO3<_Scalar>::AdjointMatrix SO3<_Scalar>::computeAd(const TangentVector &v) {
		// For SO(3), ad(v) is the hat map. Do not call ad(v) here: LieGroupBase::ad
		// dispatches back to computeAd, which would recurse forever.
		return computeHat(v);
	}

	/**
	 * @brief Right Jacobian of SO(3): Jr(w) = I - a(|w|) [w]x + b(|w|) [w]x^2
	 * (Solà et al. 2021, eq. 143).
	 */
	template<typename _Scalar>
	typename SO3<_Scalar>::AdjointMatrix SO3<_Scalar>::computeRightJacobian(const TangentVector &omega) noexcept {
		Scalar a, b, c;
		detail::so3JacobianCoefficients(omega.norm(), a, b, c);
		const Matrix W = SO3::computeHat(omega);
		return Matrix::Identity() - a * W + b * (W * W);
	}

	/**
	 * @brief Left Jacobian of SO(3): Jl(w) = Jr(-w).
	 */
	template<typename _Scalar>
	typename SO3<_Scalar>::AdjointMatrix SO3<_Scalar>::computeLeftJacobian(const TangentVector &omega) noexcept {
		return computeRightJacobian(-omega);
	}

	/**
	 * @brief Inverse right Jacobian: Jr^-1(w) = I + 1/2 [w]x + c(|w|) [w]x^2
	 * (Solà et al. 2021, eq. 144).
	 */
	template<typename _Scalar>
	typename SO3<_Scalar>::AdjointMatrix SO3<_Scalar>::computeRightJacobianInverse(const TangentVector &omega) noexcept {
		Scalar a, b, c;
		detail::so3JacobianCoefficients(omega.norm(), a, b, c);
		const Matrix W = SO3::computeHat(omega);
		return Matrix::Identity() + Scalar(0.5) * W + c * (W * W);
	}

	/**
	 * @brief Inverse left Jacobian: Jl^-1(w) = Jr^-1(-w).
	 */
	template<typename _Scalar>
	typename SO3<_Scalar>::AdjointMatrix SO3<_Scalar>::computeLeftJacobianInverse(const TangentVector &omega) noexcept {
		return computeRightJacobianInverse(-omega);
	}

} // namespace sofa::component::cosserat::liegroups
