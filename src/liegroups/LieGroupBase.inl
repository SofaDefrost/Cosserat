// This file provides improved implementations for common Lie group operations,
// such as distance, interpolation, and BCH formula, with better numerical
// stability.

// #ifndef SOFA_COMPONENT_COSSERAT_LIEGROUPS_LIEGROUPBASE_INL
// #define SOFA_COMPONENT_COSSERAT_LIEGROUPS_LIEGROUPBASE_INL
#pragma once

#include <algorithm>
#include <cmath>
#include "LieGroupBase.h"
#include "Types.h"

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
	template<typename Derived, FloatingPoint _Scalar, int _Dim, int _AlgebraDim, int _ActionDim>
		requires(_Dim > 0 && _AlgebraDim > 0 && _ActionDim > 0)
	inline typename LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::Scalar
	LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::distance(const Derived &other) const noexcept {
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
	template<typename Derived, FloatingPoint _Scalar, int _Dim, int _AlgebraDim, int _ActionDim>
		requires(_Dim > 0 && _AlgebraDim > 0 && _ActionDim > 0)
	inline typename LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::Scalar
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
	template<typename Derived, FloatingPoint _Scalar, int _Dim, int _AlgebraDim, int _ActionDim>
		requires(_Dim > 0 && _AlgebraDim > 0 && _ActionDim > 0)
	inline Derived
	LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::interpolate(const Derived &other,
																			   const Scalar &t) const noexcept {

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
	template<typename Derived, FloatingPoint _Scalar, int _Dim, int _AlgebraDim, int _ActionDim>
		requires(_Dim > 0 && _AlgebraDim > 0 && _ActionDim > 0)
	inline Derived
	LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::slerp(const Derived &other,
																		 const Scalar &t) const noexcept {

		const Scalar t_clamped = std::max(Scalar(0), std::min(Scalar(1), t));

		if (t_clamped <= Types::epsilon())
			return derived();
		if (t_clamped >= Scalar(1) - Types::epsilon())
			return other;

		// The geodesic g exp(t log(g^-1 h)) is the SLERP on any group with a bi-invariant metric,
		// and its parameter is t itself. The previous code reweighted t by
		// sin(t d) / (sin((1-t) d) + sin(t d)), which is the chord/sphere weighting for
		// interpolating *vectors* and made slerp(0.25) differ from interpolate(0.25) by 0.125 rad.
		return interpolate(other, t_clamped);
	}

	/**
	 * @brief Computes the Baker-Campbell-Hausdorff (BCH) formula for Lie algebra
	 * elements. The BCH formula provides a way to express the logarithm of the
	 * product of exponentials of two Lie algebra elements. This implementation
	 * supports up to fifth-order approximation. The result is a truncated series: it is
	 * always defined, accurate for small inputs (error ~ |.|^(order+1)) and degrading for large ones,
	 * exactly like SO3::BCH (no exception).
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
	 */
	template<typename Derived, FloatingPoint _Scalar, int _Dim, int _AlgebraDim, int _ActionDim>
		requires(_Dim > 0 && _AlgebraDim > 0 && _ActionDim > 0)
	inline typename LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::TangentVector
	LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::BCH(const TangentVector &v, const TangentVector &w,
																	   int order) {

		// Clamp order to reasonable range
		order = std::max(1, std::min(5, order));

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
						// Fifth-order terms (Dynkin / Wikipedia):
						//   -1/720 ([y,[y,[y,[y,x]]]] + [x,[x,[x,[x,y]]]])
						//   +1/360 ([x,[y,[y,[y,x]]]] + [y,[x,[x,[x,y]]]])
						//   +1/120 ([y,[x,[y,[x,y]]]] + [x,[y,[x,[y,x]]]])
						auto br = [](const AlgebraMatrix &A, const AlgebraMatrix &B) -> AlgebraMatrix {
							return A * B - B * A;
						};
						const AlgebraMatrix xy = comm_vw;
						const AlgebraMatrix yx = -comm_vw;
						const AlgebraMatrix xxy = br(vMat, xy);
						const AlgebraMatrix yyx = br(wMat, yx);
						const AlgebraMatrix xxxy = br(vMat, xxy);
						const AlgebraMatrix yyyx = br(wMat, yyx);
						const AlgebraMatrix yxy = br(wMat, xy);
						const AlgebraMatrix xyx = br(vMat, yx);
						const AlgebraMatrix t720 = br(wMat, yyyx) + br(vMat, xxxy);
						const AlgebraMatrix t360 = br(vMat, yyyx) + br(wMat, xxxy);
						const AlgebraMatrix t120 = br(wMat, br(vMat, yxy)) + br(vMat, br(wMat, xyx));
						result -= Derived::computeVee(t720) * Scalar(1.0 / 720.0);
						result += Derived::computeVee(t360) * Scalar(1.0 / 360.0);
						result += Derived::computeVee(t120) * Scalar(1.0 / 120.0);
					}
				}
			}
		}

		return result;
	}

	/**
	 * @brief Alternative implementation of the differential of the exponential map.
	 * This function calculates the Jacobian of the exponential map, which is
	 * crucial for relating velocities in the Lie algebra to velocities on the Lie
	 * group. It evaluates the series sum_k (-ad_v)^k / (k+1)! from the ad() operator, which
	 * has no singularity and no cancellation.
	 *
	 * NOTE: This is an alternative implementation to the wrapper dexp() in LieGroupBase.h
	 * that calls rightJacobian(); both return the right Jacobian (checked against each other
	 * in Tests/unit/test_LieGroupBase_Generic.cpp).
	 * 
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
	template<typename Derived, FloatingPoint _Scalar, int _Dim, int _AlgebraDim, int _ActionDim>
		requires(_Dim > 0 && _AlgebraDim > 0 && _ActionDim > 0)
	inline typename LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::AdjointMatrix
	LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::dexp_alternative(const TangentVector &v) {
		using Matrix = typename Derived::AdjointMatrix;
		using Scalar = typename Derived::Scalar;

		// Right Jacobian Jr(v) = sum_k (-ad_v)^k / (k+1)!, the same convention as dexp() ==
		// rightJacobian(). The series is valid for every group and every |v| (it converges like
		// exp), so there is no small-angle branch and no cancellation to guard against.
		const Matrix minus_ad = -Derived::ad(v);
		Matrix result = Matrix::Identity();
		Matrix term = Matrix::Identity();
		for (int k = 1; k < 64; ++k) {
			term = term * minus_ad / Scalar(k + 1);
			result += term;
			if (term.cwiseAbs().maxCoeff() <= std::numeric_limits<Scalar>::epsilon() * Scalar(1e-3)) {
				break;
			}
		}
		return result;
	}

	/**
	 * @brief Alternative implementation of the inverse differential of the exponential map.
	 * This function calculates the inverse Jacobian of the exponential
	 * map, useful for mapping velocities on the Lie group back to the Lie algebra.
	 * Computed as the inverse of dexp_alternative().
	 *
	 * NOTE: This is an alternative implementation to the wrapper dexpInv() in LieGroupBase.h
	 * that calls rightJacobianInverse(); both return the right Jacobian inverse.
	 * 
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
	template<typename Derived, FloatingPoint _Scalar, int _Dim, int _AlgebraDim, int _ActionDim>
		requires(_Dim > 0 && _AlgebraDim > 0 && _ActionDim > 0)
	inline typename LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::AdjointMatrix
	LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::dexpInv_alternative(const TangentVector &v) {
		// Inverse of the series of dexp_alternative(). The Bernoulli series
		// I + ad/2 + ad^2/12 - ad^4/720 + ... would only converge for |ad_v| < 2 pi and its
		// coefficients are easy to get wrong; inverting the always-convergent series is exact
		// up to the conditioning of Jr (singular only when exp is, i.e. at |phi| = 2 pi k).
		return dexp_alternative(v).inverse();
	}

	/**
	 * @brief Computes the differential of the logarithm map (dlog).
	 * This function calculates the Jacobian of the logarithm map, which is
	 * essential for mapping velocities on the Lie group to velocities in the Lie
	 * algebra: log(g exp(d)) = log(g) + dlog() d + O(d^2), i.e. dexpInv(log g).
	 * Groups may hide it with their own convention (SO3::dlog() is the left one).
	 * @tparam Derived The derived class implementing the specific Lie group.
	 * @tparam _Scalar The scalar type used for computations.
	 * @tparam _Dim The dimension of the group representation.
	 * @tparam _AlgebraDim The dimension of the Lie algebra.
	 * @tparam _ActionDim The dimension of vectors the group acts on.
	 * @return The adjoint matrix representing the differential of the logarithm map
	 * at the current Lie group element.
	 */
	template<typename Derived, FloatingPoint _Scalar, int _Dim, int _AlgebraDim, int _ActionDim>
		requires(_Dim > 0 && _AlgebraDim > 0 && _ActionDim > 0)
	inline typename LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::AdjointMatrix
	LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::dlog() const {
		// log(g exp(d)) = log(g) + Jr^-1(log g) d + O(d^2): the inverse of the right Jacobian
		// at log(g), i.e. dexpInv(log g), consistent with dexp()/dexpInv() of this class.
		const TangentVector v = derived().log();
		if constexpr (requires { Derived::computeRightJacobianInverse(v); }) {
			return Derived::computeRightJacobianInverse(v);
		} else {
			return dexpInv_alternative(v);
		}
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
	template<typename Derived, FloatingPoint _Scalar, int _Dim, int _AlgebraDim, int _ActionDim>
		requires(_Dim > 0 && _AlgebraDim > 0 && _ActionDim > 0)
	inline typename LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::JacobianMatrix
	LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::actionJacobian(
			const ActionVector &point) const noexcept {

		using Matrix = typename Derived::JacobianMatrix;
		using AlgVec = typename Derived::TangentVector;
		using Scalar = typename Derived::Scalar;
		using Types = typename Derived::Types;

		// Try analytical computation first (if derived class provides it)
		if constexpr (requires(const Derived &d, const ActionVector &p) { d.computeActionJacobianAnalytical(p); }) {
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
	template<typename Derived, FloatingPoint _Scalar, int _Dim, int _AlgebraDim, int _ActionDim>
		requires(_Dim > 0 && _AlgebraDim > 0 && _ActionDim > 0)
	inline typename LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::AlgebraMatrix
	LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::hat(const TangentVector &v) noexcept {
		static_assert(requires { Derived::computeHat(v); }, "Derived class must implement computeHat method");
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
	template<typename Derived, FloatingPoint _Scalar, int _Dim, int _AlgebraDim, int _ActionDim>
		requires(_Dim > 0 && _AlgebraDim > 0 && _ActionDim > 0)
	inline typename LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::TangentVector
	LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::vee(const AlgebraMatrix &X) noexcept {
		static_assert(requires { Derived::computeVee(X); }, "Derived class must implement computeVee method");
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
	template<typename Derived, FloatingPoint _Scalar, int _Dim, int _AlgebraDim, int _ActionDim>
		requires(_Dim > 0 && _AlgebraDim > 0 && _ActionDim > 0)
	inline typename LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::AdjointMatrix
	LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::ad(const TangentVector &v) {
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
	template<typename Derived, FloatingPoint _Scalar, int _Dim, int _AlgebraDim, int _ActionDim>
		requires(_Dim > 0 && _AlgebraDim > 0 && _ActionDim > 0)
	inline typename LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::AlgebraMatrix
	LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::matrixLog() const {
		// This provides the matrix logarithm representation
		const TangentVector v = derived().log();
		return Derived::computeHat(v);
	}

} // End namespace sofa::component::cosserat::liegroups

// #endif // SOFA_COMPONENT_COSSERAT_LIEGROUPS_LIEGROUPBASE_INL
