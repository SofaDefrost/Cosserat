#pragma once

#include <cmath>
#include <concepts>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <limits>
#include <random>
#include <type_traits>

namespace sofa::component::cosserat::liegroups {

	// Helper type for compile-time integer constants (for template parameters)
	template<int N>
	using IntConst = std::integral_constant<int, N>;

	/**
	 * @brief Type definitions and utilities for Lie group implementations
	 *
	 * This class provides type aliases and utility functions for different
	 * scalar types used in Lie group computations.
	 */
	template<typename Scalar_>
		requires (std::is_floating_point_v<Scalar_> || std::is_class_v<Scalar_>)
	class Types {
	public:
		using Scalar = Scalar_;

		// Eigen type aliases
		template<int Rows, int Cols = Rows>
		using Matrix = Eigen::Matrix<Scalar, Rows, Cols>;

		template<int Size>
		using Vector = Eigen::Matrix<Scalar, Size, 1>;

		// Dynamic size aliases
		using MatrixX = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
		using VectorX = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

		// Common fixed-size types
		using Matrix2 = Matrix<2, 2>;
		using Matrix3 = Matrix<3, 3>;
		using Matrix4 = Matrix<4, 4>;
		using Matrix6 = Matrix<6, 6>;

		using Vector2 = Vector<2>;
		using Vector3 = Vector<3>;
		using Vector4 = Vector<4>;
		using Vector6 = Vector<6>;

		// Quaternion type
		using Quaternion = Eigen::Quaternion<Scalar>;

		// Rotation types
		using AngleAxis = Eigen::AngleAxis<Scalar>;
		using Rotation2D = Eigen::Rotation2D<Scalar>;

		// Transform types
		using Transform2 = Eigen::Transform<Scalar, 2, Eigen::Affine>;
		using Transform3 = Eigen::Transform<Scalar, 3, Eigen::Affine>;
		using Isometry2 = Eigen::Transform<Scalar, 2, Eigen::Isometry>;
		using Isometry3 = Eigen::Transform<Scalar, 3, Eigen::Isometry>;

		// Constants
		static constexpr Scalar SMALL_ANGLE_THRESHOLD = Scalar(1e-4);
		static constexpr Scalar PI = Scalar(M_PI);
		static constexpr Scalar TWO_PI = Scalar(2 * M_PI);

		/**
		 * @brief Get machine epsilon for the scalar type
		 */
		static constexpr Scalar epsilon() noexcept { return std::numeric_limits<Scalar>::epsilon(); }

		/**
		 * @brief Get a small tolerance value for comparisons
		 */
		static constexpr Scalar tolerance() noexcept { return Scalar(100) * epsilon(); }

		/**
		 * @brief Check if a value is effectively zero within a given tolerance.
		 * @param value The scalar value to check.
		 * @param tol The tolerance for considering a value as zero. Defaults to
		 * `tolerance()`.
		 * @return True if the absolute value of `value` is less than or equal to
		 * `tol`, false otherwise.
		 */
		static constexpr bool isZero(const Scalar &value, const Scalar &tol = tolerance()) noexcept {
			return std::abs(value) <= tol;
		}

		/**
		 * @brief Check if two scalar values are approximately equal within a given
		 * tolerance.
		 * @param a The first scalar value.
		 * @param b The second scalar value.
		 * @param tol The tolerance for considering values as approximately equal.
		 * Defaults to `tolerance()`.
		 * @return True if the absolute difference between `a` and `b` is less than or
		 * equal to `tol`, false otherwise.
		 */
		static constexpr bool isApprox(const Scalar &a, const Scalar &b, const Scalar &tol = tolerance()) noexcept {
			return std::abs(a - b) <= tol;
		}

		/**
		 * @brief Computes sin(x)/x with numerical stability for small x.
		 * Uses a Taylor expansion for x near zero to avoid division by zero.
		 * @param x The input value.
		 * @return The value of sin(x)/x.
		 */
		static Scalar sinc(const Scalar &x) noexcept {
			if (isZero(x, SMALL_ANGLE_THRESHOLD)) {
				// Taylor expansion: sin(x)/x ≈ 1 - x²/6 + x⁴/120 - ...
				const Scalar x2 = x * x;
				return Scalar(1) - x2 / Scalar(6) + x2 * x2 / Scalar(120);
			}
			return std::sin(x) / x;
		}

		/**
		 * @brief Angle below which the Taylor series of (1-cos x)/x^2 and (x-sin x)/x^3 are used.
		 *
		 * Both closed forms cancel catastrophically for small x (relative error ~ eps/x^2 and eps/x^2 again),
		 * so the 6-term series (truncation < 1e-16 at 0.25 rad) replaces them up to this angle. float cancels
		 * earlier than double and the series is still exact to float precision up to 1 rad, hence 1.
		 */
		static constexpr Scalar seriesAngle() noexcept {
			if constexpr (std::is_floating_point_v<Scalar> && sizeof(Scalar) <= 4) {
				return Scalar(1);
			} else {
				return Scalar(0.25);
			}
		}

		/**
		 * @brief Computes (1 - cos(x))/x^2 = 1/2 - x^2/24 + x^4/720 - ... without cancellation.
		 * Taylor series below seriesAngle(), 2 sin^2(x/2)/x^2 above (no 1 - cos subtraction).
		 * @param x The input value (x >= 0 or not: the function is even).
		 * @return The value of (1 - cos(x))/x^2.
		 */
		static Scalar sinc2(const Scalar &x) noexcept {
			const Scalar x2 = x * x;
			if (std::abs(x) < seriesAngle()) {
				const Scalar x4 = x2 * x2;
				return Scalar(1) / Scalar(2) - x2 / Scalar(24) + x4 / Scalar(720) - x4 * x2 / Scalar(40320) +
					   x4 * x4 / Scalar(3628800) - x4 * x4 * x2 / Scalar(479001600);
			}
			const Scalar s_half = std::sin(Scalar(0.5) * x);
			return Scalar(2) * s_half * s_half / x2;
		}

		/**
		 * @brief Computes (x - sin(x))/x^3 = 1/6 - x^2/120 + x^4/5040 - ... without cancellation.
		 * Taylor series below seriesAngle(), closed form above (the cancellation there is at most eps/x^2 * x).
		 * @param x The input value.
		 * @return The value of (x - sin(x))/x^3.
		 */
		static Scalar sinc3(const Scalar &x) noexcept {
			const Scalar x2 = x * x;
			if (std::abs(x) < seriesAngle()) {
				const Scalar x4 = x2 * x2;
				return Scalar(1) / Scalar(6) - x2 / Scalar(120) + x4 / Scalar(5040) - x4 * x2 / Scalar(362880) +
					   x4 * x4 / Scalar(39916800) - x4 * x4 * x2 / Scalar(6227020800.0);
			}
			return (x - std::sin(x)) / (x2 * x);
		}

		/**
		 * @brief Computes 1 - cos(x) as 2 sin^2(x/2): no subtraction of two nearly equal numbers
		 * (1 - std::cos(x) has a relative error of eps/x^2, i.e. 8e-8 at x = 1e-5 and 100% below 1e-8).
		 * @param x The input value.
		 * @return The value of 1 - cos(x).
		 */
		static Scalar oneMinusCos(const Scalar &x) noexcept {
			const Scalar s_half = std::sin(Scalar(0.5) * x);
			return Scalar(2) * s_half * s_half;
		}

		/**
		 * @brief Checks if an angle is effectively zero.
		 * @param angle The angle to check.
		 * @param tol The tolerance.
		 * @return True if the angle is near zero.
		 */
		static bool isAngleNearZero(const Scalar &angle, const Scalar &tol = tolerance()) noexcept {
			return isZero(normalizeAngle(angle), tol);
		}

		/**
		 * @brief Checks if two angles are effectively equal (modulo 2pi).
		 * @param a The first angle.
		 * @param b The second angle.
		 * @param tol The tolerance.
		 * @return True if the angles are nearly equal.
		 */
		static bool areAnglesNearlyEqual(const Scalar &a, const Scalar &b, const Scalar &tol = tolerance()) noexcept {
			return isZero(angleDistance(a, b), tol);
		}

		/**
		 * @brief Computes the arc tangent of y/x, using the signs of both arguments
		 * to determine the correct quadrant.
		 * @param y The y-coordinate.
		 * @param x The x-coordinate.
		 * @return The angle in radians between the positive x-axis and the point (x,
		 * y).
		 */
		static Scalar atan2(const Scalar &y, const Scalar &x) noexcept { return std::atan2(y, x); }

		/**
		 * @brief Computes the safe square root of a non-negative number.
		 * Handles negative inputs gracefully by returning 0 for negative values.
		 * @param x The input value.
		 * @return The square root of x, or 0 if x is negative.
		 */
		static Scalar safeSqrt(const Scalar &x) noexcept { return std::sqrt(std::max(Scalar(0), x)); }

		/**
		 * @brief Normalizes an angle to the closed range [-pi, pi].
		 * An angle already in the range is returned unchanged (so pi stays pi and -pi stays -pi); other angles
		 * are wrapped by a multiple of 2 pi.
		 * @param angle The angle to normalize in radians.
		 * @return The normalized angle in radians.
		 */
		static Scalar normalizeAngle(const Scalar &angle) noexcept {
			if (angle >= -PI && angle <= PI) {
				return angle;
			}
			Scalar result = std::fmod(angle, TWO_PI);
			if (result > PI) {
				result -= TWO_PI;
			} else if (result < -PI) {
				result += TWO_PI;
			}
			return result;
		}

		/**
		 * @brief Computes the angle difference a - b with proper wrapping.
		 * Always returns the shortest signed angular difference, in [-pi, pi].
		 * @param a The first angle.
		 * @param b The second angle.
		 * @return The shortest angle difference in [-π, π].
		 */
		static Scalar angleDifference(Scalar a, Scalar b) { return normalizeAngle(a - b); }

		/**
		 * @brief Performs spherical linear interpolation (SLERP) between two angles.
		 * Always takes the shortest path.
		 * @param a The start angle.
		 * @param b The end angle.
		 * @param t The interpolation parameter [0, 1].
		 * @return The interpolated angle.
		 */
		static Scalar slerpAngle(Scalar a, Scalar b, Scalar t) {
			Scalar diff = angleDifference(b, a);
			return normalizeAngle(a + t * diff);
		}

		/**
		 * @brief Computes the bi-invariant distance between two angles.
		 * @param a The first angle.
		 * @param b The second angle.
		 * @return The absolute angular distance.
		 */
		static Scalar angleDistance(Scalar a, Scalar b) { return std::abs(angleDifference(a, b)); }

		/**
		 * @brief Safely normalizes a vector, handling near-zero cases.
		 * @tparam Derived The Eigen derived type.
		 * @param v The vector to normalize.
		 * @param fallback Optional fallback unit vector if normalization fails.
		 * @return The normalized vector or fallback.
		 */
		template<typename Derived>
		static auto safeNormalize(const Eigen::MatrixBase<Derived> &v, const Eigen::MatrixBase<Derived> &fallback) {
			using VectorType = typename Derived::PlainObject;
			const Scalar norm = v.norm();

			if (norm < tolerance()) {
				if (fallback.norm() > tolerance()) {
					return fallback.normalized();
				}
				// Return first standard basis vector as ultimate fallback
				VectorType result = VectorType::Zero(v.rows());
				if (result.rows() > 0) {
					result(0) = Scalar(1);
				}
				return result;
			}

			return (v / norm).eval();
		}

		/**
		 * @brief Safely normalizes a vector, handling near-zero cases (no fallback version).
		 * @tparam Derived The Eigen derived type.
		 * @param v The vector to normalize.
		 * @return The normalized vector or first basis vector if near-zero.
		 */
		template<typename Derived>
		static auto safeNormalize(const Eigen::MatrixBase<Derived> &v) {
			using VectorType = typename Derived::PlainObject;
			const Scalar norm = v.norm();

			if (norm < tolerance()) {
				// Return first standard basis vector as fallback
				VectorType result = VectorType::Zero(v.rows());
				if (result.rows() > 0) {
					result(0) = Scalar(1);
				}
				return result;
			}

			return (v / norm).eval();
		}

		/**
		 * @brief Projects a vector onto another vector safely.
		 * @tparam Derived1 The type of the vector to project.
		 * @tparam Derived2 The type of the vector to project onto.
		 * @param v The vector to project.
		 * @param onto The vector to project onto.
		 * @return The projected vector.
		 */
		template<typename Derived1, typename Derived2>
		static auto projectVector(const Eigen::MatrixBase<Derived1> &v, const Eigen::MatrixBase<Derived2> &onto) {
			using VectorType = typename Derived1::PlainObject;
			const Scalar norm_squared = onto.squaredNorm();

			if (norm_squared < tolerance()) {
				return VectorType::Zero(v.rows()).eval();
			}

			return (onto * (v.dot(onto) / norm_squared)).eval();
		}

		/**
		 * @brief Clamps a value between a minimum and maximum value.
		 * @param value The value to clamp.
		 * @param min_val The minimum allowed value.
		 * @param max_val The maximum allowed value.
		 * @return The clamped value.
		 */
		static constexpr Scalar clamp(const Scalar &value, const Scalar &min_val, const Scalar &max_val) noexcept {
			return std::max(min_val, std::min(max_val, value));
		}

		/**
		 * @brief Performs linear interpolation between two scalar values.
		 * @param a The start value.
		 * @param b The end value.
		 * @param t The interpolation parameter (typically between 0 and 1).
		 * @return The interpolated value.
		 */
		static constexpr Scalar lerp(const Scalar &a, const Scalar &b, const Scalar &t) noexcept {
			return a + t * (b - a);
		}

		/**
		 * @brief Generates a random scalar value within the range [0, 1].
		 * @tparam Generator The type of the random number generator.
		 * @param gen The random number generator.
		 * @return A random scalar value.
		 */
		template<typename Generator>
		static Scalar randomScalar(Generator &gen) noexcept {
			std::uniform_real_distribution<Scalar> dist(Scalar(0), Scalar(1));
			return dist(gen);
		}

		/**
		 * @brief Generates a random scalar value within a specified range.
		 * @tparam Generator The type of the random number generator.
		 * @param gen The random number generator.
		 * @param min_val The minimum value.
		 * @param max_val The maximum value.
		 * @return A random scalar value.
		 */
		template<typename Generator>
		static Scalar randomScalar(Generator &gen, const Scalar &min_val, const Scalar &max_val) noexcept {
			std::uniform_real_distribution<Scalar> dist(min_val, max_val);
			return dist(gen);
		}

		/**
		 * @brief Generates a random vector with components in the range [-1, 1].
		 * @tparam N The dimension of the vector.
		 * @tparam Generator The type of the random number generator.
		 * @param gen The random number generator.
		 * @return A random vector.
		 */
		template<int N, typename Generator>
		static Vector<N> randomVector(Generator &gen) noexcept {
			std::uniform_real_distribution<Scalar> dist(Scalar(-1), Scalar(1));
			Vector<N> result;
			for (int i = 0; i < N; ++i) {
				result(i) = dist(gen);
			}
			return result;
		}

		/**
		 * @brief Generates a random unit vector, uniformly distributed on the sphere S^(N-1).
		 *
		 * The components are independent standard normals (a rotation-invariant distribution) and the
		 * vector is normalised. Normalising a vector drawn in the cube [-1, 1]^N instead favours the
		 * cube diagonals (E|x| = 0.516 instead of 0.5 on S^2).
		 * @tparam N The dimension of the vector.
		 * @tparam Generator The type of the random number generator.
		 * @param gen The random number generator.
		 * @return A random unit vector.
		 */
		template<int N, typename Generator>
		static Vector<N> randomUnitVector(Generator &gen) noexcept {
			std::normal_distribution<Scalar> dist(Scalar(0), Scalar(1));
			Vector<N> v;
			for (int i = 0; i < N; ++i) {
				v(i) = dist(gen);
			}
			const Scalar norm = v.norm();
			if (norm > epsilon()) {
				v /= norm;
			} else {
				v = Vector<N>::Unit(0); // Fallback to first basis vector
			}
			return v;
		}

	private:
		// Make constructor private to prevent instantiation
		Types() = default;
	};

	// Convenience aliases for common scalar types
	using Typesf = Types<float>;
	using Typesd = Types<double>;

} // namespace sofa::component::cosserat::liegroups
