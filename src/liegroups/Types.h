#pragma once

#include <cmath>
#include <concepts>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
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
	template<typename _Scalar>
		requires std::is_floating_point_v<_Scalar>
	class Types {
	public:
		using Scalar = _Scalar;

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
		 * @brief Computes cos(x)/x with numerical stability for small x.
		 * Uses a Taylor expansion for x near zero to avoid division by zero.
		 * @param x The input value.
		 * @return The value of cos(x)/x.
		 */
		static Scalar cosc(const Scalar &x) noexcept {
			if (isZero(x, SMALL_ANGLE_THRESHOLD)) {
				// Taylor expansion: cos(x)/x ≈ 1 - x²/6 + x⁴/120 - ...
				const Scalar x2 = x * x;
				return Scalar(1) - x2 / Scalar(6) + x2 * x2 / Scalar(120);
			}
			return std::cos(x) / x;
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
		 * @brief Computes (1 - cos(x))/x^2 with numerical stability for small x.
		 * Uses a Taylor expansion for x near zero to avoid division by zero.
		 * @param x The input value.
		 * @return The value of (1 - cos(x))/x^2.
		 */
		static Scalar sinc2(const Scalar &x) noexcept {
			if (isZero(x, SMALL_ANGLE_THRESHOLD)) {
				// Taylor expansion: (1 - cos(x))/x² ≈ 1/2 - x²/24 + x⁴/720 - ...
				const Scalar x2 = x * x;
				return Scalar(0.5) - x2 / Scalar(24) + x2 * x2 / Scalar(720);
			}
			const Scalar x_sq = x * x;
			return (Scalar(1) - std::cos(x)) / x_sq;
		}

		/**
		 * @brief Computes (x - sin(x))/x^3 with numerical stability for small x.
		 * Useful for Lie group computations.
		 * @param x The input value.
		 * @return The value of (x - sin(x))/x^3.
		 */
		static Scalar sinc3(const Scalar &x) noexcept {
			if (isZero(x, SMALL_ANGLE_THRESHOLD)) {
				// Taylor expansion: (x - sin(x))/x³ ≈ 1/6 - x²/120 + x⁴/5040 - ...
				const Scalar x2 = x * x;
				return Scalar(1.0 / 6.0) - x2 / Scalar(120) + x2 * x2 / Scalar(5040);
			}
			const Scalar x_cubed = x * x * x;
			return (x - std::sin(x)) / x_cubed;
		}

		/**
		 * @brief Computes 1 - cos(x) with numerical stability.
		 * @param x The input value.
		 * @return The value of 1 - cos(x).
		 */
		static Scalar oneMinusCos(const Scalar &x) noexcept { return Scalar(1) - std::cos(x); }

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
		 * @brief Normalizes an angle to the range [-pi, pi].
		 * @param angle The angle to normalize in radians.
		 * @return The normalized angle in radians.
		 */
		static Scalar normalizeAngle(const Scalar &angle) noexcept {
			Scalar result = std::fmod(angle + PI, TWO_PI);
			if (result < Scalar(0)) {
				result += TWO_PI;
			}
			return result - PI;
		}

		/**
		 * @brief Computes the angle difference with proper wrapping.
		 * Always returns the shortest angular distance between two angles.
		 * @param a The first angle.
		 * @param b The second angle.
		 * @return The shortest angle difference in [-π, π].
		 */
		static Scalar angleDifference(Scalar a, Scalar b) {
			Scalar diff = std::fmod(a - b + PI, TWO_PI);
			if (diff < Scalar(0)) {
				diff += TWO_PI;
			}
			return diff - PI;
		}

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
		 * @brief Checks if a square matrix is approximately skew-symmetric.
		 * A matrix A is skew-symmetric if A = -A^T.
		 * @tparam N The dimension of the square matrix.
		 * @param mat The matrix to check.
		 * @param tol The tolerance for approximation. Defaults to `tolerance()`.
		 * @return True if the matrix is approximately skew-symmetric, false
		 * otherwise.
		 */
		template<int N>
		static bool isSkewSymmetric(const Matrix<N, N> &mat, const Scalar &tol = tolerance()) noexcept {
			return (mat + mat.transpose()).cwiseAbs().maxCoeff() <= tol;
		}

		/**
		 * @brief Checks if a square matrix is approximately symmetric.
		 * A matrix A is symmetric if A = A^T.
		 * @tparam N The dimension of the square matrix.
		 * @param mat The matrix to check.
		 * @param tol The tolerance for approximation. Defaults to `tolerance()`.
		 * @return True if the matrix is approximately symmetric, false otherwise.
		 */
		template<int N>
		static bool isSymmetric(const Matrix<N, N> &mat, const Scalar &tol = tolerance()) noexcept {
			return (mat - mat.transpose()).cwiseAbs().maxCoeff() <= tol;
		}

		/**
		 * @brief Checks if a square matrix is approximately orthogonal.
		 * A matrix A is orthogonal if A * A^T = I (identity matrix).
		 * @tparam N The dimension of the square matrix.
		 * @param mat The matrix to check.
		 * @param tol The tolerance for approximation. Defaults to `tolerance()`.
		 * @return True if the matrix is approximately orthogonal, false otherwise.
		 */
		template<int N>
		static bool isOrthogonal(const Matrix<N, N> &mat, const Scalar &tol = tolerance()) noexcept {
			const auto should_be_identity = mat * mat.transpose();
			return (should_be_identity - Matrix<N, N>::Identity()).cwiseAbs().maxCoeff() <= tol;
		}

		/**
		 * @brief Checks if a matrix is approximately special orthogonal (rotation matrix).
		 * A matrix is SO(n) if it's orthogonal and has determinant 1.
		 * @tparam N The dimension of the square matrix.
		 * @param mat The matrix to check.
		 * @param tol The tolerance for approximation. Defaults to `tolerance()`.
		 * @return True if the matrix is approximately in SO(n), false otherwise.
		 */
		template<int N>
		static bool isSpecialOrthogonal(const Matrix<N, N> &mat, const Scalar &tol = tolerance()) noexcept {
			return isOrthogonal(mat, tol) && isApprox(mat.determinant(), Scalar(1), tol);
		}

		/**
		 * @brief Extracts the skew-symmetric part of a square matrix.
		 * The skew-symmetric part of A is (A - A^T) / 2.
		 * @tparam N The dimension of the square matrix.
		 * @param mat The input matrix.
		 * @return The skew-symmetric part of the matrix.
		 */
		template<int N>
		static Matrix<N, N> skewPart(const Matrix<N, N> &mat) noexcept {
			return Scalar(0.5) * (mat - mat.transpose());
		}

		/**
		 * @brief Extracts the symmetric part of a square matrix.
		 * The symmetric part of A is (A + A^T) / 2.
		 * @tparam N The dimension of the square matrix.
		 * @param mat The input matrix.
		 * @return The symmetric part of the matrix.
		 */
		template<int N>
		static Matrix<N, N> symmetricPart(const Matrix<N, N> &mat) noexcept {
			return Scalar(0.5) * (mat + mat.transpose());
		}

		/**
		 * @brief Creates a 3x3 skew-symmetric matrix from a 3D vector.
		 * This is often used to represent the cross product as a matrix
		 * multiplication.
		 * @param v The 3D vector.
		 * @return The 3x3 skew-symmetric matrix.
		 */
		static Matrix3 skew3(const Vector3 &v) noexcept {
			Matrix3 result;
			result << Scalar(0), -v.z(), v.y(), v.z(), Scalar(0), -v.x(), -v.y(), v.x(), Scalar(0);
			return result;
		}

		/**
		 * @brief Extracts the 3D vector from a 3x3 skew-symmetric matrix.
		 * This is the inverse operation of `skew3`.
		 * @param mat The 3x3 skew-symmetric matrix.
		 * @return The 3D vector.
		 */
		static Vector3 unskew3(const Matrix3 &mat) noexcept { return Vector3(mat(2, 1), mat(0, 2), mat(1, 0)); }

		/**
		 * @brief Projects a matrix onto the Special Orthogonal group SO(n).
		 * Uses SVD decomposition to find the nearest rotation matrix.
		 * @tparam N The dimension of the square matrix.
		 * @param mat The input matrix.
		 * @return The nearest rotation matrix.
		 */
		template<int N>
		static Matrix<N, N> projectToSO(const Matrix<N, N> &mat) noexcept {
			Eigen::JacobiSVD<Matrix<N, N>> svd(mat, Eigen::ComputeFullU | Eigen::ComputeFullV);
			Matrix<N, N> result = svd.matrixU() * svd.matrixV().transpose();

			// Ensure determinant is 1 (not -1)
			if (result.determinant() < Scalar(0)) {
				auto U = svd.matrixU();
				U.col(N - 1) *= Scalar(-1);
				result = U * svd.matrixV().transpose();
			}

			return result;
		}

		/**
		 * @brief Computes the matrix logarithm for skew-symmetric matrices.
		 * Useful for SE(3) and SO(3) computations.
		 * @param mat The input skew-symmetric matrix.
		 * @return The matrix logarithm.
		 */
		static Matrix3 logSO3(const Matrix3 &R) noexcept {
			const Scalar trace = R.trace();
			const Scalar cos_angle = Scalar(0.5) * (trace - Scalar(1));

			if (cos_angle >= Scalar(1) - tolerance()) {
				// Near identity
				return Scalar(0.5) * (R - R.transpose());
			} else if (cos_angle <= Scalar(-1) + tolerance()) {
				// Near 180 degree rotation
				const Scalar angle = PI;
				// Find the axis (largest diagonal element)
				int i = 0;
				for (int j = 1; j < 3; ++j) {
					if (R(j, j) > R(i, i))
						i = j;
				}

				Vector3 axis = Vector3::Zero();
				axis(i) = std::sqrt((R(i, i) + Scalar(1)) / Scalar(2));

				for (int j = 0; j < 3; ++j) {
					if (j != i) {
						axis(j) = R(i, j) / (Scalar(2) * axis(i));
					}
				}

				return angle * skew3(axis);
			} else {
				const Scalar angle = std::acos(clamp(cos_angle, Scalar(-1), Scalar(1)));
				const Scalar sin_angle = std::sin(angle);

				return (angle / (Scalar(2) * sin_angle)) * (R - R.transpose());
			}
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
		 * @brief Generates a random unit vector (a vector with a norm of 1).
		 * @tparam N The dimension of the vector.
		 * @tparam Generator The type of the random number generator.
		 * @param gen The random number generator.
		 * @return A random unit vector.
		 */
		template<int N, typename Generator>
		static Vector<N> randomUnitVector(Generator &gen) noexcept {
			Vector<N> v = randomVector<N>(gen);
			const Scalar norm = v.norm();
			if (norm > epsilon()) {
				v /= norm;
			} else {
				v = Vector<N>::Unit(0); // Fallback to first basis vector
			}
			return v;
		}

		/**
		 * @brief Generates a random rotation matrix in SO(3).
		 * @tparam Generator The type of the random number generator.
		 * @param gen The random number generator.
		 * @return A random rotation matrix.
		 */
		template<typename Generator>
		static Matrix3 randomRotation(Generator &gen) noexcept {
			// Generate random axis
			Vector3 axis = randomUnitVector<3>(gen);

			// Generate random angle
			Scalar angle = randomScalar(gen, Scalar(-PI), Scalar(PI));

			// Create rotation matrix using Rodrigues' formula
			const Scalar cos_angle = std::cos(angle);
			const Scalar sin_angle = std::sin(angle);
			const Matrix3 K = skew3(axis);

			return Matrix3::Identity() + sin_angle * K + (Scalar(1) - cos_angle) * K * K;
		}

		/**
		 * @brief Performs SE(2) interpolation [angle, x, y].
		 * Uses SLERP for rotation and LERP for translation.
		 * @param start The starting SE(2) element.
		 * @param end The ending SE(2) element.
		 * @param t The interpolation parameter [0, 1].
		 * @return The interpolated SE(2) element.
		 */
		static Vector3 interpolateSE2(const Vector3 &start, const Vector3 &end, Scalar t) {
			Vector3 result;
			result(0) = slerpAngle(start(0), end(0), t);
			result(1) = start(1) + t * (end(1) - start(1));
			result(2) = start(2) + t * (end(2) - start(2));
			return result;
		}

		/**
		 * @brief Performs SE(3) interpolation using matrix representation.
		 * Uses SLERP for rotation and LERP for translation.
		 * @param T1 The starting SE(3) transformation.
		 * @param T2 The ending SE(3) transformation.
		 * @param t The interpolation parameter [0, 1].
		 * @return The interpolated SE(3) transformation.
		 */
		static Matrix4 interpolateSE3(const Matrix4 &T1, const Matrix4 &T2, Scalar t) {
			// Extract rotation and translation
			Matrix3 R1 = T1.template block<3, 3>(0, 0);
			Matrix3 R2 = T2.template block<3, 3>(0, 0);
			Vector3 t1 = T1.template block<3, 1>(0, 3);
			Vector3 t2 = T2.template block<3, 1>(0, 3);

			// Interpolate rotation using quaternion SLERP
			Quaternion q1(R1);
			Quaternion q2(R2);
			Quaternion q_interp = q1.slerp(t, q2);

			// Interpolate translation
			Vector3 t_interp = t1 + t * (t2 - t1);

			// Construct result
			Matrix4 result = Matrix4::Identity();
			result.template block<3, 3>(0, 0) = q_interp.toRotationMatrix();
			result.template block<3, 1>(0, 3) = t_interp;

			return result;
		}

		/**
		 * @brief Computes SE(3) interpolation using exponential coordinates.
		 * More mathematically principled than matrix interpolation.
		 * @param xi1 The starting SE(3) element in exponential coordinates.
		 * @param xi2 The ending SE(3) element in exponential coordinates.
		 * @param t The interpolation parameter [0, 1].
		 * @return The interpolated SE(3) element in exponential coordinates.
		 */
		static Vector6 interpolateSE3Exponential(const Vector6 &xi1, const Vector6 &xi2, Scalar t) {
			return xi1 + t * (xi2 - xi1);
		}

		/**
		 * @brief Computes the geodesic distance on SO(3) between two rotation matrices.
		 * @param R1 The first rotation matrix.
		 * @param R2 The second rotation matrix.
		 * @return The geodesic distance.
		 */
		static Scalar SO3Distance(const Matrix3 &R1, const Matrix3 &R2) {
			const Matrix3 R_diff = R1.transpose() * R2;
			const Scalar trace = R_diff.trace();
			const Scalar cos_angle = (trace - Scalar(1)) / Scalar(2);

			// Clamp to avoid numerical issues with acos
			const Scalar clamped = std::max(Scalar(-1), std::min(Scalar(1), cos_angle));
			return std::acos(clamped);
		}

		/**
		 * @brief Computes the geodesic distance on SE(3) between two transformations.
		 * @param T1 The first transformation matrix.
		 * @param T2 The second transformation matrix.
		 * @param w_rot Weight for rotation component.
		 * @param w_trans Weight for translation component.
		 * @return The weighted geodesic distance.
		 */
		static Scalar SE3Distance(const Matrix4 &T1, const Matrix4 &T2, Scalar w_rot = Scalar(1),
								  Scalar w_trans = Scalar(1)) {
			// Rotation distance
			const Matrix3 R1 = T1.template block<3, 3>(0, 0);
			const Matrix3 R2 = T2.template block<3, 3>(0, 0);
			const Scalar rot_dist = SO3Distance(R1, R2);

			// Translation distance
			const Vector3 t1 = T1.template block<3, 1>(0, 3);
			const Vector3 t2 = T2.template block<3, 1>(0, 3);
			const Scalar trans_dist = (t1 - t2).norm();

			return w_rot * rot_dist + w_trans * trans_dist;
		}

		/**
		 * @brief Generates a smooth trajectory between multiple SE(3) waypoints.
		 * @param waypoints Vector of SE(3) transformations.
		 * @param num_points Number of interpolation points per segment.
		 * @return Vector of interpolated SE(3) transformations.
		 */
		static std::vector<Matrix4> generateSE3Trajectory(const std::vector<Matrix4> &waypoints, int num_points = 10) {
			if (waypoints.size() < 2) {
				return waypoints;
			}

			std::vector<Matrix4> trajectory;
			trajectory.reserve((waypoints.size() - 1) * num_points + 1);

			for (size_t i = 0; i < waypoints.size() - 1; ++i) {
				for (int j = 0; j < num_points; ++j) {
					const Scalar t = Scalar(j) / Scalar(num_points);
					trajectory.push_back(interpolateSE3(waypoints[i], waypoints[i + 1], t));
				}
			}

			// Add final waypoint
			trajectory.push_back(waypoints.back());

			return trajectory;
		}

		/**
		 * @brief Computes the exponential map for SE(3) using exponential coordinates.
		 * @param xi The 6D exponential coordinates [rho, phi] where rho is translation, phi is rotation.
		 * @return The SE(3) transformation matrix.
		 */
		static Matrix4 SE3Exp(const Vector6 &xi) {
			const Vector3 rho = xi.template head<3>();
			const Vector3 phi = xi.template tail<3>();

			const Scalar angle = phi.norm();
			Matrix3 R;
			Matrix3 V;

			if (angle < SMALL_ANGLE_THRESHOLD) {
				// Near identity case
				R = Matrix3::Identity() + skew3(phi);
				V = Matrix3::Identity() + Scalar(0.5) * skew3(phi);
			} else {
				// General case
				const Vector3 axis = phi / angle;
				const Matrix3 K = skew3(axis);

				// Rodrigues' formula for rotation
				const Scalar sin_angle = std::sin(angle);
				const Scalar cos_angle = std::cos(angle);
				R = Matrix3::Identity() + sin_angle * K + (Scalar(1) - cos_angle) * K * K;

				// V matrix for translation
				V = Matrix3::Identity() + (Scalar(1) - cos_angle) / angle * K + (angle - sin_angle) / angle * K * K;
			}

			Matrix4 result = Matrix4::Identity();
			result.template block<3, 3>(0, 0) = R;
			result.template block<3, 1>(0, 3) = V * rho;

			return result;
		}

		/**
		 * @brief Computes the logarithm map for SE(3) to exponential coordinates.
		 * @param T The SE(3) transformation matrix.
		 * @return The 6D exponential coordinates.
		 */
		static Vector6 SE3Log(const Matrix4 &T) {
			const Matrix3 R = T.template block<3, 3>(0, 0);
			const Vector3 t = T.template block<3, 1>(0, 3);

			// Compute rotation part
			const Vector3 phi = SO3Log(R);
			const Scalar angle = phi.norm();

			Vector3 rho;
			if (angle < SMALL_ANGLE_THRESHOLD) {
				// Near identity case
				rho = t;
			} else {
				// General case
				const Vector3 axis = phi / angle;
				const Matrix3 K = skew3(axis);

				// Compute V^(-1)
				const Scalar sin_angle = std::sin(angle);
				const Scalar cos_angle = std::cos(angle);
				const Matrix3 V_inv = Matrix3::Identity() - Scalar(0.5) * K +
									  (Scalar(2) * sin_angle - angle * (Scalar(1) + cos_angle)) /
											  (Scalar(2) * angle * angle * sin_angle) * K * K;

				rho = V_inv * t;
			}

			Vector6 result;
			result.template head<3>() = rho;
			result.template tail<3>() = phi;

			return result;
		}

		/**
		 * @brief Computes the matrix logarithm for SO(3).
		 * @param R The rotation matrix.
		 * @return The axis-angle representation as a 3D vector.
		 */
		static Vector3 SO3Log(const Matrix3 &R) {
			const Scalar trace = R.trace();
			const Scalar cos_angle = (trace - Scalar(1)) / Scalar(2);

			if (cos_angle >= Scalar(1) - SMALL_ANGLE_THRESHOLD) {
				// Near identity
				return Scalar(0.5) * Vector3(R(2, 1) - R(1, 2), R(0, 2) - R(2, 0), R(1, 0) - R(0, 1));
			} else if (cos_angle <= Scalar(-1) + SMALL_ANGLE_THRESHOLD) {
				// Near 180 degree rotation
				const Scalar angle = PI;

				// Find the axis (largest diagonal element)
				int i = 0;
				for (int j = 1; j < 3; ++j) {
					if (R(j, j) > R(i, i))
						i = j;
				}

				Vector3 axis = Vector3::Zero();
				axis(i) = std::sqrt((R(i, i) + Scalar(1)) / Scalar(2));

				for (int j = 0; j < 3; ++j) {
					if (j != i) {
						axis(j) = R(i, j) / (Scalar(2) * axis(i));
					}
				}

				return angle * axis;
			} else {
				const Scalar angle = std::acos(std::max(Scalar(-1), std::min(Scalar(1), cos_angle)));
				const Scalar sin_angle = std::sin(angle);

				return (angle / (Scalar(2) * sin_angle)) *
					   Vector3(R(2, 1) - R(1, 2), R(0, 2) - R(2, 0), R(1, 0) - R(0, 1));
			}
		}

	private:
		// Make constructor private to prevent instantiation
		Types() = default;
	};

	// Convenience aliases for common scalar types
	using Typesf = Types<float>;
	using Typesd = Types<double>;

} // namespace sofa::component::cosserat::liegroups
