// Enhanced Lie Group Utilities
// This file provides specialized utility functions for Lie groups that complement
// the basic Types.h functionality with higher-level operations.

#ifndef COSSERAT_LIEGROUPS_ENHANCED_UTILS_H
#define COSSERAT_LIEGROUPS_ENHANCED_UTILS_H

#include <cmath>
#include <eigen3/Eigen/Dense>
#include <limits>
#include <type_traits>
#include <vector>

namespace sofa::component::cosserat::liegroups {

	/**
	 * @brief Advanced utility functions for Lie groups that complement Types.h
	 *
	 * This class provides higher-level operations for Lie group computations,
	 * including interpolation, path planning, and geometric utilities.
	 */
	template<typename _Scalar>
		requires std::is_floating_point_v<_Scalar>
	class LieGroupUtils {
	public:
		using Scalar = _Scalar;

		// Eigen type aliases
		template<int Rows, int Cols = Rows>
		using Matrix = Eigen::Matrix<Scalar, Rows, Cols>;
		template<int Size>
		using Vector = Eigen::Matrix<Scalar, Size, 1>;

		using Matrix2 = Matrix<2, 2>;
		using Matrix3 = Matrix<3, 3>;
		using Matrix4 = Matrix<4, 4>;

		using Vector2 = Vector<2>;
		using Vector3 = Vector<3>;
		using Vector4 = Vector<4>;
		using Vector6 = Vector<6>;

		using Quaternion = Eigen::Quaternion<Scalar>;
		using AngleAxis = Eigen::AngleAxis<Scalar>;
		using Transform3 = Eigen::Transform<Scalar, 3, Eigen::Affine>;
		using Isometry3 = Eigen::Transform<Scalar, 3, Eigen::Isometry>;

		static constexpr Scalar PI = Scalar(M_PI);
		static constexpr Scalar TWO_PI = Scalar(2 * M_PI);
		static constexpr Scalar EPSILON = std::numeric_limits<Scalar>::epsilon() * Scalar(100);

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
		 * @brief Normalizes an angle to [-π, π] range.
		 * @param angle The angle to normalize.
		 * @return The normalized angle.
		 */
		static Scalar normalizeAngle(Scalar angle) {
			Scalar result = std::fmod(angle + PI, TWO_PI);
			if (result < Scalar(0)) {
				result += TWO_PI;
			}
			return result - PI;
		}

		/**
		 * @brief Safely normalizes a vector, handling near-zero cases.
		 * @tparam Derived The Eigen derived type.
		 * @param v The vector to normalize.
		 * @param fallback Optional fallback unit vector if normalization fails.
		 * @return The normalized vector or fallback.
		 */
		template<typename Derived>
		static auto safeNormalize(const Eigen::MatrixBase<Derived> &v,
								  const Eigen::MatrixBase<Derived> &fallback = Eigen::MatrixBase<Derived>::Zero()) {
			using VectorType = typename Derived::PlainObject;
			const Scalar norm = v.norm();

			if (norm < EPSILON) {
				if (fallback.norm() > EPSILON) {
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

			if (norm_squared < EPSILON) {
				return VectorType::Zero(v.rows());
			}

			return (onto * (v.dot(onto) / norm_squared)).eval();
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

			if (angle < EPSILON) {
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
			if (angle < EPSILON) {
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

	private:
		/**
		 * @brief Creates a 3x3 skew-symmetric matrix from a 3D vector.
		 */
		static Matrix3 skew3(const Vector3 &v) {
			Matrix3 result;
			result << Scalar(0), -v.z(), v.y(), v.z(), Scalar(0), -v.x(), -v.y(), v.x(), Scalar(0);
			return result;
		}

		/**
		 * @brief Computes the matrix logarithm for SO(3).
		 */
		static Vector3 SO3Log(const Matrix3 &R) {
			const Scalar trace = R.trace();
			const Scalar cos_angle = (trace - Scalar(1)) / Scalar(2);

			if (cos_angle >= Scalar(1) - EPSILON) {
				// Near identity
				return Scalar(0.5) * Vector3(R(2, 1) - R(1, 2), R(0, 2) - R(2, 0), R(1, 0) - R(0, 1));
			} else if (cos_angle <= Scalar(-1) + EPSILON) {
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

		// Make constructor private to prevent instantiation
		LieGroupUtils() = default;
	};

	// Convenience aliases
	using LieGroupUtilsf = LieGroupUtils<float>;
	using LieGroupUtilsd = LieGroupUtils<double>;

} // namespace sofa::component::cosserat::liegroups

#endif // COSSERAT_LIEGROUPS_ENHANCED_UTILS_H
