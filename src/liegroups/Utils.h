// This file provides utility functions for Lie groups, including numerical
// stability helpers and interpolation methods.

#ifndef COSSERAT_LIEGROUPS_UTILS_H
#define COSSERAT_LIEGROUPS_UTILS_H

#include <cmath>
#include <eigen3/Eigen/Dense>
#include <limits>
#include <type_traits>

namespace sofa::component::cosserat::liegroups {

/**
 * Utility functions for Lie groups.
 */
template <typename Scalar> struct LieGroupUtils {
  using Matrix2 = Eigen::Matrix<Scalar, 2, 2>;
  using Vector2 = Eigen::Matrix<Scalar, 2, 1>;
  using Vector3 = Eigen::Matrix<Scalar, 3, 1>;

  static constexpr Scalar epsilon =
      std::numeric_limits<Scalar>::epsilon() * 100;
  static constexpr Scalar pi = M_PI;
  static constexpr Scalar two_pi = 2 * M_PI;

  /**
   * @brief Normalizes an angle to the range [-π, π].
   * @param angle The angle to normalize.
   * @return The normalized angle.
   */
  static Scalar normalizeAngle(Scalar angle) {
    // Normalize angle to [-π, π]
    angle = std::fmod(angle, two_pi);
    if (angle > pi) {
      angle -= two_pi;
    } else if (angle < -pi) {
      angle += two_pi;
    }
    return angle;
  }

  /**
   * @brief Computes sinc(x) = sin(x)/x with numerical stability for small x.
   * Uses a Taylor series approximation for small angles to avoid division by
   * zero.
   * @param x The input value.
   * @return The value of sinc(x).
   */
  static Scalar sinc(Scalar x) {
    if (std::abs(x) < epsilon) {
      // Taylor series approximation for small angles
      // sinc(x) ≈ 1 - x²/6 + x⁴/120 - ...
      return Scalar(1) - x * x / Scalar(6);
    }
    return std::sin(x) / x;
  }

  /**
   * @brief Computes the difference between two angles with proper wrapping.
   * The result is normalized to [-π, π].
   * @param a The first angle.
   * @param b The second angle.
   * @return The difference between the angles.
   */
  static Scalar angleDifference(Scalar a, Scalar b) {
    return normalizeAngle(a - b);
  }

  /**
   * @brief Performs linear interpolation between two scalars.
   * @param a The start value.
   * @param b The end value.
   * @param t The interpolation parameter (typically between 0 and 1).
   * @return The interpolated value.
   */
  static Scalar lerp(Scalar a, Scalar b, Scalar t) { return a + t * (b - a); }

  /**
   * @brief Performs spherical linear interpolation (SLERP) between two angles.
   * Ensures the shortest path is taken.
   * @param a The start angle.
   * @param b The end angle.
   * @param t The interpolation parameter (typically between 0 and 1).
   * @return The interpolated angle.
   */
  static Scalar slerpAngle(Scalar a, Scalar b, Scalar t) {
    // Ensure shortest path
    Scalar diff = angleDifference(b, a);
    return normalizeAngle(a + t * diff);
  }

  /**
   * @brief Computes the bi-invariant distance between two angles (as SO(2)
   * elements).
   * @param a The first angle.
   * @param b The second angle.
   * @return The distance between the angles.
   */
  static Scalar angleDistance(Scalar a, Scalar b) {
    return std::abs(angleDifference(a, b));
  }

  /**
   * @brief Checks if an angle is near zero within the defined epsilon.
   * @param angle The angle to check.
   * @return True if the angle is near zero, false otherwise.
   */
  static bool isAngleNearZero(Scalar angle) {
    return std::abs(angle) < epsilon;
  }

  /**
   * @brief Checks if two angles are nearly equal within the defined epsilon,
   * considering wrapping.
   * @param a The first angle.
   * @param b The second angle.
   * @return True if the angles are nearly equal, false otherwise.
   */
  static bool areAnglesNearlyEqual(Scalar a, Scalar b) {
    return angleDistance(a, b) < epsilon;
  }

  /**
   * @brief Numerically stable computation of 1 - cos(x) for small x.
   * Uses a Taylor series approximation for small angles.
   * @param x The input value.
   * @return The value of 1 - cos(x).
   */
  static Scalar oneMinusCos(Scalar x) {
    if (std::abs(x) < epsilon) {
      // Taylor series approximation for small angles
      // 1 - cos(x) ≈ x²/2 - x⁴/24 + ...
      return (x * x) / Scalar(2);
    }
    return Scalar(1) - std::cos(x);
  }

  /**
   * @brief Safely normalizes a vector, handling near-zero vectors.
   * If the norm is less than epsilon, returns a zero vector.
   * @tparam Derived The Eigen derived type of the vector.
   * @param v The vector to normalize.
   * @return The normalized vector.
   */
  template <typename Derived>
  static Eigen::Matrix<typename Derived::Scalar, Derived::RowsAtCompileTime, 1>
  safeNormalize(const Eigen::MatrixBase<Derived> &v) {
    using VectorType =
        Eigen::Matrix<typename Derived::Scalar, Derived::RowsAtCompileTime, 1>;
    typename Derived::Scalar norm = v.norm();
    if (norm < epsilon) {
      // Return zero vector or throw exception based on your preference
      return VectorType::Zero(v.rows());
    }
    return v / norm;
  }

  /**
   * @brief Projects a vector onto another vector.
   * Handles cases where the vector to project onto is near-zero.
   * @tparam Derived1 The Eigen derived type of the vector to be projected.
   * @tparam Derived2 The Eigen derived type of the vector to project onto.
   * @param v The vector to be projected.
   * @param onto The vector to project onto.
   * @return The projected vector.
   */
  template <typename Derived1, typename Derived2>
  static Eigen::Matrix<typename Derived1::Scalar, Derived1::RowsAtCompileTime,
                       1>
  projectVector(const Eigen::MatrixBase<Derived1> &v,
                const Eigen::MatrixBase<Derived2> &onto) {
    using VectorType = Eigen::Matrix<typename Derived1::Scalar,
                                     Derived1::RowsAtCompileTime, 1>;
    typename Derived2::Scalar norm_squared = onto.squaredNorm();
    if (norm_squared < epsilon) {
      return VectorType::Zero(v.rows());
    }
    return onto * (v.dot(onto) / norm_squared);
  }

  /**
   * @brief Performs path interpolation between two SE(2) elements represented
   * as [angle, x, y]. Interpolates angle using SLERP and translation using
   * LERP.
   * @param start The starting SE(2) element.
   * @param end The ending SE(2) element.
   * @param t The interpolation parameter (typically between 0 and 1).
   * @return The interpolated SE(2) element as a Vector3.
   */
  static Vector3 interpolateSE2Path(const Vector3 &start, const Vector3 &end,
                                    Scalar t) {
    Vector3 result;
    // Interpolate angle using SLERP
    result[0] = slerpAngle(start[0], end[0], t);
    // Interpolate translation using LERP
    result[1] = lerp(start[1], end[1], t);
    result[2] = lerp(start[2], end[2], t);
    return result;
  }
};

} // namespace sofa::component::cosserat::liegroups

#endif // COSSERAT_LIEGROUPS_UTILS_H
