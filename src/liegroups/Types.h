// This file provides fundamental type definitions and utility functions for Lie
// group operations, primarily using Eigen for linear algebra.

// #ifndef SOFA_COMPONENT_COSSERAT_LIEGROUPS_TYPES_H
// #define SOFA_COMPONENT_COSSERAT_LIEGROUPS_TYPES_H
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
template <int N> using IntConst = std::integral_constant<int, N>;

/**
 * @brief Type definitions and utilities for Lie group implementations
 *
 * This class provides type aliases and utility functions for different
 * scalar types used in Lie group computations.
 */
template <typename _Scalar>
  requires std::is_floating_point_v<_Scalar>
struct Types {
  using Scalar = _Scalar;

  // Eigen type aliases
  template <int Rows, int Cols = Rows>
  using Matrix = Eigen::Matrix<Scalar, Rows, Cols>;

  template <int Size> using Vector = Eigen::Matrix<Scalar, Size, 1>;

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

  /**
   * @brief Get machine epsilon for the scalar type
   */
  static constexpr Scalar epsilon() noexcept {
    return std::numeric_limits<Scalar>::epsilon();
  }

  /**
   * @brief Get a small tolerance value for comparisons
   */
  static constexpr Scalar tolerance() noexcept {
    return Scalar(100) * epsilon();
  }

  static constexpr Scalar SMALL_ANGLE_THRESHOLD = Scalar(1e-4);

  /**
   * @brief Check if a value is effectively zero within a given tolerance.
   * @param value The scalar value to check.
   * @param tol The tolerance for considering a value as zero. Defaults to
   * `tolerance()`.
   * @return True if the absolute value of `value` is less than or equal to
   * `tol`, false otherwise.
   */
  static constexpr bool isZero(const Scalar &value,
                               const Scalar &tol = tolerance()) noexcept {
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
  static constexpr bool isApprox(const Scalar &a, const Scalar &b,
                                 const Scalar &tol = tolerance()) noexcept {
    return std::abs(a - b) <= tol;
  }

  /**
   * @brief Computes cos(x)/x with numerical stability for small x.
   * Uses a Taylor expansion for x near zero to avoid division by zero.
   * @param x The input value.
   * @return The value of cos(x)/x.
   */
  static Scalar cosc(const Scalar &x) noexcept {
    if (isZero(x)) {
      return Scalar(1);
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
    if (isZero(x)) {
      return Scalar(1);
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
    if (isZero(x)) {
      return Scalar(0.5);
    }
    const Scalar x_sq = x * x;
    return (Scalar(1) - std::cos(x)) / x_sq;
  }

  /**
   * @brief Computes the arc tangent of y/x, using the signs of both arguments
   * to determine the correct quadrant.
   * @param y The y-coordinate.
   * @param x The x-coordinate.
   * @return The angle in radians between the positive x-axis and the point (x,
   * y).
   */
  static Scalar atan2(const Scalar &y, const Scalar &x) noexcept {
    return std::atan2(y, x);
  }

  /**
   * @brief Computes the safe square root of a non-negative number.
   * Handles negative inputs gracefully by returning 0 for negative values.
   * @param x The input value.
   * @return The square root of x, or 0 if x is negative.
   */
  static Scalar safeSqrt(const Scalar &x) noexcept {
    return std::sqrt(std::max(Scalar(0), x));
  }

  /**
   * @brief Normalizes an angle to the range [-pi, pi].
   * @param angle The angle to normalize in radians.
   * @return The normalized angle in radians.
   */
  static Scalar normalizeAngle(const Scalar &angle) noexcept {
    Scalar result = std::fmod(angle + Scalar(M_PI), Scalar(2 * M_PI));
    if (result < Scalar(0)) {
      result += Scalar(2 * M_PI);
    }
    return result - Scalar(M_PI);
  }

  /**
   * @brief Clamps a value between a minimum and maximum value.
   * @param value The value to clamp.
   * @param min_val The minimum allowed value.
   * @param max_val The maximum allowed value.
   * @return The clamped value.
   */
  static constexpr Scalar clamp(const Scalar &value, const Scalar &min_val,
                                const Scalar &max_val) noexcept {
    return std::max(min_val, std::min(max_val, value));
  }

  /**
   * @brief Performs linear interpolation between two scalar values.
   * @param a The start value.
   * @param b The end value.
   * @param t The interpolation parameter (typically between 0 and 1).
   * @return The interpolated value.
   */
  static constexpr Scalar lerp(const Scalar &a, const Scalar &b,
                               const Scalar &t) noexcept {
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
  template <int N>
  static bool isSkewSymmetric(const Matrix<N, N> &mat,
                              const Scalar &tol = tolerance()) noexcept {
    const auto diff = mat + mat.transpose();
    return diff.cwiseAbs().maxCoeff() <= tol;
  }

  /**
   * @brief Checks if a square matrix is approximately symmetric.
   * A matrix A is symmetric if A = A^T.
   * @tparam N The dimension of the square matrix.
   * @param mat The matrix to check.
   * @param tol The tolerance for approximation. Defaults to `tolerance()`.
   * @return True if the matrix is approximately symmetric, false otherwise.
   */
  template <int N>
  static bool isSymmetric(const Matrix<N, N> &mat,
                          const Scalar &tol = tolerance()) noexcept {
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
  template <int N>
  static bool isOrthogonal(const Matrix<N, N> &mat,
                           const Scalar &tol = tolerance()) noexcept {
    const auto should_be_identity = mat * mat.transpose();
    const auto diff = should_be_identity - Matrix<N, N>::Identity();
    return diff.cwiseAbs().maxCoeff() <= tol;
  }

  /**
   * @brief Extracts the skew-symmetric part of a square matrix.
   * The skew-symmetric part of A is (A - A^T) / 2.
   * @tparam N The dimension of the square matrix.
   * @param mat The input matrix.
   * @return The skew-symmetric part of the matrix.
   */
  template <int N>
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
  template <int N>
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
    result << Scalar(0), -v(2), v(1), v(2), Scalar(0), -v(0), -v(1), v(0),
        Scalar(0);
    return result;
  }

  /**
   * @brief Extracts the 3D vector from a 3x3 skew-symmetric matrix.
   * This is the inverse operation of `skew3`.
   * @param mat The 3x3 skew-symmetric matrix.
   * @return The 3D vector.
   */
  static Vector3 unskew3(const Matrix3 &mat) noexcept {
    return Vector3(mat(2, 1), mat(0, 2), mat(1, 0));
  }

  /**
   * @brief Generates a random scalar value within the range [0, 1].
   * @tparam Generator The type of the random number generator.
   * @param gen The random number generator.
   * @return A random scalar value.
   */
  template <typename Generator>
  static Scalar randomScalar(Generator &gen) noexcept {
    std::uniform_real_distribution<Scalar> dist(Scalar(0), Scalar(1));
    return dist(gen);
  }

  /**
   * @brief Generates a random vector with components in the range [-1, 1].
   * @tparam N The dimension of the vector.
   * @tparam Generator The type of the random number generator.
   * @param gen The random number generator.
   * @return A random vector.
   */
  template <int N, typename Generator>
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
  template <int N, typename Generator>
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
};

// Convenience aliases for common scalar types
using Typesf = Types<float>;
using Typesd = Types<double>;

} // namespace sofa::component::cosserat::liegroups

// #endif // SOFA_COMPONENT_COSSERAT_LIEGROUPS_TYPES_H
