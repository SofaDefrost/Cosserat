// #ifndef SOFA_COMPONENT_COSSERAT_LIEGROUPS_TYPES_H
// #define SOFA_COMPONENT_COSSERAT_LIEGROUPS_TYPES_H
#pragma once

#include <random>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
#include <cmath>
#include <concepts>
#include <limits>
#include <type_traits>

namespace sofa::component::cosserat::liegroups {

// Helper type for compile-time integer constants (for template parameters)
template <int N>
using IntConst = std::integral_constant<int, N>;

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

  /**
   * @brief Check if a value is effectively zero
   */
  static constexpr bool isZero(const Scalar &value,
                               const Scalar &tol = tolerance()) noexcept {
    return std::abs(value) <= tol;
  }

  /**
   * @brief Check if two values are approximately equal
   */
  static constexpr bool isApprox(const Scalar &a, const Scalar &b,
                                 const Scalar &tol = tolerance()) noexcept {
    return std::abs(a - b) <= tol;
  }

  /**
   * @brief Compute cos(x)/x with numerical stability for small x
   */
  static Scalar cosc(const Scalar &x) noexcept {
    if (isZero(x)) {
      return Scalar(1);
    }
    return std::cos(x) / x;
  }

  /**
   * @brief Compute sin(x)/x with numerical stability for small x
   */
  static Scalar sinc(const Scalar &x) noexcept {
    if (isZero(x)) {
      return Scalar(1);
    }
    return std::sin(x) / x;
  }

  /**
   * @brief Compute (1 - cos(x))/x^2 with numerical stability for small x
   */
  static Scalar sinc2(const Scalar &x) noexcept {
    if (isZero(x)) {
      return Scalar(0.5);
    }
    const Scalar x_sq = x * x;
    return (Scalar(1) - std::cos(x)) / x_sq;
  }

  /**
   * @brief Compute atan2 with better numerical properties
   */
  static Scalar atan2(const Scalar &y, const Scalar &x) noexcept {
    return std::atan2(y, x);
  }

  /**
   * @brief Safe square root that handles negative inputs gracefully
   */
  static Scalar safeSqrt(const Scalar &x) noexcept {
    return std::sqrt(std::max(Scalar(0), x));
  }

  /**
   * @brief Normalize angle to [-pi, pi]
   */
  static Scalar normalizeAngle(const Scalar &angle) noexcept {
    Scalar result = std::fmod(angle + Scalar(M_PI), Scalar(2 * M_PI));
    if (result < Scalar(0)) {
      result += Scalar(2 * M_PI);
    }
    return result - Scalar(M_PI);
  }

  /**
   * @brief Clamp value between min and max
   */
  static constexpr Scalar clamp(const Scalar &value, const Scalar &min_val,
                                const Scalar &max_val) noexcept {
    return std::max(min_val, std::min(max_val, value));
  }

  /**
   * @brief Linear interpolation
   */
  static constexpr Scalar lerp(const Scalar &a, const Scalar &b,
                               const Scalar &t) noexcept {
    return a + t * (b - a);
  }

  /**
   * @brief Check if a matrix is approximately skew-symmetric
   */
  template <int N>
  static bool isSkewSymmetric(const Matrix<N, N> &mat,
                              const Scalar &tol = tolerance()) noexcept {
    const auto diff = mat + mat.transpose();
    return diff.cwiseAbs().maxCoeff() <= tol;
  }

  /**
   * @brief Check if a matrix is approximately symmetric
   */
  template <int N>
  static bool isSymmetric(const Matrix<N, N> &mat,
                          const Scalar &tol = tolerance()) noexcept {
    const auto diff = mat - mat.transpose();
    return diff.cwiseAbs().maxCoeff() <= tol;
  }

  /**
   * @brief Check if a matrix is approximately orthogonal
   */
  template <int N>
  static bool isOrthogonal(const Matrix<N, N> &mat,
                           const Scalar &tol = tolerance()) noexcept {
    const auto should_be_identity = mat * mat.transpose();
    const auto diff = should_be_identity - Matrix<N, N>::Identity();
    return diff.cwiseAbs().maxCoeff() <= tol;
  }

  /**
   * @brief Extract the skew-symmetric part of a matrix
   */
  template <int N>
  static Matrix<N, N> skewPart(const Matrix<N, N> &mat) noexcept {
    return Scalar(0.5) * (mat - mat.transpose());
  }

  /**
   * @brief Extract the symmetric part of a matrix
   */
  template <int N>
  static Matrix<N, N> symmetricPart(const Matrix<N, N> &mat) noexcept {
    return Scalar(0.5) * (mat + mat.transpose());
  }

  /**
   * @brief Create a 3x3 skew-symmetric matrix from a 3D vector
   */
  static Matrix3 skew3(const Vector3 &v) noexcept {
    Matrix3 result;
    result << Scalar(0), -v(2), v(1), v(2), Scalar(0), -v(0), -v(1), v(0),
        Scalar(0);
    return result;
  }

  /**
   * @brief Extract 3D vector from a 3x3 skew-symmetric matrix
   */
  static Vector3 unskew3(const Matrix3 &mat) noexcept {
    return Vector3(mat(2, 1), mat(0, 2), mat(1, 0));
  }

  /**
   * @brief Generate a random scalar in [0, 1]
   */
  template <typename Generator>
  static Scalar randomScalar(Generator &gen) noexcept {
    std::uniform_real_distribution<Scalar> dist(Scalar(0), Scalar(1));
    return dist(gen);
  }

  /**
   * @brief Generate a random vector with components in [-1, 1]
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
   * @brief Generate a random unit vector
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
