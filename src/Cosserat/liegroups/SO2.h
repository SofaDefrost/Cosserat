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

// #ifndef SOFA_COMPONENT_COSSERAT_LIEGROUPS_SO2_H
// #define SOFA_COMPONENT_COSSERAT_LIEGROUPS_SO2_H
#pragma once

#include <Cosserat/liegroups/Types.h>          // Then our type system
#include <Cosserat/liegroups/LieGroupBase.h>   // Then the base class interface
#include <Cosserat/liegroups/LieGroupBase.inl> // Then the base class interface
#include <eigen3/Eigen/Geometry>

namespace sofa::component::cosserat::liegroups {

/**
 * @brief Implementation of SO(2), the Special Orthogonal group in 2D
 *
 * This class implements the group of rotations in 2D space. Elements of SO(2)
 * are represented internally using complex numbers (cos θ + i sin θ), which
 * provides an efficient way to compose rotations and compute the exponential
 * map.
 *
 * The Lie algebra so(2) consists of skew-symmetric 2×2 matrices, which can be
 * identified with real numbers (the rotation angle).
 *
 * @tparam _Scalar The scalar type (must be a floating-point type)
 */
template <typename _Scalar, int _Dim = 2>
class SO2 : public LieGroupBase<_Scalar, std::integral_constant<int, _Dim>, _Dim, _Dim>{
public:
  using Types = Types<_Scalar>;
  using Base = LieGroupBase<_Scalar, std::integral_constant<int, _Dim>, _Dim, _Dim>;
  using Scalar = typename Types::Scalar;
  using Vector = typename Types::Vector2;
  using Matrix = typename Types::Matrix2;
  using TangentVector = typename Types::TangentVector2;
  using AdjointMatrix = typename Types::AdjointMatrix2;

  static constexpr int Dim = 2;
  using Complex =
      typename Types::Vector2; // Represents complex number as 2D vector

  /**
   * @brief Default constructor creates identity rotation (angle = 0)
   */
  SO2() : m_angle(Scalar(0)) { updateComplex(); }

  /**
   * @brief Construct from angle (in radians)
   */
  explicit SO2(const Scalar &angle) : m_angle(Types::normalizeAngle(angle)) {
    updateComplex();
  }

  /**
   * @brief Copy constructor
   */
  SO2(const SO2 &other) : m_angle(other.m_angle), m_complex(other.m_complex) {}

  /**
   * @brief Assignment operator
   */
  SO2 &operator=(const SO2 &other) {
    if (this != &other) {
      m_angle = other.m_angle;
      m_complex = other.m_complex;
    }
    return *this;
  }

  /**
   * @brief Group composition (rotation composition)
   */
  SO2 operator*(const SO2 &other) const {
    // Complex multiplication: (a + bi)(c + di) = (ac - bd) + (ad + bc)i
    Scalar real_part =
        m_complex(0) * other.m_complex(0) - m_complex(1) * other.m_complex(1);
    Scalar imag_part =
        m_complex(0) * other.m_complex(1) + m_complex(1) * other.m_complex(0);

    // Convert back to angle using atan2
    Scalar result_angle = std::atan2(imag_part, real_part);
    return SO2(result_angle);
  }

  /**
   * @brief In-place composition
   */
  SO2 &operator*=(const SO2 &other) {
    *this = *this * other;
    return *this;
  }

  /**
   * @brief Inverse element (opposite rotation)
   */
  SO2 inverse() const override { return SO2(-m_angle); }

  /**
   * @brief Exponential map (angle to rotation)
   * For SO(2), this is just the angle itself as rotation
   */
  static SO2 exp(const TangentVector &algebra_element) {
    return SO2(algebra_element[0]);
  }

  /**
   * @brief Logarithmic map (rotation to angle)
   */
  TangentVector log() const override {
    TangentVector result;
    result[0] = m_angle;
    return result;
  }

  /**
   * @brief Adjoint representation
   * For SO(2), this is simply the identity matrix as the group is abelian
   */
  AdjointMatrix adjoint() const override { return AdjointMatrix::Identity(); }

  /**
   * @brief Group action on a point (rotate the point)
   */
  Vector act(const Vector &point) const override {
    Vector result;
    result(0) = m_complex(0) * point(0) - m_complex(1) * point(1);
    result(1) = m_complex(1) * point(0) + m_complex(0) * point(1);
    return result;
  }

  /**
   * @brief Left Jacobian matrix for exponential coordinates
   * For SO(2), this is simply the identity matrix
   */
  AdjointMatrix leftJacobian() const { return AdjointMatrix::Identity(); }

  /**
   * @brief Right Jacobian matrix for exponential coordinates
   * For SO(2), this is simply the identity matrix
   */
  AdjointMatrix rightJacobian() const { return AdjointMatrix::Identity(); }

  /**
   * @brief Check if approximately equal to another rotation
   */
  bool isApprox(const SO2 &other, const Scalar &eps = Types::epsilon()) const {
    return Types::isZero(Types::normalizeAngle(m_angle - other.m_angle), eps);
  }

  /**
   * @brief Check if this is approximately the identity element
   */
  bool isIdentity(const Scalar &eps = Types::epsilon()) const {
    return Types::isZero(m_angle, eps);
  }

  /**
   * @brief Get the identity element (zero rotation)
   */
  static SO2 identity() { return SO2(); }

  /**
   * @brief Generate a random rotation
   */
  static SO2 random() {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    std::uniform_real_distribution<Scalar> dis(-Types::pi(), Types::pi());
    return SO2(dis(gen));
  }

  /**
   * @brief Get the dimension of the Lie algebra (1 for SO(2))
   */
  static constexpr int algebraDimension() { return 1; }

  /**
   * @brief Get the dimension of the space the group acts on (2 for SO(2))
   */
  static constexpr int actionDimension() { return 2; }

  /**
   * @brief Get the rotation angle in radians
   */
  Scalar angle() const { return m_angle; }

  /**
   * @brief Set the rotation angle in radians
   */
  void setAngle(const Scalar &angle) {
    m_angle = Types::normalizeAngle(angle);
    updateComplex();
  }

  // Required CRTP methods:
  static SO2<Scalar> computeIdentity() noexcept;
  SO2<Scalar> computeInverse() const;

  /**
   * @brief Get the rotation matrix representation
   */
  Matrix matrix() const {
    Matrix R;
    R << m_complex(0), -m_complex(1), m_complex(1), m_complex(0);
    return R;
  }

  /**
   * @brief Get the complex number representation (cos θ, sin θ)
   */
  const Complex &complex() const { return m_complex; }

  /**
   * @brief Interpolate between two rotations using SLERP
   * @param other Target rotation
   * @param t Interpolation parameter [0,1]
   * @return Interpolated rotation
   */
  SO2 slerp(const SO2 &other, const Scalar &t) const {
    // For SO(2), SLERP reduces to linear interpolation of angles
    // with proper handling of angle wrapping
    Scalar angle_diff = Types::normalizeAngle(other.m_angle - m_angle);
    return SO2(m_angle + t * angle_diff);
  }

  /**
   * @brief Get a unit vector in the direction of the rotation
   */
  Vector direction() const { return m_complex; }

  /**
   * @brief Rotate a vector by 90 degrees counterclockwise
   */
  Vector perpendicular() const {
    Vector result;
    result(0) = -m_complex(1); // -sin θ
    result(1) = m_complex(0);  //  cos θ
    return result;
  }

  /**
   * @brief Convert to string representation
   */
  std::string toString() const {
    std::ostringstream oss;
    oss << "SO2(angle=" << m_angle << " rad, "
        << (m_angle * Scalar(180) / Types::pi()) << " deg)";
    return oss.str();
  }

  /**
   * @brief Stream output operator
   */
  friend std::ostream &operator<<(std::ostream &os, const SO2 &rotation) {
    return os << rotation.toString();
  }

private:
  Scalar m_angle;    ///< The rotation angle in radians (normalized to [-π, π])
  Complex m_complex; ///< Complex number representation (cos θ, sin θ)

  /**
   * @brief Update complex number representation from angle
   */
  void updateComplex() {
    m_complex(0) = std::cos(m_angle); // Real part
    m_complex(1) = std::sin(m_angle); // Imaginary part
  }
};

// Type aliases for common scalar types
using SO2d = SO2<double>;
using SO2f = SO2<float>;

/**
 * @brief Helper function to create SO2 from degrees
 */
template <typename Scalar> SO2<Scalar> SO2FromDegrees(const Scalar &degrees) {
  return SO2<Scalar>(degrees * Types<Scalar>::pi() / Scalar(180));
}

/**
 * @brief Helper function to convert SO2 angle to degrees
 */
template <typename Scalar> Scalar SO2ToDegrees(const SO2<Scalar> &rotation) {
  return rotation.angle() * Scalar(180) / Types<Scalar>::pi();
}

} // namespace sofa::component::cosserat::liegroups

// #endif // SOFA_COMPONENT_COSSERAT_LIEGROUPS_SO2_H
