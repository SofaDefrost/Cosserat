// This file defines the SO2 (Special Orthogonal group in 2D) class,
// representing 2D rotations.

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

#include "LieGroupBase.h"   // Then the base class interface
#include "LieGroupBase.inl" // Then the base class interface
#include "Types.h"          // Then our type system
#include <cmath>
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
template <typename _Scalar>
class SO2 : public LieGroupBase<SO2<_Scalar>, _Scalar, 2, 1, 2> {
public:
  using Base = LieGroupBase<SO2<_Scalar>, _Scalar, 2, 1, 2>;
  using Scalar = typename Base::Scalar;
  using Vector = typename Base::Vector;
  using Matrix = typename Base::Matrix;
  using TangentVector = typename Base::TangentVector;
  using AdjointMatrix = typename Base::AdjointMatrix;

  static constexpr int Dim = 2;
  using Complex =
      Eigen::Vector2<_Scalar>; // Represents complex number as 2D vector

  /**
   * @brief Default constructor creates identity rotation (angle = 0)
   */
  SO2() : m_angle(Scalar(0)) { updateComplex(); }

  /**
   * @brief Construct from angle (in radians)
   * @param angle The rotation angle in radians.
   */
  explicit SO2(const Scalar &angle)
      : m_angle(Types<_Scalar>::normalizeAngle(angle)) {
    updateComplex();
  }

  /**
   * @brief Copy constructor
   * @param other The SO2 object to copy from.
   */
  SO2(const SO2 &other) : m_angle(other.m_angle), m_complex(other.m_complex) {}

  /**
   * @brief Get the rotation angle in radians.
   * @return The rotation angle in radians.
   */
  Scalar angle() const { return m_angle; }

  /**
   * @brief Set the rotation angle in radians.
   * The angle will be normalized to [-π, π].
   * @param angle The new rotation angle in radians.
   */
  void setAngle(const Scalar &angle) {
    m_angle = Types<_Scalar>::normalizeAngle(angle);
    updateComplex();
  }

  // Required CRTP methods:
  /**
   * @brief Computes the identity element of the SO(2) group.
   * @return The identity SO(2) element (angle = 0).
   */
  static SO2<Scalar> computeIdentity() noexcept { return SO2(Scalar(0)); }
  /**
   * @brief Computes the inverse of the current SO(2) element.
   * @return The inverse SO(2) element (negative angle).
   */
  SO2<Scalar> computeInverse() const { return SO2(-m_angle); }
  /**
   * @brief Computes the exponential map from the Lie algebra so(2) to the SO(2)
   * group.
   * @param algebra_element The element from the Lie algebra (a single scalar
   * representing the angle).
   * @return The corresponding SO(2) element.
   */
  static SO2<Scalar> computeExp(const TangentVector &algebra_element) {
    return SO2(algebra_element[0]);
  }
  /**
   * @brief Computes the logarithmic map from the SO(2) group to its Lie algebra
   * so(2).
   * @return The corresponding element in the Lie algebra (a single scalar
   * representing the angle).
   */
  TangentVector computeLog() const {
    TangentVector result;
    result[0] = m_angle;
    return result;
  }
  /**
   * @brief Computes the adjoint representation of the current SO(2) element.
   * For SO(2), the adjoint matrix is the identity matrix.
   * @return The adjoint matrix.
   */
  AdjointMatrix computeAdjoint() const { return AdjointMatrix::Identity(); }
  /**
   * @brief Checks if the current SO(2) element is approximately equal to
   * another.
   * @param other The other SO2 element to compare with.
   * @param eps The tolerance for approximation.
   * @return True if the elements are approximately equal, false otherwise.
   */
  bool computeIsApprox(const SO2 &other, const Scalar &eps) const {
    return Types<_Scalar>::isZero(
        Types<_Scalar>::normalizeAngle(m_angle - other.m_angle), eps);
  }
  /**
   * @brief Applies the group action of the current SO(2) element on a 2D point.
   * @param point The 2D point to act upon.
   * @return The transformed 2D point.
   */
  typename Base::ActionVector
  computeAction(const typename Base::ActionVector &point) const {
    typename Base::ActionVector result;
    result(0) = m_complex(0) * point(0) - m_complex(1) * point(1);
    result(1) = m_complex(1) * point(0) + m_complex(0) * point(1);
    return result;
  }

  /**
   * @brief Computes the adjoint representation of a Lie algebra element for
   * SO(2). For SO(2), the adjoint matrix is the zero matrix.
   * @param v The element of the Lie algebra in vector form.
   * @return The adjoint matrix.
   */
  static AdjointMatrix computeAd(const TangentVector &v) {
    return AdjointMatrix::Zero(); // Adjoint for SO(2) is zero matrix
  }

  /**
   * @brief Generates a random SO(2) element.
   * @tparam Generator The type of the random number generator.
   * @param gen The random number generator.
   * @return A random SO(2) element.
   */
  template <typename Generator>
  static SO2<Scalar> computeRandom(Generator &gen) {
    std::uniform_real_distribution<Scalar> dis(-M_PI, M_PI);
    return SO2(dis(gen));
  }

  /**
   * @brief Prints the SO(2) element to an output stream.
   * @param os The output stream.
   * @return The output stream.
   */
  std::ostream &print(std::ostream &os) const {
    os << "SO2(angle=" << m_angle << " rad, " << (m_angle * Scalar(180) / M_PI)
       << " deg)";
    return os;
  }

  /**
   * @brief Gets the type name of the SO2 class.
   * @return A string view of the type name.
   */
  static constexpr std::string_view getTypeName() { return "SO2"; }

  /**
   * @brief Checks if the current SO(2) element is valid.
   * @return True if the angle is finite and the complex representation is
   * consistent, false otherwise.
   */
  bool computeIsValid() const {
    // Check if angle is finite and complex representation is consistent
    return std::isfinite(m_angle) && m_complex.allFinite();
  }

  /**
   * @brief Normalizes the SO(2) element.
   * Re-normalizes the angle and updates the complex representation.
   */
  void computeNormalize() {
    // Re-normalize angle and complex representation
    m_angle = Types<_Scalar>::normalizeAngle(m_angle);
    updateComplex();
  }

  /**
   * @brief Hat operator - maps ℝ¹ to 2×2 skew-symmetric matrix.
   * @param omega Single scalar (rotation rate).
   * @return 2×2 skew-symmetric matrix.
   */
  static Matrix hat(const TangentVector &omega) {
    Matrix result = Matrix::Zero();
    result(0, 1) = -omega[0];
    result(1, 0) = omega[0];
    return result;
  }

  /**
   * @brief Vee operator - inverse of hat, maps 2×2 skew-symmetric matrix to ℝ¹.
   * @param matrix 2×2 skew-symmetric matrix.
   * @return Single scalar.
   */
  static TangentVector vee(const Matrix &matrix) {
    TangentVector result;
    result[0] = matrix(1, 0);
    return result;
  }

  /**
   * @brief Get the rotation matrix representation.
   * @return The 2x2 rotation matrix.
   */
  Matrix matrix() const {
    Matrix R;
    R << m_complex(0), -m_complex(1), m_complex(1), m_complex(0);
    return R;
  }

  /**
   * @brief Get the complex number representation (cos θ, sin θ).
   * @return The complex number representation as an Eigen::Vector2.
   */
  const Complex &complex() const { return m_complex; }

  /**
   * @brief Get a unit vector in the direction of the rotation.
   * @return A unit vector representing the direction.
   */
  Vector direction() const { return m_complex; }

  /**
   * @brief Rotates a vector by 90 degrees counterclockwise.
   * @return The perpendicular vector.
   */
  Vector perpendicular() const {
    Vector result;
    result(0) = -m_complex(1); // -sin θ
    result(1) = m_complex(0);  //  cos θ
    return result;
  }

private:
  Scalar m_angle;    ///< The rotation angle in radians (normalized to [-π, π])
  Complex m_complex; ///< Complex number representation (cos θ, sin θ)

  /**
   * @brief Updates the complex number representation from the current angle.
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
 * @brief Helper function to create an SO2 object from degrees.
 * @tparam Scalar The scalar type.
 * @param degrees The angle in degrees.
 * @return An SO2 object representing the given angle.
 */
template <typename Scalar> SO2<Scalar> SO2FromDegrees(const Scalar &degrees) {
  return SO2<Scalar>(degrees * Types<Scalar>::pi() / Scalar(180));
}

/**
 * @brief Helper function to convert an SO2 angle to degrees.
 * @tparam Scalar The scalar type.
 * @param rotation The SO2 object.
 * @return The angle in degrees.
 */
template <typename Scalar> Scalar SO2ToDegrees(const SO2<Scalar> &rotation) {
  return rotation.angle() * Scalar(180) / Types<Scalar>::pi();
}

} // namespace sofa::component::cosserat::liegroups

// #endif // SOFA_COMPONENT_COSSERAT_LIEGROUPS_SO2_H