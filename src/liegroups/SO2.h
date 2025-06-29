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

#include <cmath>
#include "LieGroupBase.h"   // Then the base class interface
#include "LieGroupBase.inl" // Then the base class interface
#include "Types.h"          // Then our type system
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
   */
  explicit SO2(const Scalar &angle)
      : m_angle(Types<_Scalar>::normalizeAngle(angle)) {
    updateComplex();
  }

  /**
   * @brief Copy constructor
   */
  SO2(const SO2 &other) : m_angle(other.m_angle), m_complex(other.m_complex) {}

  /**
   * @brief Get the rotation angle in radians
   */
  Scalar angle() const { return m_angle; }

  /**
   * @brief Set the rotation angle in radians
   */
  void setAngle(const Scalar &angle) {
    m_angle = Types<_Scalar>::normalizeAngle(angle);
    updateComplex();
  }

  // Required CRTP methods:
  static SO2<Scalar> computeIdentity() noexcept { return SO2(Scalar(0)); }
  SO2<Scalar> computeInverse() const { return SO2(-m_angle); }
  static SO2<Scalar> computeExp(const TangentVector &algebra_element) {
    return SO2(algebra_element[0]);
  }
  TangentVector computeLog() const {
    TangentVector result;
    result[0] = m_angle;
    return result;
  }
  AdjointMatrix computeAdjoint() const { return AdjointMatrix::Identity(); }
  bool computeIsApprox(const SO2 &other, const Scalar &eps) const {
    return Types<_Scalar>::isZero(
        Types<_Scalar>::normalizeAngle(m_angle - other.m_angle), eps);
  }
  typename Base::ActionVector
  computeAction(const typename Base::ActionVector &point) const {
    typename Base::ActionVector result;
    result(0) = m_complex(0) * point(0) - m_complex(1) * point(1);
    result(1) = m_complex(1) * point(0) + m_complex(0) * point(1);
    return result;
  }

  static AdjointMatrix computeAd(const TangentVector &v) {
    return AdjointMatrix::Zero(); // Adjoint for SO(2) is zero matrix
  }

  template <typename Generator>
  static SO2<Scalar> computeRandom(Generator &gen) {
    std::uniform_real_distribution<Scalar> dis(-M_PI,
                                               M_PI);
    return SO2(dis(gen));
  }

  std::ostream &print(std::ostream &os) const {
    os << "SO2(angle=" << m_angle << " rad, "
       << (m_angle * Scalar(180) / M_PI) << " deg)";
    return os;
  }

  static constexpr std::string_view getTypeName() {
    return "SO2";
  }

  bool computeIsValid() const {
    // Check if angle is finite and complex representation is consistent
    return std::isfinite(m_angle) && m_complex.allFinite();
  }

  void computeNormalize() {
    // Re-normalize angle and complex representation
    m_angle = Types<_Scalar>::normalizeAngle(m_angle);
    updateComplex();
  }

  /**
   * @brief Hat operator - maps ℝ¹ to 2×2 skew-symmetric matrix
   * @param omega Single scalar (rotation rate)
   * @return 2×2 skew-symmetric matrix
   */
  static Matrix hat(const TangentVector &omega) {
    Matrix result = Matrix::Zero();
    result(0, 1) = -omega[0];
    result(1, 0) = omega[0];
    return result;
  }

  /**
   * @brief Vee operator - inverse of hat, maps 2×2 skew-symmetric matrix to ℝ¹
   * @param matrix 2×2 skew-symmetric matrix
   * @return Single scalar
   */
  static TangentVector vee(const Matrix &matrix) {
    TangentVector result;
    result[0] = matrix(1, 0);
    return result;
  }

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