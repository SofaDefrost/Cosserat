// This file defines the SE2 (Special Euclidean group in 2D) class, representing
// rigid body transformations in 2D space.

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

#pragma once

#include "LieGroupBase.h"
#include "LieGroupBase.inl"
#include "SO2.h"
#include "Types.h"
#include <cmath>
#include <random>

namespace sofa::component::cosserat::liegroups {

/**
 * @brief Implementation of SE(2), the Special Euclidean group in 2D
 *
 * This class implements the group of rigid body transformations in 2D space.
 * Elements of SE(2) are represented as a combination of:
 * - An SO(2) rotation
 * - A 2D translation vector
 *
 * The Lie algebra se(2) consists of vectors in ℝ³, where:
 * - The first two components represent the translation velocity
 * - The last component represents the angular velocity
 *
 * Mathematical representation:
 * SE(2) = {(R,t) | R ∈ SO(2), t ∈ ℝ²}
 * Matrix form: [R t; 0 1] where R is 2x2 rotation matrix, t is 2x1 translation
 *
 * @tparam _Scalar The scalar type (must be a floating-point type)
 * @tparam _Dim The dimension of the ambient space (fixed at 3 for SE(2))
 */
template <typename _Scalar>
class SE2 : public LieGroupBase<SE2<_Scalar>, _Scalar, 3, 3, 2> {

  using Base = LieGroupBase<SE2<_Scalar>, _Scalar, 3, 3, 2>;
  using Scalar = typename Base::Scalar;
  using Vector = typename Base::Vector;
  using Matrix = typename Base::Matrix;
  using TangentVector = typename Base::TangentVector;
  using AdjointMatrix = typename Base::AdjointMatrix;

  using Vector2 = Eigen::Matrix<Scalar, 2, 1>;
  using Matrix2 = Eigen::Matrix<Scalar, 2, 2>;
  using Matrix3 = Eigen::Matrix<Scalar, 3, 3>;

  static constexpr int Dim = Base::Dim;
  static constexpr int DOF = 3;       // Degrees of freedom
  static constexpr int ActionDim = 2; // Dimension of space acted upon

  // ========== Constructors ==========

  /**
   * @brief Default constructor creates identity transformation.
   * Initializes rotation to identity and translation to zero vector.
   */
  SE2() : m_rotation(), m_translation(Vector2::Zero()) {}

  /**
   * @brief Construct from rotation and translation.
   * @param rotation The SO2 rotation component.
   * @param translation The 2D translation vector.
   */
  SE2(const SO2<Scalar> &rotation, const Vector2 &translation)
      : m_rotation(rotation), m_translation(translation) {}

  /**
   * @brief Construct from angle and translation.
   * @param angle The rotation angle in radians.
   * @param translation The 2D translation vector.
   */
  SE2(const Scalar &angle, const Vector2 &translation)
      : m_rotation(angle), m_translation(translation) {}

  /**
   * @brief Construct from homogeneous transformation matrix.
   * @param matrix 3x3 homogeneous transformation matrix.
   * @throws std::invalid_argument if the matrix is not a valid SE(2)
   * transformation matrix.
   */
  explicit SE2(const Matrix3 &matrix) {
    // Validate matrix structure
    if (!isValidTransformationMatrix(matrix)) {
      throw std::invalid_argument("Invalid SE(2) transformation matrix");
    }

    m_rotation = SO2<Scalar>(matrix.template block<2, 2>(0, 0));
    m_translation = matrix.template block<2, 1>(0, 2);
  }

  /**
   * @brief Construct from Lie algebra vector.
   * @param tangent_vector 3D vector [vx, vy, ω] representing the Lie algebra
   * element.
   */
  explicit SE2(const TangentVector &tangent_vector) {
    *this = exp(tangent_vector);
  }

  /**
   * @brief Copy constructor.
   * @param other The SE2 object to copy from.
   */
  SE2(const SE2 &other) = default;

  /**
   * @brief Move constructor.
   * @param other The SE2 object to move from.
   */
  SE2(SE2 &&other) noexcept = default;

  /**
   * @brief Copy assignment operator.
   * @param other The SE2 object to assign from.
   * @return A reference to the assigned object.
   */
  SE2 &operator=(const SE2 &other) = default;

  /**
   * @brief Move assignment operator.
   * @param other The SE2 object to move from.
   * @return A reference to the assigned object.
   */
  SE2 &operator=(SE2 &&other) noexcept = default;

  // ========== Group Operations ==========

  /**
   * @brief Group composition operator.
   * Computes the composition of this transformation with another.
   * @param other The other SE2 transformation.
   * @return The composed SE2 transformation.
   */
  SE2 operator*(const SE2 &other) const {
    return SE2(m_rotation * other.m_rotation,
               m_rotation.act(other.m_translation) + m_translation);
  }

  /**
   * @brief In-place group composition operator.
   * Composes this transformation with another in-place.
   * @param other The other SE2 transformation.
   * @return A reference to this object after composition.
   */
  SE2 &operator*=(const SE2 &other) {
    m_translation = m_rotation.act(other.m_translation) + m_translation;
    m_rotation *= other.m_rotation;
    return *this;
  }

  // ========== Lie Algebra Operations ==========

  /**
   * @brief Exponential map from Lie algebra se(2) to SE(2).
   * @param algebra_element Vector in ℝ³ representing (vx, vy, ω).
   *
   * For ξ = [ρ; φ] where ρ ∈ ℝ² and φ ∈ ℝ:
   * exp(ξ) = [exp_SO2(φ), V(φ)ρ; 0, 1]
   * where V(φ) = (sin φ / φ)I + ((1-cos φ)/φ²)K with K = [0 -1; 1 0]
   */
  static SE2 computeExp(const TangentVector &algebra_element) {
    const Scalar &theta = algebra_element[2];
    const Vector2 rho(algebra_element[0], algebra_element[1]);

    // Handle small angle case for numerical stability
    if (std::abs(theta) < Types<Scalar>::epsilon()) {
      return SE2(SO2<Scalar>(theta), rho);
    }

    // Compute V matrix for translation component
    const Scalar sin_theta = std::sin(theta);
    const Scalar cos_theta = std::cos(theta);
    const Scalar theta_inv = Scalar(1.0) / theta;
    const Scalar sin_over_theta = sin_theta * theta_inv;
    const Scalar one_minus_cos_over_theta =
        (Scalar(1.0) - cos_theta) * theta_inv;

    Matrix2 V;
    V << sin_over_theta, -one_minus_cos_over_theta, one_minus_cos_over_theta,
        sin_over_theta;

    return SE2(SO2<Scalar>(theta), V * rho);
  }

  /**
   * @brief Logarithmic map from SE(2) to Lie algebra se(2).
   * @return Vector in ℝ³ representing (vx, vy, ω).
   */
  TangentVector computeLog() const {
    const Scalar theta = m_rotation.angle();
    TangentVector result;

    // Handle small angle case
    if (std::abs(theta) < Types<Scalar>::epsilon()) {
      result << m_translation, theta;
      return result;
    }

    // Check for singularity (θ = ±π)
    const Scalar abs_theta = std::abs(theta);
    if (abs_theta > M_PI - Types<Scalar>::epsilon()) {
      // Near ±π, use series expansion for numerical stability
      const Scalar theta_sq = theta * theta;
      const Scalar coeff =
          Scalar(1.0) - theta_sq / Scalar(12.0); // First correction term
      Matrix2 V_inv = coeff * Matrix2::Identity();
      V_inv(0, 1) = -theta / Scalar(2.0);
      V_inv(1, 0) = theta / Scalar(2.0);

      Vector2 rho = V_inv * m_translation;
      result << rho, theta;
      return result;
    }

    // General case
    const Scalar sin_theta = std::sin(theta);
    const Scalar cos_theta = std::cos(theta);
    const Scalar half_theta = theta * Scalar(0.5);
    const Scalar cot_half =
        cos_theta / sin_theta; // cot(θ/2) = cos(θ)/sin(θ) for θ ≠ 0

    Matrix2 V_inv;
    V_inv << half_theta * cot_half, -half_theta, half_theta,
        half_theta * cot_half;

    Vector2 rho = V_inv * m_translation;
    result << rho, theta;
    return result;
  }

  /**
   * @brief Computes the adjoint representation of the current SE(2) element.
   * @return The adjoint matrix.
   */
  AdjointMatrix computeAdjoint() const {
    AdjointMatrix Ad = AdjointMatrix::Zero();

    // Rotation block (top-left 2x2)
    Ad.template block<2, 2>(0, 0) = m_rotation.matrix();

    // Translation coupling (top-right 2x1)
    Ad(0, 2) = -m_translation.y();
    Ad(1, 2) = m_translation.x();

    // Bottom row for angular component
    Ad(2, 2) = Scalar(1.0);

    return Ad;
  }

  // ========== Group Actions ==========

  /**
   * @brief Applies the group action of the current SE(2) element on a 2D point.
   * @param point The 2D point to act upon.
   * @return The transformed 2D point.
   */
  typename Base::ActionVector
  computeAction(const typename Base::ActionVector &point) const {
    return m_rotation.computeAction(point.template head<2>()) + m_translation;
  }

  // ========== Utility Functions ==========

  /**
   * @brief Checks if the current SE2 element is approximately equal to another.
   * @param other The other SE2 element to compare with.
   * @param eps The tolerance for approximation.
   * @return True if the elements are approximately equal, false otherwise.
   */
  bool computeIsApprox(const SE2 &other, const Scalar &eps) const {
    return m_rotation.computeIsApprox(other.m_rotation, eps) &&
           m_translation.isApprox(other.m_translation, eps);
  }

  /**
   * @brief Computes the inverse of the current SE2 element.
   * @return The inverse SE2 element.
   */
  SE2 computeInverse() const {
    SO2<Scalar> inv_rot = m_rotation.computeInverse();
    return SE2(inv_rot, -(inv_rot.computeAction(m_translation)));
  }

  /**
   * @brief Computes the identity element of the SE(2) group.
   * @return The identity SE(2) element.
   */
  static SE2 computeIdentity() noexcept { return SE2(); }

  /**
   * @brief Generates a random SE(2) element.
   * @tparam Generator The type of the random number generator.
   * @param gen The random number generator.
   * @return A random SE(2) element.
   */
  template <typename Generator>
  static SE2<Scalar> computeRandom(Generator &gen) {
    std::uniform_real_distribution<Scalar> angle_dist(-M_PI, M_PI);
    std::normal_distribution<Scalar> trans_dist(0, Scalar(1.0));

    Scalar angle = angle_dist(gen);
    Vector2 translation(trans_dist(gen), trans_dist(gen));

    return SE2(angle, translation);
  }

  /**
   * @brief Prints the SE(2) element to an output stream.
   * @param os The output stream.
   * @return The output stream.
   */
  std::ostream &print(std::ostream &os) const {
    os << "SE2(angle=" << m_rotation.angle() << ", translation=("
       << m_translation.transpose() << "))";
    return os;
  }

  /**
   * @brief Gets the type name of the SE2 class.
   * @return A string view of the type name.
   */
  static constexpr std::string_view getTypeName() { return "SE2"; }

  /**
   * @brief Checks if the current SE2 element is valid.
   * @return True if both rotation and translation components are valid, false
   * otherwise.
   */
  bool computeIsValid() const {
    return m_rotation.computeIsValid() && m_translation.allFinite();
  }

  /**
   * @brief Normalizes the SE2 element.
   * Normalizes the rotation component.
   */
  void computeNormalize() { m_rotation.computeNormalize(); }

  // ========== Static Members ==========

  /**
   * @brief Hat operator - maps a 3D Lie algebra vector to a 3x3 matrix
   * representation.
   * @param v The 3D Lie algebra vector [vx, vy, ω].
   * @return The 3x3 matrix representation.
   */
  static Matrix hat(const TangentVector &v) {
    Matrix R = Matrix::Zero();
    R(0, 2) = v(0);
    R(1, 2) = v(1);
    R(0, 1) = -v(2);
    R(1, 0) = v(2);
    return R;
  }

  /**
   * @brief Vee operator - inverse of hat, maps a 3x3 matrix representation to a
   * 3D Lie algebra vector.
   * @param X The 3x3 matrix representation.
   * @return The 3D Lie algebra vector [vx, vy, ω].
   */
  static TangentVector vee(const Matrix &X) {
    TangentVector v;
    v(0) = X(0, 2);
    v(1) = X(1, 2);
    v(2) = X(1, 0);
    return v;
  }

  /**
   * @brief Computes the adjoint representation of a Lie algebra element for
   * SE(2).
   * @param v The element of the Lie algebra in vector form.
   * @return The adjoint matrix.
   */
  static AdjointMatrix computeAd(const TangentVector &v) {
    AdjointMatrix Ad = AdjointMatrix::Zero();
    Ad(0, 0) = 1.0;
    Ad(1, 1) = 1.0;
    Ad(2, 2) = 1.0;
    Ad(0, 2) = -v(1);
    Ad(1, 2) = v(0);
    return Ad;
  }

  /**
   * @brief Create SE(2) element from translation only.
   * @param translation The 2D translation vector.
   * @return An SE2 object with identity rotation and the given translation.
   */
  static SE2 fromTranslation(const Vector2 &translation) {
    return SE2(SO2<Scalar>(), translation);
  }

  /**
   * @brief Create SE(2) element from rotation only.
   * @param angle The rotation angle in radians.
   * @return An SE2 object with the given rotation and zero translation.
   */
  static SE2 fromRotation(const Scalar &angle) {
    return SE2(SO2<Scalar>(angle), Vector2::Zero());
  }

  /**
   * @brief Create SE(2) element from rotation only.
   * @param rotation The SO2 rotation component.
   * @return An SE2 object with the given rotation and zero translation.
   */
  static SE2 fromRotation(const SO2<Scalar> &rotation) {
    return SE2(rotation, Vector2::Zero());
  }

  // ========== Accessors ==========

  /**
   * @brief Access the rotation component (const version).
   * @return A const reference to the SO2 rotation component.
   */
  const SO2<Scalar> &rotation() const { return m_rotation; }
  /**
   * @brief Access the rotation component (mutable version).
   * @return A mutable reference to the SO2 rotation component.
   */
  SO2<Scalar> &rotation() { return m_rotation; }

  /**
   * @brief Access the translation component (const version).
   * @return A const reference to the 2D translation vector.
   */
  const Vector2 &translation() const { return m_translation; }
  /**
   * @brief Access the translation component (mutable version).
   * @return A mutable reference to the 2D translation vector.
   */
  Vector2 &translation() { return m_translation; }

  /**
   * @brief Get the rotation angle.
   * @return The rotation angle in radians.
   */
  Scalar angle() const { return m_rotation.angle(); }

  /**
   * @brief Get the homogeneous transformation matrix.
   * @return The 3x3 homogeneous transformation matrix.
   */
  Matrix3 matrix() const {
    Matrix3 T = Matrix3::Identity();
    T.template block<2, 2>(0, 0) = m_rotation.matrix();
    T.template block<2, 1>(0, 2) = m_translation;
    return T;
  }

  /**
   * @brief Get the inverse transformation matrix.
   * @return The 3x3 inverse homogeneous transformation matrix.
   */
  Matrix3 inverseMatrix() const { return computeInverse().matrix(); };

  // ========== Stream Operators ==========

  /**
   * @brief Overload of the stream insertion operator for SE2.
   * @param os The output stream.
   * @param se2 The SE2 object to print.
   * @return The output stream.
   */
  friend std::ostream &operator<<(std::ostream &os, const SE2 &se2) {
    os << "SE2(angle=" << se2.m_rotation.angle() << ", translation=("
       << se2.m_translation.transpose() << "))";
    return os;
  }

private:
  /**
   * @brief Validates if a 3x3 matrix is a valid SE(2) transformation matrix.
   * Checks if the bottom row is [0, 0, 1] and if the rotation part is valid.
   * @param matrix The 3x3 matrix to validate.
   * @param eps The tolerance for floating-point comparisons. Defaults to
   * `Types<Scalar>::epsilon()`.
   * @return True if the matrix is a valid SE(2) transformation matrix, false
   * otherwise.
   */
  static bool
  isValidTransformationMatrix(const Matrix3 &matrix,
                              const Scalar &eps = Types<Scalar>::epsilon()) {
    // Check bottom row is [0, 0, 1]
    if (!matrix.template block<1, 2>(2, 0).isZero(eps) ||
        std::abs(matrix(2, 2) - Scalar(1)) > eps) {
      return false;
    }

    // Check if rotation part is valid (should be done by SO2 constructor)
    Matrix2 R = matrix.template block<2, 2>(0, 0);
    return std::abs(R.determinant() - Scalar(1)) < eps &&
           (R * R.transpose()).isApprox(Matrix2::Identity(), eps);
  }

  SO2<Scalar> m_rotation; ///< Rotation component
  Vector2 m_translation;  ///< Translation component
};

// ========== Type Aliases ==========
using SE2f = SE2<float>;
using SE2d = SE2<double>;

} // namespace sofa::component::cosserat::liegroups