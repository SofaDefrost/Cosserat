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

#include <random>
#include <Cosserat/liegroups/LieGroupBase.h>
#include <Cosserat/liegroups/LieGroupBase.inl>
#include <Cosserat/liegroups/SO2.h>
#include <Cosserat/liegroups/Types.h>

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
template <typename _Scalar, int _Dim = 3>
class SE2 : public LieGroupBase<_Scalar, std::integral_constant<int, 3>, _Dim, _Dim> {
public:
  // Type safety checks
  static_assert(std::is_floating_point_v<_Scalar>, "Scalar type must be floating point");
  static_assert(_Dim == 3, "SE(2) requires ambient dimension of 3");

  using Base = LieGroupBase<_Scalar, std::integral_constant<int, 3>, _Dim, _Dim>;
  using Scalar = typename Base::Scalar;
  using Vector = typename Base::Vector;
  using Matrix = typename Base::Matrix;
  using TangentVector = typename Base::TangentVector;
  using AdjointMatrix = typename Base::AdjointMatrix;

  using Vector2 = Eigen::Matrix<Scalar, 2, 1>;
  using Matrix2 = Eigen::Matrix<Scalar, 2, 2>;
  using Matrix3 = Eigen::Matrix<Scalar, 3, 3>;
  
  static constexpr int Dim = Base::Dim;
  static constexpr int DOF = 3; // Degrees of freedom
  static constexpr int ActionDim = 2; // Dimension of space acted upon

  // ========== Constructors ==========

  /**
   * @brief Default constructor creates identity transformation
   */
  SE2() : m_rotation(), m_translation(Vector2::Zero()) {}

  /**
   * @brief Construct from rotation and translation
   */
  SE2(const SO2<Scalar>& rotation, const Vector2& translation)
      : m_rotation(rotation), m_translation(translation) {}

  /**
   * @brief Construct from angle and translation
   */
  SE2(const Scalar& angle, const Vector2& translation)
      : m_rotation(angle), m_translation(translation) {}

  /**
   * @brief Construct from homogeneous transformation matrix
   * @param matrix 3x3 homogeneous transformation matrix
   */
  explicit SE2(const Matrix3& matrix) {
    // Validate matrix structure
    if (!isValidTransformationMatrix(matrix)) {
      throw std::invalid_argument("Invalid SE(2) transformation matrix");
    }
    
    m_rotation = SO2<Scalar>(matrix.template block<2, 2>(0, 0));
    m_translation = matrix.template block<2, 1>(0, 2);
  }

  /**
   * @brief Construct from Lie algebra vector
   * @param tangent_vector 3D vector [vx, vy, ω]
   */
  explicit SE2(const TangentVector& tangent_vector) {
    *this = exp(tangent_vector);
  }

  /**
   * @brief Copy constructor
   */
  SE2(const SE2& other) = default;

  /**
   * @brief Move constructor
   */
  SE2(SE2&& other) noexcept = default;

  /**
   * @brief Copy assignment
   */
  SE2& operator=(const SE2& other) = default;

  /**
   * @brief Move assignment
   */
  SE2& operator=(SE2&& other) noexcept = default;

  // ========== Group Operations ==========

  /**
   * @brief Group composition (rigid transformation composition)
   * Implements: this * other = (R₁R₂, R₁t₂ + t₁)
   */
  SE2 operator*(const SE2& other) const {
    return SE2(m_rotation * other.m_rotation,
               m_translation + m_rotation.act(other.m_translation));
  }

  /**
   * @brief In-place composition
   */
  SE2& operator*=(const SE2& other) {
    m_translation += m_rotation.act(other.m_translation);
    m_rotation *= other.m_rotation;
    return *this;
  }

  /**
   * @brief Inverse element (inverse transformation)
   * Implements: (R,t)⁻¹ = (R^T, -R^T t)
   */
  SE2 inverse() const override {
    SO2<Scalar> inv_rot = m_rotation.inverse();
    return SE2(inv_rot, -(inv_rot.act(m_translation)));
  }

  // ========== Lie Algebra Operations ==========

  /**
   * @brief Exponential map from Lie algebra se(2) to SE(2)
   * @param algebra_element Vector in ℝ³ representing (vx, vy, ω)
   * 
   * For ξ = [ρ; φ] where ρ ∈ ℝ² and φ ∈ ℝ:
   * exp(ξ) = [exp_SO2(φ), V(φ)ρ; 0, 1]
   * where V(φ) = (sin φ / φ)I + ((1-cos φ)/φ²)K with K = [0 -1; 1 0]
   */
  SE2 exp(const TangentVector& algebra_element) const override {
    const Scalar& theta = algebra_element[2];
    const Vector2 rho(algebra_element[0], algebra_element[1]);

    // Handle small angle case for numerical stability
    if (std::abs(theta) < Types<Scalar>::epsilon()) {
      return SE2(SO2<Scalar>(theta), rho);
    }

    // Compute V matrix for translation component
    const Scalar sin_theta = std::sin(theta);
    const Scalar cos_theta = std::cos(theta);
    const Scalar theta_inv = Scalar(1) / theta;
    const Scalar sin_over_theta = sin_theta * theta_inv;
    const Scalar one_minus_cos_over_theta = (Scalar(1) - cos_theta) * theta_inv;

    Matrix2 V;
    V << sin_over_theta, -one_minus_cos_over_theta,
         one_minus_cos_over_theta, sin_over_theta;

    return SE2(SO2<Scalar>(theta), V * rho);
  }

  /**
   * @brief Logarithmic map from SE(2) to Lie algebra se(2)
   * @return Vector in ℝ³ representing (vx, vy, ω)
   */
  TangentVector log() const override {
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
      const Scalar coeff = Scalar(1) - theta_sq / Scalar(12); // First correction term
      Matrix2 V_inv = coeff * Matrix2::Identity();
      V_inv(0, 1) = -theta / Scalar(2);
      V_inv(1, 0) = theta / Scalar(2);
      
      Vector2 rho = V_inv * m_translation;
      result << rho, theta;
      return result;
    }

    // General case
    const Scalar sin_theta = std::sin(theta);
    const Scalar cos_theta = std::cos(theta);
    const Scalar half_theta = theta * Scalar(0.5);
    const Scalar cot_half = cos_theta / sin_theta; // cot(θ/2) = cos(θ)/sin(θ) for θ ≠ 0

    Matrix2 V_inv;
    V_inv << half_theta * cot_half, -half_theta,
             half_theta, half_theta * cot_half;

    Vector2 rho = V_inv * m_translation;
    result << rho, theta;
    return result;
  }

  /**
   * @brief Adjoint representation Ad_g: se(2) → se(2)
   * For g = (R,t), Ad_g = [R, [t]×; 0, 1] where [t]× is the skew-symmetric matrix
   */
  AdjointMatrix adjoint() const override {
    AdjointMatrix Ad = AdjointMatrix::Zero();
    
    // Rotation block (top-left 2x2)
    Ad.template block<2, 2>(0, 0) = m_rotation.matrix();
    
    // Translation coupling (top-right 2x1)
    Ad(0, 2) = -m_translation.y();
    Ad(1, 2) = m_translation.x();
    
    // Bottom row for angular component
    Ad(2, 2) = Scalar(1);
    
    return Ad;
  }

  // ========== Group Actions ==========

  /**
   * @brief Group action on a 2D point
   * Implements: g · p = Rp + t
   */
  Vector2 act(const Vector2& point) const {
    return m_rotation.act(point) + m_translation;
  }

  /**
   * @brief Group action on a 3D vector (treats as homogeneous coordinates)
   * For [x, y, z]^T, applies SE(2) to [x, y] and preserves z
   */
  Vector act(const Vector& point) const override {
    Vector2 transformed = act(point.template head<2>());
    Vector result;
    result << transformed, point[2];
    return result;
  }

  /**
   * @brief Batch action on multiple points
   */
  template<int N>
  Eigen::Matrix<Scalar, 2, N> act(const Eigen::Matrix<Scalar, 2, N>& points) const {
    return m_rotation.matrix() * points + m_translation.replicate(1, N);
  }

  // ========== Utility Functions ==========

  /**
   * @brief Check if approximately equal to another element
   */
  bool isApprox(const SE2& other, const Scalar& eps = Types<Scalar>::epsilon()) const {
    return m_rotation.isApprox(other.m_rotation, eps) &&
           m_translation.isApprox(other.m_translation, eps);
  }

  /**
   * @brief Equality operator
   */
  bool operator==(const SE2& other) const {
    return isApprox(other);
  }

  /**
   * @brief Inequality operator
   */
  bool operator!=(const SE2& other) const {
    return !(*this == other);
  }

  /**
   * @brief Linear interpolation between two SE(2) elements
   * @param other Target transformation
   * @param t Interpolation parameter [0,1]
   */
  SE2 interpolate(const SE2& other, const Scalar& t) const {
    // Interpolate in Lie algebra for geodesic interpolation
    TangentVector delta = (inverse() * other).log();
    return *this * exp(t * delta);
  }

  /**
   * @brief Generate random SE(2) element
   * @param gen Random number generator
   * @param translation_scale Scale for translation component
   */
  template<typename Generator>
  static SE2 random(Generator& gen, const Scalar& translation_scale = Scalar(1)) {
    std::uniform_real_distribution<Scalar> angle_dist(-M_PI, M_PI);
    std::normal_distribution<Scalar> trans_dist(0, translation_scale);
    
    Scalar angle = angle_dist(gen);
    Vector2 translation(trans_dist(gen), trans_dist(gen));
    
    return SE2(angle, translation);
  }

  // ========== Static Members ==========

  /**
   * @brief Get the identity element
   */
  static const SE2& identity() {
    static const SE2 id;
    return id;
  }

  /**
   * @brief Create SE(2) element from translation only
   */
  static SE2 fromTranslation(const Vector2& translation) {
    return SE2(SO2<Scalar>(), translation);
  }

  /**
   * @brief Create SE(2) element from rotation only
   */
  static SE2 fromRotation(const Scalar& angle) {
    return SE2(SO2<Scalar>(angle), Vector2::Zero());
  }

  static SE2 fromRotation(const SO2<Scalar>& rotation) {
    return SE2(rotation, Vector2::Zero());
  }

  // ========== Accessors ==========

  /**
   * @brief Get the dimension of the Lie algebra (3 for SE(2))
   */
  int algebraDimension() const override { return DOF; }

  /**
   * @brief Get the dimension of the space the group acts on (2 for SE(2))
   */
  int actionDimension() const override { return ActionDim; }

  /**
   * @brief Access the rotation component
   */
  const SO2<Scalar>& rotation() const { return m_rotation; }
  SO2<Scalar>& rotation() { return m_rotation; }

  /**
   * @brief Access the translation component
   */
  const Vector2& translation() const { return m_translation; }
  Vector2& translation() { return m_translation; }

  /**
   * @brief Get the rotation angle
   */
  Scalar angle() const { return m_rotation.angle(); }

  /**
   * @brief Get the homogeneous transformation matrix
   */
  Matrix3 matrix() const {
    Matrix3 T = Matrix3::Identity();
    T.template block<2, 2>(0, 0) = m_rotation.matrix();
    T.template block<2, 1>(0, 2) = m_translation;
    return T;
  }

  /**
   * @brief Get the inverse transformation matrix
   */
  Matrix3 inverseMatrix() const {
    return inverse().matrix();
  }

  // ========== Stream Operators ==========

  /**
   * @brief Output stream operator
   */
  friend std::ostream& operator<<(std::ostream& os, const SE2& se2) {
    os << "SE2(angle=" << se2.angle() 
       << ", translation=(" << se2.m_translation.transpose() << "))";
    return os;
  }

private:
  /**
   * @brief Validate if a 3x3 matrix is a valid SE(2) transformation matrix
   */
  static bool isValidTransformationMatrix(const Matrix3& matrix, 
                                        const Scalar& eps = Types<Scalar>::epsilon()) {
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

  SO2<Scalar> m_rotation;    ///< Rotation component
  Vector2 m_translation;     ///< Translation component
};

// ========== Type Aliases ==========
using SE2f = SE2<float>;
using SE2d = SE2<double>;

} // namespace sofa::component::cosserat::liegroups