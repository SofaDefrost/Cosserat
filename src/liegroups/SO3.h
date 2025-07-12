/******************************************************************************
 *                 SOFA, Simulation Open-Framework Architecture * (c) 2006
 *INRIA, USTL, UJF, CNRS, MGH                     *
 *                                                                             *
 * This program is free software; you can redistribute it and/or modify it     *
 * under the terms of the GNU Lesser General Public License as published by    *
 * the Free Software Foundation; either version 2.1 of the License, or (at     *
 * your option) any later version.                                             *
 *
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
#include "Types.h"
#include <Eigen/Geometry>

namespace sofa::component::cosserat::liegroups {

/**
 * @brief Implementation of SO(3), the Special Orthogonal group in 3D
 *
 * This class implements the group of rotations in 3D space. Elements of SO(3)
 * are represented internally using unit quaternions, which provide an efficient
 * way to compose rotations and compute the exponential map.
 *
 * The Lie algebra so(3) consists of skew-symmetric 3×3 matrices, which can be
 * identified with vectors in ℝ³ (angular velocity vectors).
 *
 * @tparam _Scalar The scalar type (must be a floating-point type)
 */
template <typename _Scalar>
class SO3 : public LieGroupBase<SO3<_Scalar>, _Scalar, 3, 3, 3> {
public:
  using Base = LieGroupBase<SO3<_Scalar>, _Scalar, 3, 3, 3>;
  using Scalar = typename Base::Scalar;
  using Vector = typename Base::Vector;
  using Matrix = typename Base::Matrix;
  using TangentVector3 = typename Base::TangentVector;
  using AdjointMatrix3 = typename Base::AdjointMatrix;

  using Quaternion = Eigen::Quaternion<Scalar>;
  static constexpr int Dim = Base::Dim;

  /**
   * @brief Default constructor creates identity rotation.
   * Initializes the quaternion to identity.
   */
  SO3() : m_quat(Quaternion::Identity()) {}

  /**
   * @brief Construct from angle-axis representation.
   * @param angle Rotation angle in radians.
   * @param axis Unit vector representing rotation axis.
   */
  SO3(const Scalar &angle, const Vector &axis)
      : m_quat(Eigen::AngleAxis<Scalar>(angle, axis.normalized())) {}

  /**
   * @brief Construct from quaternion.
   * The input quaternion will be normalized.
   * @param quat Unit quaternion.
   */
  explicit SO3(const Quaternion &quat) : m_quat(quat.normalized()) {}

  /**
   * @brief Construct from rotation matrix.
   * @param mat 3x3 rotation matrix.
   */
  explicit SO3(const Matrix &mat) : m_quat(mat) {}

  /**
   * @brief Group composition (rotation composition).
   * @param other The other SO3 rotation.
   * @return The composed SO3 rotation.
   */
  SO3 operator*(const SO3 &other) const noexcept {
    return SO3(m_quat * other.m_quat);
  }

  /**
   * @brief Composes this rotation with another.
   * @param other The other SO3 rotation.
   * @return The composed SO3 rotation.
   */
  SO3 compose(const SO3 &other) const noexcept {
    return SO3(m_quat * other.m_quat);
  }

  /**
   * @brief Computes the inverse element (opposite rotation).
   * @return The inverse SO3 rotation.
   */
  SO3 inverse() const { return SO3(m_quat.conjugate()); }

  /**
   * @brief Computes the inverse element (opposite rotation).
   * This is a CRTP-required method.
   * @return The inverse SO3 rotation.
   */
  SO3 computeInverse() const { return SO3(m_quat.conjugate()); }

  /**
   * @brief Exponential map from Lie algebra so(3) to SO(3).
   * @param omega Angular velocity vector in ℝ³.
   * @return The corresponding SO3 element.
   */
  static SO3 exp(const TangentVector &omega) noexcept {
    const Scalar theta = omega.norm();

    if (theta < Types<Scalar>::epsilon()) {
      // For small rotations, use first-order approximation
      return SO3(Quaternion(Scalar(1), omega.x() * Scalar(0.5),
                            omega.y() * Scalar(0.5), omega.z() * Scalar(0.5)));
    }

    // Use Rodrigues\' formula
    const Vector axis = omega / theta;
    const Scalar half_theta = theta * Scalar(0.5);
    const Scalar sin_half_theta = std::sin(half_theta);

    return SO3(Quaternion(std::cos(half_theta), axis.x() * sin_half_theta,
                          axis.y() * sin_half_theta,
                          axis.z() * sin_half_theta));
  }

  /**
   * @brief Computes the exponential map from Lie algebra so(3) to SO(3).
   * This is a CRTP-required method.
   * @param omega Angular velocity vector in ℝ³.
   * @return The corresponding SO3 element.
   */
  static SO3 computeExp(const TangentVector &omega) noexcept {
    const Scalar theta = omega.norm();

    if (theta < Types<Scalar>::epsilon()) {
      // For small rotations, use first-order approximation
      return SO3(Quaternion(Scalar(1), omega.x() * Scalar(0.5),
                            omega.y() * Scalar(0.5), omega.z() * Scalar(0.5)));
    }

    // Use Rodrigues\' formula
    const Vector axis = omega / theta;
    const Scalar half_theta = theta * Scalar(0.5);
    const Scalar sin_half_theta = std::sin(half_theta);

    return SO3(Quaternion(std::cos(half_theta), axis.x() * sin_half_theta,
                          axis.y() * sin_half_theta,
                          axis.z() * sin_half_theta));
  }

  /**
   * @brief Logarithmic map from SO(3) to Lie algebra so(3).
   * @return Angular velocity vector in ℝ³.
   */
  TangentVector log() const {
    // Extract angle-axis representation
    Eigen::AngleAxis<Scalar> aa(m_quat);
    const Scalar theta = aa.angle();

    if (theta < Types<Scalar>::epsilon()) {
      // For small rotations, use first-order approximation
      return Vector(m_quat.x() * Scalar(2), m_quat.y() * Scalar(2),
                    m_quat.z() * Scalar(2));
    }

    return aa.axis() * theta;
  }

  /**
   * @brief Computes the logarithmic map from SO(3) to Lie algebra so(3).
   * This is a CRTP-required method.
   * @return Angular velocity vector in ℝ³.
   */
  TangentVector computeLog() const {
    // Extract angle-axis representation
    Eigen::AngleAxis<Scalar> aa(m_quat);
    const Scalar theta = aa.angle();

    if (theta < Types<Scalar>::epsilon()) {
      // For small rotations, use first-order approximation
      return Vector(m_quat.x() * Scalar(2), m_quat.y() * Scalar(2),
                    m_quat.z() * Scalar(2));
    }

    return aa.axis() * theta;
  }

  /**
   * @brief Adjoint representation of the group element.
   * For SO(3), this is the rotation matrix itself.
   * @return The adjoint matrix representing the action on the Lie algebra.
   */
  AdjointMatrix adjoint() const noexcept { return matrix(); }

  /**
   * @brief Computes the adjoint representation of the group element.
   * This is a CRTP-required method.
   * @return The adjoint matrix representing the action on the Lie algebra.
   */
  AdjointMatrix computeAdjoint() const noexcept { return matrix(); }

  /**
   * @brief Group action on a point (rotates the point).
   * @param point The point to transform.
   * @return The transformed point.
   */
  Vector act(const Vector &point) const noexcept { return m_quat * point; }

  /**
   * @brief Computes the group action on a point (rotates the point).
   * This is a CRTP-required method.
   * @param point The point to transform.
   * @return The transformed point.
   */
  Vector computeAction(const Vector &point) const noexcept {
    return m_quat * point;
  }

  /**
   * @brief Check if approximately equal to another rotation.
   * Handles antipodal representation of the same rotation.
   * @param other Another element of the same Lie group.
   * @param eps Tolerance for comparison (optional).
   * @return true if elements are approximately equal.
   */
  bool isApprox(const SO3 &other,
                const Scalar &eps = Types<Scalar>::epsilon()) const noexcept {
    // Handle antipodal representation of same rotation
    return m_quat.coeffs().isApprox(other.m_quat.coeffs(), eps) ||
           m_quat.coeffs().isApprox(-other.m_quat.coeffs(), eps);
  }

  /**
   * @brief Computes if this element is approximately equal to another.
   * This is a CRTP-required method.
   * @param other Another element of the same Lie group.
   * @param eps Tolerance for comparison (optional).
   * @return true if elements are approximately equal.
   */
  bool
  computeIsApprox(const SO3 &other,
                  const Scalar &eps = Types<Scalar>::epsilon()) const noexcept {
    // Handle antipodal representation of same rotation
    return m_quat.coeffs().isApprox(other.m_quat.coeffs(), eps) ||
           m_quat.coeffs().isApprox(-other.m_quat.coeffs(), eps);
  }

  /**
   * @brief Get the identity element of the group.
   * @return The identity element.
   */
  static SO3 identity() noexcept { return SO3(); }

  /**
   * @brief Computes the identity element of the group.
   * This is a CRTP-required method.
   * @return The identity element.
   */
  static SO3 computeIdentity() noexcept { return SO3(); }

  /**
   * @brief Get the dimension of the Lie algebra (3 for SO(3)).
   * @return The dimension of the Lie algebra.
   */
  static constexpr int algebraDimension() { return 3; }

  /**
   * @brief Get the dimension of the space the group acts on (3 for SO(3)).
   * @return The dimension of the action space.
   */
  static constexpr int actionDimension() { return 3; }

  /**
   * @brief Compute distance between two rotations using the geodesic metric.
   * @param other Another element of the same Lie group.
   * @return A scalar representing the distance.
   */
  Scalar distance(const SO3 &other) const noexcept;

  /**
   * @brief Interpolate between two rotations using SLERP.
   * @param other Target group element.
   * @param t Interpolation parameter between 0 and 1.
   * @return Interpolated group element.
   */
  SO3 interpolate(const SO3 &other, const Scalar &t) const noexcept;

  /**
   * @brief Baker-Campbell-Hausdorff formula for so(3).
   * @param v First tangent vector.
   * @param w Second tangent vector.
   * @param order Order of approximation (1-5, default: 2).
   * @return Tangent vector approximating log(exp(v)*exp(w)).
   */
  static TangentVector BCH(const TangentVector &v, const TangentVector &w,
                           int order = 2);

  /**
   * @brief Differential of the exponential map.
   * @param v Tangent vector.
   * @return Matrix representing the differential of exp at v.
   */
  static AdjointMatrix dexp(const TangentVector &v);

  /**
   * @brief Differential of the logarithm map.
   * @return Matrix representing the differential of log at the current point.
   */
  AdjointMatrix dlog() const;

  /**
   * @brief Adjoint representation of Lie algebra element.
   * @param v Element of the Lie algebra in vector form.
   * @return Matrix representing the adjoint action.
   */
  static AdjointMatrix ad(const TangentVector &v);
  /**
   * @brief Get the rotation matrix representation.
   * @return The 3x3 rotation matrix.
   */
  Matrix matrix() const { return m_quat.toRotationMatrix(); }

  /**
   * @brief Get the quaternion representation.
   * @return A const reference to the unit quaternion representing the rotation.
   */
  const Quaternion &quaternion() const { return m_quat; }

  /**
   * @brief Convert to angle-axis representation.
   * @return The Eigen AngleAxis representation.
   */
  Eigen::AngleAxis<Scalar> angleAxis() const {
    return Eigen::AngleAxis<Scalar>(m_quat);
  }

  /**
   * @brief Convert vector to skew-symmetric matrix (hat operator).
   * @param v Vector in ℝ³.
   * @return 3x3 skew-symmetric matrix.
   */
  static Matrix hat(const TangentVector &v) noexcept {
    Matrix Omega;
    Omega << 0, -v[2], v[1], v[2], 0, -v[0], -v[1], v[0], 0;
    return Omega;
  }

  /**
   * @brief Convert skew-symmetric matrix to vector (vee operator).
   * @param Omega 3x3 skew-symmetric matrix.
   * @return Vector in ℝ³.
   */
  static TangentVector vee(const Matrix &Omega) noexcept {
    return TangentVector(Omega(2, 1), Omega(0, 2), Omega(1, 0));
  }

private:
  Quaternion m_quat; ///< Unit quaternion representing the rotation
};

} // namespace sofa::component::cosserat::liegroups

// #endif // SOFA_COMPONENT_COSSERAT_LIEGROUPS_SO3_H
