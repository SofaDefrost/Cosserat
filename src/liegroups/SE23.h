// This file defines the SE23 (extended Special Euclidean group in 3D) class,
// representing rigid body transformations with linear velocity in 3D space.

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

// #ifndef SOFA_COMPONENT_COSSERAT_LIEGROUPS_SE22_H
// #define SOFA_COMPONENT_COSSERAT_LIEGROUPS_SE23_H
#pragma once

#include "LieGroupBase.h"   // Then the base class interface
#include "LieGroupBase.inl" // Then the base class interface
#include "SE3.h"            // Then the base class interface
// #include <eigen3/Eigen/Geometry.h>

namespace sofa::component::cosserat::liegroups {

/**
 * @brief Implementation of SE_2(3), the extended Special Euclidean group in 3D
 *
 * This class implements the group of rigid body transformations with linear
 * velocity in 3D space. Elements of SE_2(3) are represented as a combination
 * of:
 * - An SE(3) transformation (rotation and position)
 * - A 3D linear velocity vector
 *
 * The Lie algebra se_2(3) consists of vectors in ℝ⁹, where:
 * - First three components represent linear velocity
 * - Middle three components represent angular velocity
 * - Last three components represent linear acceleration
 *
 * @tparam _Scalar The scalar type (must be a floating-point type)
 */
template <typename _Scalar, int _Dim = 9>
class SE23 : public LieGroupBase<SE23<_Scalar, _Dim>, _Scalar, 9, 9, 6>
//,public LieGroupOperations<SE23<_Scalar>>
{
public:
  using Base = LieGroupBase<_Scalar, std::integral_constant<int, _Dim>, 3, 3>;
  using Scalar = typename Base::Scalar;
  using Vector = typename Base::Vector;
  using Matrix = typename Base::Matrix;
  using TangentVector = typename Base::TangentVector;
  using AdjointMatrix = typename Base::AdjointMatrix;

  using Vector3 = Eigen::Matrix<Scalar, 3, 1>;
  using Matrix3 = Eigen::Matrix<Scalar, 3, 3>;
  using Matrix4 = Eigen::Matrix<Scalar, 4, 4>;
  static constexpr int Dim = Base::Dim;

  /**
   * @brief Default constructor creates identity element.
   * Initializes pose to identity and velocity to zero vector.
   */
  SE23() : m_pose(), m_velocity(Vector3::Zero()) {}

  /**
   * @brief Construct from pose and velocity.
   * @param pose The SE3 pose component.
   * @param velocity The 3D linear velocity vector.
   */
  SE23(const SE3<Scalar> &pose, const Vector3 &velocity)
      : m_pose(pose), m_velocity(velocity) {}

  /**
   * @brief Construct from rotation, position, and velocity.
   * @param rotation The SO3 rotation component.
   * @param position The 3D position vector.
   * @param velocity The 3D linear velocity vector.
   */
  SE23(const SO3<Scalar> &rotation, const Vector3 &position,
       const Vector3 &velocity)
      : m_pose(rotation, position), m_velocity(velocity) {}

  /**
   * @brief Access the pose component (const version).
   * @return A const reference to the SE3 pose component.
   */
  const SE3<Scalar> &pose() const { return m_pose; }
  /**
   * @brief Access the pose component (mutable version).
   * @return A mutable reference to the SE3 pose component.
   */
  SE3<Scalar> &pose() { return m_pose; }

  /**
   * @brief Access the velocity component (const version).
   * @return A const reference to the 3D linear velocity vector.
   */
  const Vector3 &velocity() const { return m_velocity; }
  /**
   * @brief Access the velocity component (mutable version).
   * @return A mutable reference to the 3D linear velocity vector.
   */
  Vector3 &velocity() { return m_velocity; }

  /**
   * @brief Get the extended homogeneous transformation matrix.
   * This matrix represents the SE23 element in a higher-dimensional space.
   * @return The 5x5 extended homogeneous transformation matrix.
   */
  Eigen::Matrix<Scalar, 5, 5> extendedMatrix() const {
    Eigen::Matrix<Scalar, 5, 5> T = Eigen::Matrix<Scalar, 5, 5>::Identity();
    T.template block<4, 4>(0, 0) = m_pose.matrix();
    T.template block<3, 1>(0, 4) = m_velocity;
    return T;
  }

  // Required CRTP methods:
  /**
   * @brief Computes the identity element of the SE23 group.
   * @return The identity SE23 element.
   */
  static SE23<Scalar> computeIdentity() noexcept { return SE23(); }
  /**
   * @brief Computes the inverse of the current SE23 element.
   * @return The inverse SE23 element.
   */
  SE23<Scalar> computeInverse() const {
    SE3<Scalar> inv_pose = m_pose.computeInverse();
    return SE23(inv_pose, -(inv_pose.rotation().computeAction(m_velocity)));
  }
  /**
   * @brief Computes the exponential map from the Lie algebra se_2(3) to the
   * SE23 group.
   * @param algebra_element The element from the Lie algebra (a 9D vector
   * representing linear velocity, angular velocity, and linear acceleration).
   * @return The corresponding SE23 element.
   */
  static SE23<Scalar> computeExp(const TangentVector &algebra_element) {
    Vector3 v = algebra_element.template segment<3>(0); // Linear velocity
    Vector3 w = algebra_element.template segment<3>(3); // Angular velocity
    Vector3 a = algebra_element.template segment<3>(6); // Linear acceleration

    // First compute the SE(3) part using (v, w)
    typename SE3<Scalar>::TangentVector se3_element;
    se3_element << v, w;
    SE3<Scalar> pose = SE3<Scalar>::computeExp(se3_element);

    // Compute the velocity part
    // For small rotations or zero angular velocity
    if (w.norm() < Types<Scalar>::epsilon()) {
      return SE23(pose, a);
    }

    // For finite rotations, integrate the velocity
    const Scalar theta = w.norm();
    const Vector3 w_normalized = w / theta;
    const Matrix3 w_hat = SO3<Scalar>::hat(w_normalized);

    // Integration matrix for acceleration
    Matrix3 J = Matrix3::Identity() + (Scalar(1) - std::cos(theta)) * w_hat +
                (theta - std::sin(theta)) * w_hat * w_hat;

    return SE23(pose, J * a / theta);
  }
  /**
   * @brief Computes the logarithmic map from the SE23 group to its Lie algebra
   * se_2(3).
   * @return The corresponding element in the Lie algebra (a 9D vector).
   */
  TangentVector computeLog() const {
    // First get the SE(3) part
    typename SE3<Scalar>::TangentVector se3_part = m_pose.computeLog();
    Vector3 v = se3_part.template head<3>();
    Vector3 w = se3_part.template tail<3>();

    // For small rotations or zero angular velocity
    TangentVector result;
    if (w.norm() < Types<Scalar>::epsilon()) {
      result << v, w, m_velocity;
      return result;
    }

    // For finite rotations, compute acceleration
    const Scalar theta = w.norm();
    const Vector3 w_normalized = w / theta;
    const Matrix3 w_hat = SO3<Scalar>::hat(w_normalized);

    // Integration matrix inverse
    Matrix3 J_inv =
        Matrix3::Identity() - Scalar(0.5) * w_hat +
        (Scalar(1) - theta * std::cos(theta / Scalar(2)) /
                         (Scalar(2) * std::sin(theta / Scalar(2)))) /
            (theta * theta) * (w_hat * w_hat);

    result << v, w, J_inv * m_velocity * theta;
    return result;
  }
  /**
   * @brief Computes the adjoint representation of the current SE23 element.
   * @return The adjoint matrix.
   */
  AdjointMatrix computeAdjoint() const {
    AdjointMatrix Ad = AdjointMatrix::Zero();
    Matrix3 R = m_pose.rotation().matrix();
    Matrix3 p_hat = SO3<Scalar>::hat(m_pose.translation());
    Matrix3 v_hat = SO3<Scalar>::hat(m_velocity);

    // Upper-left block: rotation
    Ad.template block<3, 3>(0, 0) = R;
    // Upper-middle block: position cross rotation
    Ad.template block<3, 3>(0, 3) = p_hat * R;
    // Upper-right block: velocity cross rotation
    Ad.template block<3, 3>(0, 6) = v_hat * R;
    // Middle-middle block: rotation
    Ad.template block<3, 3>(3, 3) = R;
    // Bottom-bottom block: rotation
    Ad.template block<3, 3>(6, 6) = R;

    return Ad;
  }
  /**
   * @brief Checks if the current SE23 element is approximately equal to
   * another.
   * @param other The other SE23 element to compare with.
   * @param eps The tolerance for approximation.
   * @return True if the elements are approximately equal, false otherwise.
   */
  bool computeIsApprox(const SE23 &other,
                       const Scalar &eps = Types<Scalar>::epsilon()) const {
    return m_pose.computeIsApprox(other.m_pose, eps) &&
           m_velocity.isApprox(other.m_velocity, eps);
  }
  /**
   * @brief Applies the group action of the current SE23 element on a
   * point-velocity pair.
   * @param point_vel The point-velocity pair (6D vector: 3D point, 3D
   * velocity).
   * @return The transformed point-velocity pair.
   */
  typename Base::ActionVector
  computeAction(const typename Base::ActionVector &point_vel) const {
    Vector3 point = point_vel.template head<3>();
    Vector3 vel = point_vel.template segment<3>(3);

    // Transform position and combine velocities
    Vector3 transformed_point = m_pose.computeAction(point);
    Vector3 transformed_vel = m_pose.rotation().computeAction(vel) + m_velocity;

    typename Base::ActionVector result;
    result.resize(6);
    result << transformed_point, transformed_vel;
    return result;
  }

  /**
   * @brief Hat operator - maps a 9D Lie algebra vector to a 5x5 matrix
   * representation. This is a placeholder, actual implementation depends on the
   * specific representation.
   * @param v The 9D Lie algebra vector.
   * @return The 5x5 matrix representation.
   */
  static Matrix hat(const TangentVector &v) {
    // For SE_2(3), the hat operator maps a 9D vector to a 5x5 matrix
    // This is a placeholder, actual implementation depends on the specific
    // representation
    Matrix result = Matrix::Zero();
    // ... implement hat operator for SE_2(3)
    return result;
  }

  /**
   * @brief Vee operator - inverse of hat, maps a 5x5 matrix representation to a
   * 9D Lie algebra vector. This is a placeholder, actual implementation depends
   * on the specific representation.
   * @param X The 5x5 matrix representation.
   * @return The 9D Lie algebra vector.
   */
  static TangentVector vee(const Matrix &X) {
    // For SE_2(3), the vee operator maps a 5x5 matrix to a 9D vector
    // This is a placeholder, actual implementation depends on the specific
    // representation
    TangentVector result = TangentVector::Zero();
    // ... implement vee operator for SE_2(3)
    return result;
  }

  /**
   * @brief Computes the adjoint representation of a Lie algebra element for
   * SE23. This is a placeholder, actual implementation depends on the specific
   * representation.
   * @param v The element of the Lie algebra in vector form.
   * @return The adjoint matrix.
   */
  static AdjointMatrix computeAd(const TangentVector &v) {
    // For SE_2(3), the adjoint operator maps a 9D vector to a 9x9 matrix
    // This is a placeholder, actual implementation depends on the specific
    // representation
    AdjointMatrix result = AdjointMatrix::Zero();
    // ... implement adjoint operator for SE_2(3)
    return result;
  }

  /**
   * @brief Generates a random SE23 element.
   * @tparam Generator The type of the random number generator.
   * @param gen The random number generator.
   * @return A random SE23 element.
   */
  template <typename Generator>
  static SE23<Scalar> computeRandom(Generator &gen) {
    return SE23(SE3<Scalar>::computeRandom(gen),
                Types<Scalar>::template randomVector<3>(gen));
  }

  /**
   * @brief Prints the SE23 element to an output stream.
   * @param os The output stream.
   * @return The output stream.
   */
  std::ostream &print(std::ostream &os) const {
    os << "SE23(pose=" << m_pose << ", velocity=" << m_velocity.transpose()
       << ")";
    return os;
  }

  /**
   * @brief Gets the type name of the SE23 class.
   * @return A string view of the type name.
   */
  static constexpr std::string_view getTypeName() { return "SE23"; }

  /**
   * @brief Checks if the current SE23 element is valid.
   * @return True if both pose and velocity components are valid, false
   * otherwise.
   */
  bool computeIsValid() const {
    return m_pose.computeIsValid() && m_velocity.allFinite();
  }

  /**
   * @brief Normalizes the SE23 element.
   * Normalizes the pose component.
   */
  void computeNormalize() { m_pose.computeNormalize(); }

private:
  SE3<Scalar> m_pose; ///< Rigid body transformation
  Vector3 m_velocity; ///< Linear velocity
};

} // namespace sofa::component::cosserat::liegroups

// #endif // SOFA_COMPONENT_COSSERAT_LIEGROUPS_SE23_H