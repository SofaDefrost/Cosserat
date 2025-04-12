/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                  *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
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
* along with this program. If not, see <http://www.gnu.org/licenses/\>.        *
******************************************************************************/

#ifndef SOFA_COMPONENT_COSSERAT_LIEGROUPS_SGAL3_H
#define SOFA_COMPONENT_COSSERAT_LIEGROUPS_SGAL3_H

#include "LieGroupBase.h"
#include "LieGroupBase.inl"
#include "SE3.h"
#include <Eigen/Geometry>

namespace sofa::component::cosserat::liegroups {

/**
 * @brief Implementation of SGal(3), the Special Galilean group in 3D
 * 
 * This class implements the group of Galilean transformations in 3D space.
 * Elements of SGal(3) are represented as a combination of:
 * - An SE(3) transformation (rotation and position)
 * - A 3D velocity vector
 * - A time parameter
 * 
 * The Lie algebra sgal(3) consists of vectors in ℝ¹⁰, where:
 * - First three components represent linear velocity
 * - Next three components represent angular velocity
 * - Next three components represent boost (velocity change)
 * - Last component represents time rate
 * 
 * @tparam _Scalar The scalar type (must be a floating-point type)
 */
template<typename _Scalar>
class SGal3 : public LieGroupBase<_Scalar, 10>,
             public LieGroupOperations<SGal3<_Scalar>> {
public:
    using Base = LieGroupBase<_Scalar, 10>;
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
     * @brief Default constructor creates identity element
     */
    SGal3() : m_pose(), m_velocity(Vector3::Zero()), m_time(0) {}

    /**
     * @brief Construct from pose, velocity, and time
     */
    SGal3(const SE3<Scalar>& pose, const Vector3& velocity, const Scalar& time)
        : m_pose(pose), m_velocity(velocity), m_time(time) {}

    /**
     * @brief Construct from rotation, position, velocity, and time
     */
    SGal3(const SO3<Scalar>& rotation, const Vector3& position,
          const Vector3& velocity, const Scalar& time)
        : m_pose(rotation, position), m_velocity(velocity), m_time(time) {}

    /**
     * @brief Group composition (Galilean transformation composition)
     */
    SGal3 operator*(const SGal3& other) const {
        return SGal3(m_pose * other.m_pose,
                    m_velocity + m_pose.rotation().act(other.m_velocity),
                    m_time + other.m_time);
    }

    /**
     * @brief Inverse element
     */
    SGal3 inverse() const override {
        SE3<Scalar> inv_pose = m_pose.inverse();
        return SGal3(inv_pose,
                    -inv_pose.rotation().act(m_velocity),
                    -m_time);
    }

    /**
     * @brief Exponential map from Lie algebra to SGal(3)
     * @param algebra_element Vector in ℝ¹⁰ representing (v, ω, β, τ)
     */
    SGal3 exp(const TangentVector& algebra_element) const override {
        Vector3 v = algebra_element.template segment<3>(0);  // Linear velocity
        Vector3 w = algebra_element.template segment<3>(3);  // Angular velocity
        Vector3 beta = algebra_element.template segment<3>(6);  // Boost
        Scalar tau = algebra_element[9];  // Time rate

        // First compute the SE(3) part using (v, w)
        typename SE3<Scalar>::TangentVector se3_element;
        se3_element << v, w;
        SE3<Scalar> pose = SE3<Scalar>().exp(se3_element);

        // For small rotations or zero angular velocity
        if (w.norm() < Types<Scalar>::epsilon()) {
            return SGal3(pose, beta, tau);
        }

        // For finite rotations, integrate the velocity with boost
        const Scalar theta = w.norm();
        const Vector3 w_normalized = w / theta;
        const Matrix3 w_hat = SO3<Scalar>::hat(w_normalized);
        
        // Integration matrix for boost
        Matrix3 J = Matrix3::Identity() + 
                   (Scalar(1) - std::cos(theta)) * w_hat +
                   (theta - std::sin(theta)) * w_hat * w_hat;
        
        return SGal3(pose, J * beta / theta, tau);
    }

    /**
     * @brief Logarithmic map from SGal(3) to Lie algebra
     * @return Vector in ℝ¹⁰ representing (v, ω, β, τ)
     */
    TangentVector log() const override {
        // First get the SE(3) part
        typename SE3<Scalar>::TangentVector se3_part = m_pose.log();
        Vector3 v = se3_part.template head<3>();
        Vector3 w = se3_part.template tail<3>();
        
        // For small rotations or zero angular velocity
        TangentVector result;
        if (w.norm() < Types<Scalar>::epsilon()) {
            result << v, w, m_velocity, m_time;
            return result;
        }

        // For finite rotations, compute boost
        const Scalar theta = w.norm();
        const Vector3 w_normalized = w / theta;
        const Matrix3 w_hat = SO3<Scalar>::hat(w_normalized);
        
        // Integration matrix inverse
        Matrix3 J_inv = Matrix3::Identity() -
                       Scalar(0.5) * w_hat +
                       (Scalar(1) - theta * std::cos(theta / Scalar(2)) / (Scalar(2) * std::sin(theta / Scalar(2)))) /
                       (theta * theta) * (w_hat * w_hat);
        
        result << v, w, J_inv * m_velocity * theta, m_time;
        return result;
    }

    /**
     * @brief Adjoint representation
     */
    AdjointMatrix adjoint() const override {
        AdjointMatrix Ad = AdjointMatrix::Zero();
        Matrix3 R = m_pose.rotation().matrix();
        Matrix3 p_hat = SO3<Scalar>::hat(m_pose.translation());
        Matrix3 v_hat = SO3<Scalar>::hat(m_velocity);
        
        // Upper-left block: rotation
        Ad.template block<3,3>(0,0) = R;
        // Upper-middle block: position cross rotation
        Ad.template block<3,3>(0,3) = p_hat * R;
        // Upper-right block: velocity cross rotation
        Ad.template block<3,3>(0,6) = v_hat * R;
        // Middle-middle block: rotation
        Ad.template block<3,3>(3,3) = R;
        // Bottom-bottom block: rotation
        Ad.template block<3,3>(6,6) = R;
        // Time row and column remain zero except diagonal
        Ad(9,9) = Scalar(1);
        
        return Ad;
    }

    /**
     * @brief Group action on a point, velocity, and time
     */
    Vector act(const Vector& point_vel_time) const override {
        Vector3 point = point_vel_time.template head<3>();
        Vector3 vel = point_vel_time.template segment<3>(3);
        Vector3 boost = point_vel_time.template segment<3>(6);
        Scalar t = point_vel_time[9];
        
        // Transform position and combine velocities with time evolution
        Vector3 transformed_point = m_pose.act(point) + m_velocity * t;
        Vector3 transformed_vel = m_pose.rotation().act(vel) + m_velocity;
        Vector3 transformed_boost = m_pose.rotation().act(boost);
        
        Vector result;
        result.resize(10);
        result << transformed_point, transformed_vel, transformed_boost, t + m_time;
        return result;
    }

    /**
     * @brief Check if approximately equal to another element
     */
    bool isApprox(const SGal3& other, 
                  const Scalar& eps = Types<Scalar>::epsilon()) const {
        return m_pose.isApprox(other.m_pose, eps) &&
               m_velocity.isApprox(other.m_velocity, eps) &&
               std::abs(m_time - other.m_time) <= eps;
    }

    /**
     * @brief Get the identity element
     */
    static const SGal3& identity() {
        static const SGal3 id;
        return id;
    }

    /**
     * @brief Get the dimension of the Lie algebra (10 for SGal(3))
     */
    int algebraDimension() const override { return 10; }

    /**
     * @brief Get the dimension of the space the group acts on (7 for SGal(3))
     */
    int actionDimension() const override { return 7; }

    /**
     * @brief Access the pose component
     */
    const SE3<Scalar>& pose() const { return m_pose; }
    SE3<Scalar>& pose() { return m_pose; }

    /**
     * @brief Access the velocity component
     */
    const Vector3& velocity() const { return m_velocity; }
    Vector3& velocity() { return m_velocity; }

    /**
     * @brief Access the time component
     */
    const Scalar& time() const { return m_time; }
    Scalar& time() { return m_time; }

    /**
     * @brief Get the extended homogeneous transformation matrix
     */
    Eigen::Matrix<Scalar, 6, 6> matrix() const {
        Eigen::Matrix<Scalar, 6, 6> T = Eigen::Matrix<Scalar, 6, 6>::Identity();
        T.template block<4,4>(0,0) = m_pose.matrix();
        T.template block<3,1>(0,4) = m_velocity;
        T(4,5) = m_time;
        return T;
    }

private:
    SE3<Scalar> m_pose;      ///< Rigid body transformation
    Vector3 m_velocity;      ///< Linear velocity
    Scalar m_time;          ///< Time parameter
};

} // namespace sofa::component::cosserat::liegroups

#include "SGal3.inl"

#endif // SOFA_COMPONENT_COSSERAT_LIEGROUPS_SGAL3_H
