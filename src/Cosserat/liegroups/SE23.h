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

#ifndef SOFA_COMPONENT_COSSERAT_LIEGROUPS_SE23_H
#define SOFA_COMPONENT_COSSERAT_LIEGROUPS_SE23_H

#include "LieGroupBase.h"
#include "LieGroupBase.inl"
#include "SE3.h"
#include <Eigen/Geometry>

namespace sofa::component::cosserat::liegroups {

/**
 * @brief Implementation of SE_2(3), the extended Special Euclidean group in 3D
 * 
 * This class implements the group of rigid body transformations with linear velocity
 * in 3D space. Elements of SE_2(3) are represented as a combination of:
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
template<typename _Scalar>
class SE23 : public LieGroupBase<_Scalar, 9>,
            public LieGroupOperations<SE23<_Scalar>> {
public:
    using Base = LieGroupBase<_Scalar, 9>;
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
    SE23() : m_pose(), m_velocity(Vector3::Zero()) {}

    /**
     * @brief Construct from pose and velocity
     */
    SE23(const SE3<Scalar>& pose, const Vector3& velocity)
        : m_pose(pose), m_velocity(velocity) {}

    /**
     * @brief Construct from rotation, position, and velocity
     */
    SE23(const SO3<Scalar>& rotation, const Vector3& position, const Vector3& velocity)
        : m_pose(rotation, position), m_velocity(velocity) {}

    /**
     * @brief Group composition (extended pose composition)
     */
    SE23 operator*(const SE23& other) const {
        return SE23(m_pose * other.m_pose,
                   m_velocity + m_pose.rotation().act(other.m_velocity));
    }

    /**
     * @brief Inverse element
     */
    SE23 inverse() const override {
        SE3<Scalar> inv_pose = m_pose.inverse();
        return SE23(inv_pose,
                   -inv_pose.rotation().act(m_velocity));
    }

    /**
     * @brief Exponential map from Lie algebra to SE_2(3)
     * @param algebra_element Vector in ℝ⁹ representing (v, ω, a)
     */
    SE23 exp(const TangentVector& algebra_element) const override {
        Vector3 v = algebra_element.template segment<3>(0);  // Linear velocity
        Vector3 w = algebra_element.template segment<3>(3);  // Angular velocity
        Vector3 a = algebra_element.template segment<3>(6);  // Linear acceleration

        // First compute the SE(3) part using (v, w)
        typename SE3<Scalar>::TangentVector se3_element;
        se3_element << v, w;
        SE3<Scalar> pose = SE3<Scalar>().exp(se3_element);

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
        Matrix3 J = Matrix3::Identity() + 
                   (Scalar(1) - std::cos(theta)) * w_hat +
                   (theta - std::sin(theta)) * w_hat * w_hat;
        
        return SE23(pose, J * a / theta);
    }

    /**
     * @brief Logarithmic map from SE_2(3) to Lie algebra
     * @return Vector in ℝ⁹ representing (v, ω, a)
     */
    TangentVector log() const override {
        // First get the SE(3) part
        typename SE3<Scalar>::TangentVector se3_part = m_pose.log();
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
        Matrix3 J_inv = Matrix3::Identity() -
                       Scalar(0.5) * w_hat +
                       (Scalar(1) - theta * std::cos(theta / Scalar(2)) / (Scalar(2) * std::sin(theta / Scalar(2)))) /
                       (theta * theta) * (w_hat * w_hat);
        
        result << v, w, J_inv * m_velocity * theta;
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
        
        return Ad;
    }

    /**
     * @brief Group action on a point and its velocity
     */
    Vector act(const Vector& point_vel) const override {
        Vector3 point = point_vel.template head<3>();
        Vector3 vel = point_vel.template segment<3>(3);
        
        // Transform position and combine velocities
        Vector3 transformed_point = m_pose.act(point);
        Vector3 transformed_vel = m_pose.rotation().act(vel) + m_velocity;
        
        Vector result;
        result.resize(9);
        result << transformed_point, transformed_vel, point_vel.template tail<3>();
        return result;
    }

    /**
     * @brief Check if approximately equal to another element
     */
    bool isApprox(const SE23& other, 
                  const Scalar& eps = Types<Scalar>::epsilon()) const {
        return m_pose.isApprox(other.m_pose, eps) &&
               m_velocity.isApprox(other.m_velocity, eps);
    }

    /**
     * @brief Get the identity element
     */
    static const SE23& identity() {
        static const SE23 id;
        return id;
    }

    /**
     * @brief Get the dimension of the Lie algebra (9 for SE_2(3))
     */
    int algebraDimension() const override { return 9; }

    /**
     * @brief Get the dimension of the space the group acts on (6 for SE_2(3))
     */
    int actionDimension() const override { return 6; }

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
     * @brief Get the extended homogeneous transformation matrix
     */
    Eigen::Matrix<Scalar, 5, 5> matrix() const {
        Eigen::Matrix<Scalar, 5, 5> T = Eigen::Matrix<Scalar, 5, 5>::Identity();
        T.template block<4,4>(0,0) = m_pose.matrix();
        T.template block<3,1>(0,4) = m_velocity;
        return T;
    }

private:
    SE3<Scalar> m_pose;      ///< Rigid body transformation
    Vector3 m_velocity;      ///< Linear velocity
};

} // namespace sofa::component::cosserat::liegroups

#include "SE23.inl"

#endif // SOFA_COMPONENT_COSSERAT_LIEGROUPS_SE23_H
