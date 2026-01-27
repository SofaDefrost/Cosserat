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

#ifndef SOFA_COMPONENT_COSSERAT_LIEGROUPS_SE3_H
#define SOFA_COMPONENT_COSSERAT_LIEGROUPS_SE3_H

#include "LieGroupBase.h"
#include "LieGroupBase.inl"
#include "SO3.h"
#include <Eigen/Geometry>

namespace sofa::component::cosserat::liegroups {

/**
 * @brief Implementation of SE(3), the Special Euclidean group in 3D
 * 
 * This class implements the group of rigid body transformations in 3D space.
 * Elements of SE(3) are represented as a combination of:
 * - An SO(3) rotation (using quaternions)
 * - A 3D translation vector
 * 
 * The Lie algebra se(3) consists of vectors in ℝ⁶, where:
 * - The first three components represent the linear velocity
 * - The last three components represent the angular velocity
 * 
 * @tparam _Scalar The scalar type (must be a floating-point type)
 */
template<typename _Scalar>
class SE3 : public LieGroupBase<_Scalar, 6>,
           public LieGroupOperations<SE3<_Scalar>> {
public:
    using Base = LieGroupBase<_Scalar, 6>;
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
     * @brief Default constructor creates identity transformation
     */
    SE3() : m_rotation(), m_translation(Vector3::Zero()) {}

    /**
     * @brief Construct from rotation and translation
     */
    SE3(const SO3<Scalar>& rotation, const Vector3& translation)
        : m_rotation(rotation), m_translation(translation) {}

    /**
     * @brief Construct from homogeneous transformation matrix
     */
    explicit SE3(const Matrix4& T)
        : m_rotation(T.template block<3,3>(0,0)),
          m_translation(T.template block<3,1>(0,3)) {}

    /**
     * @brief Group composition (rigid transformation composition)
     */
    SE3 operator*(const SE3& other) const {
        return SE3(m_rotation * other.m_rotation,
                  m_translation + m_rotation.act(other.m_translation));
    }

    /**
     * @brief Inverse element (inverse transformation)
     */
    SE3 inverse() const override {
        SO3<Scalar> inv_rot = m_rotation.inverse();
        return SE3(inv_rot, -(inv_rot.act(m_translation)));
    }

    /**
     * @brief Exponential map from Lie algebra to SE(3)
     * @param twist Vector in ℝ⁶ representing (v, ω)
     */
    SE3 exp(const TangentVector& twist) const override {
        Vector3 v = twist.template head<3>();
        Vector3 omega = twist.template tail<3>();
        const Scalar theta = omega.norm();

        SO3<Scalar> R;
        Vector3 t;

        if (theta < Types<Scalar>::epsilon()) {
            // For small rotations, use first-order approximation
            R = SO3<Scalar>().exp(omega);
            t = v;
        } else {
            // Full exponential formula
            R = SO3<Scalar>().exp(omega);
            
            // Compute translation using Rodriguez formula
            Matrix3 omega_hat = SO3<Scalar>::hat(omega);
            Matrix3 V = Matrix3::Identity() +
                       (Scalar(1) - std::cos(theta)) / (theta * theta) * omega_hat +
                       (theta - std::sin(theta)) / (theta * theta * theta) * (omega_hat * omega_hat);
            
            t = V * v;
        }

        return SE3(R, t);
    }

    /**
     * @brief Logarithmic map from SE(3) to Lie algebra
     * @return Vector in ℝ⁶ representing (v, ω)
     */
    TangentVector log() const override {
        // Extract rotation vector
        Vector3 omega = m_rotation.log();
        const Scalar theta = omega.norm();

        TangentVector result;
        Matrix3 V_inv;

        if (theta < Types<Scalar>::epsilon()) {
            // For small rotations, use first-order approximation
            V_inv = Matrix3::Identity();
        } else {
            // Full logarithm formula
            Matrix3 omega_hat = SO3<Scalar>::hat(omega);
            V_inv = Matrix3::Identity() -
                   Scalar(0.5) * omega_hat +
                   (Scalar(1) - theta * std::cos(theta / Scalar(2)) / (Scalar(2) * std::sin(theta / Scalar(2)))) /
                   (theta * theta) * (omega_hat * omega_hat);
        }

        result << V_inv * m_translation, omega;
        return result;
    }

    /**
     * @brief Adjoint representation
     */
    AdjointMatrix adjoint() const override {
        AdjointMatrix Ad = AdjointMatrix::Zero();
        Matrix3 R = m_rotation.matrix();
        Matrix3 t_hat = SO3<Scalar>::hat(m_translation);
        
        // Rotation block
        Ad.template block<3,3>(0,0) = R;
        // Translation block
        Ad.template block<3,3>(0,3) = t_hat * R;
        // Bottom-right block
        Ad.template block<3,3>(3,3) = R;
        
        return Ad;
    }

    /**
     * @brief Group action on a point
     */
    Vector3 act(const Vector3& point) const {
        return m_rotation.act(point) + m_translation;
    }

    /**
     * @brief Override of act for 6D vectors (acts on position part only)
     */
    Vector act(const Vector& point) const override {
        Vector3 pos = act(point.template head<3>());
        Vector result;
        result << pos, point.template tail<3>();
        return result;
    }

    /**
     * @brief Check if approximately equal to another element
     */
    bool isApprox(const SE3& other, 
                  const Scalar& eps = Types<Scalar>::epsilon()) const {
        return m_rotation.isApprox(other.m_rotation, eps) &&
               m_translation.isApprox(other.m_translation, eps);
    }

    /**
     * @brief Get the identity element
     */
    static const SE3& identity() {
        static const SE3 id;
        return id;
    }

    /**
     * @brief Get the dimension of the Lie algebra (6 for SE(3))
     */
    int algebraDimension() const override { return 6; }

    /**
     * @brief Get the dimension of the space the group acts on (3 for SE(3))
     */
    int actionDimension() const override { return 3; }

    /**
     * @brief Access the rotation component
     */
    const SO3<Scalar>& rotation() const { return m_rotation; }
    SO3<Scalar>& rotation() { return m_rotation; }

    /**
     * @brief Access the translation component
     */
    const Vector3& translation() const { return m_translation; }
    Vector3& translation() { return m_translation; }

    /**
     * @brief Get the homogeneous transformation matrix
     */
    Matrix4 matrix() const {
        Matrix4 T = Matrix4::Identity();
        T.template block<3,3>(0,0) = m_rotation.matrix();
        T.template block<3,1>(0,3) = m_translation;
        return T;
    }

private:
    SO3<Scalar> m_rotation;     ///< Rotation component
    Vector3 m_translation;      ///< Translation component
};

} // namespace sofa::component::cosserat::liegroups

#include "SE3.inl"

#endif // SOFA_COMPONENT_COSSERAT_LIEGROUPS_SE3_H
