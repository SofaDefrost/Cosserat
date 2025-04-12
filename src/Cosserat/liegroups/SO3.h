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

#ifndef SOFA_COMPONENT_COSSERAT_LIEGROUPS_SO3_H
#define SOFA_COMPONENT_COSSERAT_LIEGROUPS_SO3_H

#include "LieGroupBase.h"
#include "LieGroupBase.inl"
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
template<typename _Scalar>
class SO3 : public LieGroupBase<_Scalar, 3>,
           public LieGroupOperations<SO3<_Scalar>> {
public:
    using Base = LieGroupBase<_Scalar, 3>;
    using Scalar = typename Base::Scalar;
    using Vector = typename Base::Vector;
    using Matrix = typename Base::Matrix;
    using TangentVector = typename Base::TangentVector;
    using AdjointMatrix = typename Base::AdjointMatrix;
    
    using Quaternion = Eigen::Quaternion<Scalar>;
    static constexpr int Dim = Base::Dim;

    /**
     * @brief Default constructor creates identity rotation
     */
    SO3() : m_quat(Quaternion::Identity()) {}

    /**
     * @brief Construct from angle-axis representation
     * @param angle Rotation angle in radians
     * @param axis Unit vector representing rotation axis
     */
    SO3(const Scalar& angle, const Vector& axis)
        : m_quat(Eigen::AngleAxis<Scalar>(angle, axis.normalized())) {}

    /**
     * @brief Construct from quaternion
     * @param quat Unit quaternion
     */
    explicit SO3(const Quaternion& quat) : m_quat(quat.normalized()) {}

    /**
     * @brief Construct from rotation matrix
     * @param mat 3x3 rotation matrix
     */
    explicit SO3(const Matrix& mat) : m_quat(mat) {}

    /**
     * @brief Group composition (rotation composition)
     */
    SO3 operator*(const SO3& other) const {
        return SO3(m_quat * other.m_quat);
    }

    /**
     * @brief Inverse element (opposite rotation)
     */
    SO3 inverse() const override {
        return SO3(m_quat.conjugate());
    }

    /**
     * @brief Exponential map from Lie algebra to SO(3)
     * @param omega Angular velocity vector in ℝ³
     */
    SO3 exp(const TangentVector& omega) const override {
        const Scalar theta = omega.norm();
        
        if (theta < Types<Scalar>::epsilon()) {
            // For small rotations, use first-order approximation
            return SO3(Quaternion(Scalar(1), 
                                omega.x() * Scalar(0.5),
                                omega.y() * Scalar(0.5),
                                omega.z() * Scalar(0.5)));
        }

        // Use Rodrigues' formula
        const Vector axis = omega / theta;
        const Scalar half_theta = theta * Scalar(0.5);
        const Scalar sin_half_theta = std::sin(half_theta);
        
        return SO3(Quaternion(std::cos(half_theta),
                            axis.x() * sin_half_theta,
                            axis.y() * sin_half_theta,
                            axis.z() * sin_half_theta));
    }

    /**
     * @brief Logarithmic map from SO(3) to Lie algebra
     * @return Angular velocity vector in ℝ³
     */
    TangentVector log() const override {
        // Extract angle-axis representation
        Eigen::AngleAxis<Scalar> aa(m_quat);
        const Scalar theta = aa.angle();
        
        if (theta < Types<Scalar>::epsilon()) {
            // For small rotations, use first-order approximation
            return Vector(m_quat.x() * Scalar(2),
                        m_quat.y() * Scalar(2),
                        m_quat.z() * Scalar(2));
        }

        return aa.axis() * theta;
    }

    /**
     * @brief Adjoint representation
     * For SO(3), this is the rotation matrix itself
     */
    AdjointMatrix adjoint() const override {
        return matrix();
    }

    /**
     * @brief Group action on a point (rotate the point)
     */
    Vector act(const Vector& point) const override {
        return m_quat * point;
    }

    /**
     * @brief Check if approximately equal to another rotation
     */
    bool isApprox(const SO3& other, 
                  const Scalar& eps = Types<Scalar>::epsilon()) const {
        // Handle antipodal representation of same rotation
        return m_quat.coeffs().isApprox(other.m_quat.coeffs(), eps) ||
               m_quat.coeffs().isApprox(-other.m_quat.coeffs(), eps);
    }

    /**
     * @brief Get the identity element
     */
    static const SO3& identity() {
        static const SO3 id;
        return id;
    }

    /**
     * @brief Get the dimension of the Lie algebra (3 for SO(3))
     */
    int algebraDimension() const override { return 3; }

    /**
     * @brief Get the dimension of the space the group acts on (3 for SO(3))
     */
    int actionDimension() const override { return 3; }

    /**
     * @brief Get the rotation matrix representation
     */
    Matrix matrix() const {
        return m_quat.toRotationMatrix();
    }

    /**
     * @brief Get the quaternion representation
     */
    const Quaternion& quaternion() const { return m_quat; }

    /**
     * @brief Convert to angle-axis representation
     */
    Eigen::AngleAxis<Scalar> angleAxis() const {
        return Eigen::AngleAxis<Scalar>(m_quat);
    }

    /**
     * @brief Convert vector to skew-symmetric matrix
     * @param v Vector in ℝ³
     * @return 3x3 skew-symmetric matrix
     */
    static Matrix hat(const Vector& v) {
        Matrix Omega;
        Omega << 0, -v[2], v[1],
                v[2], 0, -v[0],
                -v[1], v[0], 0;
        return Omega;
    }

    /**
     * @brief Convert skew-symmetric matrix to vector
     * @param Omega 3x3 skew-symmetric matrix
     * @return Vector in ℝ³
     */
    static Vector vee(const Matrix& Omega) {
        return Vector(Omega(2,1), Omega(0,2), Omega(1,0));
    }

private:
    Quaternion m_quat;  ///< Unit quaternion representing the rotation
};

} // namespace sofa::component::cosserat::liegroups

#include "SO3.inl"

#endif // SOFA_COMPONENT_COSSERAT_LIEGROUPS_SO3_H
