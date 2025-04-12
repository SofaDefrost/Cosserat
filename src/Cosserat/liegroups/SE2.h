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

#ifndef SOFA_COMPONENT_COSSERAT_LIEGROUPS_SE2_H
#define SOFA_COMPONENT_COSSERAT_LIEGROUPS_SE2_H

#include "LieGroupBase.h"
#include "LieGroupBase.inl"
#include "SO2.h"
#include <Eigen/Geometry>

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
 * @tparam _Scalar The scalar type (must be a floating-point type)
 */
template<typename _Scalar>
class SE2 : public LieGroupBase<_Scalar, 3>,
           public LieGroupOperations<SE2<_Scalar>> {
public:
    using Base = LieGroupBase<_Scalar, 3>;
    using Scalar = typename Base::Scalar;
    using Vector = typename Base::Vector;
    using Matrix = typename Base::Matrix;
    using TangentVector = typename Base::TangentVector;
    using AdjointMatrix = typename Base::AdjointMatrix;
    
    using Vector2 = Eigen::Matrix<Scalar, 2, 1>;
    using Matrix2 = Eigen::Matrix<Scalar, 2, 2>;
    static constexpr int Dim = Base::Dim;

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
     * @brief Group composition (rigid transformation composition)
     */
    SE2 operator*(const SE2& other) const {
        return SE2(m_rotation * other.m_rotation,
                  m_translation + m_rotation.act(other.m_translation));
    }

    /**
     * @brief Inverse element (inverse transformation)
     */
    SE2 inverse() const override {
        SO2<Scalar> inv_rot = m_rotation.inverse();
        return SE2(inv_rot, -(inv_rot.act(m_translation)));
    }

    /**
     * @brief Exponential map from Lie algebra to SE(2)
     * @param algebra_element Vector in ℝ³ representing (vx, vy, ω)
     */
    SE2 exp(const TangentVector& algebra_element) const override {
        const Scalar& theta = algebra_element[2];
        const Vector2 v(algebra_element[0], algebra_element[1]);

        if (std::abs(theta) < Types<Scalar>::epsilon()) {
            // For small rotations, use first-order approximation
            return SE2(SO2<Scalar>(theta), v);
        }

        // Exact formula for finite rotations
        const Scalar sin_theta = std::sin(theta);
        const Scalar cos_theta = std::cos(theta);
        const Scalar theta_inv = Scalar(1) / theta;

        Matrix2 V;
        V << sin_theta * theta_inv, -(Scalar(1) - cos_theta) * theta_inv,
             (Scalar(1) - cos_theta) * theta_inv, sin_theta * theta_inv;

        return SE2(SO2<Scalar>(theta), V * v);
    }

    /**
     * @brief Logarithmic map from SE(2) to Lie algebra
     * @return Vector in ℝ³ representing (vx, vy, ω)
     */
    TangentVector log() const override {
        const Scalar theta = m_rotation.angle();
        TangentVector result;

        if (std::abs(theta) < Types<Scalar>::epsilon()) {
            // For small rotations, use first-order approximation
            result << m_translation, theta;
            return result;
        }

        // Exact formula for finite rotations
        const Scalar sin_theta = std::sin(theta);
        const Scalar cos_theta = std::cos(theta);
        const Scalar theta_inv = Scalar(1) / theta;
        const Scalar half_theta = theta * Scalar(0.5);

        Matrix2 V_inv;
        V_inv << half_theta * cos_theta / sin_theta, -half_theta,
                half_theta, half_theta * cos_theta / sin_theta;

        Vector2 v = V_inv * m_translation;
        result << v, theta;
        return result;
    }

    /**
     * @brief Adjoint representation
     */
    AdjointMatrix adjoint() const override {
        AdjointMatrix Ad = AdjointMatrix::Zero();
        // Rotation block
        Ad.template block<2,2>(0,0) = m_rotation.matrix();
        // Translation block
        Ad(0,2) = -m_translation.y();
        Ad(1,2) = m_translation.x();
        // Bottom-right corner is 1 for rotation
        Ad(2,2) = Scalar(1);
        return Ad;
    }

    /**
     * @brief Group action on a point
     */
    Vector2 act(const Vector2& point) const {
        return m_rotation.act(point) + m_translation;
    }

    /**
     * @brief Override of act for 3D vectors (ignores z component)
     */
    Vector act(const Vector& point) const override {
        Vector2 result = act(point.template head<2>());
        Vector return_val;
        return_val << result, point[2];
        return return_val;
    }

    /**
     * @brief Check if approximately equal to another element
     */
    bool isApprox(const SE2& other, 
                  const Scalar& eps = Types<Scalar>::epsilon()) const {
        return m_rotation.isApprox(other.m_rotation, eps) &&
               m_translation.isApprox(other.m_translation, eps);
    }

    /**
     * @brief Get the identity element
     */
    static const SE2& identity() {
        static const SE2 id;
        return id;
    }

    /**
     * @brief Get the dimension of the Lie algebra (3 for SE(2))
     */
    int algebraDimension() const override { return 3; }

    /**
     * @brief Get the dimension of the space the group acts on (2 for SE(2))
     */
    int actionDimension() const override { return 2; }

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
     * @brief Get the homogeneous transformation matrix
     */
    Eigen::Matrix<Scalar, 3, 3> matrix() const {
        Eigen::Matrix<Scalar, 3, 3> T = Eigen::Matrix<Scalar, 3, 3>::Identity();
        T.template block<2,2>(0,0) = m_rotation.matrix();
        T.template block<2,1>(0,2) = m_translation;
        return T;
    }

private:
    SO2<Scalar> m_rotation;    ///< Rotation component
    Vector2 m_translation;     ///< Translation component
};

} // namespace sofa::component::cosserat::liegroups

#include "SE2.inl"

#endif // SOFA_COMPONENT_COSSERAT_LIEGROUPS_SE2_H
