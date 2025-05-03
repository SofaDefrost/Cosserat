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

#include <Eigen/Geometry.h>  // Include Eigen first
#include <Cosserat/liegroups/Types.h>        // Then our type system
#include <Cosserat/liegroups/LieGroupBase.h> // Then the base class interface
#include <Cosserat/liegroups/SO3.h>         // Then other dependencies

// Forward declaration outside the namespace
namespace sofa::component::cosserat::liegroups {
template<typename Scalar> class SE3;
}

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
 * @tparam Scalar The scalar type (must be a floating-point type)
 * @tparam _Dim The dimension of the group representation, here 6
*/
static int Dim = 6;
template<typename Scalar>
class SE3 final : public LieGroupBase<SE3<Scalar>, Scalar, 6>
{  // Use default values for _AlgebraDim and _ActionDim
    static_assert(std::is_floating_point<Scalar>::value,
                 "Scalar type must be a floating-point type");
 public:
    using Base = LieGroupBase<SE3<Scalar>, Scalar, 6>;  // Same here
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
    SE3 inverse() const {
        SO3<Scalar> inv_rot = m_rotation.inverse();
        return SE3(inv_rot, -(inv_rot.act(m_translation)));
    }

    /**
     * @brief Exponential map from Lie algebra to SE(3)
     * @param twist Vector in ℝ⁶ representing (v, ω)
     */
    static SE3 exp(const TangentVector& twist) {
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
    TangentVector log() const {
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
    AdjointMatrix adjoint() const {
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
    Vector act(const Vector& point) const {
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
    static SE3 Identity() noexcept {
        return SE3();
    }

    /**
     * @brief Get the dimension of the Lie algebra (6 for SE(3))
     */
    int algebraDimension() const { return 6; }

    /**
     * @brief Get the dimension of the space the group acts on (3 for SE(3))
     */
    int actionDimension() const { return 6; }

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

    // Required methods to match base class interface
    SE3 compose(const SE3& other) const noexcept {
        return (*this) * other;
    }

    SE3 computeInverse() const {
        return inverse();
    }

    static SE3 computeExp(const TangentVector& v) noexcept {
        return exp(v);
    }

    TangentVector computeLog() const {
        return log();
    }

    AdjointMatrix computeAdjoint() const noexcept {
        return adjoint();
    }

    bool computeIsApprox(const SE3& other, const Scalar& eps) const noexcept {
        return isApprox(other, eps);
    }

    static SE3 computeIdentity() noexcept {
        return SE3();
    }

    typename Base::ActionVector computeAction(const typename Base::ActionVector& point) const noexcept {
        return act(point);
    }

    /**
     * @brief Baker-Campbell-Hausdorff formula for SE(3)
     * 
     * For SE(3), the BCH formula has a closed form up to second order:
     * BCH(X,Y) = X + Y + 1/2[X,Y] + higher order terms
     * where [X,Y] is the Lie bracket for se(3).
     */
    static TangentVector BCH(const TangentVector& X, const TangentVector& Y) {
        // Extract linear and angular components
        const auto& v1 = X.template head<3>();
        const auto& w1 = X.template tail<3>();
        const auto& v2 = Y.template head<3>();
        const auto& w2 = Y.template tail<3>();

        // Compute Lie bracket components
        const auto w1_cross_w2 = w1.cross(w2); // Angular x Angular
        const auto w1_cross_v2 = w1.cross(v2); // Angular x Linear
        const auto v1_cross_w2 = v1.cross(w2); // Linear x Angular
        // Combine terms for the BCH formula up to second order
        TangentVector result;
        result.template head<3>() = v1 + v2 + Scalar(0.5) * (w1_cross_v2 - v1_cross_w2);
        result.template tail<3>() = w1 + w2 + Scalar(0.5) * w1_cross_w2;

        return result;
    }

private:
    SO3<Scalar> m_rotation;     ///< Rotation component
    Vector3 m_translation;      ///< Translation component
};  // end of class SE3

} // namespace sofa::component::cosserat::liegroups


#endif // SOFA_COMPONENT_COSSERAT_LIEGROUPS_SE3_H
