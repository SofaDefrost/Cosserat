// This file defines the Sim(3) (Similarity group in 3D) class, representing
// similarity transformations in 3D space (rotation + translation + scaling).

#pragma once

#include <Eigen/Geometry>
#include <vector>
#include "LieGroupBase.h"
#include "SO3.h"
#include "Types.h"

namespace sofa::component::cosserat::liegroups {

/**
 * @brief Implementation of Sim(3), the Similarity group in 3D
 *
 * This class implements the group of similarity transformations in 3D space.
 * Elements consist of a rotation (SO(3)), a translation vector, and a scaling factor.
 * The group is 7-dimensional with Lie algebra sim(3) having elements [ω, v, s]
 * where ω ∈ ℝ³ (rotation), v ∈ ℝ³ (translation), and s ∈ ℝ (log-scale).
 *
 * @tparam _Scalar The scalar type (must be a floating-point type)
 */
template<typename _Scalar>
class Sim3 : public LieGroupBase<Sim3<_Scalar>, _Scalar, 4, 7, 3> {
public:
    using Base = LieGroupBase<Sim3<_Scalar>, _Scalar, 4, 7, 3>;
    using Scalar = typename Base::Scalar;
    using Vector4 = typename Base::Vector;

    using TangentVector = typename Base::TangentVector;
    using ActionVector = typename Base::ActionVector;

    using Matrix4 = typename Base::Matrix;
    using AdjointMatrix = typename Base::AdjointMatrix;
    using JacobianMatrix = typename Base::JacobianMatrix;
    using AlgebraMatrix = typename Base::AlgebraMatrix;

    using SO3Type = SO3<Scalar>;
    using Vector3 = typename SO3Type::Vector;
    using Matrix3 = typename SO3Type::Matrix;
    using Quaternion = typename SO3Type::Quaternion;

    // ========== Constructors ==========

    /** @brief Default constructor creates identity transformation */
    Sim3() : m_rotation(), m_translation(Vector3::Zero()), m_scale(Scalar(1)) {}

    /** @brief Construct from rotation, translation, and scale */
    Sim3(const SO3Type& rotation, const Vector3& translation, const Scalar& scale)
        : m_rotation(rotation), m_translation(translation), m_scale(scale) {}

    /** @brief Construct from quaternion, translation, and scale */
    Sim3(const Quaternion& quat, const Vector3& translation, const Scalar& scale)
        : m_rotation(quat), m_translation(translation), m_scale(scale) {}

    /** @brief Construct from 4x4 similarity matrix */
    explicit Sim3(const Matrix4& matrix) {
        // Extract scale from the matrix norm
        Matrix3 R = matrix.template block<3, 3>(0, 0);
        Scalar scale_squared = R.squaredNorm() / Scalar(3);
        m_scale = std::sqrt(scale_squared);

        // Extract rotation (normalized)
        R /= m_scale;
        m_rotation = SO3Type(R);

        // Extract translation
        m_translation = matrix.template block<3, 1>(0, 3);
    }

    // ========== CRTP Implementations ==========

    static Sim3 computeIdentity() noexcept { return Sim3(); }

    Sim3 computeInverse() const {
        SO3Type inv_rot = m_rotation.inverse();
        Scalar inv_scale = Scalar(1) / m_scale;
        return Sim3(inv_rot, -(inv_rot.act(m_translation)) * inv_scale, inv_scale);
    }

    Sim3 compose(const Sim3& other) const {
        return Sim3(m_rotation * other.m_rotation,
                   m_rotation.act(other.m_translation * m_scale) + m_translation,
                   m_scale * other.m_scale);
    }

    ActionVector computeAction(const ActionVector& point) const noexcept {
        return m_rotation.act(point * m_scale) + m_translation;
    }

    static Sim3 computeExp(const TangentVector& xi) {
        // Extract components: [ω, v, s] where s = log(scale)
        Vector3 omega = xi.template head<3>();
        Vector3 v = xi.template segment<3>(3);
        Scalar s = xi(6);

        // Compute rotation
        SO3Type R = SO3Type::exp(omega);

        // Compute scale
        Scalar scale = std::exp(s);

        // Compute translation using Rodrigues' formula for the translation part
        Scalar theta = omega.norm();
        Vector3 translation;

        if (theta < Types<Scalar>::epsilon()) {
            // Small angle approximation
            translation = v;
        } else {
            // General case
            Matrix3 omega_hat = SO3Type::hat(omega);
            Matrix3 omega_hat2 = omega_hat * omega_hat;
            Scalar theta2 = theta * theta;
            Scalar theta3 = theta2 * theta;

            // V matrix for Sim(3)
            Matrix3 V = Matrix3::Identity() +
                       ((Scalar(1) - std::cos(theta)) / theta2) * omega_hat +
                       ((theta - std::sin(theta)) / theta3) * omega_hat2;

            translation = V * v;
        }

        return Sim3(R, translation, scale);
    }

    TangentVector computeLog() const {
        // Get rotation logarithm
        Vector3 omega = m_rotation.log();

        // Compute translation logarithm
        Scalar theta = omega.norm();
        Vector3 v;

        if (theta < Types<Scalar>::epsilon()) {
            // Small angle approximation
            v = m_translation;
        } else {
            // General case - solve for v in: translation = V * v
            Matrix3 omega_hat = SO3Type::hat(omega);
            Matrix3 omega_hat2 = omega_hat * omega_hat;
            Scalar theta2 = theta * theta;
            Scalar theta3 = theta2 * theta;

            Matrix3 V = Matrix3::Identity() +
                       ((Scalar(1) - std::cos(theta)) / theta2) * omega_hat +
                       ((theta - std::sin(theta)) / theta3) * omega_hat2;

            // Solve V * v = translation
            v = V.inverse() * m_translation;
        }

        // Log scale
        Scalar s = std::log(m_scale);

        TangentVector result;
        result.template head<3>() = omega;
        result.template segment<3>(3) = v;
        result(6) = s;

        return result;
    }

    AdjointMatrix computeAdjoint() const {
        AdjointMatrix Ad = AdjointMatrix::Zero();
        Matrix3 R = m_rotation.matrix();

        // For Sim(3), adjoint is more complex due to scaling
        // Ad = [R, 0, 0; 0, R, 0; 0, 0, 1] but with scaling effects
        Ad.template block<3, 3>(0, 0) = R;
        Ad.template block<3, 3>(3, 3) = R;
        Ad.template block<3, 3>(0, 3) = SO3Type::hat(m_translation) * R;
        Ad.template block<3, 3>(0, 0) += SO3Type::hat(omega_from_rotation()) * R;

        // Scaling affects the translation part
        Ad.template block<3, 3>(3, 3) *= m_scale;
        Ad.template block<3, 3>(0, 3) *= m_scale;

        // The last row/column for scaling
        Ad(6, 6) = Scalar(1);

        return Ad;
    }

    // ========== CRTP Interface Implementations ==========

    bool computeIsApprox(const Sim3& other, const Scalar& eps) const {
        return m_rotation.isApprox(other.m_rotation, eps) &&
               m_translation.isApprox(other.m_translation, eps) &&
               Types<Scalar>::isApprox(m_scale, other.m_scale, eps);
    }

    // ========== Accessors ==========

    const SO3Type& rotation() const { return m_rotation; }
    SO3Type& rotation() { return m_rotation; }

    const Vector3& translation() const { return m_translation; }
    Vector3& translation() { return m_translation; }

    const Scalar& scale() const { return m_scale; }
    Scalar& scale() { return m_scale; }

    Matrix4 matrix() const {
        Matrix4 T = Matrix4::Identity();
        Matrix3 R_scaled = m_rotation.matrix() * m_scale;
        T.template block<3, 3>(0, 0) = R_scaled;
        T.template block<3, 1>(0, 3) = m_translation;
        return T;
    }

    /**
     * @brief Print the Sim3 element to an output stream.
     */
    std::ostream& print(std::ostream& os) const {
        os << "Sim3(R=" << m_rotation << ", t=[" << m_translation[0] << ", "
           << m_translation[1] << ", " << m_translation[2] << "], s=" << m_scale << ")";
        return os;
    }

private:
    /**
     * @brief Helper to extract rotation vector from current rotation
     */
    Vector3 omega_from_rotation() const {
        return m_rotation.log();
    }

    SO3Type m_rotation;
    Vector3 m_translation;
    Scalar m_scale;
};

// ========== Type Aliases ==========
using Sim3f = Sim3<float>;
using Sim3d = Sim3<double>;

} // namespace sofa::component::cosserat::liegroups