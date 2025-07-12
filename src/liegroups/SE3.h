// This file defines the SE3 (Special Euclidean group in 3D) class, representing
// rigid body transformations in 3D space.

#pragma once

#include "LieGroupBase.h"
#include "SO3.h"
#include "Types.h"
#include <Eigen/Geometry>
#include <vector>

namespace sofa::component::cosserat::liegroups {

/**
 * @brief Implementation of SE(3), the Special Euclidean group in 3D
 *
 * This class implements the group of rigid body transformations in 3D space.
 * Elements are represented by an SO(3) rotation and a 3D translation vector.
 *
 * The Lie algebra se(3) consists of 6D vectors (vx, vy, vz, ωx, ωy, ωz),
 * representing translational and angular velocities.
 *
 * @tparam _Scalar The scalar type (must be a floating-point type)
 */
template <typename _Scalar>
class SE3 : public LieGroupBase<SE3<_Scalar>, _Scalar, 4, 6, 3> {
public:
    using Base = LieGroupBase<SE3<_Scalar>, _Scalar, 4, 6, 3>;
    using Scalar = typename Base::Scalar;
	using Vector4 = typename Base::Vector;

    using TangentVector = typename Base::TangentVector;
	using ActionVector = typename Base::ActionVector;
	using AlgebraVector = typename Base::AlgebraVector;

	using Matrix4 = typename Base::Matrix;
    using AdjointMatrix = typename Base::AdjointMatrix;
	using JacobianMatrix = typename Base::JacobianMatrix;
	using AlgebraMatrix = typename Base::AlgebraMatrix;


    using SO3Type = SO3<Scalar>;
    using Vector3 = typename SO3Type::Vector;
    using Matrix3 = typename SO3Type::Matrix;
    using Quaternion = typename SO3Type::Quaternion;

    // ========== Constructors ==========

    SE3() : m_rotation(), m_translation(Vector3::Zero()) {}
    SE3(const SO3Type &rotation, const Vector3 &translation)
        : m_rotation(rotation), m_translation(translation) {}
    explicit SE3(const Matrix4 &matrix)
        : m_rotation(matrix.template block<3, 3>(0, 0)),
          m_translation(matrix.template block<3, 1>(0, 3)) {}

    // ========== CRTP Implementations ==========

    static SE3 computeIdentity() noexcept { return SE3(); }

    SE3 computeInverse() const {
        const SO3Type inv_rot = m_rotation.inverse();
        return SE3(inv_rot, -(inv_rot.act(m_translation)));
    }

    SE3 compose(const SE3 &other) const {
        return SE3(m_rotation * other.m_rotation,
                   m_rotation.act(other.m_translation) + m_translation);
    }

    ActionVector computeAction(const ActionVector &point) const {
        return m_rotation.act(point) + m_translation;
    }

    static SE3 computeExp(const TangentVector &xi) {
        const Vector3 rho = xi.template head<3>();
        const Vector3 phi = xi.template tail<3>();

        const Scalar angle = phi.norm();
        const SO3Type R = SO3Type::exp(phi);
        Matrix3 V;

        if (angle < Types<Scalar>::epsilon()) {
            V = Matrix3::Identity() + Scalar(0.5) * SO3Type::hat(phi);
        } else {
            const Vector3 axis = phi / angle;
            const Matrix3 K = SO3Type::hat(axis);
            const Scalar sin_angle = std::sin(angle);
            const Scalar cos_angle = std::cos(angle);
            V = Matrix3::Identity() + (Scalar(1) - cos_angle) / angle * K +
                (angle - sin_angle) / (angle*angle) * K * K;
        }
        return SE3(R, V * rho);
    }

    TangentVector computeLog() const {
        const Vector3 phi = m_rotation.log();
        const Scalar angle = phi.norm();
        Matrix3 V_inv;

        if (angle < Types<Scalar>::epsilon()) {
            V_inv = Matrix3::Identity() - Scalar(0.5) * SO3Type::hat(phi);
        } else {
            const Vector3 axis = phi / angle;
            const Matrix3 K = SO3Type::hat(axis);
            const Scalar sin_angle = std::sin(angle);
            const Scalar cos_angle = std::cos(angle);
            V_inv = Matrix3::Identity() - Scalar(0.5) * K +
                  (Scalar(2) * sin_angle - angle * (Scalar(1) + cos_angle)) /
                  (Scalar(2) * angle * angle * sin_angle) * K * K;
        }

        const Vector3 rho = V_inv * m_translation;
        TangentVector result;
        result.template head<3>() = rho;
        result.template tail<3>() = phi;
        return result;
    }

    AdjointMatrix computeAdjoint() const {
        AdjointMatrix Ad = AdjointMatrix::Zero();
        const Matrix3 R = m_rotation.matrix();
        Ad.template block<3, 3>(0, 0) = R;
        Ad.template block<3, 3>(3, 3) = R;
        Ad.template block<3, 3>(0, 3) = SO3Type::hat(m_translation) * R;
        return Ad;
    }

    bool computeIsApprox(const SE3 &other, const Scalar &eps) const {
        return m_rotation.isApprox(other.m_rotation, eps) &&
               m_translation.isApprox(other.m_translation, eps);
    }

    // ========== New Integrated Functions ==========

    SE3 interpolate(const SE3 &other, const Scalar &t) const {
        return (*this) * SE3::exp(t * (this->inverse() * other).log());
    }

    static SE3 interpolateExponential(const TangentVector &xi1, const TangentVector &xi2, Scalar t) {
        return SE3::exp(xi1 + t * (xi2 - xi1));
    }

    Scalar distance(const SE3 &other, Scalar w_rot = Scalar(1), Scalar w_trans = Scalar(1)) const {
        const Scalar rot_dist = m_rotation.distance(other.m_rotation);
        const Scalar trans_dist = (m_translation - other.m_translation).norm();
        return w_rot * rot_dist + w_trans * trans_dist;
    }

    static std::vector<SE3> generateTrajectory(const std::vector<SE3> &waypoints, int num_points = 10) {
        if (waypoints.size() < 2) {
            return waypoints;
        }

        std::vector<SE3> trajectory;
        trajectory.reserve((waypoints.size() - 1) * num_points + 1);

        for (size_t i = 0; i < waypoints.size() - 1; ++i) {
            for (int j = 0; j < num_points; ++j) {
                const Scalar t = Scalar(j) / Scalar(num_points);
                trajectory.push_back(waypoints[i].interpolate(waypoints[i + 1], t));
            }
        }
        trajectory.push_back(waypoints.back());
        return trajectory;
    }

    // ========== Accessors ==========

    const SO3Type& rotation() const { return m_rotation; }
    SO3Type& rotation() { return m_rotation; }
    const Vector3& translation() const { return m_translation; }
    Vector3& translation() { return m_translation; }

    Matrix4 matrix() const {
        Matrix4 T = Matrix4::Identity();
        T.template block<3, 3>(0, 0) = m_rotation.matrix();
        T.template block<3, 1>(0, 3) = m_translation;
        return T;
    }

private:
    SO3Type m_rotation;
    Vector3 m_translation;
};

// ========== Type Aliases ==========
using SE3f = SE3<float>;
using SE3d = SE3<double>;

} // namespace sofa::component::cosserat::liegroups
