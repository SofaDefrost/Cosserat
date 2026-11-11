// This file defines the SGal3 (Special Galilean group in 3D) class, which
// represents Galilean transformations.

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

#pragma once

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <random>
#include "LieGroup.h" // Then the base class interface
#include "LieGroup.inl" // Then the base class interface
#include "SE3.h" // Then other dependencies
#include "SO2.h" // Then other dependencies
#include "Types.h" // Then our type system

namespace sofa::component::cosserat::liegroups {
/**
 * @brief Implementation of SGal(3), the Special Galilean group in 3D
 *
 * This class implements the group of Galilean transformations in 3D space.
 * Elements of SGal(3) are represented as a combination of:
 * - An SE(3) transformation (rotation R and position p)
 * - A 3D velocity vector v
 * - A time parameter t
 * and, as a matrix group (Kelly 2023), X = [[R, v, p], [0, 1, t], [0, 0, 1]], so that
 *   compose:  (R1,v1,p1,t1)(R2,v2,p2,t2) = (R1 R2, R1 v2 + v1, R1 p2 + p1 + v1 t2, t1 + t2)
 *   inverse:  (R^T, -R^T v, -R^T (p - t v), -t)
 *
 * The Lie algebra sgal(3) consists of vectors xi = [rho ; w ; beta ; tau] in ℝ¹⁰, where:
 * - First three components (rho) represent linear velocity (position rate)
 * - Next three components (w) represent angular velocity
 * - Next three components (beta) represent the boost (velocity rate)
 * - Last component (tau) represents time rate
 * hat(xi) = [[w^, beta, rho], [0, 0, tau], [0, 0, 0]], and
 *   exp(xi) = (exp(w^), V beta, V rho + tau E beta, tau), E = sum_k [w]x^k/(k+2)!.
 *
 * The action vector is [point ; velocity ; time] ∈ ℝ⁷.
 *
 * @tparam Scalar_ The scalar type (must be a floating-point type)
 */

template<typename Scalar_, int Dim_ = 10>
class SGal3 : public LieGroup<SGal3<Scalar_, Dim_>, Scalar_, 10, 10, 7>
//,public LieGroupOperations<SGal3<Scalar_>> // the Utils may be needed here!
{

public:
    using Base = LieGroup<SGal3<Scalar_, Dim_>, Scalar_, 10, 10, 7>;
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
     * @brief Default constructor creates an identity element.
     * Initializes pose to identity, velocity to zero vector, and time to zero.
     */
    SGal3() : m_pose(), m_velocity(Vector3::Zero()), m_time(0) {}

    /**
     * @brief Construct from pose, velocity, and time.
     * @param pose The SE3 pose component.
     * @param velocity The 3D linear velocity vector.
     * @param time The time parameter.
     */
    SGal3(const SE3<Scalar> &pose, const Vector3 &velocity, const Scalar &time) :
            m_pose(pose), m_velocity(velocity), m_time(time) {}

    /**
     * @brief Construct from rotation, position, velocity, and time.
     * @param rotation The SO3 rotation component.
     * @param position The 3D position vector.
     * @param velocity The 3D linear velocity vector.
     * @param time The time parameter.
     */
    SGal3(const SO3<Scalar> &rotation, const Vector3 &position, const Vector3 &velocity, const Scalar &time) :
            m_pose(rotation, position), m_velocity(velocity), m_time(time) {}

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
     * @brief Access the time component (const version).
     * @return A const reference to the time parameter.
     */
    const Scalar &time() const { return m_time; }
    /**
     * @brief Access the time component (mutable version).
     * @return A mutable reference to the time parameter.
     */
    Scalar &time() { return m_time; }

    /**
     * @brief The 5x5 matrix [[R, v, p], [0, 1, t], [0, 0, 1]] representing this element.
     * (Previously a 6x6 matrix with the time at (4, 5), which is not a representation of the group
     * law implemented by compose(); the class could not even be instantiated before.)
     */
    Eigen::Matrix<Scalar, 5, 5> extendedMatrix() const {
        Eigen::Matrix<Scalar, 5, 5> T = Eigen::Matrix<Scalar, 5, 5>::Identity();
        T.template block<3, 3>(0, 0) = m_pose.rotation().matrix();
        T.template block<3, 1>(0, 3) = m_velocity;
        T.template block<3, 1>(0, 4) = m_pose.translation();
        T(3, 4) = m_time;
        return T;
    }

    // ---- CRTP hooks of LieGroup (public: the base class and the callers use them) ----

    static SGal3<Scalar> computeIdentity() noexcept { return SGal3(); }

    /// (R1,v1,p1,t1)(R2,v2,p2,t2) = (R1 R2, R1 v2 + v1, R1 p2 + p1 + v1 t2, t1 + t2).
    SGal3<Scalar> computeCompose(const SGal3 &other) const {
        const Matrix3 R = m_pose.rotation().matrix();
        return SGal3(SE3<Scalar>(m_pose.rotation() * other.m_pose.rotation(),
                                 R * other.m_pose.translation() + m_pose.translation() + m_velocity * other.m_time),
                     R * other.m_velocity + m_velocity, m_time + other.m_time);
    }

    SGal3<Scalar> computeInverse() const {
        const SO3<Scalar> Rinv = m_pose.rotation().inverse();
        return SGal3(SE3<Scalar>(Rinv, -(Rinv.act(m_pose.translation() - m_time * m_velocity))),
                     -(Rinv.act(m_velocity)), -m_time);
    }

    /// exp([rho ; w ; beta ; tau]) = (exp(w^), V beta, V rho + tau E beta, tau).
    static SGal3<Scalar> computeExp(const TangentVector &algebra_element) {
        const Vector3 rho = algebra_element.template segment<3>(0); // Linear velocity
        const Vector3 w = algebra_element.template segment<3>(3);   // Angular velocity
        const Vector3 beta = algebra_element.template segment<3>(6); // Boost
        const Scalar tau = algebra_element[9];                       // Time rate

        Scalar A, B, C, c1, c2, c3;
        detail::so3JacobianCoefficients(w.norm(), A, B, C);
        detail::se3QCoefficients(w.norm(), c1, c2, c3);
        const Matrix3 W = SO3<Scalar>::computeHat(w);
        const Matrix3 W2 = W * W;
        const Matrix3 V = Matrix3::Identity() + A * W + B * W2;
        const Matrix3 E = Scalar(0.5) * Matrix3::Identity() + c1 * W + c2 * W2; // sum W^k/(k+2)!
        return SGal3(SE3<Scalar>(SO3<Scalar>::exp(w), V * rho + tau * (E * beta)), V * beta, tau);
    }

    TangentVector computeLog() const {
        const Vector3 w = m_pose.rotation().log();
        Scalar A, B, C, c1, c2, c3;
        detail::so3JacobianCoefficients(w.norm(), A, B, C);
        detail::se3QCoefficients(w.norm(), c1, c2, c3);
        const Matrix3 W = SO3<Scalar>::computeHat(w);
        const Matrix3 W2 = W * W;
        const Matrix3 V_inv = Matrix3::Identity() - Scalar(0.5) * W + C * W2;
        const Matrix3 E = Scalar(0.5) * Matrix3::Identity() + c1 * W + c2 * W2;
        const Vector3 beta = V_inv * m_velocity;
        const Vector3 rho = V_inv * (m_pose.translation() - m_time * (E * beta));
        TangentVector result;
        result << rho, w, beta, m_time;
        return result;
    }

    /// Ad_g xi = vee(g hat(xi) g^-1), order [rho ; w ; beta ; tau]:
    ///   rho' = R rho + (p - t v)x R w - t R beta + v tau,  w' = R w,  beta' = R beta + vx R w,  tau' = tau.
    AdjointMatrix computeAdjoint() const {
        AdjointMatrix Ad = AdjointMatrix::Zero();
        const Matrix3 R = m_pose.rotation().matrix();
        Ad.template block<3, 3>(0, 0) = R;
        Ad.template block<3, 3>(0, 3) =
                SO3<Scalar>::computeHat(Vector3(m_pose.translation() - m_time * m_velocity)) * R;
        Ad.template block<3, 3>(0, 6) = -m_time * R;
        Ad.template block<3, 1>(0, 9) = m_velocity;
        Ad.template block<3, 3>(3, 3) = R;
        Ad.template block<3, 3>(6, 3) = SO3<Scalar>::computeHat(m_velocity) * R;
        Ad.template block<3, 3>(6, 6) = R;
        Ad(9, 9) = Scalar(1);
        return Ad;
    }

    bool computeIsApprox(const SGal3 &other, const Scalar &eps = Types<Scalar>::tolerance()) const {
        const Scalar scale = std::max({Scalar(1), m_velocity.norm(), other.m_velocity.norm()});
        return m_pose.computeIsApprox(other.m_pose, eps) && (m_velocity - other.m_velocity).norm() <= eps * scale &&
               std::abs(m_time - other.m_time) <= eps * std::max(Scalar(1), std::abs(m_time));
    }

    /// Action on [point ; velocity ; time]: p' = R p + p_g + v_g t, v' = R v + v_g, t' = t + t_g.
    typename Base::ActionVector computeAction(const typename Base::ActionVector &point_vel_time) const {
        const Vector3 point = point_vel_time.template head<3>();
        const Vector3 vel = point_vel_time.template segment<3>(3);
        const Scalar t = point_vel_time[6];
        typename Base::ActionVector result;
        result << m_pose.computeAction(point) + m_velocity * t, m_pose.rotation().computeAction(vel) + m_velocity,
                t + m_time;
        return result;
    }

    /// hat(xi) = [[w^, beta, rho], [0, 0, tau], [0, 0, 0]] (5x5) in the top-left corner of the 10x10 AlgebraMatrix.
    static Matrix computeHat(const TangentVector &xi) {
        Matrix result = Matrix::Zero();
        result.template block<3, 3>(0, 0) = SO3<Scalar>::computeHat(Vector3(xi.template segment<3>(3)));
        result.template block<3, 1>(0, 3) = xi.template segment<3>(6);
        result.template block<3, 1>(0, 4) = xi.template segment<3>(0);
        result(3, 4) = xi[9];
        return result;
    }

    static TangentVector computeVee(const Matrix &X) {
        TangentVector xi;
        xi.template segment<3>(0) = X.template block<3, 1>(0, 4);
        xi.template segment<3>(3) = SO3<Scalar>::computeVee(Matrix3(X.template block<3, 3>(0, 0)));
        xi.template segment<3>(6) = X.template block<3, 1>(0, 3);
        xi[9] = X(3, 4);
        return xi;
    }

    /// ad_xi y = [xi, y]:  rows [rho ; w ; beta ; tau]
    ///   [[w^, rho^, -tau I, beta], [0, w^, 0, 0], [0, beta^, w^, 0], [0, 0, 0, 0]]
    static AdjointMatrix computeAd(const TangentVector &xi) {
        AdjointMatrix ad = AdjointMatrix::Zero();
        const Vector3 rho = xi.template segment<3>(0), beta = xi.template segment<3>(6);
        const Matrix3 W = SO3<Scalar>::computeHat(Vector3(xi.template segment<3>(3)));
        ad.template block<3, 3>(0, 0) = W;
        ad.template block<3, 3>(0, 3) = SO3<Scalar>::computeHat(rho);
        ad.template block<3, 3>(0, 6) = -xi[9] * Matrix3::Identity();
        ad.template block<3, 1>(0, 9) = beta;
        ad.template block<3, 3>(3, 3) = W;
        ad.template block<3, 3>(6, 3) = SO3<Scalar>::computeHat(beta);
        ad.template block<3, 3>(6, 6) = W;
        return ad;
    }

    template<typename Generator>
    static SGal3<Scalar> computeRandom(Generator &gen) {
        return SGal3(SE3<Scalar>::computeRandom(gen), Types<Scalar>::template randomVector<3>(gen),
                     Types<Scalar>::randomScalar(gen));
    }

    std::ostream &print(std::ostream &os) const {
        os << "SGal3(pose=" << m_pose << ", velocity=" << m_velocity.transpose() << ", time=" << m_time << ")";
        return os;
    }

    static constexpr std::string_view getTypeName() { return "SGal3"; }

    bool computeIsValid() const { return m_pose.computeIsValid() && m_velocity.allFinite() && std::isfinite(m_time); }

    void computeNormalize() { m_pose.computeNormalize(); }

    /**
     * @brief Computes the interpolation between two SGal3 elements.
     * @param other The target SGal3 element.
     * @param t The interpolation parameter (between 0 and 1).
     * @return The interpolated SGal3 element.
     */
    SGal3<Scalar> computeInterpolate(const SGal3<Scalar> &other, const Scalar &t) const {
        const TangentVector delta = (this->computeInverse().computeCompose(other)).computeLog();
        return this->computeCompose(SGal3<Scalar>::computeExp(t * delta));
    }

private:
    SE3<Scalar> m_pose;   ///< Rigid body transformation (R, p)
    Vector3 m_velocity;   ///< Linear velocity v
    Scalar m_time;        ///< Time parameter t
};

using SGal3d = SGal3<double>;
using SGal3f = SGal3<float>;

} // namespace sofa::component::cosserat::liegroups

// Out-of-line free functions (operator/, Euler conversions, interpolate).
#include "SGal3.inl"

