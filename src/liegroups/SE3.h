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
	SE3(const Quaternion &quat, const Vector3 &translation)
		: m_rotation(quat.toRotationMatrix()), m_translation(translation) {}
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

	/**
   * @brief Exponential map from se(3) to SE(3) using Cosserat-style approach.
   *
   * For ξ = [ρ, φ]ᵀ ∈ se(3) where ρ ∈ ℝ³ (translation) and φ ∈ ℝ³ (rotation):
   *
   * Small case (‖φ‖ ≈ 0):
   *   T = I₄ + s·ξ̂
   *
   * General case:
   *   T = I₄ + s·ξ̂ + α·ξ̂² + β·ξ̂³
   *   where α = (1-cos(s‖φ‖))/‖φ‖², β = (s‖φ‖-sin(s‖φ‖))/‖φ‖³
   *
   * @param strain 6D strain vector [ρ, φ]ᵀ (linear and angular strain rates).
   * @param length Arc length parameter for integration.
   * @return The corresponding SE3 element.
   */
	static SE3 expCosserat(const TangentVector& strain, const Scalar& length) noexcept {
    	// Extract translation and rotation parts
    	const Vector3 rho = strain.template tail<3>();      // Linear strain (translation rate)
    	const Vector3 phi = strain.template head<3>();      // Angular strain (rotation rate)

    	const Scalar phi_norm = phi.norm();

    	// Handle small rotation case
    	if (phi_norm <= std::numeric_limits<Scalar>::epsilon()) {
    		return expCosseratSmall(rho, phi, length);
    	}

    	return expCosseratGeneral(rho, phi, phi_norm, length);
    }

    static SE3 computeExp(const TangentVector &xi) {
        const Vector3 rho = xi.template tail<3>();
        const Vector3 phi = xi.template head<3>();
    	std::cout << "SE3::computeExp rho: " << rho << ", phi: " << phi << std::endl;

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

	/**
  * @brief Extract position and orientation as separate components.
  * @param position Output position vector.
  * @param orientation Output quaternion.
  */
	void getPositionAndOrientation(Vector3& position, Quaternion& orientation) const {
    	position = m_translation;
    	orientation = m_rotation.quaternion();
    }

    /**
     * @brief Print the SE3 element to an output stream.
     * This method is required by the LieGroupBase CRTP interface.
     * @param os Output stream to write to.
     * @return Reference to the output stream.
     */
    std::ostream &print(std::ostream &os) const {
        os << "SE3(R=" << m_rotation << ", t=[" << m_translation[0] << ", " 
           << m_translation[1] << ", " << m_translation[2] << "])";
        return os;
    }

private:
    SO3Type m_rotation;
    Vector3 m_translation;
	/**
	 * @brief Small angle approximation for SE(3) exponential.
	 */
	static SE3 expCosseratSmall(const Vector3& rho, const Vector3& phi, const Scalar& length) noexcept {
		// For small rotations: T ≈ I + length * ξ̂
		// Build the se(3) matrix ξ̂
		// ξ̂ = [φ×  ρ]  where φ× is the skew-symmetric matrix of φ
		Matrix4 xi_hat = buildXiHat(rho,phi);

		Matrix4 g_x = Matrix4::Identity() + length * xi_hat;

		// Small angle quaternion: q ≈ (1, 0.5 * rotation_vec)
		//const Vector3 half_rotation = rotation_vec * 0.5;
		//Quaternion rotation(1.0, half_rotation.x(), half_rotation.y(), half_rotation.z());
		//rotation.normalize();

		return SE3(g_x);
	}
	/**
	 * @brief General case for SE(3) exponential with Cosserat-style approach.
	 */
	/**
	   * @brief General case SE(3) exponential using 3rd-order Taylor expansion.
	   */
	static SE3 expCosseratGeneral(const Vector3& rho, const Vector3& phi,
								 const Scalar& phi_norm, const Scalar& length) noexcept {
		const Scalar s_phi_norm = length * phi_norm;
		const Scalar phi_norm2 = phi_norm * phi_norm;
		const Scalar phi_norm3 = phi_norm2 * phi_norm;

		// Cosserat coefficients
		const Scalar cos_s_phi = std::cos(s_phi_norm);
		const Scalar sin_s_phi = std::sin(s_phi_norm);

		const Scalar alpha = (1.0 - cos_s_phi) / phi_norm2;
		const Scalar beta = (s_phi_norm - sin_s_phi) / phi_norm3;

		// Build se(3) matrix ξ̂
		const Matrix4 xi_hat = buildXiHat(rho, phi);

		// Compute powers: ξ̂², ξ̂³
		const Matrix4 xi_hat2 = xi_hat * xi_hat;
		const Matrix4 xi_hat3 = xi_hat2 * xi_hat;

		// Taylor expansion: T = I + s·ξ̂ + α·ξ̂² + β·ξ̂³
		const Matrix4 g_x = Matrix4::Identity() +
						 length * xi_hat +
						 alpha * xi_hat2 +
						 beta * xi_hat3;

		return SE3(g_x);
	}

	/**
   * @brief Build the se(3) matrix representation ξ̂ ∈ ℝ⁴ˣ⁴.
   *
   * ξ̂ = [φ×  ρ]  where φ× is the skew-symmetric matrix of φ
   *     [0   0]
   *
   * @param rho Translation part (3D).
   * @param phi Rotation part (3D).
   * @return 4x4 se(3) matrix.
   */
	static Matrix4 buildXiHat(const Vector3& rho, const Vector3& phi) noexcept {
		Matrix4 xi_hat = Matrix4::Zero();

		// Top-left 3x3: skew-symmetric matrix [φ]×
		xi_hat(0, 1) = -phi.z();
		xi_hat(0, 2) =  phi.y();
		xi_hat(1, 0) =  phi.z();
		xi_hat(1, 2) = -phi.x();
		xi_hat(2, 0) = -phi.y();
		xi_hat(2, 1) =  phi.x();

		// Top-right 3x1: translation part
		xi_hat(0, 3) =  1.0 + rho.x();
		xi_hat(1, 3) =  rho.y();
		xi_hat(2, 3) =  rho.z();

		// The bottom row is already zero
		return xi_hat;
	}
};

// ========== Type Aliases ==========
using SE3f = SE3<float>;
using SE3d = SE3<double>;

} // namespace sofa::component::cosserat::liegroups
