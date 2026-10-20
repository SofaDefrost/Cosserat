// Spherical (3 DoF) and two-axis (2 DoF) joints between two rigid bodies, SE3 product of exponentials.
//
// Same convention as RevoluteJoint.h: everything is expressed in the parent frame at q = 0,
// twists are [phi (angular) ; rho (linear)], and the child pose in the parent frame is
//
//     T_rel(q) = Exp(q_1 xi_1) ... Exp(q_n xi_n) * T_rest
//
// The five anatomical joint types map onto three classes:
//   - enarthrosis (ball-and-socket, 3 DoF)  -> SphericalJoint
//   - condylar joint (2 concurrent axes)    -> TwoAxisJoint with point1 == point2
//   - saddle joint (2 perpendicular, non concurrent axes) -> TwoAxisJoint with offset points
//   - ginglymus (hinge, 1 DoF)              -> RevoluteJoint
//   - trochoid (pivot, 1 DoF)               -> RevoluteJoint with the axis along the long axis of the bones
//     (mathematically identical to a hinge: only the axis direction differs)

#pragma once

#include <stdexcept>
#include "RevoluteJoint.h"

namespace sofa::component::cosserat::liegroups {

	/**
	 * @brief Ball-and-socket joint: child rotates about a fixed point with rotation vector w = (rx, ry, rz).
	 *
	 * T_rel(w) = Exp([w ; -w × pivot]) · T_rest, i.e. the rotation exp(w) about the fixed point @p pivot.
	 * q is a rotation vector, so the parametrisation is singular only at |w| = 2π (log branch limit π for
	 * jointAngles()).
	 */
	template<typename _Scalar>
	class SphericalJoint {
	public:
		using Scalar = _Scalar;
		using SE3Type = SE3<Scalar>;
		using SO3Type = SO3<Scalar>;
		using Vector3 = typename SE3Type::Vector3;
		using TangentVector = typename SE3Type::TangentVector;

		static constexpr int DoF = 3;
		int dof() const { return DoF; }

		SphericalJoint(const Vector3 &pivot, const SE3Type &restPose = SE3Type()) : m_pivot(pivot), m_restPose(restPose) {}

		const Vector3 &pivot() const { return m_pivot; }
		const SE3Type &restPose() const { return m_restPose; }

		/// Pose of the child frame in the parent frame for rotation vector @p w.
		SE3Type relativePose(const Vector3 &w) const {
			TangentVector xi;
			xi.template head<3>() = w;
			xi.template tail<3>() = -w.cross(m_pivot);
			return SE3Type::computeExp(xi) * m_restPose;
		}

		/// World pose of the child.
		SE3Type childPose(const SE3Type &parentPose, const Vector3 &w) const { return parentPose * relativePose(w); }

		/// Rotation vector recovered from a relative pose (rotation part only), norm in [0, π].
		Vector3 jointAngles(const SE3Type &relative) const {
			return (relative.rotation() * m_restPose.rotation().inverse()).log();
		}

		/// 6-vector measuring how far @p relative is from being reachable (zero iff reachable).
		TangentVector constraintResidual(const SE3Type &relative) const {
			return (relativePose(jointAngles(relative)).inverse() * relative).log();
		}

	private:
		Vector3 m_pivot;
		SE3Type m_restPose;
	};

	/**
	 * @brief Two successive hinges: T_rel(q1, q2) = Exp(q1 ξ1) · Exp(q2 ξ2) · T_rest.
	 *
	 * Both screws are expressed in the parent frame at q = 0 (product-of-exponentials convention), so
	 * axis 2 does NOT follow the rotation about axis 1 in the screw definition; the second hinge is
	 * carried by the first one through the composition order.
	 *  - point1 == point2 : condylar joint (concurrent axes, no axial rotation),
	 *  - offset points, perpendicular axes : saddle joint.
	 */
	template<typename _Scalar>
	class TwoAxisJoint {
	public:
		using Scalar = _Scalar;
		using SE3Type = SE3<Scalar>;
		using Vector3 = typename SE3Type::Vector3;
		using Vector2 = Eigen::Matrix<Scalar, 2, 1>;
		using TangentVector = typename SE3Type::TangentVector;

		static constexpr int DoF = 2;
		int dof() const { return DoF; }

		TwoAxisJoint(const Vector3 &axis1, const Vector3 &point1, const Vector3 &axis2, const Vector3 &point2,
					 const SE3Type &restPose = SE3Type()) :
			m_axis1(axis1.normalized()), m_point1(point1), m_axis2(axis2.normalized()), m_point2(point2),
			m_restPose(restPose) {
			if (!(axis1.norm() > Scalar(0)) || !(axis2.norm() > Scalar(0)))
				throw std::invalid_argument("TwoAxisJoint: axes must be non-zero");
			m_screw1 << m_axis1, -m_axis1.cross(m_point1);
			m_screw2 << m_axis2, -m_axis2.cross(m_point2);
		}

		const SE3Type &restPose() const { return m_restPose; }
		const TangentVector &screw1() const { return m_screw1; }
		const TangentVector &screw2() const { return m_screw2; }

		/// Pose of the child frame in the parent frame.
		SE3Type relativePose(const Vector2 &q) const {
			return SE3Type::computeExp(q[0] * m_screw1) * SE3Type::computeExp(q[1] * m_screw2) * m_restPose;
		}

		/// World pose of the child.
		SE3Type childPose(const SE3Type &parentPose, const Vector2 &q) const { return parentPose * relativePose(q); }

	private:
		Vector3 m_axis1, m_point1, m_axis2, m_point2;
		SE3Type m_restPose;
		TangentVector m_screw1, m_screw2;
	};

} // namespace sofa::component::cosserat::liegroups
