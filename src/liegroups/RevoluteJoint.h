// Revolute (hinge) joint between two rigid bodies, built on SE(3).
//
// Forward kinematics by product of exponentials (screw theory):
//
//     T_child(q) = T_parent * Exp(q * xi) * T_rest
//
// where xi = [axis ; -axis x pivot] is the unit screw of the hinge (angular part
// first, same [phi, rho] layout as SE3::computeExp) expressed in the parent frame,
// and T_rest is the pose of the child frame in the parent frame at q = 0.
// The single relative degree of freedom is the joint angle q.

#pragma once

#include <cmath>
#include <stdexcept>
#include "SE3.h"

namespace sofa::component::cosserat::liegroups {

	/**
	 * @brief Hinge between a parent body and a child body, one DoF (the angle q).
	 *
	 * Everything is expressed in the parent body frame:
	 *  - @p axis   : direction of the hinge axis (normalised internally),
	 *  - @p pivot  : any point lying on the hinge axis,
	 *  - @p restPose : pose of the child frame in the parent frame when q = 0.
	 *
	 * Twists follow the library convention ξ = [φ (angular), ρ (linear)].
	 * Body twists are expressed in the frame of the body they belong to.
	 *
	 * @tparam _Scalar Floating-point scalar type.
	 */
	template<typename _Scalar>
	class RevoluteJoint {
	public:
		using Scalar = _Scalar;
		using SE3Type = SE3<Scalar>;
		using SO3Type = SO3<Scalar>;
		using Vector3 = typename SE3Type::Vector3;
		using TangentVector = typename SE3Type::TangentVector;

		static constexpr int DoF = 1;
		int dof() const { return DoF; }

		/**
		 * @param axis     Hinge axis direction in the parent frame (non-zero).
		 * @param pivot    A point on the axis, in the parent frame.
		 * @param restPose Child frame in the parent frame at q = 0 (default: identity).
		 */
		RevoluteJoint(const Vector3 &axis, const Vector3 &pivot, const SE3Type &restPose = SE3Type()) :
			m_axis(axis.normalized()), m_pivot(pivot), m_restPose(restPose) {
			if (!(axis.norm() > Scalar(0)))
				throw std::invalid_argument("RevoluteJoint: the hinge axis must be non-zero");
			m_screw.template head<3>() = m_axis;
			m_screw.template tail<3>() = -m_axis.cross(m_pivot); // velocity of the origin when rotating about the line
		}

		const Vector3 &axis() const { return m_axis; }
		const Vector3 &pivot() const { return m_pivot; }
		const SE3Type &restPose() const { return m_restPose; }

		/// Unit screw ξ = [axis ; -axis × pivot] of the hinge, in the parent frame.
		const TangentVector &screw() const { return m_screw; }

		/// Pose of the child frame in the parent frame for joint angle @p q: Exp(q ξ) · T_rest.
		SE3Type relativePose(Scalar q) const { return SE3Type::computeExp(q * m_screw) * m_restPose; }

		/// Forward kinematics: world pose of the child given the parent's world pose and q.
		SE3Type childPose(const SE3Type &parentPose, Scalar q) const { return parentPose * relativePose(q); }

		/**
		 * @brief Inverse step: joint angle from the relative pose (child in parent frame).
		 *
		 * The rotation of Exp(qξ) is exp(q·axis), so q = axis · log(R_rel · R_rest⁻¹).
		 * Result is in (-π, π] (SO3::log returns an angle in [0, π]).
		 * Only the rotation is used; see constraintResidual() for the rest.
		 */
		Scalar jointAngle(const SE3Type &relative) const {
			const SO3Type Rq = relative.rotation() * m_restPose.rotation().inverse();
			return m_axis.dot(Rq.log());
		}

		/// Same as jointAngle(), from the two world poses.
		Scalar jointAngle(const SE3Type &parentPose, const SE3Type &childPose) const {
			return jointAngle(parentPose.inverse() * childPose);
		}

		/**
		 * @brief 6-vector measuring how far a relative pose is from satisfying the hinge.
		 *
		 * Returns log( relativePose(q*)⁻¹ · relative ) with q* = jointAngle(relative),
		 * i.e. the child-frame se(3) error. It is zero iff the relative pose is
		 * reachable through the hinge (5 independent constraints, 1 free DoF).
		 */
		TangentVector constraintResidual(const SE3Type &relative) const {
			return (relativePose(jointAngle(relative)).inverse() * relative).log();
		}

		/**
		 * @brief Body twist of the child (in the child frame) per unit q̇: Ad_{T_rel(q)⁻¹} ξ.
		 *
		 * Follows from T_rel(q+h) = Exp(hξ) T_rel(q) ⇒ T_rel(q)⁻¹ T_rel(q+h) = Exp(h Ad_{T_rel⁻¹} ξ).
		 */
		TangentVector jointTwistInChild(Scalar q) const {
			return relativePose(q).inverse().adjoint() * m_screw;
		}

		/**
		 * @brief Velocity propagation: body twist of the child from the parent's body twist.
		 *
		 * V_child = Ad_{T_rel(q)⁻¹} · V_parent + q̇ · jointTwistInChild(q)
		 *
		 * @param parentBodyTwist Parent twist expressed in the parent frame.
		 */
		TangentVector childBodyTwist(const TangentVector &parentBodyTwist, Scalar q, Scalar qdot) const {
			return relativePose(q).inverse().adjoint() * parentBodyTwist + qdot * jointTwistInChild(q);
		}

	private:
		Vector3 m_axis;
		Vector3 m_pivot;
		SE3Type m_restPose;
		TangentVector m_screw;
	};

	/**
	 * @brief Two rigid bodies connected by a revolute joint: the minimal articulated system.
	 *
	 * The parent pose (SE3, 6 DoF) and the angle q (1 DoF) are the state; the child
	 * pose is derived by group composition, so the two bodies stay exactly articulated.
	 */
	template<typename _Scalar>
	class ArticulatedRigidPair {
	public:
		using Scalar = _Scalar;
		using SE3Type = SE3<Scalar>;
		using JointType = RevoluteJoint<Scalar>;
		using TangentVector = typename SE3Type::TangentVector;
		using Vector3 = typename SE3Type::Vector3;

		explicit ArticulatedRigidPair(const JointType &joint, const SE3Type &parentPose = SE3Type(), Scalar q = Scalar(0)) :
			m_joint(joint), m_parentPose(parentPose), m_q(q) {}

		const JointType &joint() const { return m_joint; }
		const SE3Type &parentPose() const { return m_parentPose; }
		void setParentPose(const SE3Type &pose) { m_parentPose = pose; }
		Scalar angle() const { return m_q; }
		void setAngle(Scalar q) { m_q = q; }

		/// World pose of the child body (forward kinematics).
		SE3Type childPose() const { return m_joint.childPose(m_parentPose, m_q); }

		/// World position of the hinge pivot (same point seen from both bodies).
		Vector3 pivotWorld() const { return m_parentPose.act(m_joint.pivot()); }

		/// World direction of the hinge axis.
		Vector3 axisWorld() const { return m_parentPose.rotation().act(m_joint.axis()); }

	private:
		JointType m_joint;
		SE3Type m_parentPose;
		Scalar m_q;
	};

} // namespace sofa::component::cosserat::liegroups
