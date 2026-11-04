// Serial kinematic chain on SE(3), product-of-exponentials form.
//
//     T_i(q_i) = Exp(q_i xi_i) * T_rest_i        (xi_i unit screw expressed in the frame that precedes segment i)
//     T_tip(q) = T_0(q_0) * T_1(q_1) * ... * T_{n-1}(q_{n-1})
//
// Joint types: revolute (xi = [axis ; -axis x pivot]), prismatic (xi = [0 ; axis]) and fixed (no DoF).
// Twists follow the library layout [phi (angular) ; rho (linear)].
//
// Three Jacobians of the tip are provided (all 6 x nDof, columns ordered by joint):
//   - spatial : twist of the tip body seen from the base, V^s = qdot T^-1   (ang ; lin velocity of the base origin)
//   - body    : twist expressed in the tip frame, V^b = T^-1 qdot           (= Ad_{T^-1} J^s)
//   - hybrid  : angular velocity in the base frame and linear velocity of the tip origin in the base frame
//               (the usual "robotics" Jacobian that maps qdot to [omega ; pdot]).

#pragma once

#include <cmath>
#include <stdexcept>
#include <vector>

#include "SE3.h"

namespace sofa::component::cosserat::liegroups {

	template<typename _Scalar>
	class KinematicChain {
	public:
		using Scalar = _Scalar;
		using SE3Type = SE3<Scalar>;
		using Vector3 = typename SE3Type::Vector3;
		using TangentVector = typename SE3Type::TangentVector;
		using AdjointMatrix = typename SE3Type::AdjointMatrix;
		using JacobianMatrix = Eigen::Matrix<Scalar, 6, Eigen::Dynamic>;
		using VectorX = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

		enum class JointType { Fixed, Revolute, Prismatic };

		struct Segment {
			JointType type = JointType::Fixed;
			TangentVector screw = TangentVector::Zero(); ///< unit screw in the parent frame (zero when fixed)
			SE3Type rest; ///< pose of the segment's tip frame in its parent frame at q = 0
			int dofIndex = -1; ///< column index in q, -1 for a fixed segment
		};

		/// Revolute segment about the line (pivot, axis), both given in the parent frame.
		void addRevolute(const Vector3 &axis, const Vector3 &pivot, const SE3Type &rest = SE3Type()) {
			if (!(axis.norm() > Scalar(0))) throw std::invalid_argument("KinematicChain: null revolute axis");
			const Vector3 a = axis.normalized();
			Segment s;
			s.type = JointType::Revolute;
			s.screw << a, -a.cross(pivot);
			s.rest = rest;
			push(s);
		}

		/// Prismatic segment translating along @p axis (parent frame).
		void addPrismatic(const Vector3 &axis, const SE3Type &rest = SE3Type()) {
			if (!(axis.norm() > Scalar(0))) throw std::invalid_argument("KinematicChain: null prismatic axis");
			Segment s;
			s.type = JointType::Prismatic;
			s.screw << Vector3::Zero(), axis.normalized();
			s.rest = rest;
			push(s);
		}

		/// Rigid link: just a constant transform.
		void addFixed(const SE3Type &rest) {
			Segment s;
			s.type = JointType::Fixed;
			s.rest = rest;
			push(s);
		}

		int nDof() const { return m_nDof; }
		int nSegments() const { return static_cast<int>(m_segments.size()); }
		const std::vector<Segment> &segments() const { return m_segments; }

		/// Pose of segment i's tip frame in the base frame (i = nSegments()-1 gives the end effector).
		SE3Type segmentPose(const VectorX &q, int i) const {
			check(q);
			SE3Type T;
			for (int k = 0; k <= i; ++k) T = T * localPose(m_segments[k], q);
			return T;
		}

		SE3Type forward(const VectorX &q) const { return segmentPose(q, nSegments() - 1); }

		/// Spatial Jacobian: column of joint i is Ad_{prefix_i} xi_i (axis of joint i seen from the base).
		JacobianMatrix spatialJacobian(const VectorX &q) const {
			check(q);
			JacobianMatrix J = JacobianMatrix::Zero(6, m_nDof);
			SE3Type prefix; // product of the segments strictly before the current one
			for (const Segment &s : m_segments) {
				if (s.dofIndex >= 0) J.col(s.dofIndex) = prefix.computeAdjoint() * s.screw;
				prefix = prefix * localPose(s, q);
			}
			return J;
		}

		/// Body Jacobian: Ad_{T^-1} J^s.
		JacobianMatrix bodyJacobian(const VectorX &q) const {
			return forward(q).inverse().computeAdjoint() * spatialJacobian(q);
		}

		/// Hybrid Jacobian: [omega in base ; velocity of the tip origin in base].
		JacobianMatrix hybridJacobian(const VectorX &q) const {
			const JacobianMatrix Js = spatialJacobian(q);
			const Vector3 p = forward(q).translation();
			JacobianMatrix J = Js;
			for (int c = 0; c < m_nDof; ++c) J.col(c).template tail<3>() = Js.col(c).template tail<3>() + Js.col(c).template head<3>().cross(p);
			return J;
		}

	private:
		void push(Segment s) {
			if (s.type != JointType::Fixed) s.dofIndex = m_nDof++;
			m_segments.push_back(s);
		}

		void check(const VectorX &q) const {
			if (q.size() != m_nDof) throw std::invalid_argument("KinematicChain: q has the wrong size");
		}

		static SE3Type localPose(const Segment &s, const VectorX &q) {
			if (s.dofIndex < 0) return s.rest;
			return SE3Type::computeExp(q[s.dofIndex] * s.screw) * s.rest;
		}

		std::vector<Segment> m_segments;
		int m_nDof = 0;
	};

} // namespace sofa::component::cosserat::liegroups
