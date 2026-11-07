// Kinematic tree on SE(3), product-of-exponentials form (a serial chain is the special case of a path).
//
//     T_i(q_i)   = Exp(q_i xi_i) * T_rest_i          (xi_i unit screw expressed in the frame of the parent of segment i)
//     T^0_i(q)   = T^0_parent(i) * T_i(q_i)           (pose of segment i in the base frame, T^0_base = identity)
//
// Joint types: revolute (xi = [axis ; -axis x pivot]), prismatic (xi = [0 ; axis]) and fixed (no DoF).
// Twists follow the library layout [phi (angular) ; rho (linear)].
//
// Every segment names its parent (Base, or any segment added earlier); by default a segment hangs from the previous
// one, which gives a serial chain. The Jacobians of a segment are 6 x nDof, columns ordered by joint, and are zero on
// the joints that are not ancestors of that segment:
//   - spatial : twist of the segment seen from the base, V^s = Tdot T^-1   (ang ; lin velocity of the base origin)
//   - body    : twist expressed in the segment frame, V^b = T^-1 Tdot      (= Ad_{T^-1} J^s)
//   - hybrid  : angular velocity in the base frame and linear velocity of the segment origin in the base frame
//               (the usual "robotics" Jacobian, and the one SOFA's Rigid3d velocities correspond to).

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

		static constexpr int Base = -1; ///< parent index meaning "the base frame"
		static constexpr int Previous = -2; ///< parent index meaning "the segment added just before" (serial chain)

		struct Segment {
			JointType type = JointType::Fixed;
			TangentVector screw = TangentVector::Zero(); ///< unit screw in the parent frame (zero when fixed)
			SE3Type rest; ///< pose of the segment's frame in its parent frame at q = 0
			int parent = Base; ///< index of the parent segment, or Base
			int dofIndex = -1; ///< column index in q, -1 for a fixed segment
		};

		/// Revolute segment about the line (pivot, axis), both given in the parent frame.
		int addRevolute(const Vector3 &axis, const Vector3 &pivot, const SE3Type &rest = SE3Type(),
						int parent = Previous) {
			if (!(axis.norm() > Scalar(0))) throw std::invalid_argument("KinematicChain: null revolute axis");
			const Vector3 a = axis.normalized();
			Segment s;
			s.type = JointType::Revolute;
			s.screw << a, -a.cross(pivot);
			s.rest = rest;
			return push(s, parent);
		}

		/// Prismatic segment translating along @p axis (parent frame).
		int addPrismatic(const Vector3 &axis, const SE3Type &rest = SE3Type(), int parent = Previous) {
			if (!(axis.norm() > Scalar(0))) throw std::invalid_argument("KinematicChain: null prismatic axis");
			Segment s;
			s.type = JointType::Prismatic;
			s.screw << Vector3::Zero(), axis.normalized();
			s.rest = rest;
			return push(s, parent);
		}

		/// Rigid link: just a constant transform.
		int addFixed(const SE3Type &rest, int parent = Previous) {
			Segment s;
			s.type = JointType::Fixed;
			s.rest = rest;
			return push(s, parent);
		}

		int nDof() const { return m_nDof; }
		int nSegments() const { return static_cast<int>(m_segments.size()); }
		const std::vector<Segment> &segments() const { return m_segments; }

		/// Joint coordinates (column indices) that move segment i, from the base towards i.
		const std::vector<int> &ancestorDofs(int i) const { return m_ancestorDofs.at(i); }

		/// Poses of every segment in the base frame, in one pass (parents are always listed before their children).
		std::vector<SE3Type> allPoses(const VectorX &q) const {
			check(q);
			std::vector<SE3Type> poses(m_segments.size());
			for (std::size_t k = 0; k < m_segments.size(); ++k) {
				const Segment &s = m_segments[k];
				const SE3Type local = localPose(s, q);
				poses[k] = (s.parent == Base) ? local : poses[s.parent] * local;
			}
			return poses;
		}

		/// Pose of segment i in the base frame (i = -1: the last segment).
		SE3Type segmentPose(const VectorX &q, int i) const { return allPoses(q).at(resolve(i)); }

		/// Pose of the last segment (the end effector of a serial chain).
		SE3Type forward(const VectorX &q) const { return segmentPose(q, -1); }

		/**
		 * @brief Unit screw of every joint seen from the base frame, as columns (6 x nDof): S_k = Ad_{T^0_parent} xi_k.
		 * @param poses result of allPoses(q) (avoids recomputing it).
		 */
		JacobianMatrix jointScrews(const std::vector<SE3Type> &poses) const {
			JacobianMatrix S = JacobianMatrix::Zero(6, m_nDof);
			for (const Segment &s : m_segments) {
				if (s.dofIndex < 0) continue;
				S.col(s.dofIndex) = (s.parent == Base) ? s.screw : TangentVector(poses[s.parent].computeAdjoint() * s.screw);
			}
			return S;
		}

		/// Spatial Jacobian of segment i (default: the last one): the screws of its ancestor joints, others zero.
		JacobianMatrix spatialJacobian(const VectorX &q, int i = -1) const {
			const int seg = resolve(i);
			const JacobianMatrix S = jointScrews(allPoses(q));
			JacobianMatrix J = JacobianMatrix::Zero(6, m_nDof);
			for (int k : m_ancestorDofs[seg]) J.col(k) = S.col(k);
			return J;
		}

		/// Body Jacobian of segment i: Ad_{T^-1} J^s.
		JacobianMatrix bodyJacobian(const VectorX &q, int i = -1) const {
			return segmentPose(q, i).inverse().computeAdjoint() * spatialJacobian(q, i);
		}

		/// Hybrid Jacobian of segment i: [omega in base ; velocity of the segment origin in base].
		JacobianMatrix hybridJacobian(const VectorX &q, int i = -1) const {
			const JacobianMatrix Js = spatialJacobian(q, i);
			const Vector3 p = segmentPose(q, i).translation();
			JacobianMatrix J = Js;
			for (int c = 0; c < m_nDof; ++c)
				J.col(c).template tail<3>() = Js.col(c).template tail<3>() + Js.col(c).template head<3>().cross(p);
			return J;
		}

	private:
		int push(Segment s, int parent) {
			const int index = static_cast<int>(m_segments.size());
			if (parent == Previous) parent = index - 1; // -1 == Base for the first segment
			if (parent < Base || parent >= index)
				throw std::invalid_argument("KinematicChain: the parent must be Base or a segment added earlier");
			s.parent = parent;
			if (s.type != JointType::Fixed) s.dofIndex = m_nDof++;
			std::vector<int> dofs = (parent == Base) ? std::vector<int>() : m_ancestorDofs[parent];
			if (s.dofIndex >= 0) dofs.push_back(s.dofIndex);
			m_ancestorDofs.push_back(std::move(dofs));
			m_segments.push_back(s);
			return index;
		}

		int resolve(int i) const {
			if (i == -1) i = nSegments() - 1;
			if (i < 0 || i >= nSegments()) throw std::out_of_range("KinematicChain: no such segment");
			return i;
		}

		void check(const VectorX &q) const {
			if (q.size() != m_nDof) throw std::invalid_argument("KinematicChain: q has the wrong size");
		}

		static SE3Type localPose(const Segment &s, const VectorX &q) {
			if (s.dofIndex < 0) return s.rest;
			return SE3Type::computeExp(q[s.dofIndex] * s.screw) * s.rest;
		}

		std::vector<Segment> m_segments;
		std::vector<std::vector<int>> m_ancestorDofs; ///< per segment, the joints on the path from the base
		int m_nDof = 0;
	};

} // namespace sofa::component::cosserat::liegroups
