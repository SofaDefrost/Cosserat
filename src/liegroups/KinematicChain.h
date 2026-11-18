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

		// ---------------------------------------------------------------------------------------------------------
		// World-frame quantities, velocity-dependent terms and geometric stiffness (used by KinematicChainMapping).
		//
		// "World" = the frame the base pose T0 maps the base frame to. Screws are the unit joint screws S_k = Ad_{T0 T_parent} xi_k
		// (6 x nDof, [w ; r], reference point = world origin) and positions are the origins of the segments.
		// A spatial twist V = [omega ; u] has u = velocity of the (instantaneous) point at the world origin, so the velocity of
		// a point p is u + omega x p. SOFA's Rigid3d velocities are the hybrid form [omega ; velocity of the origin].
		// ---------------------------------------------------------------------------------------------------------

		struct WorldFrame {
			std::vector<SE3Type> poses; ///< T0 * T^0_i
			std::vector<Vector3> positions; ///< origin of every segment
			JacobianMatrix screws; ///< 6 x nDof
		};

		/// Poses, origins and joint screws in the world frame for the base pose @p T0.
		WorldFrame worldFrame(const VectorX &q, const SE3Type &T0 = SE3Type()) const {
			WorldFrame wf;
			const std::vector<SE3Type> local = allPoses(q);
			wf.poses.resize(local.size());
			wf.positions.resize(local.size());
			for (std::size_t i = 0; i < local.size(); ++i) {
				wf.poses[i] = T0 * local[i];
				wf.positions[i] = wf.poses[i].translation();
			}
			wf.screws = T0.computeAdjoint() * jointScrews(local);
			return wf;
		}

		/// ad(V) X for twists in the [angular ; linear] layout: [w x x_w ; w x x_v + v x x_w].
		static TangentVector adAct(const TangentVector &V, const TangentVector &X) {
			const Vector3 w = V.template head<3>(), v = V.template tail<3>();
			TangentVector r;
			r << w.cross(X.template head<3>()), w.cross(X.template tail<3>()) + v.cross(X.template head<3>());
			return r;
		}

		/// Spatial twist of every segment for the joint rates @p dq and the base spatial twist @p baseTwist.
		std::vector<TangentVector> spatialTwists(const JacobianMatrix &S, const VectorX &dq,
												 const TangentVector &baseTwist = TangentVector::Zero()) const {
			std::vector<TangentVector> V(m_segments.size());
			for (std::size_t k = 0; k < m_segments.size(); ++k) {
				const Segment &s = m_segments[k];
				TangentVector v = (s.parent == Base) ? baseTwist : V[s.parent];
				if (s.dofIndex >= 0) v += S.col(s.dofIndex) * dq[s.dofIndex];
				V[k] = v;
			}
			return V;
		}

		/**
		 * @brief Velocity-dependent part of the acceleration of every segment, [alpha ; a_origin] (hybrid, world frame):
		 *        what is left of d/dt(J dx) when the generalised accelerations are zero (so J xddot is NOT included).
		 *
		 * The joint screws move with their parent: dS_k/dt = ad(V_parent) S_k, hence dV_i/dt = dV_parent/dt + ad(V_parent) S_k dq_k.
		 * @param baseTwist  spatial twist of the base (zero for a fixed base)
		 * @param baseBias   spatial-twist derivative of the base when its own acceleration is zero: [0 ; -omega_b x v_b],
		 *                   v_b being the velocity of the base reference point (zero for a fixed base)
		 */
		std::vector<TangentVector> biasAcceleration(const JacobianMatrix &S, const std::vector<Vector3> &p, const VectorX &dq,
													const TangentVector &baseTwist = TangentVector::Zero(),
													const TangentVector &baseBias = TangentVector::Zero()) const {
			const std::vector<TangentVector> V = spatialTwists(S, dq, baseTwist);
			std::vector<TangentVector> Vd(m_segments.size()), out(m_segments.size());
			for (std::size_t k = 0; k < m_segments.size(); ++k) {
				const Segment &s = m_segments[k];
				const TangentVector &parentV = (s.parent == Base) ? baseTwist : V[s.parent];
				TangentVector vd = (s.parent == Base) ? baseBias : Vd[s.parent];
				if (s.dofIndex >= 0) vd += adAct(parentV, S.col(s.dofIndex)) * dq[s.dofIndex];
				Vd[k] = vd;
				const Vector3 w = V[k].template head<3>();
				const Vector3 vel = V[k].template tail<3>() + w.cross(p[k]);
				out[k] << vd.template head<3>(), vd.template tail<3>() + vd.template head<3>().cross(p[k]) + w.cross(vel);
			}
			return out;
		}

		struct GeometricStiffness {
			VectorX dQ; ///< change of the generalised forces J^T W
			Vector3 dMoment; ///< change of the total moment about the world origin (feeds the base torque)
		};

		/**
		 * @brief Directional derivative of the generalised forces Q = J^T W for FIXED world-frame wrenches (the geometric stiffness
		 *        K_g dx = d(J^T)/dx . dx . W), for the joint displacement @p dq and a base displacement of spatial twist @p dBaseTwist.
		 * @param wrench per segment [torque about the segment origin ; force], world frame
		 * Q_k = S_k . sum_{i in subtree(k)} [T_i + p_i x F_i ; F_i];  dS_k = ad(dV_parent(k)) S_k,  dp_i = origin velocity under (dq, dBase).
		 * The total force sum F_i is constant; the total moment changes by sum dp_i x F_i (dMoment).
		 */
		GeometricStiffness geometricStiffness(const JacobianMatrix &S, const std::vector<Vector3> &p,
											  const std::vector<TangentVector> &wrench, const VectorX &dq,
											  const TangentVector &dBaseTwist = TangentVector::Zero()) const {
			const std::size_t n = m_segments.size();
			const std::vector<TangentVector> dV = spatialTwists(S, dq, dBaseTwist);
			std::vector<TangentVector> W(n), dW(n);
			for (std::size_t i = 0; i < n; ++i) {
				const Vector3 T = wrench[i].template head<3>(), F = wrench[i].template tail<3>();
				W[i] << T + p[i].cross(F), F;
				const Vector3 dp = dV[i].template tail<3>() + dV[i].template head<3>().cross(p[i]);
				dW[i] << dp.cross(F), Vector3::Zero();
			}
			GeometricStiffness res;
			res.dQ = VectorX::Zero(m_nDof);
			res.dMoment = Vector3::Zero();
			for (std::size_t i = n; i-- > 0;) {
				const Segment &s = m_segments[i];
				if (s.dofIndex >= 0) {
					const TangentVector &parentDV = (s.parent == Base) ? dBaseTwist : dV[s.parent];
					const TangentVector Sk = S.col(s.dofIndex);
					res.dQ[s.dofIndex] = adAct(parentDV, Sk).dot(W[i]) + Sk.dot(dW[i]);
				}
				if (s.parent >= 0) {
					W[s.parent] += W[i];
					dW[s.parent] += dW[i];
				} else
					res.dMoment += dW[i].template head<3>();
			}
			return res;
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
