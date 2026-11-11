// Rigid-body dynamics of a KinematicChain (trees included), spatial-vector form, world frame, fixed base.
//
//   Spatial twist   V = [omega ; u]   (u = velocity of the point at the world origin; the library layout [phi ; rho])
//   Spatial wrench  F = [m ; f]       (force f and moment m about the world origin)
//   Spatial inertia of a segment about its own origin (frame of the segment):
//
//        I_b = [  I_o      m [c]x ]        I_o = I_c - m [c]x [c]x   (parallel axis, I_c about the centre of mass c)
//              [ -m [c]x   m 1    ]
//
//   carried to the world frame by  I_w = Ad_{T^-1}^T I_b Ad_{T^-1}  (T = pose of the segment).
//
// Newton-Euler for one body in the (inertial) world frame:      F = I_w Vdot - ad_V^T I_w V
// Joint screws S_k (KinematicChain::jointScrews) move with their parent:   Vdot_i = Vdot_parent + S_k qddot_k + ad(V_parent) S_k qdot_k
// Recursive Newton-Euler (RNEA):   forward pass for V, Vdot, F_i;  backward pass tau_k = S_k . sum_{i in subtree(k)} F_i.
// Gravity is a fictitious base acceleration Vdot_base = [0 ; -g], so tau also holds the gravity compensation torques.
//
//   tau = M(q) qddot + h(q, qdot)          M(q) : massMatrix      h : biasForces (Coriolis, centrifugal, gravity)
//   forward dynamics: qddot = M^-1 (tau - h)
//
// Reference for the tests: Lagrange's equations built from hybrid Jacobians and finite differences
// (simulations/articulated/lie_tests/test_lie_dynamics.py).

#pragma once

#include <memory>
#include <vector>

#include "KinematicChain.h"

namespace sofa::component::cosserat::liegroups {

	template<typename Scalar_>
	class KinematicChainDynamics {
	public:
		using Scalar = Scalar_;
		using Chain = KinematicChain<Scalar>;
		using SE3Type = typename Chain::SE3Type;
		using Vector3 = typename Chain::Vector3;
		using Matrix3 = Eigen::Matrix<Scalar, 3, 3>;
		using Matrix6 = Eigen::Matrix<Scalar, 6, 6>;
		using TangentVector = typename Chain::TangentVector;
		using VectorX = typename Chain::VectorX;
		using MatrixX = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;

		/// Inertial parameters of one segment, in the segment frame.
		struct BodyInertia {
			Scalar mass = Scalar(0);
			Vector3 com = Vector3::Zero(); ///< centre of mass
			Matrix3 inertiaCom = Matrix3::Zero(); ///< rotational inertia about the centre of mass, segment axes

			/// 6x6 spatial inertia about the segment origin, layout [angular ; linear].
			Matrix6 matrix() const {
				const Matrix3 cx = hat(com);
				Matrix6 M;
				M << inertiaCom - mass * cx * cx, mass * cx, -mass * cx, mass * Matrix3::Identity();
				return M;
			}
		};

		explicit KinematicChainDynamics(Chain chain) : m_chain(std::move(chain)), m_inertia(m_chain.nSegments()) {}

		const Chain &chain() const { return m_chain; }
		const BodyInertia &inertia(int segment) const { return m_inertia.at(segment); }
		void setInertia(int segment, const BodyInertia &inertia) { m_inertia.at(segment) = inertia; }

		static Matrix3 hat(const Vector3 &a) {
			Matrix3 m;
			m << 0, -a[2], a[1], a[2], 0, -a[0], -a[1], a[0], 0;
			return m;
		}

		/// 6x6 matrix of ad(V): [[omega x, 0], [u x, omega x]].
		static Matrix6 adMatrix(const TangentVector &V) {
			const Matrix3 w = hat(V.template head<3>()), u = hat(V.template tail<3>());
			Matrix6 m = Matrix6::Zero();
			m.template topLeftCorner<3, 3>() = w;
			m.template bottomLeftCorner<3, 3>() = u;
			m.template bottomRightCorner<3, 3>() = w;
			return m;
		}

		/// Spatial inertia of a segment in the world frame, I_w = Ad_{T^-1}^T I_b Ad_{T^-1}.
		Matrix6 worldInertia(int segment, const SE3Type &pose) const {
			const Matrix6 A = pose.inverse().computeAdjoint();
			return A.transpose() * m_inertia[segment].matrix() * A;
		}

		/**
		 * @brief Inverse dynamics (RNEA): the joint forces that produce @p ddq at (q, dq) under @p gravity.
		 * @param external optional per-segment external wrenches [torque about the segment origin ; force], world frame, which
		 *        are subtracted (they help the motion), so the result is the actuation needed on top of them.
		 */
		VectorX inverseDynamics(const VectorX &q, const VectorX &dq, const VectorX &ddq, const Vector3 &gravity = Vector3::Zero(),
								const std::vector<TangentVector> *external = nullptr) const {
			const std::vector<SE3Type> poses = m_chain.allPoses(q);
			const typename Chain::JacobianMatrix S = m_chain.jointScrews(poses);
			const std::size_t n = poses.size();
			std::vector<TangentVector> V(n), Vd(n), F(n);
			TangentVector baseVd = TangentVector::Zero();
			baseVd.template tail<3>() = -gravity;

			for (std::size_t k = 0; k < n; ++k) {
				const auto &s = m_chain.segments()[k];
				const TangentVector parentV = (s.parent == Chain::Base) ? TangentVector(TangentVector::Zero()) : V[s.parent];
				TangentVector v = parentV;
				TangentVector vd = (s.parent == Chain::Base) ? baseVd : Vd[s.parent];
				if (s.dofIndex >= 0) {
					v += S.col(s.dofIndex) * dq[s.dofIndex];
					vd += S.col(s.dofIndex) * ddq[s.dofIndex] + Chain::adAct(parentV, S.col(s.dofIndex)) * dq[s.dofIndex];
				}
				V[k] = v;
				Vd[k] = vd;
				const Matrix6 Iw = worldInertia(static_cast<int>(k), poses[k]);
				F[k] = Iw * vd - adMatrix(v).transpose() * (Iw * v);
				if (external) {
					const Vector3 T = (*external)[k].template head<3>(), f = (*external)[k].template tail<3>();
					TangentVector w;
					w << T + poses[k].translation().cross(f), f;
					F[k] -= w;
				}
			}
			VectorX tau = VectorX::Zero(m_chain.nDof());
			for (std::size_t i = n; i-- > 0;) {
				const auto &s = m_chain.segments()[i];
				if (s.dofIndex >= 0) tau[s.dofIndex] = TangentVector(S.col(s.dofIndex)).dot(F[i]);
				if (s.parent >= 0) F[s.parent] += F[i];
			}
			return tau;
		}

		/// Joint-space inertia matrix, one RNEA per column (no gravity, no velocity).
		MatrixX massMatrix(const VectorX &q) const {
			const int n = m_chain.nDof();
			MatrixX M(n, n);
			const VectorX zero = VectorX::Zero(n);
			for (int j = 0; j < n; ++j) {
				VectorX e = VectorX::Zero(n);
				e[j] = Scalar(1);
				M.col(j) = inverseDynamics(q, zero, e);
			}
			return M;
		}

		/// h(q, dq): Coriolis, centrifugal and gravity forces (the RNEA with zero acceleration).
		VectorX biasForces(const VectorX &q, const VectorX &dq, const Vector3 &gravity = Vector3::Zero(),
						   const std::vector<TangentVector> *external = nullptr) const {
			return inverseDynamics(q, dq, VectorX::Zero(m_chain.nDof()), gravity, external);
		}

		/// qddot = M^-1 (tau - h).
		VectorX forwardDynamics(const VectorX &q, const VectorX &dq, const VectorX &tau, const Vector3 &gravity = Vector3::Zero(),
								const std::vector<TangentVector> *external = nullptr) const {
			const MatrixX M = massMatrix(q);
			return M.ldlt().solve(tau - biasForces(q, dq, gravity, external));
		}

		Scalar kineticEnergy(const VectorX &q, const VectorX &dq) const {
			const std::vector<SE3Type> poses = m_chain.allPoses(q);
			const typename Chain::JacobianMatrix S = m_chain.jointScrews(poses);
			const std::vector<TangentVector> V = m_chain.spatialTwists(S, dq);
			Scalar e = 0;
			for (std::size_t k = 0; k < poses.size(); ++k) e += Scalar(0.5) * V[k].dot(worldInertia(static_cast<int>(k), poses[k]) * V[k]);
			return e;
		}

		/// Gravitational potential energy, U = - sum_i m_i g . (centre of mass of segment i).
		Scalar potentialEnergy(const VectorX &q, const Vector3 &gravity) const {
			const std::vector<SE3Type> poses = m_chain.allPoses(q);
			Scalar u = 0;
			for (std::size_t k = 0; k < poses.size(); ++k) {
				const Vector3 c = poses[k].translation() + poses[k].rotation().matrix() * m_inertia[k].com;
				u -= m_inertia[k].mass * gravity.dot(c);
			}
			return u;
		}

	private:
		Chain m_chain;
		std::vector<BodyInertia> m_inertia;
	};

} // namespace sofa::component::cosserat::liegroups
