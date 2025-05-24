/******************************************************************************
 *                 SOFA, Simulation Open-Framework Architecture                *
 *                 (c) 2006 INRIA, USTL, UJF, CNRS, MGH                       *
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
 * along with this program. If not, see <http://www.gnu.org/licenses/>.        *
 ******************************************************************************/
#pragma once

#include <Cosserat/controller/CosseratILQRController.h>
#include <sofa/core/objectmodel/Event.h>
#include <sofa/helper/logging/Message.h>
#include <sofa/simulation/AnimateBeginEvent.h>

#include <Eigen/Dense>

namespace Cosserat::controller {

	using sofa::helper::WriteAccessor;
	using sofa::type::vector;

	// ─────────────────────────────────────────────────────────────────────────
	// Constructor
	// ─────────────────────────────────────────────────────────────────────────

	template<class TIn1, class TIn2, class TOut>
	CosseratILQRController<TIn1, TIn2, TOut>::CosseratILQRController() :
		Inherit1(),
		l_mapping(initLink("mapping",
			"Link to the Strain2RigidCosseratMapping whose body Jacobian "
			"is used for the iLQR optimisation.")),
		d_targetPose(initData(&d_targetPose, OutCoord(),
			"targetPose",
			"Desired tip pose in world frame (Rigid3: position + quaternion).")),
		d_Q_tip(initData(&d_Q_tip, 1.0,
			"Q_tip",
			"Tip error weight Q (isotropic, applied on se(3) error).")),
		d_R_strain(initData(&d_R_strain, 1.0e-4,
			"R_strain",
			"Strain regularisation weight R (Tikhonov on delta-xi).")),
		d_stepSize(initData(&d_stepSize, 0.1,
			"stepSize",
			"Gradient descent step size alpha (mode 0 only).")),
		d_maxIterations(initData(&d_maxIterations, 5,
			"maxIterations",
			"Maximum iLQR iterations per simulation step.")),
		d_tolerance(initData(&d_tolerance, 1.0e-6,
			"tolerance",
			"Early-stop threshold on ||delta_xi|| (0 = always run all iterations).")),
		d_mode(initData(&d_mode, 1,
			"mode",
			"Control mode:\n"
			"  0 = gradient descent (delta_xi = -alpha * J^T * Q * e)\n"
			"  1 = Gauss-Newton     (solve (J^T*Q*J + R*I)*delta_xi = -J^T*Q*e)")),
		d_tipIndex(initData(&d_tipIndex, -1,
			"tipIndex",
			"Index of the output frame used as the tracked tip. "
			"-1 means last frame.")),
		d_active(initData(&d_active, true,
			"active",
			"Enable the controller (apply corrections each step).")),
		d_tipError(initData(&d_tipError, 0.0,
			"tipError",
			"[READ-ONLY] Current tip error ||log(g_ref^-1 * g_tip)|| "
			"(updated each step).")),
		d_manipulability(initData(&d_manipulability, 0.0,
			"manipulability",
			"[READ-ONLY] Yoshikawa manipulability at the tip section."))
	{
		this->f_listening.setValue(true);  // subscribe to AnimateBeginEvent
	}

	// ─────────────────────────────────────────────────────────────────────────
	// init
	// ─────────────────────────────────────────────────────────────────────────

	template<class TIn1, class TIn2, class TOut>
	void CosseratILQRController<TIn1, TIn2, TOut>::init() {
		Inherit1::init();

		if (!l_mapping.get()) {
			msg_error() << "CosseratILQRController: 'mapping' link is not set. "
			               "Please link to a Strain2RigidCosseratMapping.";
			return;
		}

		if (!l_mapping->getStrainState()) {
			msg_error() << "CosseratILQRController: mapping has no strain state.";
			return;
		}

		if (!l_mapping->getFramesState()) {
			msg_error() << "CosseratILQRController: mapping has no frames state.";
			return;
		}

		msg_info() << "CosseratILQRController: initialised. "
		           << "Mode=" << d_mode.getValue()
		           << "  Q=" << d_Q_tip.getValue()
		           << "  R=" << d_R_strain.getValue()
		           << "  maxIter=" << d_maxIterations.getValue();
	}

	// ─────────────────────────────────────────────────────────────────────────
	// handleEvent
	// ─────────────────────────────────────────────────────────────────────────

	template<class TIn1, class TIn2, class TOut>
	void CosseratILQRController<TIn1, TIn2, TOut>::handleEvent(
			sofa::core::objectmodel::Event *event) {

		if (!d_active.getValue()) return;
		if (!sofa::simulation::AnimateBeginEvent::checkEventType(event)) return;
		if (!l_mapping.get()) return;

		const std::vector<Coord1> corrections = computeControl();
		if (corrections.empty()) return;

		// Write updated strains back to the mechanical state
		sofa::core::State<TIn1> *strainState = l_mapping->getStrainState();
		if (!strainState) return;

		auto *strain_data = strainState->write(sofa::core::vec_id::write_access::position);
		if (!strain_data) return;

		sofa::VecCoord_t<TIn1> &strains = *strain_data->beginEdit();
		const int N = static_cast<int>(corrections.size());

		for (int k = 0; k < N && k < static_cast<int>(strains.size()); ++k) {
			for (int j = 0; j < NStrainDOF; ++j)
				strains[k][j] += corrections[k][j];
		}
		strain_data->endEdit();
	}

	// ─────────────────────────────────────────────────────────────────────────
	// computeControl
	// ─────────────────────────────────────────────────────────────────────────

	template<class TIn1, class TIn2, class TOut>
	std::vector<typename CosseratILQRController<TIn1, TIn2, TOut>::Coord1>
	CosseratILQRController<TIn1, TIn2, TOut>::computeControl() {
		if (!l_mapping.get()) return {};

		const BodyJacobian &bj = l_mapping->getBodyJacobian();
		const int N = bj.size();
		if (N == 0) return {};

		// ── Tip index ─────────────────────────────────────────────────────
		const sofa::core::State<TOut> *framesState = l_mapping->getFramesState();
		if (!framesState) return {};

		const auto &frames_data =
			framesState->read(sofa::core::vec_id::read_access::position)->getValue();
		const int n_frames = static_cast<int>(frames_data.size());
		if (n_frames == 0) return {};

		int tip_idx = d_tipIndex.getValue();
		if (tip_idx < 0 || tip_idx >= n_frames)
			tip_idx = n_frames - 1;

		// ── Target and current tip pose ────────────────────────────────────
		const SE3Types g_ref = rigidToSE3(d_targetPose.getValue());
		const SE3Types g_tip = rigidToSE3(frames_data[tip_idx]);

		// ── Tip error in se(3): e = log(g_ref⁻¹ · g_tip) ─────────────────
		const TangentVector e = g_ref.computeInverse().compose(g_tip).log();
		const double err_norm = e.norm();

		// Update diagnostic outputs
		d_tipError.setValue(err_norm);
		d_manipulability.setValue(bj.manipulability(N - 1));

		if (err_norm < d_tolerance.getValue() && d_maxIterations.getValue() > 0) {
			msg_info() << "CosseratILQRController: tip error " << err_norm
			           << " < tolerance, skipping update.";
			return {};
		}

		// ── Tip Jacobian J ∈ ℝ^{6 × 6N} ──────────────────────────────────
		// J is the body Jacobian at the last section = the tip Jacobian.
		// For strains with DOF < 6 (Vec3), we later extract only NStrainDOF
		// columns per block.
		const MatrixXd J_full = bj.getJacobianAtSection(N - 1);  // 6×6N

		// If NStrainDOF < 6, keep only the first NStrainDOF columns per block.
		// J_act ∈ ℝ^{6 × N*NStrainDOF} — the "actuated" sub-Jacobian.
		MatrixXd J_act(6, N * NStrainDOF);
		for (int k = 0; k < N; ++k)
			J_act.block(0, k * NStrainDOF, 6, NStrainDOF) =
				J_full.block(0, k * 6, 6, NStrainDOF);

		// ── iLQR iterations ───────────────────────────────────────────────
		const double Q = d_Q_tip.getValue();
		const double R = d_R_strain.getValue();
		const double alpha = d_stepSize.getValue();
		const int    mode  = d_mode.getValue();
		const int    maxIt = d_maxIterations.getValue();
		const double tol   = d_tolerance.getValue();

		VectorXd delta_xi = VectorXd::Zero(N * NStrainDOF);
		TangentVector e_k = e;  // error at current iterate (re-linearised per iter)

		for (int it = 0; it < maxIt; ++it) {
			VectorXd step;
			switch (mode) {
				case 0:
					step = gradientStep(J_act, e_k, alpha, Q);
					break;
				case 1:
				default:
					step = gaussNewtonStep(J_act, e_k, Q, R);
					break;
			}

			delta_xi += step;

			// Predicted residual after this step: e_{it+1} ≈ e_k + J_act·step
			// (first-order; re-evaluation would require calling apply() again)
			const TangentVector e_predicted = e_k + TangentVector(J_act * step);
			e_k = e_predicted;

			if (delta_xi.norm() > 0 && step.norm() < tol) break;
		}

		msg_info() << "CosseratILQRController: ||e||=" << err_norm
		           << "  ||delta_xi||=" << delta_xi.norm();

		return unpackDelta(delta_xi, N);
	}

	// ─────────────────────────────────────────────────────────────────────────
	// Private helpers
	// ─────────────────────────────────────────────────────────────────────────

	template<class TIn1, class TIn2, class TOut>
	typename CosseratILQRController<TIn1, TIn2, TOut>::SE3Types
	CosseratILQRController<TIn1, TIn2, TOut>::getCurrentTipPose(int tip_idx) const {
		const sofa::core::State<TOut> *fs = l_mapping->getFramesState();
		if (!fs) return SE3Types();

		const auto &frames = fs->read(sofa::core::vec_id::read_access::position)->getValue();
		if (frames.empty()) return SE3Types();

		const int idx = (tip_idx < 0 || tip_idx >= static_cast<int>(frames.size()))
			? static_cast<int>(frames.size()) - 1
			: tip_idx;
		return rigidToSE3(frames[idx]);
	}

	template<class TIn1, class TIn2, class TOut>
	typename CosseratILQRController<TIn1, TIn2, TOut>::SE3Types
	CosseratILQRController<TIn1, TIn2, TOut>::rigidToSE3(const OutCoord &coord) {
		const auto &c = coord.getCenter();
		const auto &q = coord.getOrientation();
		return SE3Types(
			SO3Type(Eigen::Quaternion<double>(q[3], q[0], q[1], q[2])),
			Eigen::Vector3d(c[0], c[1], c[2]));
	}

	template<class TIn1, class TIn2, class TOut>
	typename CosseratILQRController<TIn1, TIn2, TOut>::VectorXd
	CosseratILQRController<TIn1, TIn2, TOut>::gradientStep(
			const MatrixXd &J_tip, const TangentVector &e,
			double alpha, double Q) {
		// Δξ = -α · J^T · Q · e
		return -alpha * Q * (J_tip.transpose() * e);
	}

	template<class TIn1, class TIn2, class TOut>
	typename CosseratILQRController<TIn1, TIn2, TOut>::VectorXd
	CosseratILQRController<TIn1, TIn2, TOut>::gaussNewtonStep(
			const MatrixXd &J_tip, const TangentVector &e,
			double Q, double R) {
		// Solve (J^T·Q·J + R·I) · Δξ = -J^T·Q·e  via Cholesky (LDLT)
		const int n = static_cast<int>(J_tip.cols());
		const MatrixXd H = Q * (J_tip.transpose() * J_tip)
		                   + R * MatrixXd::Identity(n, n);
		const VectorXd g = Q * (J_tip.transpose() * e);
		return -H.ldlt().solve(g);
	}

	template<class TIn1, class TIn2, class TOut>
	std::vector<typename CosseratILQRController<TIn1, TIn2, TOut>::Coord1>
	CosseratILQRController<TIn1, TIn2, TOut>::unpackDelta(const VectorXd &delta_xi, int N) {
		std::vector<Coord1> result(static_cast<std::size_t>(N));
		for (int k = 0; k < N; ++k) {
			Coord1 dxi;
			for (int j = 0; j < NStrainDOF; ++j)
				dxi[j] = delta_xi[k * NStrainDOF + j];
			result[k] = dxi;
		}
		return result;
	}

} // namespace Cosserat::controller
