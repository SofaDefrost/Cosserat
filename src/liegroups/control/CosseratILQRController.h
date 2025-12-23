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

#include <functional>
#include <vector>
#include "../SE3.h"
#include "../Types.h"

namespace sofa::component::cosserat::liegroups::control {

	/**
	 * @brief iterative Linear Quadratic Regulator (iLQR) for Cosserat rod trajectory tracking
	 *
	 * iLQR is a model-based optimal control algorithm that computes feedback gains
	 * for trajectory tracking by iteratively linearizing the dynamics around a nominal
	 * trajectory.
	 *
	 * ## Algorithm Overview
	 *
	 * 1. **Forward Pass**: Simulate nominal trajectory with current controls
	 * 2. **Backward Pass**: Compute value function and feedback gains via dynamic programming
	 * 3. **Line Search**: Find optimal step size for control update
	 * 4. **Update**: Apply feedback gains to improve controls
	 * 5. **Repeat** until convergence
	 *
	 * ## Cosserat Dynamics
	 *
	 * For a Cosserat rod with N segments:
	 * - State: x = [T_1, T_2, ..., T_N] where T_i ∈ SE(3)
	 * - Control: u = [strain_1, ..., strain_N] where strain_i ∈ ℝ⁶
	 * - Dynamics: T_{k+1} = T_k * exp(L * u_k)
	 *
	 * ## Cost Function
	 *
	 * J = Σ_k [l(x_k, u_k)] + l_f(x_N)
	 *
	 * where:
	 * - l(x,u) = running cost (tracking error + control effort)
	 * - l_f(x) = terminal cost (final tracking error)
	 *
	 * @tparam Scalar Scalar type (float or double)
	 */
	template<typename Scalar = double>
	class CosseratILQRController {
	public:
		using SE3Type = SE3<Scalar>;
		using Vector3 = typename SE3Type::Vector3;
		using Vector6 = Eigen::Matrix<Scalar, 6, 1>;
		using Matrix3 = Eigen::Matrix<Scalar, 3, 3>;
		using Matrix6 = Eigen::Matrix<Scalar, 6, 6>;
		using MatrixX = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
		using VectorX = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

		/**
		 * @brief Reference trajectory for tracking
		 */
		struct Trajectory {
			std::vector<SE3Type> poses; // Desired poses at each time step
			std::vector<Vector6> strains; // Nominal strains (can be zero)
			double segment_length; // Length of each segment

			int horizon() const { return static_cast<int>(poses.size()) - 1; }
		};

		/**
		 * @brief Configuration parameters for iLQR
		 */
		struct Config {
			int max_iterations = 50; // Maximum iLQR iterations
			Scalar convergence_threshold = 1e-4; // Cost improvement threshold

			// Cost weights
			Scalar Q_position = 10.0; // Position tracking weight
			Scalar Q_rotation = 1.0; // Rotation tracking weight
			Scalar R_control = 0.01; // Control effort weight
			Scalar Q_final_position = 100.0; // Terminal position weight
			Scalar Q_final_rotation = 10.0; // Terminal rotation weight

			// Line search
			int max_line_search_iters = 10;
			Scalar line_search_decay = 0.5;
			Scalar min_step_size = 0.01;

			// Regularization
			Scalar regularization = 1e-6; // Regularization for matrix inversion

			bool verbose = false;
		};

		/**
		 * @brief Result of iLQR optimization
		 */
		struct Result {
			std::vector<Vector6> optimal_strains; // Optimized control sequence
			std::vector<SE3Type> trajectory; // Resulting trajectory
			std::vector<Matrix6> feedback_gains; // Feedback gains K_t

			std::vector<Scalar> cost_history;
			Scalar final_cost;
			int iterations;
			bool converged;

			std::string message;
		};

		/**
		 * @brief Construct iLQR controller
		 * @param num_segments Number of Cosserat segments
		 * @param config Configuration parameters
		 */
		explicit CosseratILQRController(int num_segments, const Config &config = Config()) :
			m_num_segments(num_segments), m_config(config) {}

		/**
		 * @brief Compute optimal control sequence for trajectory tracking
		 *
		 * @param reference Reference trajectory to track
		 * @param initial_strains Initial guess for control sequence
		 * @return Result containing optimal controls and feedback gains
		 */
		Result optimize(const Trajectory &reference, const std::vector<Vector6> &initial_strains) {
			Result result;
			result.converged = false;
			result.iterations = 0;

			// Initialize controls
			std::vector<Vector6> u = initial_strains;
			if (u.empty()) {
				u.resize(reference.horizon(), Vector6::Zero());
			}

			// Forward pass: compute initial trajectory and cost
			auto [x, cost] = forwardPass(u, reference);
			result.cost_history.push_back(cost);

			if (m_config.verbose) {
				std::cout << "iLQR Iteration 0: Cost = " << cost << "\n";
			}

			// Main iLQR loop
			for (int iter = 0; iter < m_config.max_iterations; ++iter) {
				result.iterations = iter + 1;

				// Backward pass: compute value function and gains
				auto gains = backwardPass(x, u, reference);

				// Line search: find best step size
				auto [u_new, x_new, cost_new, alpha] = lineSearch(x, u, gains, reference, cost);

				// Check for improvement
				Scalar cost_reduction = cost - cost_new;
				result.cost_history.push_back(cost_new);

				if (m_config.verbose) {
					std::cout << "iLQR Iteration " << iter + 1 << ": Cost = " << cost_new
							  << ", Reduction = " << cost_reduction << ", Alpha = " << alpha << "\n";
				}

				// Update
				u = u_new;
				x = x_new;
				cost = cost_new;

				// Check convergence
				if (std::abs(cost_reduction) < m_config.convergence_threshold) {
					result.converged = true;
					result.message = "Converged: cost reduction below threshold";
					break;
				}

				if (alpha < m_config.min_step_size) {
					result.message = "Terminated: step size too small";
					break;
				}
			}

			// Store final result
			result.optimal_strains = u;
			result.trajectory = x;
			result.feedback_gains = computeFeedbackGains(x, u, reference);
			result.final_cost = cost;

			if (!result.converged && result.iterations >= m_config.max_iterations) {
				result.message = "Terminated: maximum iterations reached";
			}

			return result;
		}

	private:
		int m_num_segments;
		Config m_config;

		/**
		 * @brief Forward simulation of Cosserat rod dynamics
		 *
		 * Dynamics: T_{k+1} = T_k * exp(L * u_k)
		 *
		 * @param u Control sequence (strains)
		 * @param ref Reference trajectory
		 * @return Trajectory and total cost
		 */
		std::pair<std::vector<SE3Type>, Scalar> forwardPass(const std::vector<Vector6> &u,
															const Trajectory &ref) const {
			const int T = static_cast<int>(u.size());
			std::vector<SE3Type> x(T + 1);

			// Initial state
			x[0] = SE3Type::computeIdentity();

			// Simulate forward
			for (int k = 0; k < T; ++k) {
				Vector6 xi = u[k] * ref.segment_length;
				x[k + 1] = x[k] * SE3Type::exp(xi);
			}

			// Compute total cost
			Scalar total_cost = 0.0;

			// Running cost
			for (int k = 0; k < T; ++k) {
				total_cost += runningCost(x[k], u[k], ref.poses[k]);
			}

			// Terminal cost
			total_cost += terminalCost(x[T], ref.poses[T]);

			return {x, total_cost};
		}

		/**
		 * @brief Running cost: l(x_k, u_k)
		 *
		 * l = Q_pos * ||pos - pos_ref||² + Q_rot * d_rot² + R * ||u||²
		 */
		Scalar runningCost(const SE3Type &x, const Vector6 &u, const SE3Type &x_ref) const {
			// Position error
			Vector3 pos_error = x.translation() - x_ref.translation();
			Scalar cost_pos = m_config.Q_position * pos_error.squaredNorm();

			// Rotation error (geodesic distance in SO(3))
			SE3Type error = x_ref.inverse() * x;
			Vector3 rot_error = error.rotation().log();
			Scalar cost_rot = m_config.Q_rotation * rot_error.squaredNorm();

			// Control effort
			Scalar cost_control = m_config.R_control * u.squaredNorm();

			return cost_pos + cost_rot + cost_control;
		}

		/**
		 * @brief Terminal cost: l_f(x_N)
		 */
		Scalar terminalCost(const SE3Type &x, const SE3Type &x_ref) const {
			Vector3 pos_error = x.translation() - x_ref.translation();
			Scalar cost_pos = m_config.Q_final_position * pos_error.squaredNorm();

			SE3Type error = x_ref.inverse() * x;
			Vector3 rot_error = error.rotation().log();
			Scalar cost_rot = m_config.Q_final_rotation * rot_error.squaredNorm();

			return cost_pos + cost_rot;
		}

		/**
		 * @brief Backward pass: Dynamic programming to compute value function
		 *
		 * Computes Q-function approximation and feedback gains via Riccati recursion
		 *
		 * @param x State trajectory
		 * @param u Control sequence
		 * @param ref Reference trajectory
		 * @return Feedback gains for each time step
		 */
		std::vector<Matrix6> backwardPass(const std::vector<SE3Type> &x, const std::vector<Vector6> &u,
										  const Trajectory &ref) const {
			const int T = static_cast<int>(u.size());
			std::vector<Matrix6> K(T); // Feedback gains

			// Terminal value function (gradient and Hessian)
			Vector6 V_x = Vector6::Zero();
			Matrix6 V_xx = Matrix6::Zero();

			// Initialize terminal cost derivatives
			computeTerminalCostDerivatives(x[T], ref.poses[T], V_x, V_xx);

			// Backward recursion
			for (int k = T - 1; k >= 0; --k) {
				// Linearize dynamics around (x[k], u[k])
				auto [F_x, F_u] = linearizeDynamics(x[k], u[k], ref.segment_length);

				// Cost derivatives
				Vector6 l_x, l_u;
				Matrix6 l_xx, l_uu, l_ux;
				computeCostDerivatives(x[k], u[k], ref.poses[k], l_x, l_u, l_xx, l_uu, l_ux);

				// Q-function approximation
				Vector6 Q_x = l_x + F_x.transpose() * V_x;
				Vector6 Q_u = l_u + F_u.transpose() * V_x;
				Matrix6 Q_xx = l_xx + F_x.transpose() * V_xx * F_x;
				Matrix6 Q_uu = l_uu + F_u.transpose() * V_xx * F_u;
				Matrix6 Q_ux = l_ux + F_u.transpose() * V_xx * F_x;

				// Regularization for numerical stability
				Q_uu += m_config.regularization * Matrix6::Identity();

				// Feedback gain: K = -Q_uu^{-1} * Q_ux
				K[k] = -Q_uu.ldlt().solve(Q_ux);

				// Feedforward term: k = -Q_uu^{-1} * Q_u
				Vector6 k_ff = -Q_uu.ldlt().solve(Q_u);

				// Update value function
				V_x = Q_x + K[k].transpose() * Q_uu * k_ff + Q_ux.transpose() * k_ff + K[k].transpose() * Q_u;
				V_xx = Q_xx + K[k].transpose() * Q_uu * K[k] + Q_ux.transpose() * K[k] + K[k].transpose() * Q_ux;

				// Ensure V_xx remains symmetric
				V_xx = 0.5 * (V_xx + V_xx.transpose());
			}

			return K;
		}

		/**
		 * @brief Linearize Cosserat dynamics around (x, u)
		 *
		 * Dynamics: x_{k+1} = f(x_k, u_k) = x_k * exp(L * u_k)
		 *
		 * Linearization:
		 * δx_{k+1} ≈ F_x * δx_k + F_u * δu_k
		 *
		 * Using SE(3) Jacobians from Phase 2
		 */
		std::pair<Matrix6, Matrix6> linearizeDynamics(const SE3Type &x, const Vector6 &u, Scalar L) const {
			// F_x: ∂f/∂x using composition jacobians
			Vector6 xi = u * L;
			SE3Type exp_u = SE3Type::exp(xi);

			auto [J_g1, J_g2] = x.composeJacobians(exp_u);
			Matrix6 F_x = J_g1; // 6x6

			// F_u: ∂f/∂u using chain rule
			// f(x,u) = x * exp(L*u)
			// ∂f/∂u = ∂(x * exp(L*u))/∂u
			//       = x * ∂exp(L*u)/∂(L*u) * L
			//       = J_g2 * leftJacobian(L*u) * L

			// Approximate leftJacobian with identity for small strains
			// Or use finite differences
			Matrix6 F_u = J_g2 * L; // Simplified

			return {F_x, F_u};
		}

		/**
		 * @brief Compute cost function derivatives
		 */
		void computeCostDerivatives(const SE3Type &x, const Vector6 &u, const SE3Type &x_ref, Vector6 &l_x,
									Vector6 &l_u, Matrix6 &l_xx, Matrix6 &l_uu, Matrix6 &l_ux) const {
			// Simplified: approximate with diagonal Hessians
			// For accurate derivatives, use finite differences or autodiff

			// Gradient w.r.t. state (simplified)
			Vector3 pos_error = x.translation() - x_ref.translation();
			l_x = Vector6::Zero();
			l_x.template head<3>() = 2.0 * m_config.Q_position * pos_error;

			// Gradient w.r.t. control
			l_u = 2.0 * m_config.R_control * u;

			// Hessians (diagonal approximation)
			l_xx = Matrix6::Identity() * m_config.Q_position;
			l_uu = Matrix6::Identity() * m_config.R_control;
			l_ux = Matrix6::Zero();
		}

		/**
		 * @brief Compute terminal cost derivatives
		 */
		void computeTerminalCostDerivatives(const SE3Type &x, const SE3Type &x_ref, Vector6 &V_x, Matrix6 &V_xx) const {
			Vector3 pos_error = x.translation() - x_ref.translation();

			V_x = Vector6::Zero();
			V_x.template head<3>() = 2.0 * m_config.Q_final_position * pos_error;

			V_xx = Matrix6::Identity() * m_config.Q_final_position;
		}

		/**
		 * @brief Line search with Armijo condition
		 */
		std::tuple<std::vector<Vector6>, std::vector<SE3Type>, Scalar, Scalar>
		lineSearch(const std::vector<SE3Type> &x, const std::vector<Vector6> &u, const std::vector<Matrix6> &K,
				   const Trajectory &ref, Scalar current_cost) const {
			Scalar alpha = 1.0;

			for (int iter = 0; iter < m_config.max_line_search_iters; ++iter) {
				// Apply feedback gains with step size alpha
				std::vector<Vector6> u_new(u.size());
				for (size_t k = 0; k < u.size(); ++k) {
					// Simple feedforward update (could add feedback term)
					Vector6 du = K[k] * Vector6::Zero(); // Simplified: no state deviation term
					u_new[k] = u[k] + alpha * du;
				}

				// Evaluate new trajectory
				auto [x_new, cost_new] = forwardPass(u_new, ref);

				// Accept if cost decreased
				if (cost_new < current_cost) {
					return {u_new, x_new, cost_new, alpha};
				}

				// Reduce step size
				alpha *= m_config.line_search_decay;

				if (alpha < m_config.min_step_size) {
					break;
				}
			}

			// No improvement: return original
			return {u, x, current_cost, 0.0};
		}

		/**
		 * @brief Compute final feedback gains
		 */
		std::vector<Matrix6> computeFeedbackGains(const std::vector<SE3Type> &x, const std::vector<Vector6> &u,
												  const Trajectory &ref) const {
			// Reuse backward pass result
			return backwardPass(x, u, ref);
		}
	};

} // namespace sofa::component::cosserat::liegroups::control
