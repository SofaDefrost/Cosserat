#pragma once

#include <Cosserat/config.h>
#include <Eigen/Dense>
#include <liegroups/GaussianOnManifold.h>
#include <liegroups/SE3.h>

namespace Cosserat::mapping {

	using namespace sofa::component::cosserat::liegroups;

	/**
	 * @brief Class for estimating the state of a beam using Kalman filtering on Lie groups
	 */
	class BeamStateEstimator {
	public:
		using SE3Type = SE3<double>;
		using TangentVector = typename SE3Type::TangentVector;
		using CovarianceMatrix = Eigen::Matrix<double, 6, 6>;
		using GaussianSE3 = GaussianOnManifold<SE3Type>;

		/**
		 * @brief Default constructor
		 */
		BeamStateEstimator() = default;

		/**
		 * @brief Initialize the estimator
		 * @param initial_pose Initial pose mean
		 * @param initial_covariance Initial covariance
		 */
		void initialize(const SE3Type &initial_pose, const CovarianceMatrix &initial_covariance) {
			pose_estimate_ = GaussianSE3(initial_pose, initial_covariance);
		}

		/**
		 * @brief Predict step (Process model)
		 *
		 * Propagates the state based on strain (control input) and process noise.
		 * X_{k+1} = X_k * Exp(strain * dt)
		 *
		 * @param strain Strain vector (control input)
		 * @param dt Time step or length segment
		 * @param process_noise Process noise covariance (Q)
		 */
		void predict(const TangentVector &strain, double dt, const CovarianceMatrix &process_noise) {
			// Control input transformation
			SE3Type control_transform = SE3Type::computeExp(strain * dt);

			// Create a Gaussian for the control input with process noise
			GaussianSE3 control_input(control_transform, process_noise);

			// Compose current estimate with control input
			// X_{k+1} = X_k * U_k
			pose_estimate_ = pose_estimate_.compose(control_input);
		}

		/**
		 * @brief Update step (Measurement model)
		 *
		 * Updates the state based on a pose measurement.
		 * Z = X * V, where V is measurement noise
		 *
		 * @param measurement Measured pose
		 * @param measurement_noise Measurement noise covariance (R)
		 */
		void update(const SE3Type &measurement, const CovarianceMatrix &measurement_noise) {
			// Standard EKF update on manifold is complex.
			// Here we implement a simplified version or "Left-Invariant EKF" style update
			// if we assume the error is defined as eta = X^-1 * X_true

			// Kalman Gain computation
			// P = current covariance
			// H = Identity (direct measurement of pose)
			// K = P * (P + R)^-1

			CovarianceMatrix P = pose_estimate_.getCovariance();
			CovarianceMatrix R = measurement_noise;
			CovarianceMatrix S = P + R; // Innovation covariance
			CovarianceMatrix K = P * S.inverse(); // Kalman Gain

			// Innovation (Residual) in tangent space
			// r = log(X^-1 * Z)
			SE3Type X = pose_estimate_.getMean();
			SE3Type Z = measurement;
			TangentVector r = (X.computeInverse() * Z).computeLog();

			// Correction
			TangentVector correction = K * r;

			// Update Mean
			// X_new = X * Exp(correction)
			SE3Type X_new = X * SE3Type::computeExp(correction);

			// Update Covariance
			// P_new = (I - K) * P
			CovarianceMatrix I = CovarianceMatrix::Identity();
			CovarianceMatrix P_new = (I - K) * P;

			pose_estimate_ = GaussianSE3(X_new, P_new);
		}

		/**
		 * @brief Get the current pose estimate
		 * @return The estimated Gaussian pose
		 */
		const GaussianSE3 &getEstimate() const { return pose_estimate_; }

		/**
		 * @brief Get the estimation confidence (trace of covariance)
		 * @return A scalar representing uncertainty (lower is better)
		 */
		double getEstimationConfidence() const { return pose_estimate_.getCovariance().trace(); }

	private:
		GaussianSE3 pose_estimate_;
	};

} // namespace Cosserat::mapping
