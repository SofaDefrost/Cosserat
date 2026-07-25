#pragma once

#include <Cosserat/config.h>
#include <Eigen/Dense>
#include <liegroups/CosseratUncertaintyPropagator.h>
#include <liegroups/GaussianOnManifold.h>
#include <liegroups/SE3.h>

namespace Cosserat::mapping {

	using namespace sofa::component::cosserat::liegroups;

	/**
	 * @brief Kalman filter for beam state estimation on the Lie group SE(3).
	 *
	 * Implements a Left-Invariant EKF (LI-EKF) on SE(3):
	 *   - predict(): propagates mean + covariance via CosseratUncertaintyPropagator
	 *     (physically correct: uses expCosserat + tangent-exp covariance transport)
	 *   - update(): Left-Invariant EKF measurement update for a direct pose observation
	 *
	 * ## Coordinate convention
	 * Error is defined right-invariantly: ε = log(X⁻¹ · X_true) ∈ se(3).
	 * Covariance Σ is in the tangent space at the current mean X.
	 */
	class BeamStateEstimator {
	public:
		using SE3Type         = SE3<double>;
		using TangentVector   = typename SE3Type::TangentVector;
		using CovarianceMatrix = Eigen::Matrix<double, 6, 6>;
		using GaussianSE3     = GaussianOnManifold<SE3Type>;
		using Propagator      = CosseratUncertaintyPropagator<double>;

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
		 * @brief Predict step — propagate mean and covariance through one rod section.
		 *
		 * Uses CosseratUncertaintyPropagator::propagateStep() for physically correct
		 * Cosserat kinematics:
		 *
		 *   μ_{k+1} = μ_k · expCosserat(strain, length)
		 *
		 *   Σ_{k+1} = Ad_{g⁻¹} · Σ_k · Ad_{g⁻¹}ᵀ   (transport of existing uncertainty)
		 *           + T(length, strain) · Σ_ξ · T(length, strain)ᵀ  (injected strain noise)
		 *
		 * where T(L, ξ) = ∫₀ᴸ exp(s · ad_ξ) ds is the tangent-exponential matrix.
		 *
		 * This replaces the previous incorrect `computeExp(strain * dt)` which:
		 *   1. Omitted the Cosserat elongation correction (+1 in ρx)
		 *   2. Did not propagate covariance through the adjoint transport
		 *
		 * @param strain       Nominal strain vector ξ ∈ se(3) for this element
		 * @param length       Arc-length of the element (section length, metres)
		 * @param strain_cov   Strain noise covariance Σ_ξ (6×6, units: [1/m]²)
		 */
		void predict(const TangentVector &strain, double length,
		             const CovarianceMatrix &strain_cov) {
			Propagator::Section sec;
			sec.strain     = strain;
			sec.length     = length;
			sec.strain_cov = strain_cov;
			pose_estimate_ = Propagator::propagateStep(pose_estimate_, sec);
		}

		/**
		 * @brief Predict along a full sequence of rod sections (batch predict).
		 *
		 * Convenience wrapper: chains N predict() calls from root to tip.
		 * Returns all N+1 intermediate Gaussian distributions.
		 *
		 * @param strains    Per-section strain vectors (size N).
		 * @param lengths    Per-section arc-lengths (size N).
		 * @param strain_cov Strain noise covariance (shared for all sections).
		 * @return           N+1 Gaussians: [base, after sec 0, …, after sec N-1].
		 */
		[[nodiscard]] std::vector<GaussianSE3> predictAlongRod(
				const std::vector<TangentVector>& strains,
				const std::vector<double>&        lengths,
				const CovarianceMatrix&           strain_cov) const
		{
			std::vector<Propagator::Section> sections(strains.size());
			for (size_t k = 0; k < strains.size(); ++k) {
				sections[k].strain     = strains[k];
				sections[k].length     = lengths[k];
				sections[k].strain_cov = strain_cov;
			}
			return Propagator::propagateAlongRod(pose_estimate_, sections);
		}

		/**
		 * @brief Update step — right-invariant EKF measurement update.
		 *
		 * Direct pose observation: Z = X_true (pose measured by a sensor).
		 * Right-invariant error: ε = log(X⁻¹ · Z) ∈ se(3).
		 *
		 *   Innovation:   r = log(X⁻¹ · Z)
		 *   Kalman gain:  K = P · (P + R)⁻¹
		 *   Correction:   X_new = X · exp(K · r)
		 *   Cov update:   P_new = (I − K) · P
		 *
		 * @param measurement        Observed pose Z ∈ SE(3).
		 * @param measurement_noise  Measurement noise covariance R (6×6).
		 */
		void update(const SE3Type &measurement, const CovarianceMatrix &measurement_noise) {
			const CovarianceMatrix &P = pose_estimate_.getCovariance();
			const CovarianceMatrix  S = P + measurement_noise;       // Innovation cov
			const CovarianceMatrix  K = P * S.inverse();              // Kalman gain

			// Innovation in tangent space: r = log(X⁻¹ · Z)
			const SE3Type &X = pose_estimate_.getMean();
			const TangentVector r = X.computeInverse().compose(measurement).computeLog();

			// Mean update: X_new = X · exp(K · r)
			const SE3Type X_new = X.compose(SE3Type::computeExp(K * r));

			// Covariance update: P_new = (I − K) · P
			const CovarianceMatrix P_new =
				(CovarianceMatrix::Identity() - K) * P;

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
