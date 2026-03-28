#pragma once

#include <Eigen/Dense>
#include <iostream>

namespace sofa::component::cosserat::liegroups {

	/**
	 * @brief Class representing a Gaussian distribution on a Lie group
	 *
	 * The distribution is defined by a mean element on the group and a covariance matrix
	 * in the tangent space of the mean.
	 *
	 * @tparam GroupType The Lie group type (e.g., SE3, SO3)
	 */
	template<typename GroupType>
	class GaussianOnManifold {
	public:
		using Scalar = typename GroupType::Scalar;
		using TangentVector = typename GroupType::TangentVector;
		using CovarianceMatrix = Eigen::Matrix<Scalar, GroupType::DoF, GroupType::DoF>;

		/**
		 * @brief Default constructor
		 */
		GaussianOnManifold() : mean_(GroupType::computeIdentity()), covariance_(CovarianceMatrix::Identity()) {}

		/**
		 * @brief Constructor with mean and covariance
		 * @param mean Mean element on the group
		 * @param covariance Covariance matrix in the tangent space
		 */
		GaussianOnManifold(const GroupType &mean, const CovarianceMatrix &covariance) :
			mean_(mean), covariance_(covariance) {}

		/**
		 * @brief Get the mean element
		 * @return The mean element
		 */
		const GroupType &getMean() const { return mean_; }

		/**
		 * @brief Set the mean element
		 * @param mean The new mean element
		 */
		void setMean(const GroupType &mean) { mean_ = mean; }

		/**
		 * @brief Get the covariance matrix
		 * @return The covariance matrix
		 */
		const CovarianceMatrix &getCovariance() const { return covariance_; }

		/**
		 * @brief Set the covariance matrix
		 * @param covariance The new covariance matrix
		 */
		void setCovariance(const CovarianceMatrix &covariance) { covariance_ = covariance; }

		/**
		 * @brief Propagate the uncertainty through a transformation
		 *
		 * Given a transformation T such that Y = T * X, where X ~ N(mean, cov),
		 * this computes the approximate distribution of Y.
		 *
		 * @param transform The transformation applied to the random variable
		 * @return The transformed Gaussian
		 */
		GaussianOnManifold transform(const GroupType &transform) const {
			// Y = T * X
			// Mean: T * mean
			// Covariance: Adjoint(T) * cov * Adjoint(T)^T

			GroupType new_mean = transform * mean_;
			typename GroupType::AdjointMatrix adj = transform.computeAdjoint();
			CovarianceMatrix new_cov = adj * covariance_ * adj.transpose();

			return GaussianOnManifold(new_mean, new_cov);
		}

		/**
		 * @brief Compose with another Gaussian (approximate)
		 *
		 * Computes the distribution of Z = X * Y, where X ~ this and Y ~ other.
		 * Assumes independence.
		 *
		 * @param other The other Gaussian distribution
		 * @return The composed Gaussian
		 */
		GaussianOnManifold compose(const GaussianOnManifold &other) const {
			// Z = X * Y
			// Mean: mean * other.mean
			// Covariance: cov + Adjoint(mean) * other.cov * Adjoint(mean)^T

			GroupType new_mean = mean_ * other.mean_;
			typename GroupType::AdjointMatrix adj = mean_.computeAdjoint();
			CovarianceMatrix new_cov = covariance_ + adj * other.covariance_ * adj.transpose();

			return GaussianOnManifold(new_mean, new_cov);
		}

		/**
		 * @brief Inverse of the Gaussian (approximate)
		 *
		 * Computes the distribution of Y = X^-1
		 *
		 * @param other The other Gaussian distribution
		 * @return The composed Gaussian
		 */
		GaussianOnManifold inverse() const {
			// Y = X^-1
			// Mean: mean^-1
			// Covariance: Adjoint(mean^-1) * cov * Adjoint(mean^-1)^T

			GroupType new_mean = mean_.computeInverse();
			typename GroupType::AdjointMatrix adj = new_mean.computeAdjoint();
			CovarianceMatrix new_cov = adj * covariance_ * adj.transpose();

			return GaussianOnManifold(new_mean, new_cov);
		}

	private:
		GroupType mean_;
		CovarianceMatrix covariance_;
	};

} // namespace sofa::component::cosserat::liegroups
