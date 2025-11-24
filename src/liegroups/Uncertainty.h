// Uncertainty propagation for Lie groups
// This file provides classes for representing and propagating uncertainty
// on Lie group manifolds, essential for state estimation and Kalman filtering.

#pragma once

#include <Eigen/Dense>
#include <functional>
#include <type_traits>
#include "Types.h"
#include "LieGroupBase.h"
#include "SO3.h"
#include "SE3.h"

namespace sofa::component::cosserat::liegroups {

/**
 * @brief Gaussian distribution on a Lie group manifold
 *
 * Represents X ~ N(μ, Σ) where μ ∈ G is the mean and Σ ∈ ℝᵐˣᵐ is the
 * covariance in the tangent space at μ.
 *
 * This class enables uncertainty propagation through Lie group operations,
 * which is essential for state estimation algorithms like ESKF, UKF, etc.
 *
 * @tparam LieGroupType The Lie group type (e.g., SO3d, SE3d)
 */
template<typename LieGroupType>
class GaussianOnManifold {
public:
    using Scalar = typename LieGroupType::Scalar;
    using TangentVector = typename LieGroupType::TangentVector;
    using Covariance = typename LieGroupType::AdjointMatrix;

    /**
     * @brief Construct a Gaussian distribution
     * @param mean The mean element μ ∈ G
     * @param covariance The covariance Σ ∈ ℝᵐˣᵐ in the tangent space at μ
     */
    GaussianOnManifold(const LieGroupType& mean, const Covariance& covariance)
        : m_mean(mean), m_covariance(covariance) {}

    /**
     * @brief Get the mean of the distribution
     * @return The mean element μ ∈ G
     */
    const LieGroupType& mean() const { return m_mean; }

    /**
     * @brief Get the covariance matrix
     * @return The covariance Σ ∈ ℝᵐˣᵐ
     */
    const Covariance& covariance() const { return m_covariance; }

    /**
     * @brief Transform covariance from local to global frame
     *
     * Using: Σ_global = Ad_X Σ_local Ad_X^T
     *
     * @return The covariance in the global tangent space
     */
    Covariance toGlobalFrame() const {
        auto Ad = m_mean.adjoint();
        return Ad * m_covariance * Ad.transpose();
    }

    /**
     * @brief Transform covariance from global to local frame
     *
     * Using: Σ_local = Ad_X⁻¹ Σ_global (Ad_X⁻¹)^T
     *
     * @param X The group element at which to compute the transformation
     * @param global_cov The covariance in global frame
     * @return The covariance in local frame at X
     */
    static Covariance toLocalFrame(const LieGroupType& X, const Covariance& global_cov) {
        auto Ad_inv = X.inverse().adjoint();
        return Ad_inv * global_cov * Ad_inv.transpose();
    }

    /**
     * @brief Propagate uncertainty through a function f: G → H
     *
     * Uses first-order approximation: Σ_Y ≈ J Σ_X J^T
     * where J = Df(X)/DX is the Jacobian of f at X.
     *
     * @tparam OutputLieGroupType The output Lie group type
     * @param input The input Gaussian distribution
     * @param func The function f: G → H
     * @param jacobian The Jacobian matrix J = Df(X)/DX
     * @return The propagated Gaussian distribution on H
     */
    template<typename OutputLieGroupType>
    static GaussianOnManifold<OutputLieGroupType> propagate(
        const GaussianOnManifold<LieGroupType>& input,
        std::function<OutputLieGroupType(const LieGroupType&)> func,
        const typename OutputLieGroupType::AdjointMatrix& jacobian) {

        auto output_mean = func(input.mean());
        auto output_cov = jacobian * input.covariance() * jacobian.transpose();

        return GaussianOnManifold<OutputLieGroupType>(output_mean, output_cov);
    }

    /**
     * @brief Propagate through composition: Y = X ⊕ δ
     *
     * This is useful for prediction steps in Kalman filters.
     *
     * @param delta The tangent space increment
     * @param delta_cov The covariance of the increment
     * @return The propagated distribution
     */
    GaussianOnManifold<LieGroupType> composeWith(
        const TangentVector& delta,
        const Covariance& delta_cov) const {

        // New mean: μ' = μ ⊕ δ
        auto new_mean = m_mean.plus(delta);

        // New covariance: Σ' = F Σ F^T + G Σ_δ G^T
        // where F = Ad_exp(δ)⁻¹, G = Jr(δ)
        auto F = LieGroupType::exp(delta).inverse().adjoint();
        auto G = LieGroupType::rightJacobian(delta);

        auto new_cov = F * m_covariance * F.transpose() + G * delta_cov * G.transpose();

        return GaussianOnManifold<LieGroupType>(new_mean, new_cov);
    }

    /**
     * @brief Compute the Mahalanobis distance to another distribution
     *
     * This measures how many standard deviations apart the means are,
     * accounting for the covariance structure.
     *
     * @param other Another Gaussian distribution
     * @return The Mahalanobis distance
     */
    Scalar mahalanobisDistance(const GaussianOnManifold<LieGroupType>& other) const {
        // Compute difference in tangent space
        TangentVector diff = m_mean.minus(other.m_mean);

        // Mahalanobis distance: sqrt(diff^T Σ⁻¹ diff)
        // For simplicity, use combined covariance
        Covariance combined_cov = (m_covariance + other.m_covariance) * 0.5;
        Eigen::LLT<Covariance> llt(combined_cov);

        if (llt.info() != Eigen::Success) {
            // Fallback to Euclidean distance if matrix is not positive definite
            return diff.norm();
        }

        TangentVector solved = llt.solve(diff);
        return std::sqrt(diff.dot(solved));
    }

private:
    LieGroupType m_mean;
    Covariance m_covariance;
};

// Convenience type aliases
template<typename Scalar>
using GaussianSO3 = GaussianOnManifold<SO3<Scalar>>;

template<typename Scalar>
using GaussianSE3 = GaussianOnManifold<SE3<Scalar>>;

/**
 * @brief Uncertainty propagation utilities
 */
class UncertaintyPropagation {
public:
    /**
     * @brief Create an isotropic Gaussian distribution
     * @tparam LieGroupType The Lie group type
     * @param mean The mean element
     * @param std_dev The standard deviation for all dimensions
     * @return A Gaussian distribution with diagonal covariance
     */
    template<typename LieGroupType>
    static GaussianOnManifold<LieGroupType> isotropic(
        const LieGroupType& mean,
        typename LieGroupType::Scalar std_dev) {

        using Covariance = typename LieGroupType::AdjointMatrix;
        Covariance cov = Covariance::Identity() * (std_dev * std_dev);
        return GaussianOnManifold<LieGroupType>(mean, cov);
    }

    /**
     * @brief Create a Gaussian with specified diagonal covariance
     * @tparam LieGroupType The Lie group type
     * @param mean The mean element
     * @param variances The variances for each dimension
     * @return A Gaussian distribution with diagonal covariance
     */
    template<typename LieGroupType>
    static GaussianOnManifold<LieGroupType> diagonal(
        const LieGroupType& mean,
        const typename LieGroupType::TangentVector& variances) {

        using Covariance = typename LieGroupType::AdjointMatrix;
        Covariance cov = variances.asDiagonal();
        return GaussianOnManifold<LieGroupType>(mean, cov);
    }
};

} // namespace sofa::component::cosserat::liegroups