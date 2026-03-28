# GaussianOnManifold Implementation

## Overview

`GaussianOnManifold<GroupType>` implements a Gaussian distribution on a Lie group manifold. This class represents probability distributions over Lie groups, which is essential for uncertainty propagation, state estimation, and sensor fusion in robotics and computer vision applications.

## Mathematical Background

A Gaussian distribution on a Lie group is defined by:
- **Mean element**: μ ∈ G (a group element)
- **Covariance matrix**: Σ ∈ ℝ^{d×d} (in the tangent space at μ, where d is the group dimension)

The distribution represents the probability density:
```
p(g) ∝ exp(-½ (log(μ⁻¹g))ᵀ Σ⁻¹ (log(μ⁻¹g)))
```

where log is the logarithmic map from the group to its Lie algebra.

## Implementation Details

The `GaussianOnManifold` class is implemented as a template with one parameter:
- `GroupType`: The Lie group type (e.g., SE3, SO3, SE2)

### Key Methods

```cpp
// Constructors
GaussianOnManifold(); // Identity mean, identity covariance
GaussianOnManifold(const GroupType& mean, const CovarianceMatrix& covariance);

// Accessors
const GroupType& getMean() const;
void setMean(const GroupType& mean);
const CovarianceMatrix& getCovariance() const;
void setCovariance(const CovarianceMatrix& covariance);

// Uncertainty propagation
GaussianOnManifold transform(const GroupType& transform) const; // Apply transformation
GaussianOnManifold compose(const GaussianOnManifold& other) const; // Compose with another Gaussian
GaussianOnManifold inverse() const; // Inverse distribution
```

## Mathematical Operations

### Transformation Propagation

Given a random variable X ~ N(μ, Σ) and a transformation Y = T * X, the resulting distribution is:
- Mean: μ' = T * μ
- Covariance: Σ' = Adjoint(T) * Σ * Adjoint(T)ᵀ

### Composition of Gaussians

For independent random variables X ~ N(μ₁, Σ₁) and Y ~ N(μ₂, Σ₂), Z = X * Y has distribution:
- Mean: μ' = μ₁ * μ₂
- Covariance: Σ' = Σ₁ + Adjoint(μ₁) * Σ₂ * Adjoint(μ₁)ᵀ

### Inverse Distribution

For Y = X⁻¹ where X ~ N(μ, Σ):
- Mean: μ' = μ⁻¹
- Covariance: Σ' = Adjoint(μ⁻¹) * Σ * Adjoint(μ⁻¹)ᵀ

## Example Usage

### Basic Usage

```cpp
#include <Cosserat/liegroups/GaussianOnManifold.h>
#include <iostream>

int main() {
    // Create a Gaussian distribution on SE(3)
    Cosserat::SE3<double> mean(
        Cosserat::SO3<double>(Eigen::AngleAxisd(0.1, Eigen::Vector3d::UnitZ())),
        Eigen::Vector3d(1.0, 2.0, 3.0)
    );

    // Create covariance matrix (6x6 for SE(3))
    Eigen::Matrix<double, 6, 6> covariance = Eigen::Matrix<double, 6, 6>::Identity() * 0.01;

    Cosserat::GaussianOnManifold<Cosserat::SE3<double>> distribution(mean, covariance);

    // Access components
    const auto& mean_element = distribution.getMean();
    const auto& cov_matrix = distribution.getCovariance();

    std::cout << "Mean pose: " << mean_element << "\n";
    std::cout << "Covariance:\n" << cov_matrix << "\n";

    return 0;
}
```

### Uncertainty Propagation

```cpp
#include <Cosserat/liegroups/GaussianOnManifold.h>

int main() {
    // Initial pose distribution
    Cosserat::SE3<double> initial_mean = Cosserat::SE3<double>::identity();
    Eigen::Matrix<double, 6, 6> initial_cov = Eigen::Matrix<double, 6, 6>::Identity() * 0.1;
    Cosserat::GaussianOnManifold<Cosserat::SE3<double>> pose_dist(initial_mean, initial_cov);

    // Apply a transformation (e.g., motion command)
    Cosserat::SE3<double> motion(
        Cosserat::SO3<double>(Eigen::AngleAxisd(M_PI/4, Eigen::Vector3d::UnitZ())),
        Eigen::Vector3d(1.0, 0.0, 0.0)
    );

    // Propagate uncertainty through the motion
    auto predicted_dist = pose_dist.transform(motion);

    std::cout << "Predicted mean: " << predicted_dist.getMean() << "\n";
    std::cout << "Predicted covariance trace: " << predicted_dist.getCovariance().trace() << "\n";

    return 0;
}
```

### Sensor Fusion

```cpp
#include <Cosserat/liegroups/GaussianOnManifold.h>

int main() {
    // Prior distribution
    Cosserat::SE3<double> prior_mean = Cosserat::SE3<double>::identity();
    Eigen::Matrix<double, 6, 6> prior_cov = Eigen::Matrix<double, 6, 6>::Identity() * 1.0;
    Cosserat::GaussianOnManifold<Cosserat::SE3<double>> prior(prior_mean, prior_cov);

    // Measurement distribution (relative to prior)
    Cosserat::SE3<double> measurement(
        Cosserat::SO3<double>(Eigen::AngleAxisd(0.05, Eigen::Vector3d::UnitZ())),
        Eigen::Vector3d(0.9, 0.1, 0.05)
    );
    Eigen::Matrix<double, 6, 6> measurement_cov = Eigen::Matrix<double, 6, 6>::Identity() * 0.01;
    Cosserat::GaussianOnManifold<Cosserat::SE3<double>> measurement_dist(measurement, measurement_cov);

    // Fuse measurements (simplified composition)
    auto posterior = prior.compose(measurement_dist);

    std::cout << "Fused mean: " << posterior.getMean() << "\n";
    std::cout << "Fused covariance trace: " << posterior.getCovariance().trace() << "\n";

    return 0;
}
```

### Kalman Filter Prediction

```cpp
#include <Cosserat/liegroups/GaussianOnManifold.h>

class LieGroupKalmanFilter {
private:
    Cosserat::GaussianOnManifold<Cosserat::SE3<double>> state_;

public:
    void predict(const Cosserat::SE3<double>& motion, const Eigen::Matrix<double, 6, 6>& motion_cov) {
        // Create motion distribution
        Cosserat::GaussianOnManifold<Cosserat::SE3<double>> motion_dist(
            motion, motion_cov
        );

        // Predict new state
        state_ = state_.compose(motion_dist);
    }

    void update(const Cosserat::SE3<double>& measurement, const Eigen::Matrix<double, 6, 6>& measurement_cov) {
        // Create measurement distribution
        Cosserat::GaussianOnManifold<Cosserat::SE3<double>> measurement_dist(
            measurement, measurement_cov
        );

        // Simplified update (would need full Kalman update in practice)
        // This is just a demonstration
        auto innovation = state_.getMean().inverse().compose(measurement);
        // ... full update logic would go here
    }

    const auto& getState() const { return state_; }
};
```

## Applications

Gaussian distributions on manifolds are essential for:

- **State estimation**: Extended Kalman filters on Lie groups
- **SLAM**: Simultaneous Localization and Mapping with uncertainty
- **Sensor fusion**: Combining measurements from multiple sensors
- **Motion planning**: Planning under uncertainty
- **Control**: Stochastic control on Lie groups
- **Computer vision**: Pose estimation with uncertainty

## Implementation Notes

### Covariance Representation

The covariance is stored in the tangent space at the mean. This choice provides:
- **Linearity**: The tangent space is a vector space, enabling linear algebra operations
- **Locality**: Captures uncertainty near the mean element
- **Computational efficiency**: Avoids complex operations on the manifold

### Numerical Stability

The implementation includes considerations for:
- **Positive definiteness**: Ensuring covariance matrices remain positive definite
- **Symmetry**: Maintaining symmetry of covariance matrices
- **Adjoint accuracy**: Using numerically stable adjoint computations

### Memory Layout

- **Mean**: Stored as the group element type
- **Covariance**: Stored as Eigen matrix in column-major order for efficiency

## Best Practices

1. **Initialize covariances appropriately** based on your sensor characteristics
2. **Check covariance positive definiteness** after operations if numerical issues are suspected
3. **Use the correct group operations** for uncertainty propagation
4. **Consider the validity region** of the tangent space approximation
5. **Normalize group elements** when setting means to avoid numerical drift
6. **Scale covariances appropriately** for the units and scales in your application

## Limitations

- **Local approximation**: The tangent space approximation is only valid near the mean
- **Linear operations**: Complex nonlinear operations may require more sophisticated uncertainty propagation
- **Computational cost**: Adjoint computations can be expensive for high-dimensional groups
- **Memory usage**: Covariance matrices scale quadratically with group dimension</content>
</xai:function_call/>
<xai:function_call name="edit">
<parameter name="filePath">src/liegroups/docs/advanced_topics.md