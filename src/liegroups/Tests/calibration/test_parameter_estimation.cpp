/******************************************************************************
 * Tests for Cosserat Parameter Estimator
 ******************************************************************************/

#include <gtest/gtest.h>
#include "../../calibration/CosseratParameterEstimator.h"

using namespace sofa::component::cosserat::liegroups;
using namespace sofa::component::cosserat::liegroups::calibration;

using Scalar = double;
using Estimator = CosseratParameterEstimator<Scalar>;
using Measurement = typename Estimator::Measurement;
using Parameters = typename Estimator::Parameters;
using Config = typename Estimator::Config;
using Vector3 = Eigen::Vector3d;
using Vector6 = Eigen::Matrix<double, 6, 1>;

/**
 * @brief Generate synthetic measurements with known parameters
 */
std::vector<Measurement> generateSyntheticData(
    const Parameters& true_params,
    int num_measurements,
    Scalar segment_length,
    Scalar noise_level = 0.0
) {
    std::vector<Measurement> measurements;
    measurements.reserve(num_measurements);
    
    // Create diverse strain configurations
    for (int i = 0; i < num_measurements; ++i) {
        Measurement m;
        m.segment_length = segment_length;
        
        // Generate random strains
        int num_segments = 3 + (i % 3);  // 3-5 segments
        for (int j = 0; j < num_segments; ++j) {
            Vector6 strain;
            strain << 0.1 * (i % 5) - 0.2,  // kappa_x: bending
                      0.1 * (i % 3) - 0.1,  // kappa_y: bending
                      0.05 * (i % 4),       // kappa_z: torsion
                      0.02 * (i % 2),       // gamma_x: shear
                      0.02 * (i % 3),       // gamma_y: shear
                      0.1 + 0.05 * (i % 5); // eps_z: elongation
            m.strains.push_back(strain);
        }
        
        // Compute "measured" position with true parameters
        SE3<Scalar> T = SE3<Scalar>::computeIdentity();
        for (const auto& strain : m.strains) {
            Vector6 scaled_strain = strain;
            scaled_strain.template head<3>() *= true_params.I_scale;
            scaled_strain.template tail<3>() *= true_params.A_scale;
            
            Vector6 xi = scaled_strain * segment_length;
            T = T * SE3<Scalar>::exp(xi);
        }
        
        m.measured_position = T.translation();
        
        // Add noise
        if (noise_level > 0.0) {
            m.measured_position += noise_level * Vector3::Random();
        }
        
        measurements.push_back(m);
    }
    
    return measurements;
}

/**
 * @brief Test parameter recovery with perfect data
 */
TEST(CosseratParameterEstimator, RecoverKnownParameters) {
    // True parameters
    Parameters true_params;
    true_params.E_scale = 1.5;
    true_params.G_scale = 0.8;
    true_params.I_scale = 1.2;
    true_params.A_scale = 1.3;
    
    // Generate synthetic data
    auto measurements = generateSyntheticData(true_params, 10, 0.1);
    
    // Configure estimator
    Config config;
    config.max_iterations = 200;
    config.learning_rate = 0.05;
    config.convergence_threshold = 1e-6;
    config.verbose = false;
    
    Estimator estimator(config);
    
    // Start from different initial guess
    Parameters initial_guess;
    initial_guess.E_scale = 1.0;
    initial_guess.G_scale = 1.0;
    initial_guess.I_scale = 1.0;
    initial_guess.A_scale = 1.0;
    
    // Estimate
    auto result = estimator.estimate(measurements, initial_guess);
    
    // Verify convergence
    EXPECT_TRUE(result.converged);
    EXPECT_LT(result.final_cost, 0.01);
    
    // Verify parameters (within 5% error)
    EXPECT_NEAR(result.parameters.I_scale, true_params.I_scale, 0.1);
    EXPECT_NEAR(result.parameters.A_scale, true_params.A_scale, 0.1);
}

/**
 * @brief Test with noisy measurements
 */
TEST(CosseratParameterEstimator, HandleNoisyMeasurements) {
    Parameters true_params;
    true_params.I_scale = 1.5;
    true_params.A_scale = 1.2;
    
    // Generate noisy data
    auto measurements = generateSyntheticData(true_params, 20, 0.1, 0.005);
    
    Config config;
    config.max_iterations = 300;
    config.learning_rate = 0.03;
    config.regularization = 0.01;
    config.verbose = false;
    
    Estimator estimator(config);
    auto result = estimator.estimate(measurements);
    
    // Should still converge, but with higher error
    EXPECT_LT(result.rmse, 0.02);  // RMSE within 2cm
    
    // Parameters should be reasonable
    EXPECT_GT(result.parameters.I_scale, 0.5);
    EXPECT_LT(result.parameters.I_scale, 2.5);
    EXPECT_GT(result.parameters.A_scale, 0.5);
    EXPECT_LT(result.parameters.A_scale, 2.5);
}

/**
 * @brief Test convergence behavior
 */
TEST(CosseratParameterEstimator, ConvergenceMonotonic) {
    Parameters true_params;
    true_params.I_scale = 1.3;
    true_params.A_scale = 1.4;
    
    auto measurements = generateSyntheticData(true_params, 15, 0.1);
    
    Config config;
    config.max_iterations = 100;
    config.learning_rate = 0.02;
    config.verbose = false;
    
    Estimator estimator(config);
    auto result = estimator.estimate(measurements);
    
    // Cost should decrease (mostly monotonic, allowing small increases)
    int num_increases = 0;
    for (size_t i = 1; i < result.cost_history.size(); ++i) {
        if (result.cost_history[i] > result.cost_history[i-1]) {
            num_increases++;
        }
    }
    
    // Allow at most 10% of iterations to increase
    EXPECT_LT(num_increases, result.cost_history.size() / 10);
    
    // Final cost should be lower than initial
    EXPECT_LT(result.final_cost, result.cost_history[0] * 0.5);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
