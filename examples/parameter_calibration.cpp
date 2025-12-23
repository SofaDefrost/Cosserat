/******************************************************************************
 * Parameter Calibration Example
 * 
 * Demonstrates using CosseratParameterEstimator to recover physical parameters
 * from synthetic measurement data.
 ******************************************************************************/

#include <iostream>
#include <iomanip>
#include "../src/liegroups/calibration/CosseratParameterEstimator.h"
#include "../src/liegroups/SE3.h"

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
 * @brief Generate synthetic measurement data with known parameters
 */
std::vector<Measurement> generateMeasurements(const Parameters& true_params) {
    std::vector<Measurement> measurements;
    
    // Simulate 20 different configurations
    for (int i = 0; i < 20; ++i) {
        Measurement m;
        m.segment_length = 0.1;  // 10 cm segments
        
        // Configuration 1: Bending in X
        if (i % 4 == 0) {
            Vector6 strain;
            strain << 0.2, 0.0, 0.0, 0.0, 0.0, 0.1;
            m.strains = {strain, strain, strain};
        }
        // Configuration 2: Bending in Y
        else if (i % 4 == 1) {
            Vector6 strain;
            strain << 0.0, 0.15, 0.0, 0.0, 0.0, 0.1;
            m.strains = {strain, strain, strain, strain};
        }
        // Configuration 3: Torsion
        else if (i % 4 == 2) {
            Vector6 strain;
            strain << 0.0, 0.0, 0.3, 0.0, 0.0, 0.1;
            m.strains = {strain, strain};
        }
        // Configuration 4: Mixed
        else {
            Vector6 strain1, strain2;
            strain1 << 0.1, 0.1, 0.05, 0.01, 0.01, 0.12;
            strain2 << -0.05, 0.08, 0.02, 0.0, 0.01, 0.11;
            m.strains = {strain1, strain2, strain1};
        }
        
        // Compute "measured" position using true parameters
        SE3<Scalar> T = SE3<Scalar>::computeIdentity();
        for (const auto& strain : m.strains) {
            Vector6 scaled_strain = strain;
            scaled_strain.template head<3>() *= true_params.I_scale;
            scaled_strain.template tail<3>() *= true_params.A_scale;
            
            Vector6 xi = scaled_strain * m.segment_length;
            T = T * SE3<Scalar>::exp(xi);
        }
        
        m.measured_position = T.translation();
        
        // Add realistic measurement noise (±1mm)
        m.measured_position += 0.001 * Vector3::Random();
        
        measurements.push_back(m);
    }
    
    return measurements;
}

int main() {
    std::cout << "=== Cosserat Parameter Calibration Example ===\n\n";
    
    // Define "true" parameters we want to recover
    Parameters true_params;
    true_params.E_scale = 1.5;   // 50% stiffer than nominal
    true_params.G_scale = 0.8;   // 20% softer than nominal
    true_params.I_scale = 1.3;   // 30% larger moment of inertia
    true_params.A_scale = 1.2;   // 20% larger cross-section
    
    std::cout << "True Parameters:\n";
    std::cout << "  E_scale = " << true_params.E_scale << "\n";
    std::cout << "  G_scale = " << true_params.G_scale << "\n";
    std::cout << "  I_scale = " << true_params.I_scale << "\n";
    std::cout << "  A_scale = " << true_params.A_scale << "\n\n";
    
    // Generate synthetic measurement data
    std::cout << "Generating synthetic measurement data...\n";
    auto measurements = generateMeasurements(true_params);
    std::cout << "  Generated " << measurements.size() << " measurements\n\n";
    
    // Configure estimator
    Config config;
    config.max_iterations = 300;
    config.learning_rate = 0.03;
    config.convergence_threshold = 1e-6;
    config.regularization = 0.005;
    config.verbose = true;
    
    Estimator estimator(config);
    
    // Initial guess (start from nominal parameters)
    Parameters initial_guess;
    initial_guess.E_scale = 1.0;
    initial_guess.G_scale = 1.0;
    initial_guess.I_scale = 1.0;
    initial_guess.A_scale = 1.0;
    
    std::cout << "Starting calibration from nominal parameters...\n";
    std::cout << "Initial guess: all scales = 1.0\n\n";
    
    // Run estimation
    auto result = estimator.estimate(measurements, initial_guess);
    
    // Print results
    std::cout << "\n=== Calibration Results ===\n";
    std::cout << "Status: " << (result.converged ? "CONVERGED" : "MAX ITERATIONS") << "\n";
    std::cout << "Message: " << result.message << "\n";
    std::cout << "Iterations: " << result.iterations << "\n";
    std::cout << "Final RMSE: " << std::fixed << std::setprecision(6) 
              << result.rmse << " m\n\n";
    
    std::cout << "Estimated Parameters:\n";
    std::cout << "  E_scale = " << std::setprecision(4) << result.parameters.E_scale;
    std::cout << "  (error: " << std::setprecision(2) 
              << 100.0 * std::abs(result.parameters.E_scale - true_params.E_scale) / true_params.E_scale 
              << "%)\n";
    
    std::cout << "  G_scale = " << std::setprecision(4) << result.parameters.G_scale;
    std::cout << "  (error: " << std::setprecision(2)
              << 100.0 * std::abs(result.parameters.G_scale - true_params.G_scale) / true_params.G_scale 
              << "%)\n";
    
    std::cout << "  I_scale = " << std::setprecision(4) << result.parameters.I_scale;
    std::cout << "  (error: " << std::setprecision(2)
              << 100.0 * std::abs(result.parameters.I_scale - true_params.I_scale) / true_params.I_scale 
              << "%)\n";
    
    std::cout << "  A_scale = " << std::setprecision(4) << result.parameters.A_scale;
    std::cout << "  (error: " << std::setprecision(2)
              << 100.0 * std::abs(result.parameters.A_scale - true_params.A_scale) / true_params.A_scale 
              << "%)\n\n";
    
    // Cost evolution
    std::cout << "Cost Evolution:\n";
    std::cout << "  Initial cost: " << std::setprecision(6) << result.cost_history[0] << "\n";
    std::cout << "  Final cost:   " << std::setprecision(6) << result.final_cost << "\n";
    std::cout << "  Reduction:    " << std::setprecision(1) 
              << 100.0 * (1.0 - result.final_cost / result.cost_history[0]) << "%\n\n";
    
    return 0;
}
