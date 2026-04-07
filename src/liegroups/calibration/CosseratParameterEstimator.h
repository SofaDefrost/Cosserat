/******************************************************************************
 * Cosserat Parameter Estimator - Calibrate physical parameters from measurements
 ******************************************************************************/

#pragma once

#include <vector>
#include <functional>
#include "../SE3.h"
#include "../Types.h"

namespace sofa::component::cosserat::liegroups::calibration {

/**
 * @brief Estimate physical parameters of Cosserat rod from measurements
 * 
 * Given measurements of (strain, position) pairs, estimates physical parameters:
 * - E: Young's modulus (rigidity for bending/elongation)
 * - G: Shear modulus (rigidity for torsion/shear)
 * - I_y, I_z: Second moments of inertia (bending resistance)
 * - A: Cross-sectional area (axial/shear resistance)
 * 
 * ## Approach
 * 
 * Uses gradient descent to minimize discrepancy between measured and predicted positions.
 * Parameters are represented as scaling factors from nominal values.
 * 
 * @tparam Scalar Scalar type
 */
template<typename Scalar = double>
class CosseratParameterEstimator {
public:
    using SE3Type = SE3<Scalar>;
    using Vector3 = typename SE3Type::Vector3;
    using Vector6 = Eigen::Matrix<Scalar, 6, 1>;
    using VectorX = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    
    /**
     * @brief Measurement data point
     */
    struct Measurement {
        std::vector<Vector6> strains;    // Applied strains
        Vector3 measured_position;       // Measured tip position
        Scalar segment_length;           // Segment length
        Scalar weight = 1.0;             // Measurement weight
    };
    
    /**
     * @brief Configuration
     */
    struct Config {
        int max_iterations = 100;
        Scalar learning_rate = 0.01;
        Scalar convergence_threshold = 1e-5;
        Scalar regularization = 0.001;    // Regularize parameters near 1.0
        
        // Parameter bounds
        Scalar min_scale = 0.1;
        Scalar max_scale = 10.0;
        
        bool verbose = false;
    };
    
    /**
     * @brief Estimated parameters
     */
    struct Parameters {
        Scalar E_scale = 1.0;     // Young's modulus scaling
        Scalar G_scale = 1.0;     // Shear modulus scaling
        Scalar I_scale = 1.0;     // Second moment scaling
        Scalar A_scale = 1.0;     // Area scaling
        
        VectorX toVector() const {
            VectorX v(4);
            v << E_scale, G_scale, I_scale, A_scale;
            return v;
        }
        
        void fromVector(const VectorX& v) {
            E_scale = v[0];
            G_scale = v[1];
            I_scale = v[2];
            A_scale = v[3];
        }
    };
    
    /**
     * @brief Estimation result
     */
    struct Result {
        Parameters parameters;
        std::vector<Scalar> cost_history;
        Scalar final_cost;
        Scalar rmse;              // Root mean squared error
        int iterations;
        bool converged;
        std::string message;
    };
    
    explicit CosseratParameterEstimator(const Config& config = Config())
        : m_config(config)
    {}
    
    /**
     * @brief Estimate parameters from measurements
     * 
     * @param measurements Training data
     * @param initial_guess Initial parameter guess
     * @return Estimated parameters
     */
    Result estimate(
        const std::vector<Measurement>& measurements,
        const Parameters& initial_guess = Parameters()
    ) {
        Result result;
        result.converged = false;
        result.iterations = 0;
        
        Parameters params = initial_guess;
        
        // Initial cost
        Scalar cost = computeCost(params, measurements);
        result.cost_history.push_back(cost);
        
        if (m_config.verbose) {
            std::cout << "Calibration Iteration 0: Cost = " << cost << "\n";
        }
        
        // Gradient descent
        for (int iter = 0; iter < m_config.max_iterations; ++iter) {
            result.iterations = iter + 1;
            
            // Compute gradient
            VectorX gradient = computeGradient(params, measurements);
            
            // Update parameters
            VectorX params_vec = params.toVector();
            params_vec -= m_config.learning_rate * gradient;
            
            // Clamp to bounds
            params_vec = params_vec.cwiseMax(m_config.min_scale).cwiseMin(m_config.max_scale);
            params.fromVector(params_vec);
            
            // Compute new cost
            Scalar new_cost = computeCost(params, measurements);
            result.cost_history.push_back(new_cost);
            
            if (m_config.verbose && iter % 10 == 0) {
                std::cout << "Calibration Iteration " << iter + 1 
                          << ": Cost = " << new_cost 
                          << ", RMSE = " << std::sqrt(new_cost / measurements.size()) << "\n";
            }
            
            // Check convergence
            Scalar cost_change = std::abs(cost - new_cost);
            if (cost_change < m_config.convergence_threshold) {
                result.converged = true;
                result.message = "Converged: cost change below threshold";
                cost = new_cost;
                break;
            }
            
            cost = new_cost;
        }
        
        result.parameters = params;
        result.final_cost = cost;
        result.rmse = std::sqrt(cost / measurements.size());
        
        if (!result.converged) {
            result.message = "Max iterations reached";
        }
        
        return result;
    }
    
private:
    Config m_config;
    
    /**
     * @brief Compute total cost (MSE + regularization)
     */
    Scalar computeCost(const Parameters& params, const std::vector<Measurement>& measurements) const {
        Scalar total_cost = 0.0;
        
        // Data term: squared error
        for (const auto& m : measurements) {
            Vector3 predicted = predictPosition(params, m.strains, m.segment_length);
            Vector3 error = predicted - m.measured_position;
            total_cost += m.weight * error.squaredNorm();
        }
        
        // Regularization: prefer parameters near 1.0
        VectorX p = params.toVector();
        VectorX deviation = p - VectorX::Ones(4);
        total_cost += m_config.regularization * deviation.squaredNorm();
        
        return total_cost;
    }
    
    /**
     * @brief Compute gradient via finite differences
     */
    VectorX computeGradient(const Parameters& params, const std::vector<Measurement>& measurements) const {
        const Scalar eps = 1e-6;
        VectorX gradient(4);
        
        Scalar cost_0 = computeCost(params, measurements);
        
        for (int i = 0; i < 4; ++i) {
            Parameters params_plus = params;
            VectorX p = params_plus.toVector();
            p[i] += eps;
            params_plus.fromVector(p);
            
            Scalar cost_plus = computeCost(params_plus, measurements);
            gradient[i] = (cost_plus - cost_0) / eps;
        }
        
        return gradient;
    }
    
    /**
     * @brief Predict tip position given parameters and strains
     * 
     * Note: parameter scaling affects strain directly in simplified model
     */
    Vector3 predictPosition(
        const Parameters& params,
        const std::vector<Vector6>& strains,
        Scalar length
    ) const {
        // Apply parameter scaling to strains (simplified model)
        // In reality, parameters affect stiffness matrix, but for estimation
        // we can approximate by scaling the strain directly
        
        SE3Type T = SE3Type::computeIdentity();
        
        for (const auto& strain : strains) {
            // Scale strains by parameters (simplified)
            Vector6 scaled_strain = strain;
            scaled_strain.template head<3>() *= params.I_scale;  // Rotation affected by I
            scaled_strain.template tail<3>() *= params.A_scale;  // Translation affected by A
            
            Vector6 xi = scaled_strain * length;
            T = T * SE3Type::exp(xi);
        }
        
        return T.translation();
    }
};

} // namespace sofa::component::cosserat::liegroups::calibration
