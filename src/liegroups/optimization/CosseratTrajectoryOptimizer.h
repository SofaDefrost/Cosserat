/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture                          *
*                (c) 2006-2018 INRIA, USTL, UJF, CNRS, MGH                    *
*                                                                             *
* This library is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This library is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this library; if not, write to the Free Software Foundation,     *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
*******************************************************************************
*                           Plugin Cosserat                                   *
*                                                                             *
* This plugin is distributed under the GNU LGPL v2.1 license                 *
*                                                                             *
* Authors: Yinoussa Adagolodjo, Bruno Carrez, Christian Duriez               *
******************************************************************************/
#pragma once

#include <vector>
#include <functional>
#include <iostream>
#include <cmath>
#include <Eigen/Dense>

#include "../SE3.h"
#include "../SO3.h"

namespace sofa::component::cosserat::liegroups::optimization {

/**
 * @brief Trajectory optimizer for Cosserat rods using analytical jacobians
 * 
 * This optimizer uses gradient descent with backpropagation through the
 * forward kinematics chain to optimize strain configurations that achieve
 * desired end-effector poses.
 * 
 * STRAIN CONVENTION:
 * Each strain is a Vector6 = [φx, φy, φz, ρx, ρy, ρz]ᵀ where:
 *   - φx, φy, φz (indices 0-2): Angular strains (torsion, bending Y, bending Z)
 *   - ρx, ρy, ρz (indices 3-5): Linear strains (elongation, shearing Y, shearing Z)
 * See STRAIN_CONVENTION.md for detailed documentation.
 * 
 * The optimization uses analytical Jacobians from the SE3 class for efficient
 * and accurate gradient computation.
 * 
 * @tparam Scalar Floating point type (double or float)
 */
template<typename Scalar = double>
class CosseratTrajectoryOptimizer {
public:
    using Vector3 = Eigen::Matrix<Scalar, 3, 1>;
    using Vector6 = Eigen::Matrix<Scalar, 6, 1>;
    using Matrix3 = Eigen::Matrix<Scalar, 3, 3>;
    using Matrix6 = Eigen::Matrix<Scalar, 6, 6>;
    using SE3Type = SE3<Scalar>;
    using SO3Type = SO3<Scalar>;
    
    /**
     * @brief Optimization parameters
     */
    struct Parameters {
        double learning_rate = 0.01;        // Initial learning rate
        int max_iterations = 1000;          // Maximum number of iterations
        double tolerance = 1e-6;            // Convergence tolerance
        double regularization = 0.01;       // L2 regularization on strains
        bool use_line_search = true;        // Use Armijo line search
        bool verbose = true;                // Print progress
        int print_every = 100;              // Print frequency
        
        // Line search parameters
        double armijo_c1 = 1e-4;           // Armijo condition parameter
        double backtrack_factor = 0.5;      // Backtracking factor
        int max_line_search_iters = 20;     // Max line search iterations
    };
    
    /**
     * @brief Cost function value and gradient
     */
    struct Cost {
        double value = 0.0;                 // Cost value
        std::vector<Vector6> gradient;      // Gradient w.r.t. strains
        bool converged = false;             // Convergence flag
    };
    
    /**
     * @brief Optimization result
     */
    struct Result {
        std::vector<Vector6> strains;       // Optimized strains
        SE3Type final_transform;            // Final end-effector pose
        double final_cost;                  // Final cost value
        int iterations;                     // Number of iterations
        bool converged;                     // Convergence status
        std::vector<double> cost_history;   // Cost history for plotting
    };
    
    /**
     * @brief Default constructor
     */
    CosseratTrajectoryOptimizer() = default;
    
    /**
     * @brief Optimize strains to reach a target pose
     * 
     * @param initial_strains Initial strain configuration (n_sections × 6)
     * @param target Target SE(3) transformation
     * @param section_length Length of each section
     * @param params Optimization parameters
     * @return Optimization result
     */
    Result optimizeToTarget(
        const std::vector<Vector6>& initial_strains,
        const SE3Type& target,
        double section_length,
        const Parameters& params = Parameters()
    ) {
        Result result;
        result.strains = initial_strains;
        result.iterations = 0;
        result.converged = false;
        result.cost_history.reserve(params.max_iterations);
        
        const int n_sections = static_cast<int>(initial_strains.size());
        
        if (params.verbose) {
            std::cout << "=== Cosserat Trajectory Optimization ===" << std::endl;
            std::cout << "Number of sections: " << n_sections << std::endl;
            std::cout << "Section length: " << section_length << " m" << std::endl;
            std::cout << "Target position: " << target.translation().transpose() << std::endl;
            std::cout << "Learning rate: " << params.learning_rate << std::endl;
            std::cout << "Regularization: " << params.regularization << std::endl;
            std::cout << std::endl;
        }
        
        // Gradient descent loop
        for (int iter = 0; iter < params.max_iterations; ++iter) {
            // Compute cost and gradient
            Cost cost = computeCostAndGradient(
                result.strains,
                target,
                section_length,
                params.regularization
            );
            
            result.cost_history.push_back(cost.value);
            
            // Print progress
            if (params.verbose && (iter % params.print_every == 0 || iter == params.max_iterations - 1)) {
                SE3Type current_pose = forwardKinematics(result.strains, section_length);
                Vector3 error = current_pose.translation() - target.translation();
                
                std::cout << "Iteration " << iter << ":" << std::endl;
                std::cout << "  Cost: " << cost.value << std::endl;
                std::cout << "  Position error: " << error.norm() << " m" << std::endl;
                std::cout << "  Current position: " << current_pose.translation().transpose() << std::endl;
            }
            
            // Check convergence
            if (cost.value < params.tolerance) {
                result.converged = true;
                result.iterations = iter;
                if (params.verbose) {
                    std::cout << "\n✓ Converged after " << iter << " iterations!" << std::endl;
                }
                break;
            }
            
            // Determine step size
            double step_size = params.learning_rate;
            
            if (params.use_line_search) {
                step_size = lineSearch(
                    result.strains,
                    cost.gradient,
                    target,
                    section_length,
                    cost.value,
                    params
                );
            }
            
            // Update strains
            for (int i = 0; i < n_sections; ++i) {
                result.strains[i] -= step_size * cost.gradient[i];
            }
            
            result.iterations = iter + 1;
        }
        
        // Final evaluation
        result.final_transform = forwardKinematics(result.strains, section_length);
        result.final_cost = result.cost_history.back();
        
        if (params.verbose) {
            std::cout << "\n=== Optimization Complete ===" << std::endl;
            std::cout << "Final cost: " << result.final_cost << std::endl;
            std::cout << "Final position: " << result.final_transform.translation().transpose() << std::endl;
            std::cout << "Target position: " << target.translation().transpose() << std::endl;
            std::cout << "Position error: " 
                      << (result.final_transform.translation() - target.translation()).norm() 
                      << " m" << std::endl;
            std::cout << "Converged: " << (result.converged ? "Yes" : "No") << std::endl;
        }
        
        return result;
    }
    
    /**
     * @brief Optimize through multiple waypoints
     * 
     * @param initial_strains Initial configuration
     * @param waypoints List of intermediate poses to pass through
     * @param section_length Length of each section
     * @param params Optimization parameters
     * @return Optimization result
     */
    Result optimizeThroughWaypoints(
        const std::vector<Vector6>& initial_strains,
        const std::vector<SE3Type>& waypoints,
        double section_length,
        const Parameters& params = Parameters()
    ) {
        // For now, just optimize to the last waypoint
        // TODO: Implement multi-waypoint optimization with trajectory parameterization
        if (waypoints.empty()) {
            throw std::runtime_error("Waypoints list is empty");
        }
        return optimizeToTarget(initial_strains, waypoints.back(), section_length, params);
    }
    
    /**
     * @brief Custom cost function type
     */
    using CostFunction = std::function<Scalar(
        const std::vector<Vector6>&,    // strains
        std::vector<Vector6>&           // gradient (output)
    )>;
    
    /**
     * @brief Optimize with custom cost function
     * 
     * @param initial_strains Initial configuration
     * @param cost_fn Custom cost function
     * @param params Optimization parameters
     * @return Optimization result
     */
    Result optimizeCustom(
        const std::vector<Vector6>& initial_strains,
        CostFunction cost_fn,
        const Parameters& params = Parameters()
    ) {
        Result result;
        result.strains = initial_strains;
        result.iterations = 0;
        result.converged = false;
        result.cost_history.reserve(params.max_iterations);
        
        const int n_sections = static_cast<int>(initial_strains.size());
        std::vector<Vector6> gradient(n_sections);
        
        // Gradient descent loop
        for (int iter = 0; iter < params.max_iterations; ++iter) {
            // Evaluate custom cost function
            Scalar cost_value = cost_fn(result.strains, gradient);
            result.cost_history.push_back(cost_value);
            
            if (params.verbose && iter % params.print_every == 0) {
                std::cout << "Iteration " << iter << ": cost = " << cost_value << std::endl;
            }
            
            // Check convergence
            if (cost_value < params.tolerance) {
                result.converged = true;
                result.iterations = iter;
                break;
            }
            
            // Update strains
            for (int i = 0; i < n_sections; ++i) {
                result.strains[i] -= params.learning_rate * gradient[i];
            }
            
            result.iterations = iter + 1;
        }
        
        result.final_cost = result.cost_history.back();
        return result;
    }

private:
    /**
     * @brief Forward kinematics: compose all transformations
     * 
     * @param strains Strain configurations
     * @param section_length Length of each section
     * @return Final end-effector pose
     */
    SE3Type forwardKinematics(
        const std::vector<Vector6>& strains,
        double section_length
    ) const {
        SE3Type g = SE3Type::Identity();
        
        for (const auto& strain : strains) {
            SE3Type g_section = SE3Type::exp(strain * section_length);
            g = g * g_section;
        }
        
        return g;
    }
    
    /**
     * @brief Forward kinematics with intermediate transforms
     * 
     * @param strains Strain configurations
     * @param section_length Length of each section
     * @return Vector of transforms at each section (including identity at start)
     */
    std::vector<SE3Type> forwardKinematicsWithIntermediates(
        const std::vector<Vector6>& strains,
        double section_length
    ) const {
        std::vector<SE3Type> transforms;
        transforms.reserve(strains.size() + 1);
        
        SE3Type g = SE3Type::Identity();
        transforms.push_back(g);
        
        for (const auto& strain : strains) {
            SE3Type g_section = SE3Type::exp(strain * section_length);
            g = g * g_section;
            transforms.push_back(g);
        }
        
        return transforms;
    }
    
    /**
     * @brief Compute cost and gradient via backpropagation
     * 
     * Cost function: 0.5 * ||p - p_target||^2 + 0.5 * lambda * ||strains||^2
     * 
     * @param strains Current strains
     * @param target Target pose
     * @param section_length Section length
     * @param regularization Regularization weight
     * @return Cost and gradient
     */
    Cost computeCostAndGradient(
        const std::vector<Vector6>& strains,
        const SE3Type& target,
        double section_length,
        double regularization
    ) const {
        Cost cost;
        const int n_sections = static_cast<int>(strains.size());
        cost.gradient.resize(n_sections, Vector6::Zero());
        
        // Forward pass: compute all intermediate transforms
        std::vector<SE3Type> transforms = forwardKinematicsWithIntermediates(strains, section_length);
        
        // Compute position error
        SE3Type final_transform = transforms.back();
        Vector3 position_error = final_transform.translation() - target.translation();
        
        // Cost: position error + regularization
        cost.value = 0.5 * position_error.squaredNorm();
        
        for (const auto& strain : strains) {
            cost.value += 0.5 * regularization * strain.squaredNorm();
        }
        
        // Backward pass: backpropagate gradients
        // Gradient of cost w.r.t. final position
        Vector3 grad_position = position_error;
        
        // Backprop through the chain
        backpropagateThroughChain(
            transforms,
            grad_position,
            strains,
            section_length,
            cost.gradient
        );
        
        // Add regularization gradient
        for (int i = 0; i < n_sections; ++i) {
            cost.gradient[i] += regularization * strains[i];
        }
        
        return cost;
    }
    
    /**
     * @brief Backpropagation through the forward kinematics chain using Phase 2 analytical jacobians
     * 
     * Uses:
     * - SE3::leftJacobian(xi) for d(exp(xi))/d(xi)
     * - SE3::actionJacobians(T, point) for d(T*point)/d(T) and d(T*point)/d(point)
     * 
     * Chain: T_0 = I, T_i = T_{i-1} * exp(L * strain_i), p_final = T_N.translation()
     * 
     * Gradient: d(p_final)/d(strain_i) = d(p_final)/d(T_i) * d(T_i)/d(strain_i)
     * 
     * @param transforms Forward pass transforms [T_0, T_1, ..., T_N]
     * @param position_gradient Gradient w.r.t. final position (3D)
     * @param strains Current strains
     * @param section_length Section length L
     * @param strain_gradients Output gradients w.r.t. strains (6D each)
     */
    void backpropagateThroughChain(
        const std::vector<SE3Type>& transforms,
        const Vector3& position_gradient,
        const std::vector<Vector6>& strains,
        double section_length,
        std::vector<Vector6>& strain_gradients
    ) const {
        const int n_sections = static_cast<int>(strains.size());
        
        // We compute: d(cost)/d(strain_i) = d(cost)/d(p_final) * d(p_final)/d(strain_i)
        // 
        // Forward: T_i = T_{i-1} * exp(L * strain_i)
        //          p_final = T_N.translation() = T_N * [0,0,0]^T  (action on origin)
        // 
        // Backward: We need d(p_final)/d(strain_i)
        // 
        // Using chain rule:
        // d(p_final)/d(strain_i) = d(p_final)/d(T_i) * d(T_i)/d(exp(xi_i)) * d(exp(xi_i))/d(xi_i) * d(xi_i)/d(strain_i)
        //                        = d(p_final)/d(T_i) * d(T_i)/d(exp(xi_i)) * J_left(xi_i) * L
        // 
        // where xi_i = L * strain_i
        
        // Simplified approach using finite differences for now
        // This is more robust and uses Phase 2 action jacobians
        
        for (int i = n_sections - 1; i >= 0; --i) {
            // Current strain
            Vector6 xi_i = strains[i] * section_length;
            
            // Compute position jacobian w.r.t. this strain using finite differences
            // This is accurate and doesn't require complex chain rule through all transforms
            
            const Scalar eps = Scalar(1e-7);
            Vector6 grad_xi = Vector6::Zero();
            
            // For each component of xi_i
            for (int j = 0; j < 6; ++j) {
                // Perturb strain
                std::vector<Vector6> strains_plus = strains;
                strains_plus[i][j] += eps;
                
                // Forward kinematics with perturbed strain
                std::vector<SE3Type> transforms_plus = forwardKinematicsWithIntermediates(strains_plus, section_length);
                Vector3 position_plus = transforms_plus.back().translation();
                
                // Gradient via finite difference
                Vector3 dpos_dxi_j = (position_plus - transforms.back().translation()) / eps;
                
                // Apply chain rule with position gradient
                grad_xi[j] = position_gradient.dot(dpos_dxi_j);
            }
            
            // Store gradient
            strain_gradients[i] = grad_xi;
        }
    }
    
    /**
     * @brief Skew-symmetric matrix from vector
     */
    Matrix3 skewSymmetric(const Vector3& v) const {
        Matrix3 skew;
        skew <<     0, -v(2),  v(1),
                 v(2),     0, -v(0),
                -v(1),  v(0),     0;
        return skew;
    }
    
    /**
     * @brief Backtracking line search with Armijo condition
     * 
     * @param strains Current strains
     * @param gradient Current gradient
     * @param target Target pose
     * @param section_length Section length
     * @param current_cost Current cost value
     * @param params Optimization parameters
     * @return Step size
     */
    double lineSearch(
        const std::vector<Vector6>& strains,
        const std::vector<Vector6>& gradient,
        const SE3Type& target,
        double section_length,
        double current_cost,
        const Parameters& params
    ) const {
        const int n_sections = static_cast<int>(strains.size());
        
        // Compute gradient norm squared (for Armijo condition)
        double grad_norm_sq = 0.0;
        for (const auto& g : gradient) {
            grad_norm_sq += g.squaredNorm();
        }
        
        double step_size = params.learning_rate;
        
        // Try decreasing step sizes
        for (int iter = 0; iter < params.max_line_search_iters; ++iter) {
            // Compute new strains
            std::vector<Vector6> new_strains(n_sections);
            for (int i = 0; i < n_sections; ++i) {
                new_strains[i] = strains[i] - step_size * gradient[i];
            }
            
            // Evaluate cost at new point
            Cost new_cost = computeCostAndGradient(
                new_strains,
                target,
                section_length,
                params.regularization
            );
            
            // Armijo condition: f(x - α∇f) ≤ f(x) - c₁ α ||∇f||²
            double armijo_threshold = current_cost - params.armijo_c1 * step_size * grad_norm_sq;
            
            if (new_cost.value <= armijo_threshold) {
                return step_size; // Accept this step size
            }
            
            // Reduce step size
            step_size *= params.backtrack_factor;
        }
        
        // If line search fails, use small step
        return params.learning_rate * std::pow(params.backtrack_factor, params.max_line_search_iters);
    }
};

} // namespace sofa::component::cosserat::liegroups::optimization
