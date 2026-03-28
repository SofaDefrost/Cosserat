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

/**
 * @file simple_trajectory_optimization.cpp
 * @brief Simple example of trajectory optimization for Cosserat rods
 * 
 * This example demonstrates how to use the CosseratTrajectoryOptimizer
 * to find strain configurations that move the end-effector to a target position.
 * 
 * The optimizer uses gradient descent with backpropagation through the
 * forward kinematics chain.
 */

#include <iostream>
#include <vector>
#include <Eigen/Dense>

#include "../src/liegroups/SE3.h"
#include "../src/liegroups/optimization/CosseratTrajectoryOptimizer.h"

using namespace sofa::component::cosserat::liegroups;
using namespace sofa::component::cosserat::liegroups::optimization;

int main() {
    std::cout << "\n";
    std::cout << "╔══════════════════════════════════════════════════════════╗\n";
    std::cout << "║   Cosserat Trajectory Optimization - Simple Example     ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════╝\n";
    std::cout << "\n";
    
    // ======== CONFIGURATION ========
    const int n_sections = 10;                  // Number of sections
    const double section_length = 0.1;          // 10cm per section (total 1m)
    
    // Target: 80cm in X direction, 20cm up in Z
    Eigen::Vector3d target_position(0.8, 0.0, 0.2);
    
    std::cout << "Configuration:" << std::endl;
    std::cout << "  - Number of sections: " << n_sections << std::endl;
    std::cout << "  - Section length: " << section_length << " m" << std::endl;
    std::cout << "  - Total length: " << n_sections * section_length << " m" << std::endl;
    std::cout << "  - Target position: [" 
              << target_position.x() << ", "
              << target_position.y() << ", "
              << target_position.z() << "] m" << std::endl;
    std::cout << "\n";
    
    // ======== INITIALIZATION ========
    // Start with straight beam (zero strains)
    std::vector<Eigen::Matrix<double, 6, 1>> initial_strains(n_sections);
    for (auto& strain : initial_strains) {
        strain.setZero();
    }
    
    // Compute initial position
    SE3d initial_transform = SE3d::Identity();
    for (const auto& strain : initial_strains) {
        initial_transform = initial_transform * SE3d::exp(strain * section_length);
    }
    
    std::cout << "Initial configuration:" << std::endl;
    std::cout << "  - Position: " << initial_transform.translation().transpose() << " m" << std::endl;
    std::cout << "  - Distance to target: " 
              << (initial_transform.translation() - target_position).norm() 
              << " m" << std::endl;
    std::cout << "\n";
    
    // ======== OPTIMIZATION ========
    // Create target SE(3) transform (position only, identity rotation)
    SE3d target = SE3d::Identity();
    target.translation() = target_position;
    
    // Setup optimizer parameters
    CosseratTrajectoryOptimizer<double>::Parameters params;
    params.learning_rate = 0.05;               // Learning rate
    params.max_iterations = 500;               // Max iterations
    params.tolerance = 1e-6;                   // Convergence tolerance
    params.regularization = 0.001;             // Small regularization
    params.use_line_search = true;             // Use line search
    params.verbose = true;                     // Print progress
    params.print_every = 50;                   // Print every 50 iterations
    
    // Create optimizer
    CosseratTrajectoryOptimizer<double> optimizer;
    
    // Run optimization
    std::cout << "Starting optimization...\n" << std::endl;
    auto result = optimizer.optimizeToTarget(
        initial_strains,
        target,
        section_length,
        params
    );
    
    // ======== RESULTS ========
    std::cout << "\n";
    std::cout << "╔══════════════════════════════════════════════════════════╗\n";
    std::cout << "║                     Results Summary                      ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════╝\n";
    std::cout << "\n";
    
    std::cout << "Optimization completed!" << std::endl;
    std::cout << "  - Converged: " << (result.converged ? "✓ Yes" : "✗ No") << std::endl;
    std::cout << "  - Iterations: " << result.iterations << std::endl;
    std::cout << "  - Final cost: " << result.final_cost << std::endl;
    std::cout << "\n";
    
    std::cout << "Final configuration:" << std::endl;
    std::cout << "  - Position: " << result.final_transform.translation().transpose() << " m" << std::endl;
    std::cout << "  - Target:   " << target_position.transpose() << " m" << std::endl;
    
    Eigen::Vector3d final_error = result.final_transform.translation() - target_position;
    std::cout << "  - Error:    " << final_error.norm() << " m" << std::endl;
    std::cout << "\n";
    
    // Print optimized strains for first few sections
    std::cout << "Optimized strains (first 3 sections):" << std::endl;
    for (int i = 0; i < std::min(3, n_sections); ++i) {
        std::cout << "  Section " << i << ": [";
        for (int j = 0; j < 6; ++j) {
            std::cout << result.strains[i](j);
            if (j < 5) std::cout << ", ";
        }
        std::cout << "]" << std::endl;
    }
    std::cout << "  ..." << std::endl;
    std::cout << "\n";
    
    // ======== CONVERGENCE ANALYSIS ========
    if (result.cost_history.size() > 1) {
        std::cout << "Cost history (every 50 iterations):" << std::endl;
        for (size_t i = 0; i < result.cost_history.size(); i += 50) {
            std::cout << "  Iter " << i << ": " << result.cost_history[i] << std::endl;
        }
        if ((result.cost_history.size() - 1) % 50 != 0) {
            std::cout << "  Iter " << (result.cost_history.size() - 1) 
                      << ": " << result.cost_history.back() << std::endl;
        }
    }
    std::cout << "\n";
    
    // ======== VALIDATION ========
    // Verify forward kinematics
    SE3d verification = SE3d::Identity();
    for (const auto& strain : result.strains) {
        verification = verification * SE3d::exp(strain * section_length);
    }
    
    double verification_error = (verification.translation() - result.final_transform.translation()).norm();
    std::cout << "Verification:" << std::endl;
    std::cout << "  - FK recomputation error: " << verification_error;
    if (verification_error < 1e-10) {
        std::cout << " ✓ (numerical precision)" << std::endl;
    } else {
        std::cout << " ✗ (something is wrong!)" << std::endl;
    }
    std::cout << "\n";
    
    // ======== SUCCESS CRITERIA ========
    bool success = result.converged && (final_error.norm() < 0.01); // Less than 1cm error
    
    if (success) {
        std::cout << "╔══════════════════════════════════════════════════════════╗\n";
        std::cout << "║                  ✓ SUCCESS                              ║\n";
        std::cout << "║  The optimizer successfully reached the target!          ║\n";
        std::cout << "╚══════════════════════════════════════════════════════════╝\n";
        return 0;
    } else {
        std::cout << "╔══════════════════════════════════════════════════════════╗\n";
        std::cout << "║                  ⚠ PARTIAL SUCCESS                       ║\n";
        std::cout << "║  The optimizer made progress but didn't fully converge.  ║\n";
        std::cout << "║  Try adjusting the parameters or increasing iterations.  ║\n";
        std::cout << "╚══════════════════════════════════════════════════════════╝\n";
        return 1;
    }
}
