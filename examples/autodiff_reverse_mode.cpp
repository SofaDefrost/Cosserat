/******************************************************************************
 *                 SOFA, Simulation Open-Framework Architecture * (c) 2006
 *INRIA, USTL, UJF, CNRS, MGH                     *
 *                                                                             *
 * This program is free software; you can redistribute it and/or modify it     *
 * under the terms of the GNU Lesser General Public License as published by    *
 * the Free Software Foundation; either version 2.1 of the License, or (at     *
 * your option) any later version.                                             *
 *                                                                             *
 * This program is distributed in the hope that it will be useful, but WITHOUT *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
 * for more details.                                                           *
 *                                                                             *
 * You should have received a copy of the GNU Lesser General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/\>. *
 ******************************************************************************/

/**
 * @file autodiff_reverse_mode.cpp
 * @brief Example demonstrating reverse mode automatic differentiation for Cosserat rods
 * 
 * This example shows how to:
 * - Use autodiff::var for reverse mode differentiation
 * - Compute gradients through multi-segment Cosserat forward kinematics
 * - Optimize rod strains to reach target positions
 * 
 * Reverse mode is efficient when you have:
 * - Many inputs (e.g., 60 strain parameters)
 * - One scalar output (e.g., cost function)
 * 
 * For N segments (6N parameters), reverse mode is ~N times faster than forward mode!
 * 
 * Compile with:
 *   cmake -DCOSSERAT_BUILD_EXAMPLES=ON -DCOSSERAT_WITH_AUTODIFF=ON ..
 *   make
 */

#ifdef COSSERAT_WITH_AUTODIFF

#include <iostream>
#include <iomanip>
#include <vector>
#include <autodiff/reverse/var.hpp>
#include "../src/liegroups/SE3.h"
#include "../src/liegroups/AutodiffSupport.h"

using namespace autodiff;
using namespace sofa::component::cosserat::liegroups;

// ============================================================================
// Example 1: Single Segment - Distance to Target
// ============================================================================

void example1_single_segment() {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "Example 1: Single Segment - Gradient w.r.t. Strain\n";
    std::cout << std::string(70, '=') << "\n\n";

    using SE3var = SE3<var>;

    // Strain parameters: [φx, φy, φz, ρx, ρy, ρz]
    Eigen::Matrix<double, 6, 1> strain_val;
    strain_val << 0.0, 0.2, 0.0, 0.0, 0.0, 0.0;  // Bending in Y

    // Convert to var
    Eigen::Matrix<var, 6, 1> strain;
    for (int i = 0; i < 6; i++) {
        strain[i] = strain_val[i];
    }

    // Segment length
    double L = 1.0;

    // Forward kinematics: T = exp(L * strain)
    Eigen::Matrix<var, 6, 1> xi = L * strain;
    SE3var T = SE3var::exp(xi);
    auto pos = T.translation();

    // Target: straight rod at (1, 0, 0)
    Eigen::Vector3d target(1.0, 0.0, 0.0);

    // Cost: squared distance to target
    var dx = pos[0] - target[0];
    var dy = pos[1] - target[1];
    var dz = pos[2] - target[2];
    var cost = dx*dx + dy*dy + dz*dz;

    // Compute ALL gradients in ONE reverse pass
    Derivatives dcost = derivatives(cost);

    Eigen::Matrix<double, 6, 1> gradient;
    for (int i = 0; i < 6; i++) {
        gradient[i] = dcost(strain[i]);
    }

    // Display results
    std::cout << "Segment length:    " << L << " m\n";
    std::cout << "Initial strain:    [" << strain_val.transpose() << "]\n";
    std::cout << "                   (bending in Y direction)\n\n";
    
    std::cout << "Tip position:      [" << toDouble(pos).transpose() << "]\n";
    std::cout << "Target position:   [" << target.transpose() << "]\n";
    std::cout << "Distance cost:     " << val(cost) << "\n\n";

    std::cout << "Gradient ∂cost/∂strain:\n";
    std::cout << "  φx (torsion):    " << std::setw(12) << gradient[0] << "\n";
    std::cout << "  φy (bending Y):  " << std::setw(12) << gradient[1] << " ← most sensitive\n";
    std::cout << "  φz (bending Z):  " << std::setw(12) << gradient[2] << "\n";
    std::cout << "  ρx (elongation): " << std::setw(12) << gradient[3] << "\n";
    std::cout << "  ρy (shear Y):    " << std::setw(12) << gradient[4] << "\n";
    std::cout << "  ρz (shear Z):    " << std::setw(12) << gradient[5] << "\n";

    std::cout << "\n✓ Reverse mode computed all 6 gradients in ONE pass!\n";
}

// ============================================================================
// Example 2: Multi-Segment Cosserat Rod
// ============================================================================

void example2_multi_segment() {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "Example 2: Multi-Segment Cosserat Rod (10 segments)\n";
    std::cout << std::string(70, '=') << "\n\n";

    using SE3var = SE3<var>;

    const int N = 10;  // 10 segments = 60 parameters
    const double L = 0.15;  // 15 cm per segment

    // Strain for each segment
    std::vector<Eigen::Matrix<var, 6, 1>> strains(N);
    std::vector<Eigen::Matrix<double, 6, 1>> strain_vals(N);

    // Initialize with gradually increasing bending
    for (int i = 0; i < N; i++) {
        strain_vals[i] << 0.0, 0.05 * (i + 1), 0.0, 0.0, 0.0, 0.0;
        for (int j = 0; j < 6; j++) {
            strains[i][j] = strain_vals[i][j];
        }
    }

    // Forward kinematics: T = T1 * T2 * ... * T10
    SE3var T = SE3var::identity();
    for (int i = 0; i < N; i++) {
        Eigen::Matrix<var, 6, 1> xi = L * strains[i];
        SE3var Ti = SE3var::exp(xi);
        T = T * Ti;
    }

    // Target: rod should reach (1.5, 0, 0)
    Eigen::Vector3d target(1.5, 0.0, 0.0);
    auto pos = T.translation();

    // Cost: distance to target
    var dx = pos[0] - target[0];
    var dy = pos[1] - target[1];
    var dz = pos[2] - target[2];
    var cost = dx*dx + dy*dy + dz*dz;

    // ONE reverse pass → gradients for ALL 60 parameters!
    Derivatives dcost = derivatives(cost);

    std::vector<Eigen::Matrix<double, 6, 1>> gradients(N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < 6; j++) {
            gradients[i][j] = dcost(strains[i][j]);
        }
    }

    // Display results
    std::cout << "Configuration:\n";
    std::cout << "  Number of segments:  " << N << "\n";
    std::cout << "  Segment length:      " << L << " m\n";
    std::cout << "  Total length:        " << N * L << " m\n";
    std::cout << "  Total parameters:    " << N * 6 << " (6 per segment)\n\n";

    std::cout << "Results:\n";
    std::cout << "  Tip position:        [" << toDouble(pos).transpose() << "]\n";
    std::cout << "  Target position:     [" << target.transpose() << "]\n";
    std::cout << "  Distance cost:       " << val(cost) << "\n\n";

    std::cout << "Gradient magnitudes per segment:\n";
    for (int i = 0; i < N; i++) {
        std::cout << "  Segment " << std::setw(2) << i + 1 << ": "
                  << "||∇||₂ = " << std::setw(10) << std::setprecision(6) 
                  << gradients[i].norm() << "\n";
    }

    // Show gradient for first and last segment
    std::cout << "\nGradient for Segment 1 (base):\n";
    std::cout << "  [" << gradients[0].transpose() << "]\n";
    std::cout << "\nGradient for Segment " << N << " (tip):\n";
    std::cout << "  [" << gradients[N-1].transpose() << "]\n";

    std::cout << "\n✓ Reverse mode computed ALL 60 gradients in ONE pass!\n";
    std::cout << "  (Forward mode would need 60 passes)\n";
    std::cout << "  Speedup: ~60× faster!\n";
}

// ============================================================================
// Example 3: Gradient Descent Optimization
// ============================================================================

void example3_optimization() {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "Example 3: Gradient Descent Optimization\n";
    std::cout << std::string(70, '=') << "\n\n";

    using SE3var = SE3<var>;

    const int N = 3;
    const double L = 0.5;
    const double learning_rate = 0.05;
    const int max_iterations = 50;

    // Initial guess: zero strains
    std::vector<Eigen::Matrix<double, 6, 1>> strains(N);
    for (int i = 0; i < N; i++) {
        strains[i].setZero();
    }

    // Target
    Eigen::Vector3d target(1.2, 0.3, 0.0);

    std::cout << "Configuration:\n";
    std::cout << "  Segments:        " << N << "\n";
    std::cout << "  Length per seg:  " << L << " m\n";
    std::cout << "  Learning rate:   " << learning_rate << "\n";
    std::cout << "  Target:          [" << target.transpose() << "]\n\n";

    std::cout << "Optimization progress:\n";
    std::cout << "Iter    Cost         Position\n";
    std::cout << "----  ----------  ------------------------\n";

    for (int iter = 0; iter < max_iterations; iter++) {
        // Convert to var
        std::vector<Eigen::Matrix<var, 6, 1>> strains_var(N);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < 6; j++) {
                strains_var[i][j] = strains[i][j];
            }
        }

        // Forward kinematics
        SE3var T = SE3var::identity();
        for (int i = 0; i < N; i++) {
            SE3var Ti = SE3var::exp(L * strains_var[i]);
            T = T * Ti;
        }

        auto pos = T.translation();

        // Cost
        var dx = pos[0] - target[0];
        var dy = pos[1] - target[1];
        var dz = pos[2] - target[2];
        var cost = dx*dx + dy*dy + dz*dz;

        // Compute gradients
        Derivatives dcost = derivatives(cost);

        std::vector<Eigen::Matrix<double, 6, 1>> gradients(N);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < 6; j++) {
                gradients[i][j] = dcost(strains_var[i][j]);
            }
        }

        // Print every 5 iterations
        if (iter % 5 == 0) {
            auto pos_val = toDouble(pos);
            std::cout << std::setw(4) << iter << "  "
                      << std::setw(10) << std::setprecision(6) << val(cost) << "  "
                      << "[" << std::setw(6) << std::setprecision(3) << pos_val.transpose() << "]\n";
        }

        // Gradient descent update
        for (int i = 0; i < N; i++) {
            strains[i] -= learning_rate * gradients[i];
        }

        // Check convergence
        if (val(cost) < 1e-6) {
            std::cout << "\n✓ Converged at iteration " << iter << "!\n";
            break;
        }
    }

    std::cout << "\n✓ Optimization completed using reverse mode autodiff!\n";
}

// ============================================================================
// Example 4: Performance Comparison
// ============================================================================

void example4_performance() {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "Example 4: Performance Analysis\n";
    std::cout << std::string(70, '=') << "\n\n";

    std::cout << "Reverse mode autodiff is extremely efficient for optimization!\n\n";

    std::cout << "Problem: N segments (6N parameters) → 1 cost value\n\n";

    std::cout << "Method                    Complexity    N=10    N=50    N=100\n";
    std::cout << "------------------------  ----------  ------  ------  -------\n";
    std::cout << "Finite Differences        O(12N)       120     600     1200\n";
    std::cout << "Forward Mode (naive)      O(6N)         60     300      600\n";
    std::cout << "Reverse Mode              O(1)           1       1        1\n";
    std::cout << "------------------------  ----------  ------  ------  -------\n";
    std::cout << "Speedup (vs finite diff)               120×    600×    1200×\n";
    std::cout << "Speedup (vs forward mode)               60×    300×     600×\n\n";

    std::cout << "Key insights:\n";
    std::cout << "  • Reverse mode cost is CONSTANT regardless of parameter count\n";
    std::cout << "  • For optimization problems (many params → 1 cost), always use reverse mode\n";
    std::cout << "  • Forward mode is better for few params → many outputs\n";

    std::cout << "\n✓ Reverse mode is the optimal choice for Cosserat optimization!\n";
}

// ============================================================================
// Main Program
// ============================================================================

int main() {
    std::cout << "\n";
    std::cout << "╔══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║                                                                  ║\n";
    std::cout << "║    Reverse Mode Automatic Differentiation for Cosserat Rods     ║\n";
    std::cout << "║                                                                  ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════════╝\n";

    std::cout << "\nThis example demonstrates reverse mode autodiff using autodiff::var.\n";
    std::cout << "Reverse mode is efficient for: many inputs → one output (optimization!).\n";

    try {
        example1_single_segment();
        example2_multi_segment();
        example3_optimization();
        example4_performance();

        std::cout << "\n" << std::string(70, '=') << "\n";
        std::cout << "All examples completed successfully!\n";
        std::cout << std::string(70, '=') << "\n\n";

        std::cout << "Next steps:\n";
        std::cout << "  • See CosseratTrajectoryOptimizer.h for full implementation\n";
        std::cout << "  • Run simple_trajectory_optimization example\n";
        std::cout << "  • Read DIFFERENTIATION.md for complete guide\n\n";

    } catch (const std::exception& e) {
        std::cerr << "\n❌ Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}

#else

// Fallback when autodiff is not enabled
int main() {
    std::cerr << "❌ This example requires autodiff support.\n";
    std::cerr << "   Recompile with: cmake -DCOSSERAT_WITH_AUTODIFF=ON ..\n";
    return 1;
}

#endif // COSSERAT_WITH_AUTODIFF
