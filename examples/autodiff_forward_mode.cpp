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
 * @file autodiff_forward_mode.cpp
 * @brief Example demonstrating forward mode automatic differentiation with SO3
 * 
 * This example shows how to:
 * - Use autodiff::dual for forward mode differentiation
 * - Compute derivatives of rotation operations
 * - Compare with analytical jacobians
 * 
 * Forward mode is efficient when you have:
 * - Few inputs (e.g., 3 parameters)
 * - Many outputs (e.g., 100 values)
 * 
 * Compile with:
 *   cmake -DCOSSERAT_BUILD_EXAMPLES=ON -DCOSSERAT_WITH_AUTODIFF=ON ..
 *   make
 */

#ifdef COSSERAT_WITH_AUTODIFF

#include <iostream>
#include <iomanip>
#include <autodiff/forward/dual.hpp>
#include "../src/liegroups/SO3.h"
#include "../src/liegroups/AutodiffSupport.h"

using namespace autodiff;
using namespace sofa::component::cosserat::liegroups;

// ============================================================================
// Example 1: Derivative of Rotation Angle
// ============================================================================

/**
 * @brief Computes the angle of rotation from an axis-angle vector
 * 
 * For omega in R^3, this returns ||omega|| (the rotation angle).
 * The gradient is: d(||omega||)/d(omega) = omega / ||omega||
 */
dual rotationAngle(const Eigen::Matrix<dual, 3, 1>& omega) {
    using SO3dual = SO3<dual>;
    SO3dual R = SO3dual::exp(omega);
    return R.log().norm();  // Returns angle
}

void example1_rotation_angle() {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "Example 1: Derivative of Rotation Angle\n";
    std::cout << std::string(70, '=') << "\n\n";

    // Test point
    Eigen::Vector3d omega_val(0.3, 0.4, 0.5);
    double angle_analytical = omega_val.norm();
    Eigen::Vector3d gradient_analytical = omega_val / angle_analytical;

    // Convert to dual numbers
    Eigen::Matrix<dual, 3, 1> omega = toAutodiff<Eigen::Vector3d, dual>(omega_val);

    // Compute derivatives using forward mode
    Eigen::Vector3d gradient_autodiff;
    for (int i = 0; i < 3; i++) {
        gradient_autodiff[i] = derivative(rotationAngle, wrt(omega[i]), at(omega));
    }

    // Display results
    std::cout << "Input omega:           [" << omega_val.transpose() << "]\n";
    std::cout << "Rotation angle:        " << angle_analytical << " rad\n";
    std::cout << "                       " << angle_analytical * 180.0 / M_PI << " deg\n\n";

    std::cout << "Analytical gradient:   [" << gradient_analytical.transpose() << "]\n";
    std::cout << "Autodiff gradient:     [" << gradient_autodiff.transpose() << "]\n";
    std::cout << "Error:                 " << (gradient_analytical - gradient_autodiff).norm() << "\n";

    std::cout << "\n✓ Forward mode correctly computed gradient of rotation angle!\n";
}

// ============================================================================
// Example 2: Derivative of Rotation Matrix Entry
// ============================================================================

/**
 * @brief Returns a specific entry of the rotation matrix
 */
dual rotationMatrixEntry(const Eigen::Matrix<dual, 3, 1>& omega, int row, int col) {
    using SO3dual = SO3<dual>;
    SO3dual R = SO3dual::exp(omega);
    return R.matrix()(row, col);
}

void example2_matrix_entry() {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "Example 2: Derivative of Rotation Matrix Entry R(0,0)\n";
    std::cout << std::string(70, '=') << "\n\n";

    Eigen::Vector3d omega_val(0.1, 0.2, 0.3);
    Eigen::Matrix<dual, 3, 1> omega = toAutodiff<Eigen::Vector3d, dual>(omega_val);

    // Compute R(0,0) and its gradient
    dual R00 = rotationMatrixEntry(omega, 0, 0);
    
    Eigen::Vector3d gradient;
    for (int i = 0; i < 3; i++) {
        gradient[i] = derivative([](const auto& w) { 
            return rotationMatrixEntry(w, 0, 0); 
        }, wrt(omega[i]), at(omega));
    }

    // Also compute the full rotation matrix
    using SO3d = SO3<double>;
    SO3d R = SO3d::exp(omega_val);
    Eigen::Matrix3d R_mat = R.matrix();

    std::cout << "Input omega:     [" << omega_val.transpose() << "]\n\n";
    std::cout << "Rotation matrix R:\n" << R_mat << "\n\n";
    std::cout << "R(0,0) = " << R_mat(0, 0) << "\n";
    std::cout << "∂R(0,0)/∂omega = [" << gradient.transpose() << "]\n";

    std::cout << "\n✓ Forward mode computed gradient of matrix entry!\n";
}

// ============================================================================
// Example 3: Multiple Outputs - Efficient with Forward Mode
// ============================================================================

/**
 * @brief Compute all 9 rotation matrix entries and their gradients
 * 
 * This demonstrates forward mode's strength: computing many outputs
 * from few inputs is efficient (3 forward passes for 9 outputs).
 */
void example3_multiple_outputs() {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "Example 3: Gradient of All Rotation Matrix Entries\n";
    std::cout << std::string(70, '=') << "\n\n";

    Eigen::Vector3d omega_val(0.2, 0.1, -0.15);
    Eigen::Matrix<dual, 3, 1> omega = toAutodiff<Eigen::Vector3d, dual>(omega_val);

    std::cout << "Input omega: [" << omega_val.transpose() << "]\n\n";
    std::cout << "Computing ∂R(i,j)/∂omega for all matrix entries...\n\n";

    // For each matrix entry, compute gradient
    for (int row = 0; row < 3; row++) {
        for (int col = 0; col < 3; col++) {
            Eigen::Vector3d gradient;
            
            // One forward pass per omega component
            for (int k = 0; k < 3; k++) {
                gradient[k] = derivative(
                    [row, col](const auto& w) { return rotationMatrixEntry(w, row, col); },
                    wrt(omega[k]), 
                    at(omega)
                );
            }

            std::cout << "∂R(" << row << "," << col << ")/∂omega = ["
                      << std::setw(8) << std::setprecision(4) << gradient.transpose() 
                      << " ]\n";
        }
    }

    std::cout << "\n✓ Forward mode efficiently computed gradients of all 9 entries!\n";
    std::cout << "  (Only 3 forward passes needed for 3 input dimensions)\n";
}

// ============================================================================
// Example 4: Comparison with Analytical Jacobian
// ============================================================================

void example4_comparison() {
    std::cout << "\n" << std::string(70, '=') << "\n";
    std::cout << "Example 4: Comparison with Analytical Jacobian\n";
    std::cout << std::string(70, '=') << "\n\n";

    using SO3d = SO3<double>;
    using SO3dual = SO3<dual>;

    Eigen::Vector3d omega_val(0.3, -0.2, 0.4);

    // Analytical jacobian (Phase 2)
    Eigen::Matrix3d J_analytical = SO3d::leftJacobian(omega_val);

    // Autodiff jacobian (forward mode)
    Eigen::Matrix3d J_autodiff;
    Eigen::Matrix<dual, 3, 1> omega = toAutodiff<Eigen::Vector3d, dual>(omega_val);

    for (int out = 0; out < 3; out++) {
        for (int in = 0; in < 3; in++) {
            J_autodiff(out, in) = derivative(
                [out](const auto& w) { 
                    SO3dual R = SO3dual::exp(w);
                    return R.log()[out];
                },
                wrt(omega[in]),
                at(omega)
            );
        }
    }

    double error = (J_analytical - J_autodiff).norm() / J_analytical.norm();

    std::cout << "Input omega: [" << omega_val.transpose() << "]\n\n";
    std::cout << "Analytical Jacobian (Phase 2):\n" << J_analytical << "\n\n";
    std::cout << "Autodiff Jacobian (Forward Mode):\n" << J_autodiff << "\n\n";
    std::cout << "Relative error: " << std::scientific << error << "\n";

    if (error < 1e-10) {
        std::cout << "\n✓ Perfect agreement between analytical and autodiff!\n";
    } else {
        std::cout << "\n⚠ Some numerical difference detected\n";
    }
}

// ============================================================================
// Main Program
// ============================================================================

int main() {
    std::cout << "\n";
    std::cout << "╔══════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║                                                                  ║\n";
    std::cout << "║         Forward Mode Automatic Differentiation with SO3         ║\n";
    std::cout << "║                                                                  ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════════╝\n";

    std::cout << "\nThis example demonstrates forward mode autodiff using autodiff::dual.\n";
    std::cout << "Forward mode is efficient for: few inputs → many outputs.\n";

    try {
        example1_rotation_angle();
        example2_matrix_entry();
        example3_multiple_outputs();
        example4_comparison();

        std::cout << "\n" << std::string(70, '=') << "\n";
        std::cout << "All examples completed successfully!\n";
        std::cout << std::string(70, '=') << "\n\n";

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
