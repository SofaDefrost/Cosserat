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
 * @file test_autodiff_integration.cpp
 * @brief Tests for automatic differentiation integration with Lie groups
 * 
 * This file contains tests demonstrating the use of the autodiff library
 * (https://autodiff.github.io/) with SO3 and SE3 Lie groups.
 * 
 * Tests cover:
 * - Forward mode differentiation with autodiff::dual
 * - Reverse mode differentiation with autodiff::var
 * - Gradient computation for optimization problems
 * - Comparison with analytical jacobians
 */

#include <gtest/gtest.h>

#ifdef COSSERAT_WITH_AUTODIFF

#include <autodiff/forward/dual.hpp>
#include <autodiff/reverse/var.hpp>
#include "../../SO3.h"
#include "../../SE3.h"
#include "../../AutodiffSupport.h"
#include "../DifferentiationTestUtils.h"

using namespace sofa::component::cosserat::liegroups;
using namespace autodiff;

// ============================================================================
// Forward Mode Tests (autodiff::dual)
// ============================================================================

/**
 * @brief Test forward mode differentiation of SO3 exponential map
 * 
 * Computes d(exp(omega))/d(omega) using forward mode autodiff.
 */
TEST(AutodiffIntegration, ForwardMode_SO3_Exponential) {
	using SO3d = SO3<double>;
	using SO3dual = SO3<dual>;
	
	// Test point in so(3)
	Eigen::Vector3d omega_val(0.1, 0.2, 0.3);
	
	// Convert to dual numbers
	Eigen::Matrix<dual, 3, 1> omega = toAutodiff<Eigen::Vector3d, dual>(omega_val);
	
	// Function to differentiate: rotation matrix entry [0,0]
	auto func = [](const Eigen::Matrix<dual, 3, 1>& w) -> dual {
		SO3dual R = SO3dual::exp(w);
		return R.matrix()(0, 0);
	};
	
	// Compute derivatives using forward mode
	Eigen::Vector3d gradient;
	for (int i = 0; i < 3; i++) {
		gradient[i] = derivative(func, wrt(omega[i]), at(omega));
	}
	
	// Verify derivatives are non-zero
	EXPECT_GT(gradient.norm(), 1e-10);
	
	std::cout << "Forward mode gradient of R(0,0) w.r.t. omega: " 
	          << gradient.transpose() << std::endl;
}

/**
 * @brief Test forward mode differentiation of SE3 exponential map
 * 
 * Computes d(exp(xi).translation)/d(xi) using forward mode autodiff.
 */
TEST(AutodiffIntegration, ForwardMode_SE3_Translation) {
	using SE3d = SE3<double>;
	using SE3dual = SE3<dual>;
	
	// Test point in se(3): [rotation | translation]
	Eigen::Matrix<double, 6, 1> xi_val;
	xi_val << 0.1, 0.2, 0.3, 0.5, 0.0, 0.0;
	
	// Convert to dual numbers
	Eigen::Matrix<dual, 6, 1> xi = toAutodiff<Eigen::Matrix<double, 6, 1>, dual>(xi_val);
	
	// Function to differentiate: x-component of translation
	auto func = [](const Eigen::Matrix<dual, 6, 1>& x) -> dual {
		SE3dual T = SE3dual::exp(x);
		return T.translation()[0];
	};
	
	// Compute derivatives using forward mode
	Eigen::Matrix<double, 6, 1> gradient;
	for (int i = 0; i < 6; i++) {
		gradient[i] = derivative(func, wrt(xi[i]), at(xi));
	}
	
	// Translation x should be most sensitive to xi[3] (ρx)
	EXPECT_GT(std::abs(gradient[3]), std::abs(gradient[0]));
	
	std::cout << "Forward mode gradient of T.translation.x w.r.t. xi: " 
	          << gradient.transpose() << std::endl;
}

/**
 * @brief Test forward mode differentiation of rotation angle
 * 
 * Computes d(||log(exp(omega))||)/d(omega) = d(||omega||)/d(omega)
 * Analytical result: gradient = omega / ||omega||
 */
TEST(AutodiffIntegration, ForwardMode_RotationAngle) {
	using SO3dual = SO3<dual>;
	
	Eigen::Vector3d omega_val(0.3, 0.4, 0.5);
	double angle_analytical = omega_val.norm();
	Eigen::Vector3d gradient_analytical = omega_val / angle_analytical;
	
	// Convert to dual
	Eigen::Matrix<dual, 3, 1> omega = toAutodiff<Eigen::Vector3d, dual>(omega_val);
	
	// Function: rotation angle
	auto func = [](const Eigen::Matrix<dual, 3, 1>& w) -> dual {
		SO3dual R = SO3dual::exp(w);
		return R.log().norm();
	};
	
	// Compute derivatives
	Eigen::Vector3d gradient_autodiff;
	for (int i = 0; i < 3; i++) {
		gradient_autodiff[i] = derivative(func, wrt(omega[i]), at(omega));
	}
	
	// Compare with analytical
	EXPECT_NEAR((gradient_autodiff - gradient_analytical).norm(), 0.0, 1e-8);
	
	std::cout << "Forward mode gradient (autodiff): " << gradient_autodiff.transpose() << std::endl;
	std::cout << "Analytical gradient:              " << gradient_analytical.transpose() << std::endl;
}

// ============================================================================
// Reverse Mode Tests (autodiff::var)
// ============================================================================

/**
 * @brief Test reverse mode differentiation of SE3 distance cost
 * 
 * Given a target position, compute gradient of distance cost w.r.t. se(3) parameters.
 * This is efficient with reverse mode since we have many inputs (6) and one output.
 */
TEST(AutodiffIntegration, ReverseMode_SE3_DistanceCost) {
	using SE3d = SE3<double>;
	using SE3var = SE3<var>;
	
	// Initial parameters in se(3)
	Eigen::Matrix<double, 6, 1> xi_val;
	xi_val << 0.0, 0.0, 0.0, 0.5, 0.1, 0.2;
	
	// Convert to var
	Eigen::Matrix<var, 6, 1> xi;
	for (int i = 0; i < 6; i++) {
		xi[i] = xi_val[i];
	}
	
	// Target position
	Eigen::Vector3d target(1.0, 0.0, 0.0);
	
	// Cost function: squared distance to target
	SE3var T = SE3var::exp(xi);
	auto pos = T.translation();
	var dx = pos[0] - target[0];
	var dy = pos[1] - target[1];
	var dz = pos[2] - target[2];
	var cost = dx*dx + dy*dy + dz*dz;
	
	// Compute all gradients in ONE reverse pass (efficient!)
	Derivatives dcost = derivatives(cost);
	
	Eigen::Matrix<double, 6, 1> gradient;
	for (int i = 0; i < 6; i++) {
		gradient[i] = dcost(xi[i]);
	}
	
	// Gradient should be non-zero
	EXPECT_GT(gradient.norm(), 1e-10);
	
	// Gradient should point towards decreasing cost
	// Since we're too far in x (0.5 < 1.0), gradient[3] (ρx) should be positive
	EXPECT_GT(gradient[3], 0.0);
	
	std::cout << "Reverse mode gradient of distance cost w.r.t. xi: " 
	          << gradient.transpose() << std::endl;
	std::cout << "Initial position: " << toDouble(pos).transpose() << std::endl;
	std::cout << "Target position:  " << target.transpose() << std::endl;
	std::cout << "Cost: " << val(cost) << std::endl;
}

/**
 * @brief Test reverse mode for Cosserat strain optimization
 * 
 * Simulate a simple Cosserat rod with one segment and compute the gradient
 * of a tip position cost w.r.t. strain parameters.
 */
TEST(AutodiffIntegration, ReverseMode_CosseratStrain) {
	using SE3var = SE3<var>;
	
	// Strain parameters: [φx, φy, φz, ρx, ρy, ρz]
	Eigen::Matrix<double, 6, 1> strain_val;
	strain_val << 0.0, 0.1, 0.0, 0.0, 0.0, 0.0;  // Bending in Y
	
	// Convert to var
	Eigen::Matrix<var, 6, 1> strain;
	for (int i = 0; i < 6; i++) {
		strain[i] = strain_val[i];
	}
	
	// Length of segment
	double L = 1.0;
	
	// Compute tip position: T = exp(L * strain)
	Eigen::Matrix<var, 6, 1> xi = L * strain;
	SE3var T = SE3var::exp(xi);
	auto pos = T.translation();
	
	// Target: tip should be at (1, 0, 0) - straight rod
	Eigen::Vector3d target(1.0, 0.0, 0.0);
	
	// Cost: distance to target
	var dx = pos[0] - target[0];
	var dy = pos[1] - target[1];
	var dz = pos[2] - target[2];
	var cost = dx*dx + dy*dy + dz*dz;
	
	// Compute gradient w.r.t. strain
	Derivatives dcost = derivatives(cost);
	
	Eigen::Matrix<double, 6, 1> gradient;
	for (int i = 0; i < 6; i++) {
		gradient[i] = dcost(strain[i]);
	}
	
	// Since we have bending in Y (φy > 0), gradient[1] should be non-zero
	EXPECT_GT(std::abs(gradient[1]), 1e-6);
	
	std::cout << "Strain gradient for tip positioning: " << gradient.transpose() << std::endl;
	std::cout << "Tip position: " << toDouble(pos).transpose() << std::endl;
	std::cout << "Target:       " << target.transpose() << std::endl;
	std::cout << "Cost:         " << val(cost) << std::endl;
}

/**
 * @brief Test reverse mode for multi-segment Cosserat rod
 * 
 * Chain multiple SE3 transformations and compute gradient through the entire chain.
 */
TEST(AutodiffIntegration, ReverseMode_CosseratChain) {
	using SE3var = SE3<var>;
	
	const int N = 3;  // 3 segments
	const double L = 0.5;  // Length per segment
	
	// Strain for each segment
	std::vector<Eigen::Matrix<var, 6, 1>> strains(N);
	std::vector<Eigen::Matrix<double, 6, 1>> strain_vals(N);
	
	// Initialize with small perturbations
	for (int i = 0; i < N; i++) {
		strain_vals[i] << 0.0, 0.1 * i, 0.0, 0.0, 0.0, 0.0;
		for (int j = 0; j < 6; j++) {
			strains[i][j] = strain_vals[i][j];
		}
	}
	
	// Forward kinematics: T = T1 * T2 * T3
	SE3var T = SE3var::identity();
	for (int i = 0; i < N; i++) {
		Eigen::Matrix<var, 6, 1> xi = L * strains[i];
		SE3var Ti = SE3var::exp(xi);
		T = T * Ti;
	}
	
	// Cost: tip should reach target
	Eigen::Vector3d target(1.5, 0.0, 0.0);
	auto pos = T.translation();
	var dx = pos[0] - target[0];
	var dy = pos[1] - target[1];
	var dz = pos[2] - target[2];
	var cost = dx*dx + dy*dy + dz*dz;
	
	// Compute gradients for ALL strains in one pass (this is the power of reverse mode!)
	Derivatives dcost = derivatives(cost);
	
	std::vector<Eigen::Matrix<double, 6, 1>> gradients(N);
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < 6; j++) {
			gradients[i][j] = dcost(strains[i][j]);
		}
	}
	
	// Display results
	std::cout << "\n=== Multi-Segment Cosserat Chain ===" << std::endl;
	std::cout << "Tip position: " << toDouble(pos).transpose() << std::endl;
	std::cout << "Target:       " << target.transpose() << std::endl;
	std::cout << "Cost:         " << val(cost) << std::endl;
	
	for (int i = 0; i < N; i++) {
		std::cout << "Gradient[" << i << "]: " << gradients[i].transpose() << std::endl;
	}
	
	// At least one gradient should be significant
	bool has_significant_gradient = false;
	for (int i = 0; i < N; i++) {
		if (gradients[i].norm() > 1e-6) {
			has_significant_gradient = true;
			break;
		}
	}
	EXPECT_TRUE(has_significant_gradient);
}

// ============================================================================
// Comparison with Analytical Jacobians
// ============================================================================

/**
 * @brief Compare autodiff gradients with analytical jacobians
 * 
 * For SE3, compare:
 * - Autodiff: d(exp(xi).translation)/d(xi) via reverse mode
 * - Analytical: Using Phase 2 jacobians
 */
TEST(AutodiffIntegration, CompareWithAnalytical_SE3) {
	using SE3d = SE3<double>;
	using SE3var = SE3<var>;
	
	Eigen::Matrix<double, 6, 1> xi_val;
	xi_val << 0.1, 0.2, 0.3, 0.4, 0.5, 0.6;
	
	// === Analytical jacobian (from Phase 2) ===
	SE3d T = SE3d::exp(xi_val);
	Eigen::Matrix<double, 6, 6> J_left = SE3d::leftJacobian(xi_val);
	
	// Translation part: d(translation)/d(xi) = J_left.bottomRows<3>()
	Eigen::Matrix<double, 3, 6> J_analytical = J_left.bottomRows<3>();
	
	// === Autodiff jacobian ===
	Eigen::Matrix<var, 6, 1> xi;
	for (int i = 0; i < 6; i++) {
		xi[i] = xi_val[i];
	}
	
	SE3var T_var = SE3var::exp(xi);
	auto pos = T_var.translation();
	
	Eigen::Matrix<double, 3, 6> J_autodiff;
	for (int row = 0; row < 3; row++) {
		Derivatives dpos = derivatives(pos[row]);
		for (int col = 0; col < 6; col++) {
			J_autodiff(row, col) = dpos(xi[col]);
		}
	}
	
	// Compare
	double error = (J_analytical - J_autodiff).norm() / J_analytical.norm();
	
	std::cout << "\n=== Analytical vs Autodiff Jacobian ===" << std::endl;
	std::cout << "Relative error: " << error << std::endl;
	std::cout << "Analytical:\n" << J_analytical << std::endl;
	std::cout << "Autodiff:\n" << J_autodiff << std::endl;
	
	EXPECT_LT(error, 1e-6);
}

#else

// Dummy test when autodiff is not enabled
TEST(AutodiffIntegration, NotEnabled) {
	GTEST_SKIP() << "Autodiff support not enabled. Compile with -DCOSSERAT_WITH_AUTODIFF=ON";
}

#endif // COSSERAT_WITH_AUTODIFF
