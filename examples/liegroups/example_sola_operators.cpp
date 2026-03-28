/**
 * @file example_sola_operators.cpp
 * @brief Example demonstrating the new plus/minus operators and Jacobians
 * 
 * Based on "A micro Lie theory for state estimation in robotics" (Solà et al., 2021)
 * 
 * This example shows:
 * 1. Using ⊕/⊖ operators for intuitive state updates
 * 2. Computing right and left Jacobians
 * 3. Uncertainty propagation through composition
 */

#include <iostream>
#include <iomanip>
#include <Eigen/Dense>
#include "../../src/liegroups/SO3.h"
#include "../../src/liegroups/SE3.h"

using namespace sofa::component::cosserat::liegroups;

void printSeparator(const std::string& title) {
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "  " << title << "\n";
    std::cout << std::string(60, '=') << "\n\n";
}

void example1_PlusMinusOperators() {
    printSeparator("Example 1: Plus/Minus Operators");
    
    // Initial rotation
    Eigen::Vector3d omega_init(0.3, 0.2, 0.1);
    SO3<double> R = SO3<double>::exp(omega_init);
    
    std::cout << "Initial rotation (angle-axis): " 
              << omega_init.transpose() << "\n\n";
    
    // Apply an increment using the new plus operator
    Eigen::Vector3d delta(0.1, 0.05, 0.02);
    
    std::cout << "Applying increment delta = " << delta.transpose() << "\n\n";
    
    // Old way (still works)
    SO3<double> R_new_old = R.compose(SO3<double>::exp(delta));
    
    // New way (more intuitive!)
    SO3<double> R_new = R + delta;  // or R.plus(delta)
    
    std::cout << "Old syntax: R.compose(SO3::exp(delta))\n";
    std::cout << "New syntax: R + delta\n";
    std::cout << "Results match: " << (R_new.isApprox(R_new_old) ? "YES ✓" : "NO ✗") << "\n\n";
    
    // Compute the difference back
    Eigen::Vector3d delta_recovered = R_new - R;  // or R_new.minus(R)
    
    std::cout << "Recovered delta: " << delta_recovered.transpose() << "\n";
    std::cout << "Original delta:  " << delta.transpose() << "\n";
    std::cout << "Match: " << (delta.isApprox(delta_recovered, 1e-10) ? "YES ✓" : "NO ✗") << "\n";
}

void example2_RightJacobian() {
    printSeparator("Example 2: Right Jacobian Jr(τ)");
    
    Eigen::Vector3d omega(0.5, 0.3, 0.2);
    
    std::cout << "Computing right Jacobian for ω = " << omega.transpose() << "\n\n";
    
    // Compute right Jacobian
    auto Jr = SO3<double>::rightJacobian(omega);
    auto Jr_inv = SO3<double>::rightJacobianInverse(omega);
    
    std::cout << "Jr(ω) =\n" << Jr << "\n\n";
    std::cout << "Jr⁻¹(ω) =\n" << Jr_inv << "\n\n";
    
    // Verify Jr * Jr⁻¹ = I
    auto product = Jr * Jr_inv;
    std::cout << "Jr(ω) * Jr⁻¹(ω) =\n" << product << "\n\n";
    std::cout << "Is identity: " 
              << (product.isApprox(Eigen::Matrix3d::Identity(), 1e-10) ? "YES ✓" : "NO ✗") << "\n";
}

void example3_LeftJacobian() {
    printSeparator("Example 3: Left Jacobian Jl(τ)");
    
    Eigen::Vector3d omega(0.4, -0.2, 0.3);
    
    std::cout << "Computing left Jacobian for ω = " << omega.transpose() << "\n\n";
    
    // Compute left Jacobian
    auto Jl = SO3<double>::leftJacobian(omega);
    auto Jr = SO3<double>::rightJacobian(omega);
    
    std::cout << "Jl(ω) =\n" << Jl << "\n\n";
    std::cout << "Jr(ω) =\n" << Jr << "\n\n";
    
    // Verify Jl(ω) = Jr(-ω)
    auto Jr_minus = SO3<double>::rightJacobian(-omega);
    std::cout << "Jr(-ω) =\n" << Jr_minus << "\n\n";
    std::cout << "Jl(ω) = Jr(-ω): " 
              << (Jl.isApprox(Jr_minus, 1e-10) ? "YES ✓" : "NO ✗") << "\n";
}

void example4_UncertaintyPropagation() {
    printSeparator("Example 4: Uncertainty Propagation");
    
    std::cout << "Demonstrating covariance propagation through composition\n\n";
    
    // Initial state with uncertainty
    Eigen::Vector3d omega_mean(0.2, 0.1, 0.15);
    SO3<double> R_mean = SO3<double>::exp(omega_mean);
    Eigen::Matrix3d Sigma_R = 0.01 * Eigen::Matrix3d::Identity();  // Covariance
    
    std::cout << "Initial state:\n";
    std::cout << "  Mean (angle-axis): " << omega_mean.transpose() << "\n";
    std::cout << "  Covariance trace: " << Sigma_R.trace() << "\n\n";
    
    // Apply increment with its own uncertainty
    Eigen::Vector3d delta_mean(0.05, 0.02, 0.01);
    Eigen::Matrix3d Sigma_delta = 0.005 * Eigen::Matrix3d::Identity();
    
    std::cout << "Increment:\n";
    std::cout << "  Mean: " << delta_mean.transpose() << "\n";
    std::cout << "  Covariance trace: " << Sigma_delta.trace() << "\n\n";
    
    // Update mean
    SO3<double> R_new = R_mean + delta_mean;
    
    // Propagate covariance: Σ_new = F * Σ_R * F^T + G * Σ_delta * G^T
    // where F = Ad_Exp(delta)^-1 and G = Jr(delta)
    auto F = SO3<double>::exp(delta_mean).inverse().adjoint();
    auto G = SO3<double>::rightJacobian(delta_mean);
    
    Eigen::Matrix3d Sigma_new = F * Sigma_R * F.transpose() + G * Sigma_delta * G.transpose();
    
    std::cout << "After composition:\n";
    std::cout << "  New mean (angle-axis): " << R_new.log().transpose() << "\n";
    std::cout << "  New covariance trace: " << Sigma_new.trace() << "\n\n";
    
    std::cout << "Jacobian F (wrt state) =\n" << F << "\n\n";
    std::cout << "Jacobian G (wrt increment) =\n" << G << "\n\n";
    
    std::cout << "Note: This is the foundation for ESKF (Error-State Kalman Filter)\n";
}

void example5_SmallAngleCase() {
    printSeparator("Example 5: Small Angle Approximations");
    
    std::cout << "Testing small angle behavior\n\n";
    
    Eigen::Vector3d omega_small(1e-6, 2e-6, 5e-7);
    
    std::cout << "Very small angle: " << omega_small.transpose() << "\n";
    std::cout << "Norm: " << omega_small.norm() << " (< ε = " 
              << Types<double>::epsilon() << ")\n\n";
    
    // Compute Jacobians (should use approximations internally)
    auto Jr = SO3<double>::rightJacobian(omega_small);
    auto Jr_inv = SO3<double>::rightJacobianInverse(omega_small);
    
    std::cout << "Jr(ω) ≈ I - ½[ω]× for small ω\n";
    std::cout << "Jr(ω) =\n" << Jr << "\n\n";
    
    std::cout << "Jr⁻¹(ω) ≈ I + ½[ω]× for small ω\n";
    std::cout << "Jr⁻¹(ω) =\n" << Jr_inv << "\n\n";
    
    // Verify they're still inverses
    auto product = Jr * Jr_inv;
    std::cout << "Jr * Jr⁻¹ still equals I: "
              << (product.isApprox(Eigen::Matrix3d::Identity(), 1e-10) ? "YES ✓" : "NO ✗") << "\n";
}

void example6_ESKFPredict() {
    printSeparator("Example 6: ESKF Prediction Step");
    
    std::cout << "Simulating an ESKF prediction step for attitude estimation\n\n";
    
    // Current state
    Eigen::Vector3d omega_current(0.3, 0.1, 0.2);
    SO3<double> R_current = SO3<double>::exp(omega_current);
    Eigen::Matrix3d P_current = Eigen::Matrix3d::Identity() * 0.01;
    
    std::cout << "Current state:\n";
    std::cout << "  Attitude: " << omega_current.transpose() << "\n";
    std::cout << "  Covariance trace: " << P_current.trace() << "\n\n";
    
    // Control input (measured angular velocity * dt)
    double dt = 0.01;  // 10ms
    Eigen::Vector3d omega_meas(0.5, -0.2, 0.3);
    Eigen::Vector3d u = omega_meas * dt;
    Eigen::Matrix3d Q = Eigen::Matrix3d::Identity() * 0.001;  // Process noise
    
    std::cout << "Control input:\n";
    std::cout << "  Measured ω: " << omega_meas.transpose() << " rad/s\n";
    std::cout << "  dt: " << dt << " s\n";
    std::cout << "  u = ω·dt: " << u.transpose() << "\n\n";
    
    // ESKF Prediction
    std::cout << "ESKF Prediction:\n\n";
    
    // State propagation
    SO3<double> R_predicted = R_current + u;  // Using new plus operator!
    
    std::cout << "  X̂ₖ = X̂ₖ₋₁ ⊕ u  (intuitive notation!)\n\n";
    
    // Covariance propagation
    auto F = SO3<double>::exp(u).inverse().adjoint();  // Ad_Exp(u)^-1
    auto G = SO3<double>::rightJacobian(u);             // Jr(u)
    
    Eigen::Matrix3d P_predicted = F * P_current * F.transpose() + G * Q * G.transpose();
    
    std::cout << "  Pₖ = F·Pₖ₋₁·Fᵀ + G·Q·Gᵀ\n";
    std::cout << "  where F = Ad_Exp(u)⁻¹, G = Jr(u)\n\n";
    
    std::cout << "Predicted state:\n";
    std::cout << "  Attitude: " << R_predicted.log().transpose() << "\n";
    std::cout << "  Covariance trace: " << P_predicted.trace() << "\n\n";
    
    std::cout << "Comparison with naive approach:\n";
    SO3<double> R_naive = R_current.compose(SO3<double>::exp(u));
    std::cout << "  States match: " << (R_predicted.isApprox(R_naive) ? "YES ✓" : "NO ✗") << "\n";
    std::cout << "  But new syntax is much clearer!\n";
}

int main() {
    std::cout << std::fixed << std::setprecision(6);
    
    std::cout << "\n";
    std::cout << "╔════════════════════════════════════════════════════════════╗\n";
    std::cout << "║  Lie Theory Operators Demo                                ║\n";
    std::cout << "║  Based on Solà et al. (2021)                              ║\n";
    std::cout << "╚════════════════════════════════════════════════════════════╝\n";
    
    try {
        example1_PlusMinusOperators();
        example2_RightJacobian();
        example3_LeftJacobian();
        example4_UncertaintyPropagation();
        example5_SmallAngleCase();
        example6_ESKFPredict();
        
        printSeparator("Summary");
        std::cout << "All examples completed successfully! ✓\n\n";
        std::cout << "Key takeaways:\n";
        std::cout << "1. ⊕/⊖ operators make code more readable and intuitive\n";
        std::cout << "2. Jr/Jl Jacobians are essential for uncertainty propagation\n";
        std::cout << "3. Small angle approximations ensure numerical stability\n";
        std::cout << "4. ESKF implementation becomes straightforward\n\n";
        std::cout << "Next steps:\n";
        std::cout << "- Implement SE(3) Jacobians (see AMELIORATIONS_PROPOSEES.md)\n";
        std::cout << "- Add Jacobians for elementary operations\n";
        std::cout << "- Create GaussianOnManifold class\n\n";
        
        return 0;
    }
    catch (const std::exception& e) {
        std::cerr << "\n❌ Error: " << e.what() << "\n";
        return 1;
    }
}
