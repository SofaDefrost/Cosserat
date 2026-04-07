/**
 * Test rapide des gradients de l'optimiseur
 */

#include <iostream>
#include <vector>
#include <Eigen/Dense>

#include "../src/liegroups/SE3.h"
#include "../src/liegroups/optimization/CosseratTrajectoryOptimizer.h"

using namespace sofa::component::cosserat::liegroups;
using namespace sofa::component::cosserat::liegroups::optimization;

int main() {
    // Configuration simple
    const int n_sections = 3;
    const double length = 0.1;
    
    // Convention strains Cosserat: [φx, φy, φz, ρx, ρy, ρz]
    // φx=torsion, φy=bending_Y, φz=bending_Z
    // ρx=elongation, ρy=shearing_Y, ρz=shearing_Z
    std::vector<Eigen::Matrix<double, 6, 1>> strains(n_sections);
    strains[0] << 0.1, 0.0, 0.0, 0.0, 0.0, 0.0;  // Torsion (rotation autour X)
    strains[1] << 0.0, 0.0, 0.0, 0.1, 0.0, 0.0;  // Elongation en X
    strains[2] << 0.05, 0.0, 0.0, 0.0, 0.0, 0.0; // Torsion
    
    // Cible
    SE3d target = SE3d::Identity();
    target.translation() = Eigen::Vector3d(0.3, 0.0, 0.0);
    
    // Créer optimiseur
    CosseratTrajectoryOptimizer<double> optimizer;
    CosseratTrajectoryOptimizer<double>::Parameters params;
    params.regularization = 0.0;  // Pas de régularisation pour ce test
    
    // Forward kinematics pour calculer position actuelle
    SE3d g = SE3d::Identity();
    for (const auto& strain : strains) {
        g = g * SE3d::exp(strain * length);
    }
    
    std::cout << "Position actuelle: " << g.translation().transpose() << std::endl;
    std::cout << "Cible: " << target.translation().transpose() << std::endl;
    
    // Calculer le coût et gradient analytique
    auto cost_gradient = [&](const std::vector<Eigen::Matrix<double, 6, 1>>& s) {
        SE3d g_local = SE3d::Identity();
        for (const auto& strain : s) {
            g_local = g_local * SE3d::exp(strain * length);
        }
        Eigen::Vector3d error = g_local.translation() - target.translation();
        return 0.5 * error.squaredNorm();
    };
    
    double cost = cost_gradient(strains);
    std::cout << "\nCoût initial: " << cost << std::endl;
    
    // Test: gradient par différences finies
    std::cout << "\n=== Test des Gradients ===" << std::endl;
    const double h = 1e-7;
    
    for (int i = 0; i < n_sections; ++i) {
        std::cout << "\nSection " << i << ":" << std::endl;
        for (int j = 0; j < 6; ++j) {
            // Perturbation +h
            auto strains_plus = strains;
            strains_plus[i](j) += h;
            double cost_plus = cost_gradient(strains_plus);
            
            // Perturbation -h
            auto strains_minus = strains;
            strains_minus[i](j) -= h;
            double cost_minus = cost_gradient(strains_minus);
            
            // Gradient numérique
            double grad_num = (cost_plus - cost_minus) / (2.0 * h);
            
            if (std::abs(grad_num) > 1e-10) {
                std::cout << "  strain[" << j << "]: grad = " << grad_num << std::endl;
            }
        }
    }
    
    return 0;
}
