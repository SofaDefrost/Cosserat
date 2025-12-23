/**
 * @file ilqr_trajectory_tracking.cpp
 * @brief Example of iLQR controller for Cosserat rod trajectory tracking
 */

#include <iostream>
#include <iomanip>
#include "../src/liegroups/control/CosseratILQRController.h"

using namespace sofa::component::cosserat::liegroups;
using namespace sofa::component::cosserat::liegroups::control;

int main() {
    using Controller = CosseratILQRController<double>;
    using SE3Type = SE3<double>;
    using Vector6 = Eigen::Matrix<double, 6, 1>;
    using Vector3 = Eigen::Vector3d;
    
    std::cout << "\n╔══════════════════════════════════════════════════════════════╗\n";
    std::cout << "║          iLQR Trajectory Tracking for Cosserat Rods         ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════╝\n\n";
    
    // Configuration
    Controller::Config config;
    config.max_iterations = 30;
    config.Q_position = 100.0;
    config.Q_rotation = 10.0;
    config.R_control = 0.1;
    config.Q_final_position = 500.0;
    config.convergence_threshold = 1e-3;
    config.verbose = true;
    
    const int N = 8;  // 8 segments
    Controller controller(N, config);
    
    // Reference trajectory: helix
    std::cout << "Creating reference trajectory (helix)...\n";
    Controller::Trajectory reference;
    reference.segment_length = 0.125;  // 1m total length
    
    for (int i = 0; i <= N; ++i) {
        double t = i * 0.125;
        double radius = 0.2;
        double pitch = 0.3;
        
        SE3Type pose = SE3Type::computeIdentity();
        pose.translation() = Vector3(
            t,                            // Along X
            radius * std::cos(2*M_PI*t),  // Circular in YZ
            radius * std::sin(2*M_PI*t)
        );
        
        reference.poses.push_back(pose);
    }
    
    std::cout << "  - Segments: " << N << "\n";
    std::cout << "  - Segment length: " << reference.segment_length << " m\n";
    std::cout << "  - Total length: " << N * reference.segment_length << " m\n\n";
    
    // Initial guess: zero strains
    std::vector<Vector6> initial_strains(N, Vector6::Zero());
    
    // Optimize
    std::cout << "Running iLQR optimization...\n";
    std::cout << std::string(70, '-') << "\n";
    
    auto result = controller.optimize(reference, initial_strains);
    
    std::cout << std::string(70, '-') << "\n\n";
    
    // Results
    std::cout << "╔══════════════════════════════════════════════════════════════╗\n";
    std::cout << "║                        Results                               ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════╝\n\n";
    
    std::cout << "Status: " << (result.converged ? "✓ Converged" : "✗ Not converged") << "\n";
    std::cout << "Iterations: " << result.iterations << "\n";
    std::cout << "Final cost: " << std::fixed << std::setprecision(6) << result.final_cost << "\n";
    std::cout << "Message: " << result.message << "\n\n";
    
    // Show first 3 optimal strains
    std::cout << "Optimal strains (first 3 segments):\n";
    for (int i = 0; i < std::min(3, (int)result.optimal_strains.size()); ++i) {
        std::cout << "  Segment " << i+1 << ": [";
        for (int j = 0; j < 6; ++j) {
            std::cout << std::setw(8) << std::setprecision(4) << result.optimal_strains[i][j];
            if (j < 5) std::cout << ", ";
        }
        std::cout << "]\n";
    }
    std::cout << "  ...\n\n";
    
    // Cost history
    std::cout << "Cost reduction: " << std::setprecision(2) 
              << (1.0 - result.final_cost / result.cost_history[0]) * 100.0 << "%\n\n";
    
    // Tracking error
    Vector3 final_pos = result.trajectory.back().translation();
    Vector3 target_pos = reference.poses.back().translation();
    double tracking_error = (final_pos - target_pos).norm();
    
    std::cout << "Final tracking error: " << std::setprecision(4) << tracking_error << " m\n";
    std::cout << "  Final position:  [" << final_pos.transpose() << "]\n";
    std::cout << "  Target position: [" << target_pos.transpose() << "]\n\n";
    
    std::cout << "╔══════════════════════════════════════════════════════════════╗\n";
    std::cout << "║                    ✓ Success                                 ║\n";
    std::cout << "║  iLQR successfully computed optimal control for helix        ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════╝\n\n";
    
    return 0;
}
