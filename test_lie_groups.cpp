// Test program for Lie group implementations
#include <iostream>
#include <Eigen/Dense>
#include "src/liegroups/Types.h"
#include "src/liegroups/LieGroupBase.h"
#include "src/liegroups/SO3.h"
#include "src/liegroups/SE3.h"
#include "src/liegroups/Sim3.h"
#include "src/liegroups/Uncertainty.h"

int main() {
    using namespace sofa::component::cosserat::liegroups;

    std::cout << "Testing Lie Group Implementations" << std::endl;
    std::cout << "=================================" << std::endl;

    // Test SO3
    std::cout << "\nTesting SO3:" << std::endl;
    Eigen::Vector3d omega(0.1, 0.2, 0.3);
    SO3<double> R = SO3<double>::exp(omega);
    std::cout << "SO3 exp([0.1, 0.2, 0.3]): " << R << std::endl;

    // Test SE3
    std::cout << "\nTesting SE3:" << std::endl;
    Eigen::Matrix<double, 6, 1> xi;
    xi << 0.1, 0.2, 0.3, 1.0, 2.0, 3.0; // rotation + translation
    SE3<double> T = SE3<double>::exp(xi);
    std::cout << "SE3 exp([0.1, 0.2, 0.3, 1.0, 2.0, 3.0]):" << std::endl;
    std::cout << T.matrix() << std::endl;

    // Test SE3 Jacobians
    std::cout << "\nTesting SE3 Jacobians:" << std::endl;
    auto Jr = SE3<double>::rightJacobian(xi);
    std::cout << "Right Jacobian shape: " << Jr.rows() << "x" << Jr.cols() << std::endl;

    // Test Sim3
    std::cout << "\nTesting Sim3:" << std::endl;
    Eigen::Matrix<double, 7, 1> xi_sim3;
    xi_sim3 << 0.1, 0.2, 0.3, 1.0, 2.0, 3.0, 0.1; // rotation + translation + log(scale)
    Sim3<double> S = Sim3<double>::exp(xi_sim3);
    std::cout << "Sim3 exp result - scale: " << S.scale() << std::endl;

    // Test Uncertainty propagation
    std::cout << "\nTesting Uncertainty Propagation:" << std::endl;
    using SE3d = SE3<double>;
    SE3d mean = SE3d::Identity();
    Eigen::Matrix<double, 6, 6> cov = Eigen::Matrix<double, 6, 6>::Identity() * 0.01;
    GaussianOnManifold<SE3d> gaussian(mean, cov);
    std::cout << "Created Gaussian distribution on SE3" << std::endl;

    std::cout << "\nAll tests completed successfully!" << std::endl;
    return 0;
}