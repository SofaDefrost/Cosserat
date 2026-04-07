#pragma once

#include <cmath>
#include <functional>
#include <iostream>
#include <Eigen/Dense>
#include "../../Types.h"

namespace sofa::component::cosserat::liegroups::testing {

/**
 * @brief Utility class for testing differentiability of Lie group operations
 * 
 * This class provides methods to compute numerical derivatives using finite differences
 * and compare them with analytical Jacobians.
 */
template<typename Scalar = double>
class DifferentiationTestUtils {
public:
    using Types = liegroups::Types<Scalar>;
    
    /**
     * @brief Compute numerical gradient using forward finite differences
     * 
     * @param f Function to differentiate: R^n -> R
     * @param x Point at which to compute gradient
     * @param h Step size (default: sqrt(epsilon))
     * @return Gradient vector
     */
    template<int N>
    static Eigen::Matrix<Scalar, N, 1> 
    finiteDifferenceGradient(
        const std::function<Scalar(const Eigen::Matrix<Scalar, N, 1>&)>& f,
        const Eigen::Matrix<Scalar, N, 1>& x,
        Scalar h = std::sqrt(Types::epsilon())
    ) {
        Eigen::Matrix<Scalar, N, 1> grad;
        const Scalar f_x = f(x);
        
        for (int i = 0; i < N; ++i) {
            Eigen::Matrix<Scalar, N, 1> x_plus = x;
            x_plus(i) += h;
            grad(i) = (f(x_plus) - f_x) / h;
        }
        
        return grad;
    }
    
    /**
     * @brief Compute numerical gradient using central finite differences (more accurate)
     * 
     * @param f Function to differentiate: R^n -> R
     * @param x Point at which to compute gradient
     * @param h Step size (default: cbrt(epsilon))
     * @return Gradient vector
     */
    template<int N>
    static Eigen::Matrix<Scalar, N, 1> 
    centralDifferenceGradient(
        const std::function<Scalar(const Eigen::Matrix<Scalar, N, 1>&)>& f,
        const Eigen::Matrix<Scalar, N, 1>& x,
        Scalar h = std::cbrt(Types::epsilon())
    ) {
        Eigen::Matrix<Scalar, N, 1> grad;
        
        for (int i = 0; i < N; ++i) {
            Eigen::Matrix<Scalar, N, 1> x_plus = x;
            Eigen::Matrix<Scalar, N, 1> x_minus = x;
            x_plus(i) += h;
            x_minus(i) -= h;
            grad(i) = (f(x_plus) - f(x_minus)) / (Scalar(2) * h);
        }
        
        return grad;
    }
    
    /**
     * @brief Compute numerical Jacobian using forward finite differences
     * 
     * @param f Function to differentiate: R^n -> R^m
     * @param x Point at which to compute Jacobian
     * @param h Step size
     * @return Jacobian matrix (m x n)
     */
    template<int M, int N>
    static Eigen::Matrix<Scalar, M, N> 
    finiteDifferenceJacobian(
        const std::function<Eigen::Matrix<Scalar, M, 1>(const Eigen::Matrix<Scalar, N, 1>&)>& f,
        const Eigen::Matrix<Scalar, N, 1>& x,
        Scalar h = std::sqrt(Types::epsilon())
    ) {
        Eigen::Matrix<Scalar, M, N> jacobian;
        const Eigen::Matrix<Scalar, M, 1> f_x = f(x);
        
        for (int j = 0; j < N; ++j) {
            Eigen::Matrix<Scalar, N, 1> x_plus = x;
            x_plus(j) += h;
            jacobian.col(j) = (f(x_plus) - f_x) / h;
        }
        
        return jacobian;
    }
    
    /**
     * @brief Compute numerical Jacobian using central finite differences
     * 
     * @param f Function to differentiate: R^n -> R^m
     * @param x Point at which to compute Jacobian
     * @param h Step size
     * @return Jacobian matrix (m x n)
     */
    template<int M, int N>
    static Eigen::Matrix<Scalar, M, N> 
    centralDifferenceJacobian(
        const std::function<Eigen::Matrix<Scalar, M, 1>(const Eigen::Matrix<Scalar, N, 1>&)>& f,
        const Eigen::Matrix<Scalar, N, 1>& x,
        Scalar h = std::cbrt(Types::epsilon())
    ) {
        Eigen::Matrix<Scalar, M, N> jacobian;
        
        for (int j = 0; j < N; ++j) {
            Eigen::Matrix<Scalar, N, 1> x_plus = x;
            Eigen::Matrix<Scalar, N, 1> x_minus = x;
            x_plus(j) += h;
            x_minus(j) -= h;
            jacobian.col(j) = (f(x_plus) - f(x_minus)) / (Scalar(2) * h);
        }
        
        return jacobian;
    }
    
    /**
     * @brief Compare two matrices and check if they are approximately equal
     * 
     * @param analytical Analytically computed matrix
     * @param numerical Numerically computed matrix
     * @param tolerance Relative tolerance
     * @param verbose Print detailed comparison
     * @return true if matrices are approximately equal
     */
    template<int M, int N>
    static bool compareMatrices(
        const Eigen::Matrix<Scalar, M, N>& analytical,
        const Eigen::Matrix<Scalar, M, N>& numerical,
        Scalar tolerance = Scalar(1e-5),
        bool verbose = false
    ) {
        const Scalar max_diff = (analytical - numerical).cwiseAbs().maxCoeff();
        const Scalar analytical_norm = analytical.norm();
        const Scalar relative_error = analytical_norm > Types::epsilon() 
            ? max_diff / analytical_norm 
            : max_diff;
        
        const bool passed = relative_error < tolerance;
        
        if (verbose || !passed) {
            std::cout << "=== Matrix Comparison ===" << std::endl;
            std::cout << "Analytical:\n" << analytical << std::endl;
            std::cout << "Numerical:\n" << numerical << std::endl;
            std::cout << "Difference:\n" << (analytical - numerical) << std::endl;
            std::cout << "Max absolute error: " << max_diff << std::endl;
            std::cout << "Relative error: " << relative_error << std::endl;
            std::cout << "Tolerance: " << tolerance << std::endl;
            std::cout << "Status: " << (passed ? "PASSED ✓" : "FAILED ✗") << std::endl;
            std::cout << "========================" << std::endl;
        }
        
        return passed;
    }
    
    /**
     * @brief Test if a Jacobian is correct using finite differences
     * 
     * @param f Function R^n -> R^m
     * @param analytical_jacobian Analytically computed Jacobian
     * @param x Point at which to test
     * @param tolerance Relative tolerance
     * @param use_central Use central differences (more accurate but slower)
     * @return true if Jacobian is correct
     */
    template<int M, int N>
    static bool testJacobian(
        const std::function<Eigen::Matrix<Scalar, M, 1>(const Eigen::Matrix<Scalar, N, 1>&)>& f,
        const Eigen::Matrix<Scalar, M, N>& analytical_jacobian,
        const Eigen::Matrix<Scalar, N, 1>& x,
        Scalar tolerance = Scalar(1e-5),
        bool use_central = true,
        bool verbose = false
    ) {
        Eigen::Matrix<Scalar, M, N> numerical_jacobian;
        
        if (use_central) {
            numerical_jacobian = centralDifferenceJacobian<M, N>(f, x);
        } else {
            numerical_jacobian = finiteDifferenceJacobian<M, N>(f, x);
        }
        
        if (verbose) {
            std::cout << "\n=== Jacobian Test ===" << std::endl;
            std::cout << "Testing at point: " << x.transpose() << std::endl;
            std::cout << "Method: " << (use_central ? "Central" : "Forward") << " differences" << std::endl;
        }
        
        return compareMatrices<M, N>(analytical_jacobian, numerical_jacobian, tolerance, verbose);
    }
    
    /**
     * @brief Dynamic version of Jacobian testing for runtime-sized matrices
     */
    static bool testJacobianDynamic(
        const std::function<Eigen::VectorXd(const Eigen::VectorXd&)>& f,
        const Eigen::MatrixXd& analytical_jacobian,
        const Eigen::VectorXd& x,
        Scalar tolerance = Scalar(1e-5),
        bool use_central = true,
        bool verbose = false
    ) {
        const int m = analytical_jacobian.rows();
        const int n = analytical_jacobian.cols();
        
        Eigen::MatrixXd numerical_jacobian(m, n);
        const Scalar h = use_central ? std::cbrt(Types::epsilon()) : std::sqrt(Types::epsilon());
        
        if (use_central) {
            for (int j = 0; j < n; ++j) {
                Eigen::VectorXd x_plus = x;
                Eigen::VectorXd x_minus = x;
                x_plus(j) += h;
                x_minus(j) -= h;
                numerical_jacobian.col(j) = (f(x_plus) - f(x_minus)) / (Scalar(2) * h);
            }
        } else {
            const Eigen::VectorXd f_x = f(x);
            for (int j = 0; j < n; ++j) {
                Eigen::VectorXd x_plus = x;
                x_plus(j) += h;
                numerical_jacobian.col(j) = (f(x_plus) - f_x) / h;
            }
        }
        
        const Scalar max_diff = (analytical_jacobian - numerical_jacobian).cwiseAbs().maxCoeff();
        const Scalar analytical_norm = analytical_jacobian.norm();
        const Scalar relative_error = analytical_norm > Types::epsilon() 
            ? max_diff / analytical_norm 
            : max_diff;
        
        const bool passed = relative_error < tolerance;
        
        if (verbose || !passed) {
            std::cout << "\n=== Jacobian Test (Dynamic) ===" << std::endl;
            std::cout << "Dimensions: " << m << " x " << n << std::endl;
            std::cout << "Testing at point: " << x.transpose() << std::endl;
            std::cout << "Method: " << (use_central ? "Central" : "Forward") << " differences" << std::endl;
            std::cout << "Max absolute error: " << max_diff << std::endl;
            std::cout << "Relative error: " << relative_error << std::endl;
            std::cout << "Tolerance: " << tolerance << std::endl;
            std::cout << "Status: " << (passed ? "PASSED ✓" : "FAILED ✗") << std::endl;
            std::cout << "==============================" << std::endl;
        }
        
        return passed;
    }
    
    /**
     * @brief Compute relative error between two values
     */
    static Scalar relativeError(Scalar a, Scalar b) {
        const Scalar abs_a = std::abs(a);
        const Scalar abs_b = std::abs(b);
        const Scalar max_val = std::max(abs_a, abs_b);
        
        if (max_val < Types::epsilon()) {
            return std::abs(a - b);
        }
        
        return std::abs(a - b) / max_val;
    }
    
    /**
     * @brief Check if two scalars are approximately equal
     */
    static bool isApprox(Scalar a, Scalar b, Scalar tolerance = Types::tolerance()) {
        return relativeError(a, b) < tolerance;
    }
};

// Convenience aliases
using DifferentiationTestUtilsf = DifferentiationTestUtils<float>;
using DifferentiationTestUtilsd = DifferentiationTestUtils<double>;

} // namespace sofa::component::cosserat::liegroups::testing
