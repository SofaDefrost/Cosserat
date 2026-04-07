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
#pragma once

#include <type_traits>

/**
 * @file AutodiffSupport.h
 * @brief Integration layer for automatic differentiation libraries
 * 
 * This header provides utilities to use Lie groups with automatic differentiation
 * libraries like autodiff (https://autodiff.github.io/).
 * 
 * ## Usage with autodiff library
 * 
 * The autodiff library provides two main types for automatic differentiation:
 * - `autodiff::dual` for forward mode (efficient for f: ℝⁿ → ℝ with small n)
 * - `autodiff::var` for reverse mode (efficient for f: ℝⁿ → ℝ with large n)
 * 
 * ### Forward Mode Example:
 * @code
 * #include <autodiff/forward/dual.hpp>
 * #include "SO3.h"
 * #include "AutodiffSupport.h"
 * 
 * using namespace autodiff;
 * using namespace sofa::component::cosserat::liegroups;
 * 
 * // Define SO3 with dual numbers
 * using SO3dual = SO3<dual>;
 * 
 * // Function to differentiate: angle of rotation
 * dual rotationAngle(const Eigen::Matrix<dual, 3, 1>& omega) {
 *     SO3dual R = SO3dual::exp(omega);
 *     return R.log().norm();  // Returns angle
 * }
 * 
 * // Compute derivative
 * Eigen::Vector3d omega_val(0.1, 0.2, 0.3);
 * Eigen::Matrix<dual, 3, 1> omega;
 * omega << dual(omega_val[0]), dual(omega_val[1]), dual(omega_val[2]);
 * 
 * dual angle;
 * Eigen::Vector3d dangle_domega;
 * 
 * // Compute all derivatives in one pass
 * for (int i = 0; i < 3; i++) {
 *     angle = derivative(rotationAngle, wrt(omega[i]), at(omega));
 *     dangle_domega[i] = angle;
 * }
 * @endcode
 * 
 * ### Reverse Mode Example:
 * @code
 * #include <autodiff/reverse/var.hpp>
 * #include "SE3.h"
 * #include "AutodiffSupport.h"
 * 
 * using namespace autodiff;
 * using namespace sofa::component::cosserat::liegroups;
 * 
 * // Define SE3 with var for reverse mode
 * using SE3var = SE3<var>;
 * 
 * // Cost function: distance to target position
 * var distanceToTarget(const Eigen::Matrix<var, 6, 1>& xi,
 *                      const Eigen::Matrix<double, 3, 1>& target) {
 *     SE3var T = SE3var::exp(xi);
 *     auto pos = T.translation();
 *     var dx = pos[0] - target[0];
 *     var dy = pos[1] - target[1];
 *     var dz = pos[2] - target[2];
 *     return dx*dx + dy*dy + dz*dz;
 * }
 * 
 * // Compute all derivatives at once (efficient!)
 * Eigen::Matrix<var, 6, 1> xi;
 * // ... initialize xi ...
 * 
 * Eigen::Vector3d target(1.0, 0.0, 0.0);
 * var cost = distanceToTarget(xi, target);
 * 
 * // Get all gradients in one reverse pass
 * Derivatives dcost = derivatives(cost);
 * Eigen::Matrix<double, 6, 1> gradient;
 * for (int i = 0; i < 6; i++) {
 *     gradient[i] = dcost(xi[i]);
 * }
 * @endcode
 * 
 * ## Compilation
 * 
 * To enable autodiff support, compile with:
 * @code{.sh}
 * cmake -DCOSSERAT_WITH_AUTODIFF=ON ..
 * @endcode
 * 
 * The CMake configuration will automatically detect the autodiff library if
 * it's located at `../autodiff` relative to the plugin root.
 */

#ifdef COSSERAT_WITH_AUTODIFF
#include <autodiff/forward/dual.hpp>
#include <autodiff/reverse/var.hpp>

namespace sofa::component::cosserat::liegroups {

	/**
	 * @brief Type trait to detect if a type is an autodiff type
	 */
	template<typename T>
	struct is_autodiff_type : std::false_type {};

	template<>
	struct is_autodiff_type<autodiff::dual> : std::true_type {};

	template<>
	struct is_autodiff_type<autodiff::var> : std::true_type {};

	template<typename T>
	inline constexpr bool is_autodiff_type_v = is_autodiff_type<T>::value;

	/**
	 * @brief Extract the underlying scalar type from potentially autodiff types
	 */
	template<typename T>
	struct scalar_type {
		using type = T;
	};

	template<>
	struct scalar_type<autodiff::dual> {
		using type = double;
	};

	template<>
	struct scalar_type<autodiff::var> {
		using type = double;
	};

	template<typename T>
	using scalar_type_t = typename scalar_type<T>::type;

	/**
	 * @brief Convert autodiff types to their value (strip derivatives)
	 */
	template<typename T>
	inline auto value(const T& x) -> std::enable_if_t<!is_autodiff_type_v<T>, T> {
		return x;
	}

	inline double value(const autodiff::dual& x) {
		return autodiff::val(x);
	}

	inline double value(const autodiff::var& x) {
		return autodiff::val(x);
	}

	/**
	 * @brief Helper to convert Eigen vectors/matrices between autodiff and regular types
	 */
	template<typename Derived, typename Scalar>
	auto toAutodiff(const Eigen::MatrixBase<Derived>& mat) 
		-> Eigen::Matrix<Scalar, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime> {
		
		Eigen::Matrix<Scalar, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime> result;
		result.resize(mat.rows(), mat.cols());
		
		for (int i = 0; i < mat.rows(); i++) {
			for (int j = 0; j < mat.cols(); j++) {
				result(i, j) = Scalar(mat(i, j));
			}
		}
		return result;
	}

	/**
	 * @brief Convert autodiff vectors/matrices to double
	 */
	template<typename Derived>
	auto toDouble(const Eigen::MatrixBase<Derived>& mat)
		-> Eigen::Matrix<double, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime> {
		
		using Scalar = typename Derived::Scalar;
		Eigen::Matrix<double, Derived::RowsAtCompileTime, Derived::ColsAtCompileTime> result;
		result.resize(mat.rows(), mat.cols());
		
		for (int i = 0; i < mat.rows(); i++) {
			for (int j = 0; j < mat.cols(); j++) {
				result(i, j) = value(mat(i, j));
			}
		}
		return result;
	}

} // namespace sofa::component::cosserat::liegroups

#else
// Autodiff not enabled - provide empty stubs

namespace sofa::component::cosserat::liegroups {

	template<typename T>
	struct is_autodiff_type : std::false_type {};

	template<typename T>
	inline constexpr bool is_autodiff_type_v = false;

	template<typename T>
	struct scalar_type {
		using type = T;
	};

	template<typename T>
	using scalar_type_t = typename scalar_type<T>::type;

	template<typename T>
	inline T value(const T& x) {
		return x;
	}

} // namespace sofa::component::cosserat::liegroups

#endif // COSSERAT_WITH_AUTODIFF
