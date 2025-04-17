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
 *                           Plugin Cosserat    v1.0                           *
 *				                                                              *
 * This plugin is also distributed under the GNU LGPL (Lesser General          *
 * Public License) license with the same conditions than SOFA.                 *
 *                                                                             *
 * Contributors: Defrost team  (INRIA, University of Lille, CNRS,              *
 *               Ecole Centrale de Lille)                                      *
 *                                                                             *
 * Contact information: https://project.inria.fr/softrobot/contact/            *
 *                                                                             *
 ******************************************************************************/
#pragma once

#include <Eigen/Dense>
#include <limits>
#include <cmath>
#include <type_traits>

namespace sofa::component::cosserat::liegroups {

/**
 * @brief Types container for Lie groups
 *
 * This template class provides type definitions for various
 * data types used in Lie group implementations.
 *
 * @tparam _Scalar The scalar type used for computations (must be a floating-point type)
 */
template<typename _Scalar>
struct Types {
    static_assert(std::is_floating_point<_Scalar>::value,
                 "Scalar type must be a floating-point type");

    using Scalar = _Scalar;

    // Common Eigen types
    template<int Rows, int Cols = Rows>
    using Matrix = Eigen::Matrix<Scalar, Rows, Cols>;

    template<int Rows>
    using Vector = Eigen::Matrix<Scalar, Rows, 1>;

    using Vector2 = Vector<2>;
    using Vector3 = Vector<3>;
    using Vector4 = Vector<4>;
    using Vector6 = Vector<6>;
    using Vector7 = Vector<7>;

    using Matrix2 = Matrix<2, 2>;
    using Matrix3 = Matrix<3, 3>;
    using Matrix4 = Matrix<4, 4>;
    using Matrix6 = Matrix<6, 6>;
    using Matrix7 = Matrix<7, 7>;

    using Quaternion = Eigen::Quaternion<Scalar>;
    using AngleAxis = Eigen::AngleAxis<Scalar>;
    using Translation = Eigen::Translation<Scalar, 3>;

    // Common matrices for transformation groups
    using RotationMatrix = Matrix3;
    using HomogeneousMatrix = Matrix4;

    // Common types for tangent spaces
    using TangentVector2 = Vector<2>;
    using TangentVector3 = Vector<3>;
    using TangentVector6 = Vector<6>;

    // Common types for adjoint representations
    using AdjointMatrix2 = Matrix<2, 2>;
    using AdjointMatrix3 = Matrix<3, 3>;
    using AdjointMatrix6 = Matrix<6, 6>;

    // Dynamic-sized matrices
    using MatrixX = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    using VectorX = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;

    // Array types for element-wise operations
    template<int Rows, int Cols = Rows>
    using Array = Eigen::Array<Scalar, Rows, Cols>;

    /**
     * @brief Get the machine epsilon for the scalar type
     * @return The smallest representable positive number epsilon such that 1 + epsilon > 1
     */
    static constexpr Scalar epsilon() {
        return std::numeric_limits<Scalar>::epsilon();
    }

    /**
     * @brief Get the smallest positive value representable by the scalar type
     * @return The minimum positive value
     */
    static constexpr Scalar minPositive() {
        return std::numeric_limits<Scalar>::min();
    }

    /**
     * @brief Get the maximum value representable by the scalar type
     * @return The maximum value
     */
    static constexpr Scalar maxValue() {
        return std::numeric_limits<Scalar>::max();
    }

    /**
     * @brief Get pi constant value
     * @return Pi with the precision of the scalar type
     */
    static constexpr Scalar pi() {
        return Scalar(M_PI);
    }

    /**
     * @brief Get 2*pi constant value
     * @return 2*Pi with the precision of the scalar type
     */
    static constexpr Scalar twoPi() {
        return Scalar(2 * M_PI);
    }

    /**
     * @brief Get pi/2 constant value
     * @return Pi/2 with the precision of the scalar type
     */
    static constexpr Scalar halfPi() {
        return Scalar(M_PI_2);
    }

    /**
     * @brief Check if a value is approximately zero
     * @param value The value to check
     * @param eps The tolerance (optional)
     * @return true if the value is approximately zero
     */
    static bool isZero(const Scalar& value, const Scalar& eps = epsilon()) {
        return std::abs(value) <= eps;
    }

    /**
     * @brief Check if two values are approximately equal
     * @param a First value
     * @param b Second value
     * @param eps The tolerance (optional)
     * @return true if the values are approximately equal
     */
    static bool isApprox(const Scalar& a, const Scalar& b, const Scalar& eps = epsilon()) {
        if (a == b) return true;
        const Scalar absA = std::abs(a);
        const Scalar absB = std::abs(b);
        const Scalar diff = std::abs(a - b);
        if (a == 0 || b == 0 || (absA + absB < minPositive())) {
            return diff < eps * minPositive();
        } else {
            return diff < eps * std::max(absA, absB);
        }
    }

    /**
     * @brief Normalize an angle to the range [-pi, pi]
     * @param angle The angle to normalize (in radians)
     * @return The normalized angle
     */
    static Scalar normalizeAngle(Scalar angle) {
        angle = std::fmod(angle + pi(), twoPi());
        if (angle < 0)
            angle += twoPi();
        return angle - pi();
    }

    /**
     * @brief Compute a safe inverse for small denominators
     * @param value The value to invert
     * @param eps The minimum absolute value below which to apply special handling
     * @return The inverted value, or regularized value if near zero
     */
    static Scalar safeInverse(const Scalar& value, const Scalar& eps = epsilon()) {
        if (std::abs(value) < eps) {
            // Regularized inverse: sign(value) * 1/eps
            return (value >= 0 ? Scalar(1) : Scalar(-1)) / eps;
        }
        return Scalar(1) / value;
    }

    /**
     * @brief Compute sinc(x) = sin(x)/x with proper limiting behavior at x=0
     * @param x The input value
     * @return sin(x)/x for x≠0, or 1 for x=0
     */
    static Scalar sinc(const Scalar& x) {
        if (isZero(x)) {
            return Scalar(1);
        }
        return std::sin(x) / x;
    }

    /**
     * @brief Compute cosc(x) = (1-cos(x))/x² with proper limiting behavior at x=0
     * @param x The input value
     * @return (1-cos(x))/x² for x≠0, or 0.5 for x=0
     */
    static Scalar cosc(const Scalar& x) {
        if (isZero(x)) {
            return Scalar(0.5);
        }
        return (Scalar(1) - std::cos(x)) / (x * x);
    }

    /**
     * @brief Normalize a vector to unit length
     * @param v Vector to normalize
     * @param eps Threshold for zero-length check
     * @return Normalized vector, or zero vector if input length < eps
     */
    template<int Dim>
    static Vector<Dim> normalize(const Vector<Dim>& v, const Scalar& eps = epsilon()) {
        const Scalar norm = v.norm();
        if (norm < eps) {
            return Vector<Dim>::Zero();
        }
        return v / norm;
    }

    /**
     * @brief Linear interpolation between two scalars
     * @param a Starting value
     * @param b Ending value
     * @param t Interpolation parameter [0,1]
     * @return Interpolated value: a*(1-t) + b*t
     */
    static Scalar lerp(const Scalar& a, const Scalar& b, const Scalar& t) {
        return a + t * (b - a);
    }

    /**
     * @brief Linear interpolation between two vectors
     * @param a Starting vector
     * @param b Ending vector
     * @param t Interpolation parameter [0,1]
     * @return Interpolated vector: a*(1-t) + b*t
     */
    template<int Dim>
    static Vector<Dim> lerp(const Vector<Dim>& a, const Vector<Dim>& b, const Scalar& t) {
        return a + t * (b - a);
    }

    /**
     * @brief Clamp a value between minimum and maximum
     * @param val Value to clamp
     * @param min Minimum allowed value
     * @param max Maximum allowed value
     * @return Clamped value
     */
    static Scalar clamp(const Scalar& val, const Scalar& min, const Scalar& max) {
        return std::max(min, std::min(max, val));
    }

    /**
     * @brief Squared norm of a vector
     * @param v Input vector
     * @return Squared norm (avoids computing square root)
     */
    template<int Dim>
    static Scalar squaredNorm(const Vector<Dim>& v) {
        return v.squaredNorm();
    }

    /**
     * @brief Distance between two points
     * @param a First point
     * @param b Second point
     * @return Euclidean distance
     */
    template<int Dim>
    static Scalar distance(const Vector<Dim>& a, const Vector<Dim>& b) {
        return (a - b).norm();
    }

    /**
     * @brief Squared distance between two points (more efficient)
     * @param a First point
     * @param b Second point
     * @return Squared Euclidean distance
     */
    template<int Dim>
    static Scalar squaredDistance(const Vector<Dim>& a, const Vector<Dim>& b) {
        return (a - b).squaredNorm();
    }

    /**
     * @brief Create a skew-symmetric 3x3 matrix from a 3D vector (hat operator)
     * @param v 3D vector [x, y, z]
     * @return Skew-symmetric matrix [0, -z, y; z, 0, -x; -y, x, 0]
     */
    static Matrix3 skew(const Vector3& v) {
        Matrix3 result = Matrix3::Zero();
        result(0, 1) = -v(2);
        result(0, 2) = v(1);
        result(1, 0) = v(2);
        result(1, 2) = -v(0);
        result(2, 0) = -v(1);
        result(2, 1) = v(0);
        return result;
    }

    /**
     * @brief Extract vector from skew-symmetric matrix (vee operator)
     * @param S Skew-symmetric 3x3 matrix
     * @return 3D vector [S(2,1), S(0,2), S(1,0)]
     */
    static Vector3 unskew(const Matrix3& S) {
        Vector3 result;
        result(0) = S(2, 1);
        result(1) = S(0, 2);
        result(2) = S(1, 0);
        return result;
    }

    /**
     * @brief Create a rotation matrix from axis-angle representation
     * @param axis Rotation axis (will be normalized)
     * @param angle Rotation angle in radians
     * @return 3x3 rotation matrix
     */
    static Matrix3 rotationMatrix(const Vector3& axis, const Scalar& angle) {
        Vector3 n = normalize(axis);
        Scalar c = std::cos(angle);
        Scalar s = std::sin(angle);
        Scalar oneminusc = Scalar(1) - c;

        Matrix3 R;
        R(0, 0) = n[0] * n[0] * oneminusc + c;
        R(0, 1) = n[0] * n[1] * oneminusc - n[2] * s;
        R(0, 2) = n[0] * n[2] * oneminusc + n[1] * s;
        R(1, 0) = n[1] * n[0] * oneminusc + n[2] * s;
        R(1, 1) = n[1] * n[1] * oneminusc + c;
        R(1, 2) = n[1] * n[2] * oneminusc - n[0] * s;
        R(2, 0) = n[2] * n[0] * oneminusc - n[1] * s;
        R(2, 1) = n[2] * n[1] * oneminusc + n[0] * s;
        R(2, 2) = n[2] * n[2] * oneminusc + c;

        return R;
    }

    /**
     * @brief Fourth-order Runge-Kutta integration method
     *
     * @tparam F Function type for the derivative calculation
     * @tparam State Type of state vector
     * @param f Function that computes derivatives (should take state and return derivative)
     * @param y0 Initial state
     * @param dt Time step
     * @return Updated state after integration step
     */
    template <typename F, typename State>
    static State rungeKutta4(F f, const State& y0, const Scalar& dt) {
        State k1 = f(y0);
        State k2 = f(y0 + dt * k1 / Scalar(2));
        State k3 = f(y0 + dt * k2 / Scalar(2));
        State k4 = f(y0 + dt * k3);

        return y0 + dt * (k1 + Scalar(2) * k2 + Scalar(2) * k3 + k4) / Scalar(6);
    }
};

// Common type aliases for convenience
using Typesd = Types<double>;
using Typesf = Types<float>;

} // namespace sofa::component::cosserat::liegroups