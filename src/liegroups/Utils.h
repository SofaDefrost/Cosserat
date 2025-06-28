#ifndef COSSERAT_LIEGROUPS_UTILS_H
#define COSSERAT_LIEGROUPS_UTILS_H

#include <Eigen/Core>
#include <cmath>
#include <limits>
#include <type_traits>

namespace sofa::component::cosserat::liegroups {

/**
 * Utility functions for Lie groups.
 */
template <typename Scalar>
struct LieGroupUtils {
    using Matrix2 = Eigen::Matrix<Scalar, 2, 2>;
    using Vector2 = Eigen::Matrix<Scalar, 2, 1>;
    using Vector3 = Eigen::Matrix<Scalar, 3, 1>;

    static constexpr Scalar epsilon = std::numeric_limits<Scalar>::epsilon() * 100;
    static constexpr Scalar pi = M_PI;
    static constexpr Scalar two_pi = 2 * M_PI;

    /**
     * Normalize angle to [-π, π]
     */
    static Scalar normalizeAngle(Scalar angle) {
        // Normalize angle to [-π, π]
        angle = std::fmod(angle, two_pi);
        if (angle > pi) {
            angle -= two_pi;
        } else if (angle < -pi) {
            angle += two_pi;
        }
        return angle;
    }

    /**
     * Compute sinc(x) = sin(x)/x with numerical stability for small x
     */
    static Scalar sinc(Scalar x) {
        if (std::abs(x) < epsilon) {
            // Taylor series approximation for small angles
            // sinc(x) ≈ 1 - x²/6 + x⁴/120 - ...
            return Scalar(1) - x * x / Scalar(6);
        }
        return std::sin(x) / x;
    }

    /**
     * Compute the difference between two angles with proper wrapping
     */
    static Scalar angleDifference(Scalar a, Scalar b) {
        return normalizeAngle(a - b);
    }

    /**
     * Linear interpolation between two scalars
     */
    static Scalar lerp(Scalar a, Scalar b, Scalar t) {
        return a + t * (b - a);
    }

    /**
     * Spherical linear interpolation (SLERP) between two angles
     */
    static Scalar slerpAngle(Scalar a, Scalar b, Scalar t) {
        // Ensure shortest path
        Scalar diff = angleDifference(b, a);
        return normalizeAngle(a + t * diff);
    }

    /**
     * Bi-invariant distance between two angles (as SO(2) elements)
     */
    static Scalar angleDistance(Scalar a, Scalar b) {
        return std::abs(angleDifference(a, b));
    }

    /**
     * Check if an angle is near zero (within epsilon)
     */
    static bool isAngleNearZero(Scalar angle) {
        return std::abs(angle) < epsilon;
    }

    /**
     * Check if two angles are nearly equal (within epsilon, considering wrapping)
     */
    static bool areAnglesNearlyEqual(Scalar a, Scalar b) {
        return angleDistance(a, b) < epsilon;
    }

    /**
     * Numerically stable computation of 1 - cos(x) for small x
     */
    static Scalar oneMinusCos(Scalar x) {
        if (std::abs(x) < epsilon) {
            // Taylor series approximation for small angles
            // 1 - cos(x) ≈ x²/2 - x⁴/24 + ...
            return (x * x) / Scalar(2);
        }
        return Scalar(1) - std::cos(x);
    }

    /**
     * Safe vector normalization that handles near-zero vectors
     */
    template <typename Derived>
    static Eigen::Matrix<typename Derived::Scalar, Derived::RowsAtCompileTime, 1>
    safeNormalize(const Eigen::MatrixBase<Derived>& v) {
        using VectorType = Eigen::Matrix<typename Derived::Scalar, Derived::RowsAtCompileTime, 1>;
        typename Derived::Scalar norm = v.norm();
        if (norm < epsilon) {
            // Return zero vector or throw exception based on your preference
            return VectorType::Zero(v.rows());
        }
        return v / norm;
    }

    /**
     * Project a vector onto another vector
     */
    template <typename Derived1, typename Derived2>
    static Eigen::Matrix<typename Derived1::Scalar, Derived1::RowsAtCompileTime, 1>
    projectVector(const Eigen::MatrixBase<Derived1>& v, const Eigen::MatrixBase<Derived2>& onto) {
        using VectorType = Eigen::Matrix<typename Derived1::Scalar, Derived1::RowsAtCompileTime, 1>;
        typename Derived2::Scalar norm_squared = onto.squaredNorm();
        if (norm_squared < epsilon) {
            return VectorType::Zero(v.rows());
        }
        return onto * (v.dot(onto) / norm_squared);
    }

    /**
     * Path interpolation between two SE(2) elements represented as [angle, x, y]
     */
    static Vector3 interpolateSE2Path(const Vector3& start, const Vector3& end, Scalar t) {
        Vector3 result;
        // Interpolate angle using SLERP
        result[0] = slerpAngle(start[0], end[0], t);
        // Interpolate translation using LERP
        result[1] = lerp(start[1], end[1], t);
        result[2] = lerp(start[2], end[2], t);
        return result;
    }
};

} // namespace sofa::component::cosserat::liegroups

#endif // COSSERAT_LIEGROUPS_UTILS_H

