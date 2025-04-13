/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                  *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
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
* along with this program. If not, see <http://www.gnu.org/licenses/\>.        *
******************************************************************************/

#ifndef SOFA_COMPONENT_COSSERAT_LIEGROUPS_SO2_H
#define SOFA_COMPONENT_COSSERAT_LIEGROUPS_SO2_H

#include "LieGroupBase.h"
#include "LieGroupBase.inl"
#include <Eigen/Geometry>

namespace sofa::component::cosserat::liegroups {

/**
 * @brief Implementation of SO(2), the Special Orthogonal group in 2D
 * 
 * This class implements the group of rotations in 2D space. Elements of SO(2)
 * are represented internally using complex numbers (cos θ + i sin θ), which
 * provides an efficient way to compose rotations and compute the exponential map.
 * 
 * The Lie algebra so(2) consists of skew-symmetric 2×2 matrices, which can be
 * identified with real numbers (the rotation angle).
 * 
 * @tparam _Scalar The scalar type (must be a floating-point type)
 */
template<typename _Scalar>
class SO2 : public LieGroupBase<_Scalar, 2>,
           public LieGroupOperations<SO2<_Scalar>> {
public:
    using Base = LieGroupBase<_Scalar, 2>;
    using Scalar = typename Base::Scalar;
    using Vector = typename Base::Vector;
    using Matrix = typename Base::Matrix;
    using TangentVector = typename Base::TangentVector;
    using AdjointMatrix = typename Base::AdjointMatrix;
    
    static constexpr int Dim = Base::Dim;
    using Complex = Eigen::Matrix<Scalar, 2, 1>;  // Represents complex number as 2D vector

    /**
     * @brief Default constructor creates identity rotation (angle = 0)
     */
    SO2() : m_angle(0) {
        updateComplex();
    }

    /**
     * @brief Construct from angle (in radians)
     */
    explicit SO2(const Scalar& angle) : m_angle(angle) {
        updateComplex();
    }

    /**
     * @brief Group composition (rotation composition)
     */
    SO2 operator*(const SO2& other) const {
        // Complex multiplication
        Complex result;
        result(0) = m_complex(0) * other.m_complex(0) - m_complex(1) * other.m_complex(1);
        result(1) = m_complex(0) * other.m_complex(1) + m_complex(1) * other.m_complex(0);
        return SO2(std::atan2(result(1), result(0)));
    }

    /**
     * @brief Inverse element (opposite rotation)
     */
    SO2 inverse() const override {
        return SO2(-m_angle);
    }

    /**
     * @brief Exponential map (angle to rotation)
     * For SO(2), this is just the angle itself as rotation
     */
    SO2 exp(const TangentVector& algebra_element) const override {
        return SO2(algebra_element[0]);
    }

    /**
     * @brief Logarithmic map (rotation to angle)
     */
    TangentVector log() const override {
        return TangentVector::Constant(m_angle);
    }

    /**
     * @brief Adjoint representation
     * For SO(2), this is simply the identity matrix as the group is abelian
     */
    AdjointMatrix adjoint() const override {
        return AdjointMatrix::Identity();
    }

    /**
     * @brief Group action on a point (rotate the point)
     */
    Vector act(const Vector& point) const override {
        Vector result;
        result(0) = m_complex(0) * point(0) - m_complex(1) * point(1);
        result(1) = m_complex(1) * point(0) + m_complex(0) * point(1);
        return result;
    }

    /**
     * @brief Check if approximately equal to another rotation
     */
    bool isApprox(const SO2& other, 
                  const Scalar& eps = Types<Scalar>::epsilon()) const {
        return std::abs(normalizeAngle(m_angle - other.m_angle)) <= eps;
    }

    /**
     * @brief Get the identity element (zero rotation)
     */
    static const SO2& identity() {
        static const SO2 id;
        return id;
    }

    /**
     * @brief Get the dimension of the Lie algebra (1 for SO(2))
     */
    int algebraDimension() const override { return 1; }

    /**
     * @brief Get the dimension of the space the group acts on (2 for SO(2))
     */
    int actionDimension() const override { return 2; }

    /**
     * @brief Get the rotation angle in radians
     */
    Scalar angle() const { return m_angle; }

    /**
     * @brief Get the rotation matrix representation
     */
    Matrix matrix() const {
        Matrix R;
        R << m_complex(0), -m_complex(1),
             m_complex(1),  m_complex(0);
        return R;
    }

private:
    Scalar m_angle;    ///< The rotation angle in radians
    Complex m_complex; ///< Complex number representation (cos θ, sin θ)

    /**
     * @brief Update complex number representation from angle
     */
    void updateComplex() {
        m_complex << std::cos(m_angle), std::sin(m_angle);
    }

    /**
     * @brief Normalize angle to [-π, π]
     */
    static Scalar normalizeAngle(Scalar angle) {
        const Scalar two_pi = 2 * Types<Scalar>::pi();
        angle = std::fmod(angle + Types<Scalar>::pi(), two_pi);
        if (angle < 0) angle += two_pi;
        return angle - Types<Scalar>::pi();
    }
};

} // namespace sofa::component::cosserat::liegroups

#include "SO2.inl"

#endif // SOFA_COMPONENT_COSSERAT_LIEGROUPS_SO2_H
