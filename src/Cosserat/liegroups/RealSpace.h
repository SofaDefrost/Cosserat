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

#ifndef SOFA_COMPONENT_COSSERAT_LIEGROUPS_REALSPACE_H
#define SOFA_COMPONENT_COSSERAT_LIEGROUPS_REALSPACE_H

#include "LieGroupBase.h"
#include "LieGroupBase.inl"

namespace sofa::component::cosserat::liegroups {

/**
 * @brief Implementation of the Real space ℝ(n) as a Lie group
 * 
 * This class implements the vector space ℝ(n) as a Lie group where:
 * - The group operation is vector addition
 * - The inverse operation is vector negation
 * - The Lie algebra is the same space with the same operations
 * - The exponential and logarithm maps are identity functions
 * - The adjoint representation is the identity matrix
 * 
 * @tparam _Scalar The scalar type (must be a floating-point type)
 * @tparam _Dim The dimension of the space
 */
template<typename _Scalar, int _Dim>
class RealSpace : public LieGroupBase<_Scalar, _Dim>,
                 public LieGroupOperations<RealSpace<_Scalar, _Dim>> {
public:
    using Base = LieGroupBase<_Scalar, _Dim>;
    using Scalar = typename Base::Scalar;
    using Vector = typename Base::Vector;
    using Matrix = typename Base::Matrix;
    using TangentVector = typename Base::TangentVector;
    using AdjointMatrix = typename Base::AdjointMatrix;
    
    static constexpr int Dim = Base::Dim;

    /**
     * @brief Default constructor initializes to zero (identity element)
     */
    RealSpace() : m_data(Vector::Zero()) {}

    /**
     * @brief Construct from a vector
     */
    explicit RealSpace(const Vector& v) : m_data(v) {}

    /**
     * @brief Group composition (vector addition)
     */
    RealSpace operator*(const RealSpace& other) const {
        return RealSpace(m_data + other.m_data);
    }

    /**
     * @brief Inverse element (negation)
     */
    RealSpace inverse() const override {
        return RealSpace(-m_data);
    }

    /**
     * @brief Exponential map (identity map for ℝ(n))
     */
    RealSpace exp(const TangentVector& algebra_element) const override {
        return RealSpace(algebra_element);
    }

    /**
     * @brief Logarithmic map (identity map for ℝ(n))
     */
    TangentVector log() const override {
        return m_data;
    }

    /**
     * @brief Adjoint representation (identity matrix for ℝ(n))
     */
    AdjointMatrix adjoint() const override {
        return AdjointMatrix::Identity();
    }

    /**
     * @brief Group action on a point (translation)
     */
    Vector act(const Vector& point) const override {
        return point + m_data;
    }

    /**
     * @brief Check if approximately equal to another element
     */
    bool isApprox(const RealSpace& other, 
                  const Scalar& eps = Types<Scalar>::epsilon()) const {
        return m_data.isApprox(other.m_data, eps);
    }

    /**
     * @brief Get the identity element (zero vector)
     */
    static const RealSpace& identity() {
        static const RealSpace id;
        return id;
    }

    /**
     * @brief Get the dimension of the Lie algebra
     */
    int algebraDimension() const override { return Dim; }

    /**
     * @brief Get the dimension of the space the group acts on
     */
    int actionDimension() const override { return Dim; }

    /**
     * @brief Access the underlying data
     */
    const Vector& data() const { return m_data; }
    Vector& data() { return m_data; }

private:
    Vector m_data;  ///< The underlying vector data
};

} // namespace sofa::component::cosserat::liegroups

#include "RealSpace.inl"

#endif // SOFA_COMPONENT_COSSERAT_LIEGROUPS_REALSPACE_H
