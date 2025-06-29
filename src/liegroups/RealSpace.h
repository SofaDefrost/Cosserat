/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                  *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* iThis program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/\>.        *
******************************************************************************/

#pragma once

#include "LieGroupBase.h"
#include "LieGroupBase.inl"
#include "Types.h"
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
class RealSpace : public LieGroupBase<RealSpace<_Scalar, _Dim>, _Scalar, _Dim, _Dim, _Dim>
                 //,public LieGroupOperations<RealSpace<_Scalar, _Dim>> 
                 {
public:
    using Base = LieGroupBase<_Scalar, std::integral_constant<int, _Dim>, _Dim, _Dim>;
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
    // Implementations of LieGroupBase pure virtual methods
    RealSpace computeInverse() const {
        return RealSpace(-m_data);
    }

    static RealSpace computeExp(const TangentVector& algebra_element) {
        return RealSpace(algebra_element);
    }

    TangentVector computeLog() const {
        return m_data;
    }

    AdjointMatrix computeAdjoint() const {
        return AdjointMatrix::Identity();
    }

    Vector computeAction(const Vector& point) const {
        return point + m_data;
    }

    bool computeIsApprox(const RealSpace& other, 
                  const Scalar& eps = Types<Scalar>::epsilon()) const {
        return m_data.isApprox(other.m_data, eps);
    }

    static RealSpace computeIdentity() {
        return RealSpace(Vector::Zero());
    }

    static Matrix computeHat(const TangentVector& v) {
        Matrix result = Matrix::Zero();
        for (int i = 0; i < Dim; ++i) {
            result(i, i) = v(i);
        }
        return result;
    }

    static TangentVector computeVee(const Matrix& X) {
        TangentVector result;
        for (int i = 0; i < Dim; ++i) {
            result(i) = X(i, i);
        }
        return result;
    }

    static AdjointMatrix computeAd(const TangentVector& v) {
        return AdjointMatrix::Zero(); // Adjoint for R^n is zero matrix
    }

    template <typename Generator>
    static RealSpace computeRandom(Generator& gen) {
        return RealSpace(Types<Scalar>::template randomVector<Dim>(gen));
    }

    std::ostream& print(std::ostream& os) const {
        os << m_data.transpose();
        return os;
    }

    static constexpr std::string_view getTypeName() {
        return "RealSpace";
    }

    bool computeIsValid() const {
        return m_data.allFinite(); // Check if all elements are finite
    }

    void computeNormalize() {
        // No normalization needed for RealSpace
    }

    Scalar squaredDistance(const RealSpace& other) const {
        return (m_data - other.m_data).squaredNorm();
    }

private:
    Vector m_data;  ///< The underlying vector data
};

} // namespace sofa::component::cosserat::liegroups