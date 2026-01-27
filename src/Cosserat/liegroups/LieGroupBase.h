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

#ifndef SOFA_COMPONENT_COSSERAT_LIEGROUPS_LIEGROUPBASE_H
#define SOFA_COMPONENT_COSSERAT_LIEGROUPS_LIEGROUPBASE_H

#include "Types.h"
#include <type_traits>

namespace sofa::component::cosserat::liegroups {

/**
 * @brief Base class template for all Lie group implementations
 * 
 * This class defines the interface that all Lie group implementations must satisfy.
 * It provides pure virtual methods for the fundamental operations of Lie groups:
 * composition, inversion, exponential map, logarithmic map, adjoint representation,
 * and group action on points.
 * 
 * @tparam _Scalar The scalar type used for computations (must be a floating-point type)
 * @tparam _Dim The dimension of the group representation
 */
template<typename _Scalar, int _Dim>
class LieGroupBase {
public:
    static_assert(std::is_floating_point<_Scalar>::value,
                 "Scalar type must be a floating-point type");
    static_assert(_Dim > 0, "Dimension must be positive");

    using Scalar = _Scalar;
    using Types = Types<Scalar>;
    static constexpr int Dim = _Dim;

    // Define commonly used types
    using Vector = typename Types::template Vector<Dim>;
    using Matrix = typename Types::template Matrix<Dim, Dim>;
    using TangentVector = typename Types::template Vector<Dim>;
    using AdjointMatrix = typename Types::template Matrix<Dim, Dim>;

    virtual ~LieGroupBase() = default;

    /**
     * @brief Group composition operation
     * @param other Another element of the same Lie group
     * @return The composition this * other
     */
    virtual LieGroupBase& operator*(const LieGroupBase& other) const = 0;

    /**
     * @brief Compute the inverse element
     * @return The inverse element such that this * inverse() = identity()
     */
    virtual LieGroupBase inverse() const = 0;

    /**
     * @brief Exponential map from Lie algebra to Lie group
     * @param algebra_element Element of the Lie algebra (tangent space at identity)
     * @return The corresponding element in the Lie group
     */
    virtual LieGroupBase exp(const TangentVector& algebra_element) const = 0;

    /**
     * @brief Logarithmic map from Lie group to Lie algebra
     * @return The corresponding element in the Lie algebra
     */
    virtual TangentVector log() const = 0;

    /**
     * @brief Adjoint representation of the group element
     * @return The adjoint matrix representing the action on the Lie algebra
     */
    virtual AdjointMatrix adjoint() const = 0;

    /**
     * @brief Group action on a point
     * @param point The point to transform
     * @return The transformed point
     */
    virtual Vector act(const Vector& point) const = 0;

    /**
     * @brief Check if this element is approximately equal to another
     * @param other Another element of the same Lie group
     * @param eps Tolerance for comparison (optional)
     * @return true if elements are approximately equal
     */
    virtual bool isApprox(const LieGroupBase& other, 
                         const Scalar& eps = Types::epsilon()) const = 0;

    /**
     * @brief Get the identity element of the group
     * @return Reference to the identity element
     */
    virtual const LieGroupBase& identity() const = 0;

    /**
     * @brief Get the dimension of the Lie algebra (tangent space)
     * @return The dimension of the Lie algebra
     */
    virtual int algebraDimension() const = 0;

    /**
     * @brief Get the dimension of the space the group acts on
     * @return The dimension of the space the group acts on
     */
    virtual int actionDimension() const = 0;
};

} // namespace sofa::component::cosserat::liegroups

#endif // SOFA_COMPONENT_COSSERAT_LIEGROUPS_LIEGROUPBASE_H
