/******************************************************************************
 *                 SOFA, Simulation Open-Framework Architecture * (c) 2006
 *INRIA, USTL, UJF, CNRS, MGH                     *
 *                                                                             *
 * This program is free software; you can redistribute it and/or modify it     *
 * under the terms of the GNU Lesser General Public License as published by    *
 * the Free Software Foundation; either version 2.1 of the License, or (at     *
 * your option) any later version.                                             *
 *                                                                             *
 * iThis program is distributed in the hope that it will be useful, but WITHOUT
 ** ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or *
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
 * for more details.                                                           *
 *                                                                             *
 * You should have received a copy of the GNU Lesser General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/\>. *
 ******************************************************************************/

#pragma once

#include <iostream>
#include "LieGroup.h"
#include "Types.h"

namespace sofa::component::cosserat::liegroups {

/**
 * @brief Implementation of the Real space ℝ(n) as a Lie group
 *  * This class implements the vector space ℝ(n) as a Lie group where:
 * - The group operation is vector addition
 * - The inverse operation is vector negation
 * - The Lie algebra is the same space with the same operations
 * - The exponential and logarithm maps are identity functions
 * - The adjoint representation is the identity matrix
 *  * @tparam Scalar_ The scalar type (must be a floating-point type)
 * @tparam Dim_ The dimension of the space
 */
template<typename Scalar_, int Dim_>
class Rn : public LieGroup<GroupProperties<Rn<Scalar_, Dim_>, Scalar_, Dim_, Dim_, Dim_>>
{
   public:
    using Base = LieGroup<GroupProperties<Rn<Scalar_, Dim_>, Scalar_, Dim_, Dim_, Dim_>>;
    using Scalar = typename Base::Scalar;
    using Vector = typename Base::Vector;
    using Matrix = typename Base::Matrix;
    using TangentVector = typename Base::TangentVector;
    using AdjointMatrix = typename Base::AdjointMatrix;

    static constexpr int Dim = Base::Dim;

    /**
     * @brief Default constructor initializes to zero (identity element)
     */
    Rn() : m_data(Vector::Zero()) {}

    /**
     * @brief Construct from a vector
     * @param v The vector to construct from.
     */
    explicit Rn(const Vector &v) : m_data(v) {}

    /**
     * @brief Computes the inverse of the current Rn element.
     * For Rn, the inverse is simply the negation of the vector.
     * @return The inverse Rn element.
     */
    Rn computeInverse() const { return Rn(-m_data); }

    /**
     * @brief Composes two Rn elements (group operation).
     * For Rn, composition is vector addition.
     * @param other The other Rn element to compose with.
     * @return The composed Rn element.
     */
    Rn computeCompose(const Rn &other) const { return Rn(m_data + other.m_data); }

    /**
     * @brief Computes the exponential map from the Lie algebra to the Rn
     * group. For Rn, the exponential map is the identity function.
     * @param algebra_element The element from the Lie algebra.
     * @return The corresponding Rn element.
     */
    static Rn computeExp(const TangentVector &algebra_element) { return Rn(algebra_element); }

    /**
     * @brief Computes the logarithmic map from the Rn group to its Lie
     * algebra. For Rn, the logarithmic map is the identity function.
     * @return The corresponding element in the Lie algebra.
     */
    TangentVector computeLog() const { return m_data; }

    /**
     * @brief Computes the adjoint representation of the current Rn
     * element. For Rn, the adjoint matrix is the identity matrix.
     * @return The adjoint matrix.
     */
    AdjointMatrix computeAdjoint() const { return AdjointMatrix::Identity(); }

    /**
     * @brief Applies the group action of the current Rn element on a
     * point. For Rn, the action is vector addition.
     * @param point The point to act upon.
     * @return The transformed point.
     */
    Vector computeAction(const Vector &point) const { return point + m_data; }

    /**
     * @brief Checks if the current Rn element is approximately equal to
     * another.
     * @param other The other Rn element to compare with.
     * @param eps The tolerance for approximation. Defaults to
     * `Types<Scalar>::epsilon()`.
     * @return True if the elements are approximately equal, false otherwise.
     */
    bool computeIsApprox(const Rn &other, const Scalar &eps = Types<Scalar>::epsilon()) const {
        return m_data.isApprox(other.m_data, eps);
    }

    /**
     * @brief Computes the identity element of the Rn group.
     * For Rn, the identity element is the zero vector.
     * @return The identity Rn element.
     */
    static Rn computeIdentity() { return Rn(Vector::Zero()); }

    /**
     * @brief Computes the hat operator for Rn.
     * This maps a tangent vector to its matrix representation in the Lie algebra.
     * For Rn, this is a diagonal matrix with the vector elements on the
     * diagonal.
     * @param v The tangent vector.
     * @return The matrix representation in the Lie algebra.
     */
    static Matrix computeHat(const TangentVector &v) {
        Matrix result = Matrix::Zero();
        for (int i = 0; i < Dim; ++i) {
            result(i, i) = v(i);
        }
        return result;
    }

    /**
     * @brief Computes the vee operator for Rn.
     * This maps a matrix representation in the Lie algebra back to its tangent
     * vector form. This is the inverse of the hat operator.
     * @param X The matrix representation in the Lie algebra.
     * @return The tangent vector.
     */
    static TangentVector computeVee(const Matrix &X) {
        TangentVector result;
        for (int i = 0; i < Dim; ++i) {
            result(i) = X(i, i);
        }
        return result;
    }

    /**
     * @brief Computes the adjoint representation of a Lie algebra element for
     * Rn. For Rn, the adjoint matrix is the zero matrix.
     * @param v The element of the Lie algebra in vector form.
     * @return The adjoint matrix.
     */
    static AdjointMatrix computeAd([[maybe_unused]] const TangentVector &v) {
        return AdjointMatrix::Zero(); // Adjoint for R^n is zero matrix
    }

    /**
     * @brief Generates a random Rn element.
     * @tparam Generator The type of the random number generator.
     * @param gen The random number generator.
     * @return A random Rn element.
     */
    template<typename Generator>
    static Rn computeRandom(Generator &gen) {
        return Rn(Types<Scalar>::template randomVector<Dim>(gen));
    }

    /**
     * @brief Prints the Rn element to an output stream.
     * @param os The output stream.
     * @return The output stream.
     */
    std::ostream &print(std::ostream &os) const {
        os << m_data.transpose();
        return os;
    }

    /**
     * @brief Gets the type name of the Rn class.
     * @return A string view of the type name.
     */
    static constexpr std::string_view getTypeName() { return "Rn"; }

    /**
     * @brief Checks if the current Rn element is valid.
     * @return True if all elements of the underlying vector are finite, false
     * otherwise.
     */
    bool computeIsValid() const {
        return m_data.allFinite(); // Check if all elements are finite
    }

    /**
     * @brief Normalizes the Rn element.
     * For Rn, no normalization is needed.
     */
    void computeNormalize() {
        // No normalization needed for Rn
    }

    /**
     * @brief Computes the squared distance between the current Rn element
     * and another.
     * @param other The other Rn element.
     * @return The squared Euclidean distance between the underlying vectors.
     */
    Scalar squaredDistance(const Rn &other) const { return (m_data - other.m_data).squaredNorm(); }

   private:
    Vector m_data; ///< The underlying vector data
};

} // namespace sofa::component::cosserat::liegroups
