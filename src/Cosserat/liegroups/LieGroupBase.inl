#ifndef SOFA_COMPONENT_COSSERAT_LIEGROUPS_LIEGROUPBASE_INL
#define SOFA_COMPONENT_COSSERAT_LIEGROUPS_LIEGROUPBASE_INL

#include "LieGroupBase.h"
#include "Types.h"

namespace sofa::component::cosserat::liegroups {

/**
 * @brief Default implementations for common operations
 *
 * This file provides default implementations of common Lie group operations
 * in terms of the core operations that derived classes must define.
 */


// Default implementation of distance using the logarithm
template<typename Derived, typename _Scalar, int _Dim, int _AlgebraDim, int _ActionDim>
inline typename LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::Scalar
LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::distance(const Derived& other) const noexcept {
    try {
        // Default implementation: norm of the difference in the Lie algebra
        // This uses g^(-1)*h which maps to the tangent space at identity
        Derived rel = derived().inverse().compose(other);
        return rel.log().norm();
    } catch (const std::exception&) {
        // Fallback if inverse or log fails
        if (derived().computeIsApprox(other)) {
            return Scalar(0);
        }
        return std::numeric_limits<Scalar>::max();
    }
}

// Default implementation of interpolation using exp/log
template<typename Derived, typename _Scalar, int _Dim, int _AlgebraDim, int _ActionDim>
inline Derived
LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::interpolate(
    const Derived& other, const Scalar& t) const noexcept {

    if (t <= Scalar(0)) return derived();
    if (t >= Scalar(1)) return other;

    try {
        // Default implementation using the Lie group geodesic
        Derived rel = derived().inverse().compose(other);
        TangentVector delta = rel.log();
        return derived().compose(Derived::exp(t * delta));
    } catch (const std::exception&) {
        // Fallback simple linear interpolation (not generally valid for Lie groups)
        return derived();
    }
}

// @Todo: the BCH implementation could benefit from a note about convergence
// radius and a warning when the approximation might be less accurate for large
// inputs.
// Default implementation of Baker-Campbell-Hausdorff formula
template<typename Derived, typename _Scalar, int _Dim, int _AlgebraDim, int _ActionDim>
inline typename LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::TangentVector
LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::BCH(
    const TangentVector& v, const TangentVector& w, int order) {

    using Adjoint = typename Derived::AdjointMatrix;

    // First-order approximation is just v + w
    TangentVector result = v + w;

    // Include higher order terms if requested
    if (order >= 2) {
        // Add second-order term: 1/2 * [v, w]
        AlgebraMatrix vMat = Derived::computeHat(v);
        AlgebraMatrix wMat = Derived::computeHat(w);
        AlgebraMatrix commutator = vMat * wMat - wMat * vMat;
        result += Derived::computeVee(commutator) * Scalar(0.5);

        if (order >= 3) {
            // Add third-order term: 1/12 * [v, [v, w]] + 1/12 * [w, [w, v]]
            AlgebraMatrix vvw = vMat * commutator - commutator * vMat;
            AlgebraMatrix wwv = wMat * (commutator * Scalar(-1)) - (commutator * Scalar(-1)) * wMat;

            result += Derived::computeVee(vvw) * Scalar(1.0/12.0);
            result += Derived::computeVee(wwv) * Scalar(1.0/12.0);
        }
    }

    return result;
}

// Default implementation for differential of exponential
template<typename Derived, typename _Scalar, int _Dim, int _AlgebraDim, int _ActionDim>
inline typename LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::AdjointMatrix
LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::dexp(const TangentVector& v) {
    using Matrix = typename Derived::AdjointMatrix;
    using Scalar = typename Derived::Scalar;
    using Types = typename Derived::Types;

    // Compute using series approximation
    const Scalar theta = v.norm();

    if (Types::isZero(theta)) {
        return Matrix::Identity();
    }

    // Create the adjoint representation of v
    Matrix ad_v = Derived::ad(v);

    // Series approximation of dexp
    Matrix result = Matrix::Identity();
    result += ad_v * Scalar(0.5);

    const Scalar theta_sq = theta * theta;
    Scalar factor = Scalar(1);

    // Add higher-order terms (Bernoulli numbers)
    if (!Types::isZero(theta_sq)) {
        const Scalar inv_theta_sq = Scalar(1) / theta_sq;

        // (1-cos(θ))/θ²
        factor = Types::cosc(theta);
        result += ad_v * ad_v * factor * Scalar(1.0/6.0);

        // (θ-sin(θ))/θ³
        factor = (theta - std::sin(theta)) * inv_theta_sq / theta;
        result += ad_v * ad_v * ad_v * factor * Scalar(1.0/24.0);
    }

    return result;
}

// Default implementation for differential of logarithm
template<typename Derived, typename _Scalar, int _Dim, int _AlgebraDim, int _ActionDim>
inline typename LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::AdjointMatrix
LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::dlog() const {
    using Matrix = typename Derived::AdjointMatrix;
    using Vector = typename Derived::TangentVector;
    using Scalar = typename Derived::Scalar;
    using Types = typename Derived::Types;

    // Get the logarithm of the current element
    Vector v = derived().log();
    const Scalar theta = v.norm();

    if (Types::isZero(theta)) {
        return Matrix::Identity();
    }

    // Create the adjoint representation of v
    Matrix ad_v = Derived::ad(v);

    // Series approximation of dlog
    Matrix result = Matrix::Identity();

    // Add adjustment factors
    const Scalar half_theta = theta * Scalar(0.5);
    Scalar factor = Scalar(0.5) * half_theta * std::cos(half_theta) / std::sin(half_theta);
    result -= ad_v * Scalar(0.5);

    if (!Types::isZero(theta * theta)) {
        // Add higher-order terms
        result += ad_v * ad_v * (Scalar(1.0/12.0) - factor * Scalar(1.0/6.0));
    }

    return result;
}

// Default implementation of action Jacobian
template<typename Derived, typename _Scalar, int _Dim, int _AlgebraDim, int _ActionDim>
inline typename LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::JacobianMatrix
LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::actionJacobian(
    const ActionVector& point) const noexcept {

    using Matrix = typename Derived::JacobianMatrix;
    using AlgVec = typename Derived::TangentVector;
    using Scalar = typename Derived::Scalar;
    using Types = typename Derived::Types;

    const Scalar eps = std::sqrt(Types::epsilon());
    Matrix J = Matrix::Zero();

    // Numerical differentiation approach
    for (int i = 0; i < AlgebraDim; ++i) {
        AlgVec delta = AlgVec::Zero();
        delta[i] = eps;

        // Forward difference
        Derived perturbed = derived().compose(Derived::exp(delta));
        ActionVector point_plus = perturbed.act(point);

        // Compute column of Jacobian
        J.col(i) = (point_plus - derived().act(point)) / eps;
    }

    return J;
}

// Default hat operator (to be specialized for different groups)
template<typename Derived, typename _Scalar, int _Dim, int _AlgebraDim, int _ActionDim>
inline typename LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::AlgebraMatrix
LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::hat(const TangentVector& v) noexcept {
    return Derived::computeHat(v);
}

// Default vee operator (to be specialized for different groups)
template<typename Derived, typename _Scalar, int _Dim, int _AlgebraDim, int _ActionDim>
inline typename LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::TangentVector
LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::vee(const AlgebraMatrix& X) noexcept {
    return Derived::computeVee(X);
}

// Helper for computing the adjoint representation of a Lie algebra element
template<typename Derived, typename _Scalar, int _Dim, int _AlgebraDim, int _ActionDim>
inline typename LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::AdjointMatrix
LieGroupBase<Derived, _Scalar, _Dim, _AlgebraDim, _ActionDim>::ad(const TangentVector& v) {
    using Matrix = typename Derived::AdjointMatrix;
    using AlgMat = typename Derived::AlgebraMatrix;

    // This is a default implementation for matrix Lie groups
    // Derived classes should specialize this if needed
    AlgMat X = Derived::computeHat(v);

    // For matrix Lie algebras, ad_X(Y) = [X, Y] = XY - YX
    // This builds the matrix representation of this operator
    Matrix result = Matrix::Zero();

    for (int j = 0; j < AlgebraDim; ++j) {
        // Create basis vector and convert to algebra matrix
        TangentVector e_j = TangentVector::Zero();
        e_j[j] = Scalar(1);
        AlgMat E_j = Derived::computeHat(e_j);

        // Compute [X, E_j]
        AlgMat commutator = X * E_j - E_j * X;

        // Extract vector form and set as column
        result.col(j) = Derived::computeVee(commutator);
    }

    return result;
}

} // End namespace sofa::component::cosserat::liegroups

#endif // SOFA_COMPONENT_COSSERAT_LIEGROUPS_LIEGROUPBASE_INL
