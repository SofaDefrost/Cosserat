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


#include <tuple>
#include <type_traits>
#include <array>
#include <iostream>
#include <random>
#include <Cosserat/liegroups/Types.h>
#include <Cosserat/liegroups/LieGroupBase.h>
#include <Cosserat/liegroups/LieGroupBase.inl>
#include <Cosserat/liegroups/SO2.h>
#include <Cosserat/liegroups/RealSpace.h>
#include <Cosserat/liegroups/SE2.h>
#include <Cosserat/liegroups/SE3.h>


namespace sofa::component::cosserat::liegroups {

namespace detail {

// Helper to compute total dimension of all groups
template<typename... Groups>
struct TotalDimension;

template<>
struct TotalDimension<> {
    static constexpr int value = 0;
};

template<typename Group, typename... Rest>
struct TotalDimension<Group, Rest...> {
    static constexpr int value = Group::Dim + TotalDimension<Rest...>::value;
};

// Helper to compute total action dimension
template<typename... Groups>
struct TotalActionDimension;

template<>
struct TotalActionDimension<> {
    static constexpr int value = 0;
};

template<typename Group, typename... Rest>
struct TotalActionDimension<Group, Rest...> {
    // This assumes we know actionDimension at compile time
    // For runtime computation, we'll use a different approach
    static constexpr int value = -1; // Indicates runtime computation needed
};

// Helper to check if all types are LieGroupBase derivatives with same scalar type
template<typename FirstScalar, typename... Groups>
struct AllAreLieGroups;

template<typename FirstScalar>
struct AllAreLieGroups<FirstScalar> : std::true_type {};

template<typename FirstScalar, typename Group, typename... Rest>
struct AllAreLieGroups<FirstScalar, Group, Rest...> {
    static constexpr bool value = 
        std::is_base_of_v<LieGroupBase<typename Group::Scalar, 
                                     std::integral_constant<int, Group::Dim>, 
                                     Group::Dim, Group::Dim>, Group> &&
        std::is_same_v<FirstScalar, typename Group::Scalar> &&
        AllAreLieGroups<FirstScalar, Rest...>::value;
};

// Compile-time offset computation
template<int Index, typename... Groups>
struct OffsetAt;

template<int Index>
struct OffsetAt<Index> {
    static constexpr int value = 0;
};

template<int Index, typename Group, typename... Rest>
struct OffsetAt<Index, Group, Rest...> {
    static constexpr int value = (Index == 0) ? 0 : 
                                Group::Dim + OffsetAt<Index - 1, Rest...>::value;
};

// Runtime offset computation for action dimensions
template<typename... Groups>
class ActionOffsets {
public:
    explicit ActionOffsets(const std::tuple<Groups...>& groups) {
        computeOffsets(groups, std::index_sequence_for<Groups...>());
    }
    
    int operator[](std::size_t i) const { return offsets_[i]; }
    int total() const { return total_; }

private:
    std::array<int, sizeof...(Groups)> offsets_;
    int total_ = 0;
    
    template<std::size_t... Is>
    void computeOffsets(const std::tuple<Groups...>& groups, std::index_sequence<Is...>) {
        int offset = 0;
        ((offsets_[Is] = offset, 
          offset += std::get<Is>(groups).actionDimension()), ...);
        total_ = offset;
    }
};

} // namespace detail

/**
 * @brief Implementation of product manifold bundle of Lie groups
 * 
 * This class implements a Cartesian product of multiple Lie groups, allowing them to be
 * treated as a single composite Lie group. The bundle maintains the product structure 
 * while providing all necessary group operations.
 * 
 * Mathematical Background:
 * For Lie groups G₁, G₂, ..., Gₙ, the product manifold G = G₁ × G₂ × ... × Gₙ
 * is also a Lie group with:
 * - Group operation: (g₁, g₂, ..., gₙ) * (h₁, h₂, ..., hₙ) = (g₁h₁, g₂h₂, ..., gₙhₙ)
 * - Identity: (e₁, e₂, ..., eₙ) where eᵢ is identity of Gᵢ
 * - Inverse: (g₁, g₂, ..., gₙ)⁻¹ = (g₁⁻¹, g₂⁻¹, ..., gₙ⁻¹)
 * - Lie algebra: 𝔤 = 𝔤₁ ⊕ 𝔤₂ ⊕ ... ⊕ 𝔤ₙ (direct sum)
 * - Adjoint: block diagonal with Adᵢ on diagonal blocks
 * 
 * Usage Examples:
 * ```cpp
 * // Rigid body pose with velocity
 * using RigidBodyState = Bundle<SE3<double>, RealSpace<double, 6>>;
 * 
 * // 2D robot with multiple joints
 * using Robot2D = Bundle<SE2<double>, SO2<double>, SO2<double>>;
 * 
 * // Create and manipulate
 * auto state1 = RigidBodyState(SE3<double>::identity(), 
 *                             RealSpace<double, 6>::zero());
 * auto state2 = RigidBodyState(pose, velocity);
 * auto combined = state1 * state2;
 * ```
 * 
 * @tparam Groups The Lie groups to bundle together (must have same scalar type)
 */
template<typename... Groups>
class Bundle : public LieGroupBase<typename std::tuple_element_t<0, std::tuple<Groups...>>::Scalar,
                                std::integral_constant<int, detail::TotalDimension<Groups...>::value>,
                                detail::TotalDimension<Groups...>::value,
                                detail::TotalDimension<Groups...>::value> {
    
    // Compile-time validation
    static_assert(sizeof...(Groups) > 0, "Bundle must contain at least one group");
    
    using FirstGroup = std::tuple_element_t<0, std::tuple<Groups...>>;
    using FirstScalar = typename FirstGroup::Scalar;
    
    static_assert(detail::AllAreLieGroups<FirstScalar, Groups...>::value,
                 "All template parameters must be Lie groups with the same scalar type");

public:
    using Base = LieGroupBase<FirstScalar,
                            std::integral_constant<int, detail::TotalDimension<Groups...>::value>,
                            detail::TotalDimension<Groups...>::value,
                            detail::TotalDimension<Groups...>::value>;
    using Scalar = typename Base::Scalar;
    using Vector = typename Base::Vector;
    using Matrix = typename Base::Matrix;
    using TangentVector = typename Base::TangentVector;
    using AdjointMatrix = typename Base::AdjointMatrix;
    
    static constexpr int Dim = Base::Dim;
    static constexpr std::size_t NumGroups = sizeof...(Groups);
    
    using GroupTuple = std::tuple<Groups...>;
    
    // Compile-time offset table for algebra elements
    template<std::size_t I>
    static constexpr int AlgebraOffset = detail::OffsetAt<I, Groups...>::value;
    
    template<std::size_t I>
    using GroupType = std::tuple_element_t<I, GroupTuple>;

    // ========== Constructors ==========

    /**
     * @brief Default constructor creates identity bundle
     */
    Bundle() : m_groups(), m_action_offsets(m_groups) {}

    /**
     * @brief Construct from individual group elements
     */
    explicit Bundle(const Groups&... groups) 
        : m_groups(groups...), m_action_offsets(m_groups) {}

    /**
     * @brief Construct from tuple of groups
     */
    explicit Bundle(const GroupTuple& groups) 
        : m_groups(groups), m_action_offsets(m_groups) {}

    /**
     * @brief Construct from Lie algebra vector
     */
    explicit Bundle(const TangentVector& algebra_element) 
        : Bundle() {
        *this = exp(algebra_element);
    }

    /**
     * @brief Copy constructor
     */
    Bundle(const Bundle& other) = default;

    /**
     * @brief Move constructor  
     */
    Bundle(Bundle&& other) noexcept = default;

    /**
     * @brief Copy assignment
     */
    Bundle& operator=(const Bundle& other) = default;

    /**
     * @brief Move assignment
     */
    Bundle& operator=(Bundle&& other) noexcept = default;

    // ========== Group Operations ==========

    /**
     * @brief Group composition (component-wise)
     * Implements: (g₁, ..., gₙ) * (h₁, ..., hₙ) = (g₁h₁, ..., gₙhₙ)
     */
    Bundle operator*(const Bundle& other) const {
        return multiply_impl(other, std::index_sequence_for<Groups...>());
    }

    /**
     * @brief In-place group composition
     */
    Bundle& operator*=(const Bundle& other) {
        multiply_assign_impl(other, std::index_sequence_for<Groups...>());
        return *this;
    }

    /**
     * @brief Inverse element (component-wise)
     * Implements: (g₁, ..., gₙ)⁻¹ = (g₁⁻¹, ..., gₙ⁻¹)
     */
    Bundle inverse() const override {
        return inverse_impl(std::index_sequence_for<Groups...>());
    }

    // ========== Lie Algebra Operations ==========

    /**
     * @brief Exponential map from Lie algebra to bundle
     * The Lie algebra of the product is the direct sum of individual algebras
     */
    Bundle exp(const TangentVector& algebra_element) const override {
        validateAlgebraElement(algebra_element);
        return exp_impl(algebra_element, std::index_sequence_for<Groups...>());
    }

    /**
     * @brief Logarithmic map from bundle to Lie algebra
     * Maps to the direct sum of individual Lie algebras
     */
    TangentVector log() const override {
        return log_impl(std::index_sequence_for<Groups...>());
    }

    /**
     * @brief Adjoint representation (block diagonal structure)
     * Ad_{(g₁,...,gₙ)} = diag(Ad_{g₁}, ..., Ad_{gₙ})
     */
    AdjointMatrix adjoint() const override {
        return adjoint_impl(std::index_sequence_for<Groups...>());
    }

    // ========== Group Actions ==========

    /**
     * @brief Group action on a point (component-wise on appropriate subspaces)
     * Each group acts on its corresponding portion of the input vector
     */
    Vector act(const Vector& point) const override {
        validateActionInput(point);
        return act_impl(point, std::index_sequence_for<Groups...>());
    }

    /**
     * @brief Batch group action on multiple points
     */
    template<int N>
    Eigen::Matrix<Scalar, Eigen::Dynamic, N> act(
        const Eigen::Matrix<Scalar, Eigen::Dynamic, N>& points) const {
        
        if (points.rows() != actionDimension()) {
            throw std::invalid_argument("Point matrix has wrong dimension");
        }
        
        Eigen::Matrix<Scalar, Eigen::Dynamic, N> result(actionDimension(), N);
        
        for (int i = 0; i < N; ++i) {
            result.col(i) = act(points.col(i));
        }
        
        return result;
    }

    // ========== Utility Functions ==========

    /**
     * @brief Check if approximately equal to another bundle
     */
    bool isApprox(const Bundle& other, 
                  const Scalar& eps = Types<Scalar>::epsilon()) const {
        return isApprox_impl(other, eps, std::index_sequence_for<Groups...>());
    }

    /**
     * @brief Equality operator
     */
    bool operator==(const Bundle& other) const {
        return isApprox(other);
    }

    /**
     * @brief Inequality operator
     */
    bool operator!=(const Bundle& other) const {
        return !(*this == other);
    }

    /**
     * @brief Linear interpolation between two bundles (geodesic in product space)
     * @param other Target bundle
     * @param t Interpolation parameter [0,1]
     */
    Bundle interpolate(const Bundle& other, const Scalar& t) const {
        if (t < 0 || t > 1) {
            throw std::invalid_argument("Interpolation parameter must be in [0,1]");
        }
        
        TangentVector delta = (inverse() * other).log();
        return *this * exp(t * delta);
    }

    /**
     * @brief Generate random bundle element
     */
    template<typename Generator>
    static Bundle random(Generator& gen, const Scalar& scale = Scalar(1)) {
        return random_impl(gen, scale, std::index_sequence_for<Groups...>());
    }

    // ========== Accessors ==========

    /**
     * @brief Get the identity element
     */
    static const Bundle& identity() {
        static const Bundle id;
        return id;
    }

    /**
     * @brief Get the dimension of the Lie algebra
     */
    int algebraDimension() const override { return Dim; }

    /**
     * @brief Get the dimension of the space the group acts on (computed at runtime)
     */
    int actionDimension() const override {
        return m_action_offsets.total();
    }

    /**
     * @brief Access individual group elements (const)
     */
    template<std::size_t I>
    const auto& get() const {
        static_assert(I < NumGroups, "Index out of bounds");
        return std::get<I>(m_groups);
    }

    /**
     * @brief Access individual group elements (mutable)
     */
    template<std::size_t I>
    auto& get() {
        static_assert(I < NumGroups, "Index out of bounds");
        // Need to recompute action offsets if groups are modified
        auto& result = std::get<I>(m_groups);
        // In practice, you might want to make this const and provide setters
        return result;
    }

    /**
     * @brief Set individual group element
     */
    template<std::size_t I>
    void set(const GroupType<I>& group) {
        std::get<I>(m_groups) = group;
        m_action_offsets = detail::ActionOffsets<Groups...>(m_groups);
    }

    /**
     * @brief Get the underlying tuple
     */
    const GroupTuple& groups() const { return m_groups; }

    /**
     * @brief Get algebra element for specific group
     */
    template<std::size_t I>
    typename GroupType<I>::TangentVector getAlgebraElement() const {
        return std::get<I>(m_groups).log();
    }

    /**
     * @brief Set from algebra element for specific group  
     */
    template<std::size_t I>
    void setFromAlgebra(const typename GroupType<I>::TangentVector& algebra) {
        std::get<I>(m_groups) = GroupType<I>().exp(algebra);
        m_action_offsets = detail::ActionOffsets<Groups...>(m_groups);
    }

    // ========== Stream Output ==========

    /**
     * @brief Output stream operator
     */
    friend std::ostream& operator<<(std::ostream& os, const Bundle& bundle) {
        os << "Bundle<" << NumGroups << ">(";
        bundle.print_impl(os, std::index_sequence_for<Groups...>());
        os << ")";
        return os;
    }

private:
    GroupTuple m_groups;  ///< Tuple of group elements
    detail::ActionOffsets<Groups...> m_action_offsets;  ///< Cached action dimension offsets

    // ========== Implementation Helpers ==========

    // Validation helpers
    void validateAlgebraElement(const TangentVector& element) const {
        if (element.size() != Dim) {
            throw std::invalid_argument("Algebra element has wrong dimension");
        }
    }

    void validateActionInput(const Vector& point) const {
        if (point.size() != actionDimension()) {
            throw std::invalid_argument("Action input has wrong dimension");
        }
    }

    // Multiplication implementation
    template<std::size_t... Is>
    Bundle multiply_impl(const Bundle& other, std::index_sequence<Is...>) const {
        return Bundle((std::get<Is>(m_groups) * std::get<Is>(other.m_groups))...);
    }

    template<std::size_t... Is>
    void multiply_assign_impl(const Bundle& other, std::index_sequence<Is...>) {
        ((std::get<Is>(m_groups) *= std::get<Is>(other.m_groups)), ...);
    }

    // Inverse implementation
    template<std::size_t... Is>
    Bundle inverse_impl(std::index_sequence<Is...>) const {
        return Bundle((std::get<Is>(m_groups).inverse())...);
    }

    // Exponential map implementation
    template<std::size_t... Is>
    Bundle exp_impl(const TangentVector& algebra_element, std::index_sequence<Is...>) const {
        return Bundle((GroupType<Is>().exp(
            algebra_element.template segment<GroupType<Is>::Dim>(AlgebraOffset<Is>)
        ))...);
    }

    // Logarithmic map implementation
    template<std::size_t... Is>
    TangentVector log_impl(std::index_sequence<Is...>) const {
        TangentVector result;
        ((result.template segment<GroupType<Is>::Dim>(AlgebraOffset<Is>) = 
          std::get<Is>(m_groups).log()), ...);
        return result;
    }

    // Adjoint implementation (block diagonal)
    template<std::size_t... Is>
    AdjointMatrix adjoint_impl(std::index_sequence<Is...>) const {
        AdjointMatrix result = AdjointMatrix::Zero();
        ((result.template block<GroupType<Is>::Dim, GroupType<Is>::Dim>
          (AlgebraOffset<Is>, AlgebraOffset<Is>) = std::get<Is>(m_groups).adjoint()), ...);
        return result;
    }

    // Group action implementation
    template<std::size_t... Is>
    Vector act_impl(const Vector& point, std::index_sequence<Is...>) const {
        Vector result(actionDimension());
        
        // Apply each group's action to its corresponding subspace
        ((applyGroupAction<Is>(result, point)), ...);
        
        return result;
    }

    template<std::size_t I>
    void applyGroupAction(Vector& result, const Vector& point) const {
        const auto& group = std::get<I>(m_groups);
        int in_offset = m_action_offsets[I];
        int out_offset = in_offset; // Same offset for output
        int dim = group.actionDimension();
        
        auto input_segment = point.segment(in_offset, dim);
        auto output_segment = result.segment(out_offset, dim);
        output_segment = group.act(input_segment);
    }

    // Approximate equality implementation
    template<std::size_t... Is>
    bool isApprox_impl(const Bundle& other, const Scalar& eps, 
                      std::index_sequence<Is...>) const {
        return (std::get<Is>(m_groups).isApprox(std::get<Is>(other.m_groups), eps) && ...);
    }

    // Random generation implementation
    template<typename Generator, std::size_t... Is>
    static Bundle random_impl(Generator& gen, const Scalar& scale, 
                             std::index_sequence<Is...>) {
        return Bundle((GroupType<Is>::random(gen, scale))...);
    }

    // Stream output implementation
    template<std::size_t... Is>
    void print_impl(std::ostream& os, std::index_sequence<Is...>) const {
        bool first = true;
        ((os << (first ? (first = false, "") : ", ") << std::get<Is>(m_groups)), ...);
    }
};

// ========== Type Aliases ==========

// Common bundles for robotics applications
template<typename Scalar>
using SE3_Velocity = Bundle<SE3<Scalar>, RealSpace<Scalar, 6>>;

template<typename Scalar>
using SE2_Velocity = Bundle<SE2<Scalar>, RealSpace<Scalar, 3>>;

template<typename Scalar, int N>
using SE3_Joints = Bundle<SE3<Scalar>, RealSpace<Scalar, N>>;

// Convenience aliases
template<typename... Groups>
using Bundlef = Bundle<Groups...>;

template<typename... Groups>
using Bundled = Bundle<Groups...>;

} // namespace sofa::component::cosserat::liegroups