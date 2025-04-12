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

#ifndef SOFA_COMPONENT_COSSERAT_LIEGROUPS_BUNDLE_H
#define SOFA_COMPONENT_COSSERAT_LIEGROUPS_BUNDLE_H

#include "LieGroupBase.h"
#include "LieGroupBase.inl"
#include <tuple>
#include <type_traits>

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

// Helper to check if all types are LieGroupBase derivatives
template<typename... Groups>
struct AllAreLieGroups;

template<>
struct AllAreLieGroups<> : std::true_type {};

template<typename Group, typename... Rest>
struct AllAreLieGroups<Group, Rest...> {
    static constexpr bool value = 
        std::is_base_of_v<LieGroupBase<typename Group::Scalar, Group::Dim>, Group> &&
        AllAreLieGroups<Rest...>::value;
};

// Helper to split vector into segments
template<typename Vector, int offset, typename Group>
void setSegment(Vector& vec, const Group& group) {
    vec.template segment<Group::Dim>(offset) = group.log();
}

template<typename Vector, int offset, typename Group, typename... Rest>
void setSegment(Vector& vec, const Group& group, const Rest&... rest) {
    vec.template segment<Group::Dim>(offset) = group.log();
    setSegment<Vector, offset + Group::Dim>(vec, rest...);
}

// Helper to extract segments and create groups
template<typename Vector, int offset, typename Group>
Group getGroup(const Vector& vec) {
    return Group().exp(vec.template segment<Group::Dim>(offset));
}

} // namespace detail

/**
 * @brief Implementation of product manifold bundle of Lie groups
 * 
 * This class implements a product of multiple Lie groups, allowing them to be
 * treated as a single composite Lie group. The bundle maintains the product
 * structure while providing all the necessary group operations.
 * 
 * Example usage:
 * using PoseVel = Bundle<SE3<double>, RealSpace<double, 3>>;
 * 
 * @tparam Groups The Lie groups to bundle together
 */
template<typename... Groups>
class Bundle : public LieGroupBase<typename std::tuple_element_t<0, std::tuple<Groups...>>::Scalar,
                                detail::TotalDimension<Groups...>::value>,
              public LieGroupOperations<Bundle<Groups...>> {
    static_assert(sizeof...(Groups) > 0, "Bundle must contain at least one group");
    static_assert(detail::AllAreLieGroups<Groups...>::value,
                 "All template parameters must be Lie groups");

public:
    using FirstGroup = std::tuple_element_t<0, std::tuple<Groups...>>;
    using Base = LieGroupBase<typename FirstGroup::Scalar,
                            detail::TotalDimension<Groups...>::value>;
    using Scalar = typename Base::Scalar;
    using Vector = typename Base::Vector;
    using Matrix = typename Base::Matrix;
    using TangentVector = typename Base::TangentVector;
    using AdjointMatrix = typename Base::AdjointMatrix;
    
    static constexpr int Dim = Base::Dim;
    using GroupTuple = std::tuple<Groups...>;

    /**
     * @brief Default constructor creates identity bundle
     */
    Bundle() = default;

    /**
     * @brief Construct from individual group elements
     */
    Bundle(const Groups&... groups) : m_groups(groups...) {}

    /**
     * @brief Group composition (component-wise)
     */
    Bundle operator*(const Bundle& other) const {
        return multiply(other, std::index_sequence_for<Groups...>());
    }

    /**
     * @brief Inverse element (component-wise)
     */
    Bundle inverse() const override {
        return inverse_impl(std::index_sequence_for<Groups...>());
    }

    /**
     * @brief Exponential map from Lie algebra to bundle
     */
    Bundle exp(const TangentVector& algebra_element) const override {
        return exp_impl(algebra_element, std::index_sequence_for<Groups...>());
    }

    /**
     * @brief Logarithmic map from bundle to Lie algebra
     */
    TangentVector log() const override {
        TangentVector result;
        detail::setSegment<TangentVector, 0>(result, std::get<Groups>(m_groups)...);
        return result;
    }

    /**
     * @brief Adjoint representation (block diagonal)
     */
    AdjointMatrix adjoint() const override {
        return adjoint_impl(std::index_sequence_for<Groups...>());
    }

    /**
     * @brief Group action on a point (component-wise)
     */
    Vector act(const Vector& point) const override {
        return act_impl(point, std::index_sequence_for<Groups...>());
    }

    /**
     * @brief Check if approximately equal to another bundle
     */
    bool isApprox(const Bundle& other,
                  const Scalar& eps = Types<Scalar>::epsilon()) const {
        return isApprox_impl(other, eps, std::index_sequence_for<Groups...>());
    }

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
     * @brief Get the dimension of the space the group acts on
     */
    int actionDimension() const override {
        return (std::get<Groups>(m_groups).actionDimension() + ...);
    }

    /**
     * @brief Access individual group elements
     */
    template<std::size_t I>
    const auto& get() const {
        return std::get<I>(m_groups);
    }

    template<std::size_t I>
    auto& get() {
        return std::get<I>(m_groups);
    }

private:
    GroupTuple m_groups;  ///< Tuple of group elements

    // Helper for multiplication
    template<std::size_t... Is>
    Bundle multiply(const Bundle& other, std::index_sequence<Is...>) const {
        return Bundle((std::get<Is>(m_groups) * std::get<Is>(other.m_groups))...);
    }

    // Helper for inverse
    template<std::size_t... Is>
    Bundle inverse_impl(std::index_sequence<Is...>) const {
        return Bundle((std::get<Is>(m_groups).inverse())...);
    }

    // Helper for exponential map
    template<std::size_t... Is>
    Bundle exp_impl(const TangentVector& algebra_element, std::index_sequence<Is...>) const {
        int offset = 0;
        return Bundle((detail::getGroup<TangentVector, 
                      (offset = (Is == 0 ? 0 : 
                               offset + std::tuple_element_t<Is-1, GroupTuple>::Dim))>
                      (algebra_element))...);
    }

    // Helper for adjoint
    template<std::size_t... Is>
    AdjointMatrix adjoint_impl(std::index_sequence<Is...>) const {
        AdjointMatrix result = AdjointMatrix::Zero();
        int offset = 0;
        ((result.template block<std::tuple_element_t<Is, GroupTuple>::Dim,
                              std::tuple_element_t<Is, GroupTuple>::Dim>
          (offset, offset) = std::get<Is>(m_groups).adjoint(),
          offset += std::tuple_element_t<Is, GroupTuple>::Dim), ...);
        return result;
    }

    // Helper for group action
    template<std::size_t... Is>
    Vector act_impl(const Vector& point, std::index_sequence<Is...>) const {
        Vector result;
        int inOffset = 0, outOffset = 0;
        ((result.template segment<std::get<Is>(m_groups).actionDimension()>(outOffset) = 
          std::get<Is>(m_groups).act(point.template segment<std::get<Is>(m_groups).actionDimension()>(inOffset)),
          inOffset += std::get<Is>(m_groups).actionDimension(),
          outOffset += std::get<Is>(m_groups).actionDimension()), ...);
        return result;
    }

    // Helper for approximate equality
    template<std::size_t... Is>
    bool isApprox_impl(const Bundle& other, const Scalar& eps, std::index_sequence<Is...>) const {
        return (std::get<Is>(m_groups).isApprox(std::get<Is>(other.m_groups), eps) && ...);
    }
};

} // namespace sofa::component::cosserat::liegroups

#include "Bundle.inl"

#endif // SOFA_COMPONENT_COSSERAT_LIEGROUPS_BUNDLE_H
