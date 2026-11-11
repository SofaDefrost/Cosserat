// This file defines the Bundle class, which implements a product manifold of
// multiple Lie groups.

/******************************************************************************
 *                 SOFA, Simulation Open-Framework Architecture * (c) 2006
 *INRIA, USTL, UJF, CNRS, MGH                     *
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
 * along with this program. If not, see <http://www.gnu.org/licenses/\>. *
 ******************************************************************************/

#pragma once

#include <array>
#include <iostream>
#include <random>
#include <tuple>
#include <type_traits>
#include "LieGroup.h"
#include "LieGroup.inl"
#include "Rn.h"
#include "SE2.h"
#include "SE3.h"
#include "SO2.h"
#include "Types.h"
// Rn.h is already included above

namespace sofa::component::cosserat::liegroups {

	namespace detail {

		/// Sum of Group::<member> over the pack, evaluated at compile time. Dim (size of the matrix
		/// representation), AlgebraDim (tangent coordinates) and ActionDim (size of the vectors acted on)
		/// differ for SE3 (4, 6, 3), SO2 (2, 1, 2), SE2 (3, 3, 2): each one needs its own offsets.
		template<typename... Groups>
		struct Sums {
			static constexpr int dim = (0 + ... + Groups::Dim);
			static constexpr int algebra = (0 + ... + Groups::AlgebraDim);
			static constexpr int action = (0 + ... + Groups::ActionDim);
		};

		/// Offset of group I in the concatenation of Groups, for the member selected by Which
		/// (0: Dim, 1: AlgebraDim, 2: ActionDim).
		template<int Which, std::size_t I, typename... Groups>
		struct OffsetOf {
			static constexpr int value = [] {
				constexpr int sizes[] = {(Which == 0 ? Groups::Dim : Which == 1 ? Groups::AlgebraDim : Groups::ActionDim)...};
				int off = 0;
				for (std::size_t k = 0; k < I; ++k)
					off += sizes[k];
				return off;
			}();
		};

		template<typename FirstScalar, typename... Groups>
		struct AllAreLieGroups {
			static constexpr bool value =
					(std::is_base_of_v<LieGroup<GroupProperties<Groups, typename Groups::Scalar, Groups::Dim,
																Groups::AlgebraDim, Groups::ActionDim>>,
									   Groups> &&
					 ...) &&
					(std::is_same_v<FirstScalar, typename Groups::Scalar> && ...);
		};

	} // namespace detail

	/**
	 * @brief Implementation of product manifold bundle of Lie groups
	 *
	 * This class implements a Cartesian product of multiple Lie groups, allowing
	 * them to be treated as a single composite Lie group. The bundle maintains the
	 * product structure while providing all necessary group operations.
	 *
	 * Mathematical Background:
	 * For Lie groups G₁, G₂, ..., Gₙ, the product manifold G = G₁ × G₂ × ... × Gₙ
	 * is also a Lie group with:
	 * - Group operation: (g₁, g₂, ..., gₙ) * (h₁, h₂, ..., hₙ) = (g₁h₁, g₂h₂, ...,
	 * gₙhₙ)
	 * - Identity: (e₁, e₂, ..., eₙ) where eᵢ is identity of Gᵢ
	 * - Inverse: (g₁, g₂, ..., gₙ)⁻¹ = (g₁⁻¹, g₂⁻¹, ..., gₙ⁻¹)
	 * - Lie algebra: 𝔤 = 𝔤₁ ⊕ 𝔤₂ ⊕ ... ⊕ 𝔤ₙ (direct sum)
	 * - Adjoint: block diagonal with Adᵢ on diagonal blocks
	 *
	 * Usage Examples:
	 * ```cpp
	 * // Rigid body pose with velocity
	 * using RigidBodyState = Bundle<SE3<double>, Rn<double, 6>>;
	 *
	 * // 2D robot with multiple joints
	 * using Robot2D = Bundle<SE2<double>, SO2<double>, SO2<double>>;
	 *
	 * // Create and manipulate
	 * auto state1 = RigidBodyState(SE3<double>::identity(),
	 *                             Rn<double, 6>::zero());
	 * auto state2 = RigidBodyState(pose, velocity);
	 * auto combined = state1 * state2;
	 * ```
	 *
	 * @tparam Groups The Lie groups to bundle together (must have same scalar type)
	 */
	template<typename... Groups>
	class Bundle
		: public LieGroup<GroupProperties<Bundle<Groups...>, typename std::tuple_element_t<0, std::tuple<Groups...>>::Scalar,
										  detail::Sums<Groups...>::dim, detail::Sums<Groups...>::algebra,
										  detail::Sums<Groups...>::action>> {
		static_assert(sizeof...(Groups) > 0, "Bundle must contain at least one group");
		using FirstGroup = std::tuple_element_t<0, std::tuple<Groups...>>;
		using FirstScalar = typename FirstGroup::Scalar;
		static_assert(detail::AllAreLieGroups<FirstScalar, Groups...>::value,
					  "All template parameters must be Lie groups with the same scalar type");

	public:
		using Base = LieGroup<GroupProperties<Bundle<Groups...>, FirstScalar, detail::Sums<Groups...>::dim,
											  detail::Sums<Groups...>::algebra, detail::Sums<Groups...>::action>>;
		using Scalar = typename Base::Scalar;
		using Vector = typename Base::Vector;
		using Matrix = typename Base::Matrix;
		using TangentVector = typename Base::TangentVector;
		using AdjointMatrix = typename Base::AdjointMatrix;
		using ActionVector = typename Base::ActionVector;
		static constexpr int Dim = Base::Dim;
		static constexpr int AlgebraDim = Base::AlgebraDim;
		static constexpr int ActionDim = Base::ActionDim;
		static constexpr std::size_t NumGroups = sizeof...(Groups);

		using GroupTuple = std::tuple<Groups...>;

		/// Offset of group I among the tangent coordinates.
		template<std::size_t I>
		static constexpr int AlgebraOffset = detail::OffsetOf<1, I, Groups...>::value;
		/// Offset of group I among the coordinates of the vectors the bundle acts on.
		template<std::size_t I>
		static constexpr int ActionOffset = detail::OffsetOf<2, I, Groups...>::value;
		/// Offset of group I in the (block diagonal) matrix representation.
		template<std::size_t I>
		static constexpr int DimOffset = detail::OffsetOf<0, I, Groups...>::value;

		template<std::size_t I>
		using GroupType = std::tuple_element_t<I, GroupTuple>;

		Bundle() : m_groups() {}
		explicit Bundle(const Groups &...groups) : m_groups(groups...) {}
		explicit Bundle(const GroupTuple &groups) : m_groups(groups) {}
		explicit Bundle(const TangentVector &algebra_element) : Bundle() { *this = computeExp(algebra_element); }

		Bundle(const Bundle &other) = default;
		Bundle(Bundle &&other) noexcept = default;
		Bundle &operator=(const Bundle &other) = default;
		Bundle &operator=(Bundle &&other) noexcept = default;

		// ---- group operations (also the CRTP hooks of LieGroup) ----

		Bundle operator*(const Bundle &other) const { return computeCompose(other); }
		Bundle computeCompose(const Bundle &other) const { return multiply_impl(other, std::index_sequence_for<Groups...>()); }
		Bundle &operator*=(const Bundle &other) {
			multiply_assign_impl(other, std::index_sequence_for<Groups...>());
			return *this;
		}

		Bundle computeInverse() const { return inverse_impl(std::index_sequence_for<Groups...>()); }

		static Bundle computeExp(const TangentVector &algebra_element) {
			return exp_impl(algebra_element, std::index_sequence_for<Groups...>());
		}

		TangentVector computeLog() const { return log_impl(std::index_sequence_for<Groups...>()); }

		AdjointMatrix computeAdjoint() const { return adjoint_impl(std::index_sequence_for<Groups...>()); }

		ActionVector act(const ActionVector &point) const { return computeAction(point); }
		ActionVector computeAction(const ActionVector &point) const {
			return act_impl(point, std::index_sequence_for<Groups...>());
		}

		template<int N>
		Eigen::Matrix<Scalar, Eigen::Dynamic, N> act(const Eigen::Matrix<Scalar, Eigen::Dynamic, N> &points) const {
			if (points.rows() != ActionDim) {
				throw std::invalid_argument("Point matrix has wrong dimension");
			}
			Eigen::Matrix<Scalar, Eigen::Dynamic, N> result(ActionDim, N);
			for (int i = 0; i < N; ++i) {
				result.col(i) = act(ActionVector(points.col(i)));
			}
			return result;
		}

		bool computeIsApprox(const Bundle &other, const Scalar &eps = Types<Scalar>::tolerance()) const {
			return isApprox_impl(other, eps, std::index_sequence_for<Groups...>());
		}
		bool operator==(const Bundle &other) const { return computeIsApprox(other); }
		bool operator!=(const Bundle &other) const { return !(*this == other); }

		Bundle interpolate(const Bundle &other, const Scalar &t) const {
			if (t < 0 || t > 1) {
				throw std::invalid_argument("Interpolation parameter must be in [0,1]");
			}
			const TangentVector delta = (computeInverse() * other).computeLog();
			return *this * computeExp(t * delta);
		}

		/// Random element, each group drawn with its own Random(gen); `scale` is kept for API compatibility
		/// and ignored (the groups have no scale parameter).
		template<typename Generator>
		static Bundle random(Generator &gen, [[maybe_unused]] const Scalar &scale = Scalar(1)) {
			return computeRandom(gen);
		}
		template<typename Generator>
		static Bundle computeRandom(Generator &gen) {
			return Bundle(Groups::computeRandom(gen)...);
		}

		static const Bundle &identity() {
			static const Bundle id;
			return id;
		}
		static Bundle computeIdentity() noexcept { return Bundle(); }

		int algebraDimension() const { return AlgebraDim; }
		int actionDimension() const { return ActionDim; }

		// ---- Lie algebra hooks: block diagonal, so that the generic ad / Jacobians / BCH work ----

		/// Block-diagonal embedding of the groups' own matrix representations.
		static Matrix computeHat(const TangentVector &v) {
			Matrix result = Matrix::Zero();
			hat_impl(result, v, std::index_sequence_for<Groups...>());
			return result;
		}
		static TangentVector computeVee(const Matrix &X) {
			TangentVector result;
			vee_impl(result, X, std::index_sequence_for<Groups...>());
			return result;
		}
		static AdjointMatrix computeAd(const TangentVector &v) {
			AdjointMatrix result = AdjointMatrix::Zero();
			ad_impl(result, v, std::index_sequence_for<Groups...>());
			return result;
		}

		bool computeIsValid() const { return isValid_impl(std::index_sequence_for<Groups...>()); }
		void computeNormalize() { normalize_impl(std::index_sequence_for<Groups...>()); }
		static constexpr std::string_view getTypeName() { return "Bundle"; }
		std::ostream &print(std::ostream &os) const {
			os << "Bundle<" << NumGroups << ">(";
			print_impl(os, std::index_sequence_for<Groups...>());
			os << ")";
			return os;
		}

		template<std::size_t I>
		const auto &get() const {
			static_assert(I < NumGroups, "Index out of bounds");
			return std::get<I>(m_groups);
		}
		template<std::size_t I>
		auto &get() {
			static_assert(I < NumGroups, "Index out of bounds");
			return std::get<I>(m_groups);
		}
		template<std::size_t I>
		void set(const GroupType<I> &group) {
			std::get<I>(m_groups) = group;
		}
		const GroupTuple &groups() const { return m_groups; }

		template<std::size_t I>
		typename GroupType<I>::TangentVector getAlgebraElement() const {
			return std::get<I>(m_groups).log();
		}
		template<std::size_t I>
		void setFromAlgebra(const typename GroupType<I>::TangentVector &algebra) {
			std::get<I>(m_groups) = GroupType<I>::exp(algebra);
		}

		friend std::ostream &operator<<(std::ostream &os, const Bundle &bundle) { return bundle.print(os); }

	private:
		GroupTuple m_groups; ///< Tuple of group elements

		template<std::size_t... Is>
		Bundle multiply_impl(const Bundle &other, std::index_sequence<Is...>) const {
			return Bundle((std::get<Is>(m_groups) * std::get<Is>(other.m_groups))...);
		}
		template<std::size_t... Is>
		void multiply_assign_impl(const Bundle &other, std::index_sequence<Is...>) {
			((std::get<Is>(m_groups) = std::get<Is>(m_groups) * std::get<Is>(other.m_groups)), ...);
		}
		template<std::size_t... Is>
		Bundle inverse_impl(std::index_sequence<Is...>) const {
			return Bundle((std::get<Is>(m_groups).inverse())...);
		}
		template<std::size_t... Is>
		static Bundle exp_impl(const TangentVector &x, std::index_sequence<Is...>) {
			return Bundle((GroupType<Is>::exp(typename GroupType<Is>::TangentVector(
								  x.template segment<GroupType<Is>::AlgebraDim>(AlgebraOffset<Is>))))...);
		}
		template<std::size_t... Is>
		TangentVector log_impl(std::index_sequence<Is...>) const {
			TangentVector result;
			((result.template segment<GroupType<Is>::AlgebraDim>(AlgebraOffset<Is>) = std::get<Is>(m_groups).log()),
			 ...);
			return result;
		}
		template<std::size_t... Is>
		AdjointMatrix adjoint_impl(std::index_sequence<Is...>) const {
			AdjointMatrix result = AdjointMatrix::Zero();
			((result.template block<GroupType<Is>::AlgebraDim, GroupType<Is>::AlgebraDim>(
					  AlgebraOffset<Is>, AlgebraOffset<Is>) = std::get<Is>(m_groups).adjoint()),
			 ...);
			return result;
		}
		template<std::size_t... Is>
		ActionVector act_impl(const ActionVector &point, std::index_sequence<Is...>) const {
			ActionVector result;
			((result.template segment<GroupType<Is>::ActionDim>(ActionOffset<Is>) = std::get<Is>(m_groups).act(
					  typename GroupType<Is>::ActionVector(point.template segment<GroupType<Is>::ActionDim>(ActionOffset<Is>)))),
			 ...);
			return result;
		}
		template<std::size_t... Is>
		bool isApprox_impl(const Bundle &other, const Scalar &eps, std::index_sequence<Is...>) const {
			return (std::get<Is>(m_groups).isApprox(std::get<Is>(other.m_groups), eps) && ...);
		}
		template<std::size_t... Is>
		static void hat_impl(Matrix &out, const TangentVector &v, std::index_sequence<Is...>) {
			((out.template block<GroupType<Is>::Dim, GroupType<Is>::Dim>(DimOffset<Is>, DimOffset<Is>) =
					  GroupType<Is>::computeHat(typename GroupType<Is>::TangentVector(
							  v.template segment<GroupType<Is>::AlgebraDim>(AlgebraOffset<Is>)))),
			 ...);
		}
		template<std::size_t... Is>
		static void vee_impl(TangentVector &out, const Matrix &X, std::index_sequence<Is...>) {
			((out.template segment<GroupType<Is>::AlgebraDim>(AlgebraOffset<Is>) = GroupType<Is>::computeVee(
					  typename GroupType<Is>::AlgebraMatrix(
							  X.template block<GroupType<Is>::Dim, GroupType<Is>::Dim>(DimOffset<Is>, DimOffset<Is>)))),
			 ...);
		}
		template<std::size_t... Is>
		static void ad_impl(AdjointMatrix &out, const TangentVector &v, std::index_sequence<Is...>) {
			((out.template block<GroupType<Is>::AlgebraDim, GroupType<Is>::AlgebraDim>(
					  AlgebraOffset<Is>, AlgebraOffset<Is>) =
					  GroupType<Is>::computeAd(typename GroupType<Is>::TangentVector(
							  v.template segment<GroupType<Is>::AlgebraDim>(AlgebraOffset<Is>)))),
			 ...);
		}
		template<std::size_t... Is>
		bool isValid_impl(std::index_sequence<Is...>) const {
			return (std::get<Is>(m_groups).isValid() && ...);
		}
		template<std::size_t... Is>
		void normalize_impl(std::index_sequence<Is...>) {
			(std::get<Is>(m_groups).normalize(), ...);
		}
		template<std::size_t... Is>
		void print_impl(std::ostream &os, std::index_sequence<Is...>) const {
			bool first = true;
			((os << (first ? (first = false, "") : ", ") << std::get<Is>(m_groups)), ...);
		}
	};

	template<typename Scalar>
	using SE3_Velocity = Bundle<SE3<Scalar>, Rn<Scalar, 6>>;

	template<typename Scalar>
	using SE2_Velocity = Bundle<SE2<Scalar>, Rn<Scalar, 3>>;

	template<typename Scalar, int N>
	using SE3_Joints = Bundle<SE3<Scalar>, Rn<Scalar, N>>;

	// Convenience aliases
	template<typename... Groups>
	using Bundlef = Bundle<Groups...>;

	template<typename... Groups>
	using Bundled = Bundle<Groups...>;

} // namespace sofa::component::cosserat::liegroups
