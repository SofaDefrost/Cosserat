/******************************************************************************
 *       SOFA, Simulation Open-Framework Architecture, development version     *
 *                (c) 2006-2019 INRIA, USTL, UJF, CNRS, MGH                    *
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
 * along with this program. If not, see <http://www.gnu.org/licenses/>.        *
 ******************************************************************************/
#pragma once

#include <Cosserat/config.h>
#include <liegroups/LinearStrainElement.h>
#include <sofa/core/Multi2Mapping.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/type/RGBAColor.h>
#include <sofa/type/fixed_array.h>
#include <sofa/type/vector.h>

#include <Eigen/StdVector>
#include <type_traits>
#include <vector>

namespace Cosserat::mapping {

	/**
	 * @brief Linear Strain Element (LSE) mapping: nodal poses + strain slopes -> mean strains.
	 *
	 * Implements the geometrically explicit Cosserat-rod element of Xun & Tamadazte,
	 * "Geometrically explicit Cosserat-rod modeling with piecewise linear strain for complex rod
	 * systems", Computers & Structures 325 (2026) 108181 (equation numbers below refer to it).
	 *
	 * The rod is a graph: `elements` lists node pairs (a, b), so branched rods, closed loops and
	 * lattices are assembled by sharing nodes, exactly like nodal finite elements.
	 *
	 *   parent 1 : Rigid3 nodal poses g_i            (carry the mass, contacts, boundary conditions)
	 *   parent 2 : Vec6 strain slopes β_e            (one per element)
	 *   child    : Vec6 mean strains ξ̄_e            (carry the elastic stiffness)
	 *
	 * Per element e = (a, b) of length h:
	 *
	 *     ξ̄_e = A⁻¹ Log(g_a⁻¹ g_b),    A = h·I − (h³/12) ad_β                    (Eq. 6)
	 *
	 * (the strain varies linearly, ξ(s) = ξ̄ + (s − h/2) β, and g_b = g_a Exp(A ξ̄) is the
	 * fourth-order Magnus end pose, Eq. 4-5). With β = 0 it reduces to the constant-strain
	 * element ξ̄ = Log(g_a⁻¹ g_b)/h used by Frames2StrainCosseratMapping.
	 *
	 * The output is the *full* se(3) strain ξ̄ = [κ ; ε] (a straight unstretched rod has
	 * ε = (1, 0, 0)), the same convention as Frames2StrainCosseratMapping.
	 *
	 * Jacobians (Eq. 7-10) are expressed in the body frame of each node; SOFA's world-frame
	 * `[v ; ω]` velocities are converted with the projector P(R) = [[0 R],[R 0]], as in the other
	 * Cosserat mappings. The kernel is liegroups/LinearStrainElement.h.
	 *
	 * Tangent: `applyDJT` is deliberately not provided. The paper's tangent is the Gauss–Newton
	 * one (Eq. 30), Jᵀ(hK)J with Jacobians frozen, which is what the assembled mapping pull-back
	 * J^T K J yields. It is exact when the strain equals its rest value.
	 *
	 * @tparam TIn1 nodal poses  (Rigid3Types)
	 * @tparam TIn2 slopes       (Vec6Types)
	 * @tparam TOut mean strains (Vec6Types)
	 */
	template<class TIn1, class TIn2, class TOut>
	class LinearStrainElementMapping : public sofa::core::Multi2Mapping<TIn1, TIn2, TOut> {
	public:
		SOFA_CLASS(SOFA_TEMPLATE3(LinearStrainElementMapping, TIn1, TIn2, TOut),
				   SOFA_TEMPLATE3(sofa::core::Multi2Mapping, TIn1, TIn2, TOut));

		using Inherit = sofa::core::Multi2Mapping<TIn1, TIn2, TOut>;
		using In1 = TIn1;
		using In2 = TIn2;
		using Out = TOut;

		using LSE = sofa::component::cosserat::liegroups::LinearStrainElement<double>;
		using Vec6 = LSE::Vec6;
		using Mat6 = LSE::Mat6;
		using Matrix3 = Eigen::Matrix3d;
		using Edge = sofa::type::fixed_array<unsigned int, 2>;

		static_assert(std::is_same_v<TIn1, sofa::defaulttype::Rigid3Types>,
					  "LinearStrainElementMapping: the first input must be Rigid3Types");

		// ── Data ────────────────────────────────────────────────────────────
		/// Node pairs (a, b), indices in the pose (first input) state.
		sofa::Data<sofa::type::vector<Edge>> d_elements;
		/// Reference length h_e of each element. If empty, the chord length of the initial poses
		/// is used (exact for straight elements; give the arc length for initially curved ones).
		sofa::Data<sofa::type::vector<SReal>> d_lengths;
		/// Colour used to draw the elements.
		sofa::Data<sofa::type::RGBAColor> d_color;

	public:
		void init() override;
		void draw(const sofa::core::visual::VisualParams *vparams) override;

		void apply(const sofa::core::MechanicalParams *mparams,
				   const sofa::type::vector<sofa::DataVecCoord_t<Out> *> &dataVecOutPos,
				   const sofa::type::vector<const sofa::DataVecCoord_t<In1> *> &dataVecIn1Pos,
				   const sofa::type::vector<const sofa::DataVecCoord_t<In2> *> &dataVecIn2Pos) override;

		void applyJ(const sofa::core::MechanicalParams *mparams,
					const sofa::type::vector<sofa::DataVecDeriv_t<Out> *> &dataVecOutVel,
					const sofa::type::vector<const sofa::DataVecDeriv_t<In1> *> &dataVecIn1Vel,
					const sofa::type::vector<const sofa::DataVecDeriv_t<In2> *> &dataVecIn2Vel) override;

		void applyJT(const sofa::core::MechanicalParams *mparams,
					 const sofa::type::vector<sofa::DataVecDeriv_t<In1> *> &dataVecOut1Force,
					 const sofa::type::vector<sofa::DataVecDeriv_t<In2> *> &dataVecOut2Force,
					 const sofa::type::vector<const sofa::DataVecDeriv_t<Out> *> &dataVecInForce) override;

		void applyJT(const sofa::core::ConstraintParams *cparams,
					 const sofa::type::vector<sofa::DataMatrixDeriv_t<In1> *> &dataMatOut1Const,
					 const sofa::type::vector<sofa::DataMatrixDeriv_t<In2> *> &dataMatOut2Const,
					 const sofa::type::vector<const sofa::DataMatrixDeriv_t<Out> *> &dataMatInConst) override;

		/// Geometric stiffness ∂(Jᵀf)/∂q: intentionally a no-op. The tangent of the reference
		/// formulation is Gauss–Newton (Eq. 30): Jacobians frozen, Jᵀ(hK)J only (see class doc).
		void applyDJT(const sofa::core::MechanicalParams *mparams, sofa::core::MultiVecDerivId inForce,
					  sofa::core::ConstMultiVecDerivId outForce) override;

		/// Number of elements.
		[[nodiscard]] std::size_t getNbElements() const { return d_elements.getValue().size(); }

	protected:
		LinearStrainElementMapping();
		~LinearStrainElementMapping() override = default;

		/// Kinematics of the last `apply`, reused by applyJ / applyJT (like most SOFA mappings).
		using KinematicsVector = std::vector<LSE::Kinematics, Eigen::aligned_allocator<LSE::Kinematics>>;
		KinematicsVector m_kinematics;
		/// Rotation matrix of every node at the last `apply` (world-frame <-> body-frame projector).
		std::vector<Matrix3> m_nodeRotations;
		/// Element lengths actually used (d_lengths, or chord lengths computed at init).
		std::vector<double> m_lengths;

		/// True when the element list, the lengths and the cached kinematics are consistent.
		bool isReady(std::size_t nbNodes, const char *caller) const;
	};

#if !defined(SOFA_COSSERAT_CPP_LinearStrainElementMapping)
	extern template class SOFA_COSSERAT_API LinearStrainElementMapping<
			sofa::defaulttype::Rigid3Types, sofa::defaulttype::Vec6Types, sofa::defaulttype::Vec6Types>;
#endif

} // namespace Cosserat::mapping
