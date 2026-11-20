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

#include <Cosserat/mapping/LinearStrainElementMapping.h>
#include <Cosserat/mapping/SofaLieGroupsUtils.h>
#include <sofa/core/Multi2Mapping.inl>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/logging/Message.h>
#include <sofa/helper/visual/DrawTool.h>

#include <cmath>

namespace Cosserat::mapping {

	namespace lse_detail {
		using Vec6 = Eigen::Matrix<double, 6, 1>;
		using Matrix3 = Eigen::Matrix3d;

		/// SOFA world-frame rigid velocity [v ; ω] -> body twist [ω_b ; v_b] = P(R)ᵀ [v ; ω],
		/// with P(R) = [[0 R],[R 0]] (same projector as the other Cosserat mappings).
		template<class Deriv>
		inline Vec6 worldToBodyTwist(const Matrix3 &R, const Deriv &d) {
			const Eigen::Vector3d v(d[0], d[1], d[2]);
			const Eigen::Vector3d w(d[3], d[4], d[5]);
			Vec6 twist;
			twist.template head<3>() = R.transpose() * w;
			twist.template tail<3>() = R.transpose() * v;
			return twist;
		}

		/// Body-frame generalised force [a ; b] -> SOFA world-frame [R b ; R a] = P(R) [a ; b].
		inline Vec6 bodyToWorldForce(const Matrix3 &R, const Vec6 &x) {
			Vec6 f;
			f.template head<3>() = R * x.template tail<3>();
			f.template tail<3>() = R * x.template head<3>();
			return f;
		}

		template<class T>
		inline Vec6 toEigen(const T &v) {
			Vec6 e;
			for (int k = 0; k < 6; ++k)
				e[k] = v[k];
			return e;
		}
	} // namespace lse_detail

	// ──────────────────────────────────────────────────────────────────────────
	// Constructor
	// ──────────────────────────────────────────────────────────────────────────
	template<class TIn1, class TIn2, class TOut>
	LinearStrainElementMapping<TIn1, TIn2, TOut>::LinearStrainElementMapping()
		: d_elements(initData(&d_elements, "elements",
							  "Node pairs (a, b) of each element, indices in the pose (first input) state."))
		, d_lengths(initData(&d_lengths, "lengths",
							 "Reference length of each element. If empty, the chord length of the initial "
							 "poses is used (give the arc length for initially curved elements)."))
		, d_color(initData(&d_color, sofa::type::RGBAColor(40 / 255.0, 104 / 255.0, 137 / 255.0, 0.8), "color",
						   "Colour used to draw the elements.")) {}

	// ──────────────────────────────────────────────────────────────────────────
	// init
	// ──────────────────────────────────────────────────────────────────────────
	template<class TIn1, class TIn2, class TOut>
	void LinearStrainElementMapping<TIn1, TIn2, TOut>::init() {
		using sofa::core::objectmodel::ComponentState;
		this->d_componentState.setValue(ComponentState::Invalid);

		if (this->fromModels1.empty() || this->fromModels2.empty() || this->toModels.empty()) {
			msg_error() << "Three states are required: nodal poses (input1), strain slopes (input2) "
						   "and mean strains (output).";
			return;
		}

		const auto &poses = this->fromModels1[0]->read(sofa::core::vec_id::read_access::position)->getValue();
		const auto &elements = d_elements.getValue();
		const std::size_t nbElements = elements.size();

		if (nbElements == 0) {
			msg_error() << "'elements' is empty: give the node pairs of the rod elements.";
			return;
		}
		for (std::size_t e = 0; e < nbElements; ++e) {
			const unsigned int a = elements[e][0], b = elements[e][1];
			if (a >= poses.size() || b >= poses.size() || a == b) {
				msg_error() << "Element " << e << " = (" << a << ", " << b << ") is invalid for " << poses.size()
							<< " nodal poses.";
				return;
			}
		}
		if (this->fromModels2[0]->getSize() != nbElements) {
			msg_error() << "The slope state must hold one Vec6 per element: " << nbElements << " expected, got "
						<< this->fromModels2[0]->getSize() << ".";
			return;
		}

		// Element lengths: given, or chord length of the initial poses.
		m_lengths.assign(nbElements, 0.0);
		const auto &lengths = d_lengths.getValue();
		if (!lengths.empty() && lengths.size() != nbElements) {
			msg_error() << "'lengths' must be empty or hold one value per element (" << nbElements << "), got "
						<< lengths.size() << ".";
			return;
		}
		for (std::size_t e = 0; e < nbElements; ++e) {
			const double h = lengths.empty() ? (poses[elements[e][1]].getCenter() - poses[elements[e][0]].getCenter()).norm()
											 : static_cast<double>(lengths[e]);
			if (!(h > 1e-12)) {
				msg_error() << "Element " << e << " has a non-positive length (" << h << ").";
				return;
			}
			m_lengths[e] = h;
		}

		this->toModels[0]->resize(nbElements);

		// Multi2Mapping::init() calls apply() and applyJ(): the component must already be usable.
		this->d_componentState.setValue(ComponentState::Valid);
		Inherit::init();
		this->d_componentState.setValue(ComponentState::Valid);
	}

	template<class TIn1, class TIn2, class TOut>
	bool LinearStrainElementMapping<TIn1, TIn2, TOut>::isReady(std::size_t nbNodes, const char *caller) const {
		if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
			return false;
		const std::size_t n = getNbElements();
		if (m_lengths.size() != n || m_kinematics.size() != n || m_nodeRotations.size() != nbNodes) {
			msg_error() << caller << ": the cached kinematics do not match the current state sizes "
						<< "(apply() must run first).";
			return false;
		}
		return true;
	}

	// ──────────────────────────────────────────────────────────────────────────
	// apply — ξ̄_e = A⁻¹ Log(g_a⁻¹ g_b)
	// ──────────────────────────────────────────────────────────────────────────
	template<class TIn1, class TIn2, class TOut>
	void LinearStrainElementMapping<TIn1, TIn2, TOut>::apply(
			const sofa::core::MechanicalParams * /*mparams*/,
			const sofa::type::vector<sofa::DataVecCoord_t<Out> *> &dataVecOutPos,
			const sofa::type::vector<const sofa::DataVecCoord_t<In1> *> &dataVecIn1Pos,
			const sofa::type::vector<const sofa::DataVecCoord_t<In2> *> &dataVecIn2Pos) {
		if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
			return;
		if (dataVecOutPos.empty() || dataVecIn1Pos.empty() || dataVecIn2Pos.empty())
			return;

		const auto &poses = dataVecIn1Pos[0]->getValue();
		const auto &slopes = dataVecIn2Pos[0]->getValue();
		const auto &elements = d_elements.getValue();
		const std::size_t n = elements.size();
		if (slopes.size() != n || m_lengths.size() != n)
			return;

		std::vector<SE3Lie> g(poses.size());
		m_nodeRotations.resize(poses.size());
		for (std::size_t i = 0; i < poses.size(); ++i) {
			g[i] = rigidCoordToSE3(poses[i]);
			m_nodeRotations[i] = g[i].rotation().matrix();
		}

		m_kinematics.resize(n);
		auto &strains = *dataVecOutPos[0]->beginEdit();
		strains.resize(n);
		for (std::size_t e = 0; e < n; ++e) {
			const Vec6 beta = lse_detail::toEigen(slopes[e]);
			m_kinematics[e] = LSE::evaluate(g[elements[e][0]], g[elements[e][1]], beta, m_lengths[e]);
			for (int k = 0; k < 6; ++k)
				strains[e][k] = m_kinematics[e].xiBar[k];
		}
		dataVecOutPos[0]->endEdit();
	}

	// ──────────────────────────────────────────────────────────────────────────
	// applyJ — ξ̄̇_e = J₁ η_a + J₂ η_b + J₃ β̇_e
	// ──────────────────────────────────────────────────────────────────────────
	template<class TIn1, class TIn2, class TOut>
	void LinearStrainElementMapping<TIn1, TIn2, TOut>::applyJ(
			const sofa::core::MechanicalParams * /*mparams*/,
			const sofa::type::vector<sofa::DataVecDeriv_t<Out> *> &dataVecOutVel,
			const sofa::type::vector<const sofa::DataVecDeriv_t<In1> *> &dataVecIn1Vel,
			const sofa::type::vector<const sofa::DataVecDeriv_t<In2> *> &dataVecIn2Vel) {
		if (dataVecOutVel.empty() || dataVecIn1Vel.empty() || dataVecIn2Vel.empty())
			return;

		const auto &v1 = dataVecIn1Vel[0]->getValue();
		const auto &v2 = dataVecIn2Vel[0]->getValue();
		if (!isReady(v1.size(), "applyJ"))
			return;

		const auto &elements = d_elements.getValue();
		const std::size_t n = elements.size();
		if (v2.size() != n)
			return;

		auto &out = *dataVecOutVel[0]->beginEdit();
		out.resize(n);
		for (std::size_t e = 0; e < n; ++e) {
			const unsigned int a = elements[e][0], b = elements[e][1];
			const LSE::Kinematics &k = m_kinematics[e];
			const Vec6 etaA = lse_detail::worldToBodyTwist(m_nodeRotations[a], v1[a]);
			const Vec6 etaB = lse_detail::worldToBodyTwist(m_nodeRotations[b], v1[b]);
			const Vec6 betaDot = lse_detail::toEigen(v2[e]);
			const Vec6 xiDot = k.J1 * etaA + k.J2 * etaB + k.J3 * betaDot;
			for (int j = 0; j < 6; ++j)
				out[e][j] = xiDot[j];
		}
		dataVecOutVel[0]->endEdit();
	}

	// ──────────────────────────────────────────────────────────────────────────
	// applyJT (forces) — f_a += P J₁ᵀ λ, f_b += P J₂ᵀ λ, f_β += J₃ᵀ λ
	// ──────────────────────────────────────────────────────────────────────────
	template<class TIn1, class TIn2, class TOut>
	void LinearStrainElementMapping<TIn1, TIn2, TOut>::applyJT(
			const sofa::core::MechanicalParams * /*mparams*/,
			const sofa::type::vector<sofa::DataVecDeriv_t<In1> *> &dataVecOut1Force,
			const sofa::type::vector<sofa::DataVecDeriv_t<In2> *> &dataVecOut2Force,
			const sofa::type::vector<const sofa::DataVecDeriv_t<Out> *> &dataVecInForce) {
		if (dataVecOut1Force.empty() || dataVecOut2Force.empty() || dataVecInForce.empty())
			return;

		auto &f1 = *dataVecOut1Force[0]->beginEdit();
		auto &f2 = *dataVecOut2Force[0]->beginEdit();
		const auto &lambdas = dataVecInForce[0]->getValue();

		const auto &elements = d_elements.getValue();
		const std::size_t n = elements.size();
		if (isReady(f1.size(), "applyJT") && f2.size() == n && lambdas.size() == n) {
			for (std::size_t e = 0; e < n; ++e) {
				const unsigned int a = elements[e][0], b = elements[e][1];
				const LSE::Kinematics &k = m_kinematics[e];
				const Vec6 lambda = lse_detail::toEigen(lambdas[e]);
				const Vec6 fa = lse_detail::bodyToWorldForce(m_nodeRotations[a], k.J1.transpose() * lambda);
				const Vec6 fb = lse_detail::bodyToWorldForce(m_nodeRotations[b], k.J2.transpose() * lambda);
				const Vec6 fBeta = k.J3.transpose() * lambda;
				for (int j = 0; j < 6; ++j) {
					f1[a][j] += fa[j];
					f1[b][j] += fb[j];
					f2[e][j] += fBeta[j];
				}
			}
		}
		dataVecOut1Force[0]->endEdit();
		dataVecOut2Force[0]->endEdit();
	}

	// ──────────────────────────────────────────────────────────────────────────
	// applyJT (constraints)
	// ──────────────────────────────────────────────────────────────────────────
	template<class TIn1, class TIn2, class TOut>
	void LinearStrainElementMapping<TIn1, TIn2, TOut>::applyJT(
			const sofa::core::ConstraintParams * /*cparams*/,
			const sofa::type::vector<sofa::DataMatrixDeriv_t<In1> *> &dataMatOut1Const,
			const sofa::type::vector<sofa::DataMatrixDeriv_t<In2> *> &dataMatOut2Const,
			const sofa::type::vector<const sofa::DataMatrixDeriv_t<Out> *> &dataMatInConst) {
		if (dataMatOut1Const.empty() || dataMatOut2Const.empty() || dataMatInConst.empty())
			return;
		if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
			return;

		sofa::MatrixDeriv_t<In1> &out1 = *dataMatOut1Const[0]->beginEdit();
		sofa::MatrixDeriv_t<In2> &out2 = *dataMatOut2Const[0]->beginEdit();
		const sofa::MatrixDeriv_t<Out> &in = dataMatInConst[0]->getValue();

		const auto &elements = d_elements.getValue();
		if (isReady(m_nodeRotations.size(), "applyJT(constraint)")) {
			for (auto rowIt = in.begin(); rowIt != in.end(); ++rowIt) {
				auto colIt = rowIt.begin();
				if (colIt == rowIt.end())
					continue;

				typename sofa::MatrixDeriv_t<In1>::RowIterator o1 = out1.writeLine(rowIt.index());
				typename sofa::MatrixDeriv_t<In2>::RowIterator o2 = out2.writeLine(rowIt.index());

				for (; colIt != rowIt.end(); ++colIt) {
					const std::size_t e = colIt.index();
					if (e >= elements.size())
						continue;
					const unsigned int a = elements[e][0], b = elements[e][1];
					const LSE::Kinematics &k = m_kinematics[e];
					const Vec6 lambda = lse_detail::toEigen(colIt.val());

					const Vec6 fa = lse_detail::bodyToWorldForce(m_nodeRotations[a], k.J1.transpose() * lambda);
					const Vec6 fb = lse_detail::bodyToWorldForce(m_nodeRotations[b], k.J2.transpose() * lambda);
					const Vec6 fBeta = k.J3.transpose() * lambda;

					sofa::type::Vec<6, double> faVec, fbVec, fBetaVec;
					for (int j = 0; j < 6; ++j) {
						faVec[j] = fa[j];
						fbVec[j] = fb[j];
						fBetaVec[j] = fBeta[j];
					}
					o1.addCol(a, faVec);
					o1.addCol(b, fbVec);
					o2.addCol(e, fBetaVec);
				}
			}
		}
		dataMatOut1Const[0]->endEdit();
		dataMatOut2Const[0]->endEdit();
	}

	// ──────────────────────────────────────────────────────────────────────────
	// applyDJT — no-op by design (Gauss–Newton tangent, Eq. 30 of the reference)
	// ──────────────────────────────────────────────────────────────────────────
	template<class TIn1, class TIn2, class TOut>
	void LinearStrainElementMapping<TIn1, TIn2, TOut>::applyDJT(const sofa::core::MechanicalParams * /*mparams*/,
																sofa::core::MultiVecDerivId /*inForce*/,
																sofa::core::ConstMultiVecDerivId /*outForce*/) {
		// The reference formulation keeps the Jacobians J₁, J₂, J₃ constant in the tangent
		// (Gauss–Newton). The mapping term ∂(Jᵀf)/∂q is proportional to the strain residual and
		// vanishes at the rest strain; it is the missing exact-Hessian part the paper mentions.
	}

	// ──────────────────────────────────────────────────────────────────────────
	// draw
	// ──────────────────────────────────────────────────────────────────────────
	template<class TIn1, class TIn2, class TOut>
	void LinearStrainElementMapping<TIn1, TIn2, TOut>::draw(const sofa::core::visual::VisualParams *vparams) {
		if (!vparams->displayFlags().getShowMechanicalMappings() || this->fromModels1.empty())
			return;
		if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
			return;

		const auto stateLifeCycle = vparams->drawTool()->makeStateLifeCycle();
		const auto &poses = this->fromModels1[0]->read(sofa::core::vec_id::read_access::position)->getValue();
		const sofa::type::RGBAColor color = d_color.getValue();
		for (const auto &edge : d_elements.getValue()) {
			if (edge[0] < poses.size() && edge[1] < poses.size())
				vparams->drawTool()->drawLine(poses[edge[0]].getCenter(), poses[edge[1]].getCenter(), color);
		}
	}

} // namespace Cosserat::mapping
