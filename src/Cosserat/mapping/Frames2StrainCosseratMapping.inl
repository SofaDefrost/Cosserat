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

#include <Cosserat/mapping/Frames2StrainCosseratMapping.h>
#include <Cosserat/mapping/SofaLieGroupsUtils.h>
#include <sofa/core/Mapping.inl>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/gl/template.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/helper/logging/Message.h>
#include <sofa/helper/visual/DrawTool.h>
#include <sofa/type/Quat.h>

#include <cassert>
#include <string>

namespace Cosserat::mapping {

	using sofa::core::objectmodel::BaseContext;
	using sofa::helper::AdvancedTimer;
	using sofa::helper::WriteAccessor;
	using sofa::type::RGBAColor;
	using sofa::type::vector;
	using namespace sofa::component::cosserat::liegroups;


	// ──────────────────────────────────────────────────────────────────────────
	// Constructor
	// ──────────────────────────────────────────────────────────────────────────
	template<class TIn, class TOut>
	Frames2StrainCosseratMapping<TIn, TOut>::Frames2StrainCosseratMapping()
		: sofa::core::Mapping<TIn, TOut>()
		, CosseratBeamGeometry<TOut>()
		, d_deformationAxis(initData(&d_deformationAxis, (int)1, "deformationAxis",
									 "Axis to colour-map for the strain visualisation."))
		, d_max(initData(&d_max, (SReal)1.0e-2, "max", "Maximum strain value for the colormap."))
		, d_min(initData(&d_min, (SReal)0.0, "min", "Minimum strain value for the colormap."))
		, d_radius(initData(&d_radius, (SReal)0.05, "radius", "Cylinder radius for the beam visualisation."))
		, d_drawMapBeam(initData(&d_drawMapBeam, true, "nonColored",
								 "If false, the beam is rendered using the strain colormap."))
		, d_color(initData(&d_color,
						   sofa::type::RGBAColor(40 / 255.0, 104 / 255.0, 137 / 255.0, 0.8),
						   "color", "Default beam colour."))
		, d_index(initData(&d_index, "index",
						   "Optional list of section indices to highlight."))
	{
		// Register the Data members owned by the CosseratBeamGeometry mixin
		// (the standard initData() pattern doesn't apply here — see
		// CosseratGeometryMapping.inl for the rationale).
		this->d_curv_abs_section.setName("curv_abs_input");
		this->d_curv_abs_section.setHelp("Curvilinear abscissa of the input sections along the rod");
		this->addData(&this->d_curv_abs_section);

		this->d_curv_abs_frames.setName("curv_abs_output");
		this->d_curv_abs_frames.setHelp("Curvilinear abscissa of the output frames along the rod");
		this->addData(&this->d_curv_abs_frames);

		this->d_debug.setName("debug");
		this->d_debug.setHelp("Enable debug output");
		this->d_debug.setValue(false);
		this->addData(&this->d_debug);

		// Update callback for runtime changes to the discretisation
		this->addUpdateCallback("updateFrames",
			{&this->d_curv_abs_section, &this->d_curv_abs_frames, &this->d_debug},
			[this](const sofa::core::DataTracker &t) {
				SOFA_UNUSED(t);
				msg_info() << "Frames2StrainCosseratMapping updateFrames callback called";
				this->updateGeometryInfo();
				msg_info_when(this->d_debug.getValue()) << "====> Update Callback <====";
				return sofa::core::objectmodel::ComponentState::Valid;
			},
			{});
	}

	// ──────────────────────────────────────────────────────────────────────────
	// init() — bind mechanical states then build geometry
	// ──────────────────────────────────────────────────────────────────────────
	template<class TIn, class TOut>
	void Frames2StrainCosseratMapping<TIn, TOut>::init() {
		msg_info() << "Frames2StrainCosseratMapping init() called";

		m_input  = this->fromModel.get();
		m_output = this->toModel.get();

		if (m_input == nullptr) {
			msg_error() << "Input (frames) mechanical state not found";
			return;
		}
		if (m_output == nullptr) {
			msg_error() << "Output (strain) mechanical state not found";
			return;
		}

		// IMPORTANT ordering note:
		//   updateGeometryInfo() writes into m_frame_properties[i], so we must
		//   seed it with empty entries BEFORE calling it. The legacy code did
		//   this via the derived-class `initialization()` hook on the parent
		//   CosseratGeometryMapping; the new single-input flow seeds inline.
		const auto xfromData = m_input->read(sofa::core::vec_id::read_access::position);
		const auto &frames   = xfromData->getValue();
		this->m_frame_properties.clear();
		this->m_frame_properties.reserve(frames.size());
		for (size_t i = 0; i < frames.size(); ++i) {
			this->m_frame_properties.emplace_back();
		}

		// Build the per-frame topology (frame ↔ section mapping)
		this->updateGeometryInfo();

		// Initialise section properties with zero strain (real values will be
		// computed in apply()).
		typename Out::VecCoord zero_strain;
		const auto &curv_abs_section = this->d_curv_abs_section.getValue();
		if (curv_abs_section.size() >= 2) {
			zero_strain.resize(curv_abs_section.size() - 1);
		}
		this->initializeSectionProperties(zero_strain);

		// Initialise colourmap for the draw() output.
		// (SOFA deleted setColorScheme/reinit; build the palette via the
		//  string constructor instead.)
		m_colorMap = sofa::helper::ColorMap(256, "Blue to Red");

		// Parent init
		Inherit::init();
	}

	// ──────────────────────────────────────────────────────────────────────────
	// apply() — ξ_k = log(g_a⁻¹ · g_b) / h_k
	// ──────────────────────────────────────────────────────────────────────────
	template<class TIn, class TOut>
	void Frames2StrainCosseratMapping<TIn, TOut>::apply(
			const sofa::core::MechanicalParams * /*mparams*/,
			sofa::DataVecCoord_t<Out> &dataVecOut,
			const sofa::DataVecCoord_t<In> &dataVecIn) {

		msg_info_when(this->d_debug.getValue()) << " ########## Apply Function ########";

		if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
			return;

		const sofa::VecCoord_t<In> &frames = dataVecIn.getValue();

		const auto nbSections = this->m_section_properties.size() - 1;
		sofa::VecCoord_t<Out> &strains = *dataVecOut.beginEdit();
		strains.resize(nbSections);
		for (auto &s : strains)
			s.clear();

		// Convert each frame to SE3 (single conversion site via helper)
		std::vector<SE3Types> g_frames(frames.size());
		for (unsigned int i = 0; i < frames.size(); ++i) {
			g_frames[i] = rigidCoordToSE3(frames[i]);
		}

		const auto &curv_abs_sections = this->d_curv_abs_section.getValue();
		const auto &curv_abs_frames   = this->d_curv_abs_frames.getValue();

		for (unsigned int i = 0; i < nbSections; ++i) {
			double section_start = curv_abs_sections[i];
			double section_end   = curv_abs_sections[i + 1];

			int left_frame_idx = -1, right_frame_idx = -1;
			for (size_t j = 0; j < curv_abs_frames.size(); ++j) {
				if (std::abs(curv_abs_frames[j] - section_start) < 1e-12)
					left_frame_idx = j;
				if (std::abs(curv_abs_frames[j] - section_end) < 1e-12)
					right_frame_idx = j;
			}

			if (left_frame_idx >= 0 && right_frame_idx >= 0) {
				double dx = section_end - section_start;
				SE3Types g_rel = g_frames[left_frame_idx].computeInverse() * g_frames[right_frame_idx];
				TangentVector xi = g_rel.computeLog() / dx;
				for (int j = 0; j < 6; ++j)
					strains[i][j] = xi[j];
			} else {
				msg_warning() << "Could not find boundary frames for section " << i;
			}
		}

		dataVecOut.endEdit();
	}

	// ──────────────────────────────────────────────────────────────────────────
	// applyJ — ξ̇_k = J₁·η_a + J₂·η_b
	// ──────────────────────────────────────────────────────────────────────────
	template<class TIn, class TOut>
	void Frames2StrainCosseratMapping<TIn, TOut>::applyJ(
			const sofa::core::MechanicalParams * /*mparams*/,
			sofa::DataVecDeriv_t<Out> &dataVecOut,
			const sofa::DataVecDeriv_t<In> &dataVecIn) {

		if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
			return;

		msg_info_when(this->d_debug.getValue())
			<< " ########## Frames2StrainCosseratMapping ApplyJ Function ########";

		const sofa::VecDeriv_t<In> &frame_vel = dataVecIn.getValue();
		sofa::VecDeriv_t<Out> &strain_vel     = *dataVecOut.beginEdit();

		// Current strain values are stored in the output state (we wrote them in apply())
		const sofa::VecCoord_t<Out> &strain =
				m_output->read(sofa::core::vec_id::read_access::position)->getValue();

		// Frame coords (= the input "position" — they're absolute SE3 poses)
		const sofa::VecCoord_t<In> framePositions =
				m_input->read(sofa::core::vec_id::read_access::position)->getValue();

		const auto section_count = this->d_curv_abs_section.getValue().size() - 1;
		strain_vel.resize(section_count);
		for (auto &vel : strain_vel)
			vel.clear();

		// Convert each frame to SE3 once
		std::vector<SE3Types> g_frames(framePositions.size());
		for (unsigned int i = 0; i < framePositions.size(); ++i) {
			g_frames[i] = rigidCoordToSE3(framePositions[i]);
		}

		for (size_t i = 0; i < section_count; ++i) {
			const auto &section = this->m_section_properties[i + 1];
			double dx = section.getLength();

			TangentVector strain_i = TangentVector::Zero();
			for (int j = 0; j < 6; ++j)
				strain_i[j] = strain[i][j];

			TangentVector Omega_i  = dx * strain_i;
			// dexp_inv = Jr⁻¹(Ω) = Jl⁻¹(-Ω): the right-trivialized inverse Jacobian.
			AdjointMatrix dexp_inv = SE3Types::computeRightJacobianInverse(Omega_i);
			AdjointMatrix J2       = (1. / dx) * dexp_inv;

			SE3Types g_inv       = SE3Types::computeExp(-Omega_i); // = exp(-Ω)
			AdjointMatrix Adg_inv = g_inv.computeAdjoint();         // = Ad_{exp(-Ω)}
			AdjointMatrix J1     = -(1. / dx) * dexp_inv * Adg_inv;

			// Project frame velocities (global → local body twist).
			// P(R) is anti-diagonal ⇒ P(R)ᵀ = P(R⁻¹), so we transpose the
			// projector of the (non-inverted) frame instead of inverting it.
			const SE3Types &ga            = g_frames[i];
			AdjointMatrix a_projector     = ga.buildProjectionMatrix(ga.rotation().matrix());
			TangentVector vela_global(frame_vel[i][0], frame_vel[i][1], frame_vel[i][2],
									  frame_vel[i][3], frame_vel[i][4], frame_vel[i][5]);
			TangentVector eta_a           = a_projector.transpose() * vela_global;

			const SE3Types &gb            = g_frames[i + 1];
			AdjointMatrix b_projector     = gb.buildProjectionMatrix(gb.rotation().matrix());
			TangentVector velb_global(frame_vel[i + 1][0], frame_vel[i + 1][1], frame_vel[i + 1][2],
									  frame_vel[i + 1][3], frame_vel[i + 1][4], frame_vel[i + 1][5]);
			TangentVector eta_b           = b_projector.transpose() * velb_global;

			TangentVector output_vel      = J1 * eta_a + J2 * eta_b;

			for (int k = 0; k < 6; ++k)
				strain_vel[i][k] = output_vel[k];
		}

		dataVecOut.endEdit();
	}

	// ──────────────────────────────────────────────────────────────────────────
	// applyJT — propagate strain forces back to frame forces
	// ──────────────────────────────────────────────────────────────────────────
	template<class TIn, class TOut>
	void Frames2StrainCosseratMapping<TIn, TOut>::applyJT(
			const sofa::core::MechanicalParams * /*mparams*/,
			sofa::DataVecDeriv_t<In> &dataVecOut,
			const sofa::DataVecDeriv_t<Out> &dataVecIn) {

		if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
			return;

		msg_info_when(this->d_debug.getValue())
			<< " ########## Frames2StrainCosseratMapping ApplyJT Force Function ########";

		const sofa::VecDeriv_t<Out> &strainForces = dataVecIn.getValue();
		sofa::VecDeriv_t<In> &frameForces         = *dataVecOut.beginEdit();

		const sofa::VecCoord_t<Out> &strainState =
				m_output->read(sofa::core::vec_id::read_access::position)->getValue();

		const sofa::VecCoord_t<In> framePositions =
				m_input->read(sofa::core::vec_id::read_access::position)->getValue();

		frameForces.resize(framePositions.size());

		// Convert each frame to SE3 once
		std::vector<SE3Types> g_frames(framePositions.size());
		for (unsigned int i = 0; i < framePositions.size(); ++i) {
			g_frames[i] = rigidCoordToSE3(framePositions[i]);
		}

		const auto section_count = this->d_curv_abs_section.getValue().size() - 1;

		for (unsigned int i = 0; i < section_count; ++i) {
			const auto &section = this->m_section_properties[i + 1];
			double dx = section.getLength();

			TangentVector strain_i = TangentVector::Zero();
			TangentVector lambda   = TangentVector::Zero();
			for (int j = 0; j < 6; ++j) {
				strain_i[j] = strainState[i][j];
				lambda[j]   = strainForces[i][j];
			}

			TangentVector Omega_i  = dx * strain_i;
			// dexp_inv = Jr⁻¹(Ω) = Jl⁻¹(-Ω): the right-trivialized inverse Jacobian.
			AdjointMatrix dexp_inv = SE3Types::computeRightJacobianInverse(Omega_i);
			AdjointMatrix J2       = (1. / dx) * dexp_inv;

			SE3Types g_inv        = SE3Types::computeExp(-Omega_i);
			AdjointMatrix Adg_inv = g_inv.computeAdjoint();
			AdjointMatrix J1      = -(1. / dx) * dexp_inv * Adg_inv;

			// NB: the dx factor is already baked into J1/J2 (via -1/dx and 1/dx)
			TangentVector fa_local = J1.transpose() * lambda;
			TangentVector fb_local = J2.transpose() * lambda;

			// Project (local → global). P(R) is orthogonal ⇒ P⁻ᵀ = P, so the
			// projector applies directly (no transpose/inverse needed).
			const SE3Types &ga          = g_frames[i];
			AdjointMatrix a_projector   = ga.buildProjectionMatrix(ga.rotation().matrix());
			TangentVector fa_global     = a_projector * fa_local;

			const SE3Types &gb          = g_frames[i + 1];
			AdjointMatrix b_projector   = gb.buildProjectionMatrix(gb.rotation().matrix());
			TangentVector fb_global     = b_projector * fb_local;

			for (int k = 0; k < 6; ++k) {
				frameForces[i][k]     += fa_global[k];
				frameForces[i + 1][k] += fb_global[k];
			}
		}

		dataVecOut.endEdit();
	}

	// ──────────────────────────────────────────────────────────────────────────
	// applyJT — constraint version
	// ──────────────────────────────────────────────────────────────────────────
	template<class TIn, class TOut>
	void Frames2StrainCosseratMapping<TIn, TOut>::applyJT(
			const sofa::core::ConstraintParams * /*cparams*/,
			sofa::DataMatrixDeriv_t<In> &dataMatOut,
			const sofa::DataMatrixDeriv_t<Out> &dataMatIn) {

		if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
			return;

		msg_info_when(this->d_debug.getValue())
			<< " ########## Frames2StrainCosseratMapping ApplyJT Constraint Function ########";

		sofa::MatrixDeriv_t<In> &out         = *dataMatOut.beginEdit();
		const sofa::MatrixDeriv_t<Out> &in   = dataMatIn.getValue();

		const sofa::VecCoord_t<In> framePositions =
				m_input->read(sofa::core::vec_id::read_access::position)->getValue();

		const sofa::VecCoord_t<Out> &strainState =
				m_output->read(sofa::core::vec_id::read_access::position)->getValue();

		// Convert each frame to SE3 once
		std::vector<SE3Types> g_frames(framePositions.size());
		for (unsigned int i = 0; i < framePositions.size(); ++i) {
			g_frames[i] = rigidCoordToSE3(framePositions[i]);
		}

		// Process constraints
		for (auto rowIt = in.begin(); rowIt != in.end(); ++rowIt) {
			auto colIt = rowIt.begin();
			if (colIt == rowIt.end())
				continue;

			typename sofa::MatrixDeriv_t<In>::RowIterator o = out.writeLine(rowIt.index());

			while (colIt != rowIt.end()) {
				int strainIndex = colIt.index();
				const auto &section = this->m_section_properties[strainIndex + 1];

				TangentVector strain          = TangentVector::Zero();
				TangentVector constraintValue = TangentVector::Zero();

				const sofa::Deriv_t<Out> val = colIt.val();
				for (unsigned int j = 0; j < 6; ++j) {
					strain[j]          = strainState[strainIndex][j];
					constraintValue[j] = val[j];
				}

				double dx = section.getLength();
				TangentVector Omega    = dx * strain;
				// dexp_inv = Jr⁻¹(Ω) = Jl⁻¹(-Ω): the right-trivialized inverse Jacobian.
				AdjointMatrix dexp_inv = SE3Types::computeRightJacobianInverse(Omega);
				AdjointMatrix J2       = (1. / dx) * dexp_inv;

				SE3Types g_inv        = SE3Types::computeExp(-Omega);
				AdjointMatrix Adg_inv = g_inv.computeAdjoint();
				AdjointMatrix J1      = -(1. / dx) * dexp_inv * Adg_inv;

				TangentVector fa_local = J1.transpose() * constraintValue;
				TangentVector fb_local = J2.transpose() * constraintValue;

				// Project (local → global). P(R) orthogonal ⇒ P⁻ᵀ = P.
				const SE3Types &ga          = g_frames[strainIndex];
				AdjointMatrix a_projector   = ga.buildProjectionMatrix(ga.rotation().matrix());
				TangentVector fa_global     = a_projector * fa_local;

				const SE3Types &gb          = g_frames[strainIndex + 1];
				AdjointMatrix b_projector   = gb.buildProjectionMatrix(gb.rotation().matrix());
				TangentVector fb_global     = b_projector * fb_local;

				sofa::type::Vec<6, double> fa_vec, fb_vec;
				for (int k = 0; k < 6; ++k) {
					fa_vec[k] = fa_global[k];
					fb_vec[k] = fb_global[k];
				}

				o.addCol(strainIndex,     fa_vec);
				o.addCol(strainIndex + 1, fb_vec);

				colIt++;
			}
		}

		dataMatOut.endEdit();
	}

	// ──────────────────────────────────────────────────────────────────────────
	// applyDJT — geometric stiffness (not yet implemented)
	// ──────────────────────────────────────────────────────────────────────────
	template<class TIn, class TOut>
	void Frames2StrainCosseratMapping<TIn, TOut>::applyDJT(
			const sofa::core::MechanicalParams * /*mparams*/,
			sofa::core::MultiVecDerivId          /*inForce*/,
			sofa::core::ConstMultiVecDerivId     /*outForce*/) {
		static thread_local bool s_warned = false;
		if (!s_warned) {
			s_warned = true;
			msg_warning() << "applyDJT() is not implemented for Frames2StrainCosseratMapping. "
							 "Geometric stiffness is missing — implicit solvers may converge poorly "
							 "on large deformations.";
		}
	}

	// ──────────────────────────────────────────────────────────────────────────
	// computeBBox / draw
	// ──────────────────────────────────────────────────────────────────────────
	template<class TIn, class TOut>
	void Frames2StrainCosseratMapping<TIn, TOut>::computeBBox(
			const sofa::core::ExecParams *params, bool onlyVisible) {
		Inherit::computeBBox(params, onlyVisible);
	}

	template<class TIn, class TOut>
	void Frames2StrainCosseratMapping<TIn, TOut>::draw(const sofa::core::visual::VisualParams *vparams) {
		if (!vparams->displayFlags().getShowMechanicalMappings())
			return;

		const auto stateLifeCycle = vparams->drawTool()->makeStateLifeCycle();

		// Strain output (for colourmap)
		const auto strainData = m_output->read(sofa::core::vec_id::read_access::position);
		const auto xPos       = strainData->getValue();

		// Frames input (positions to draw the rod)
		const auto framesData = m_input->read(sofa::core::vec_id::read_access::position);
		const auto xData      = framesData->getValue();

		vector<sofa::type::Vec3> positions;
		positions.reserve(xData.size());
		for (unsigned int i = 0; i < xData.size(); ++i)
			positions.push_back(xData[i].getCenter());

		RGBAColor drawColor = d_color.getValue();
		for (unsigned int i = 0; i + 1 < positions.size(); ++i)
			vparams->drawTool()->drawCylinder(positions[i], positions[i + 1], d_radius.getValue(), drawColor);

		SReal min = d_min.getValue();
		SReal max = d_max.getValue();
		sofa::helper::ColorMap::evaluator<SReal> _eval = m_colorMap.getEvaluator(min, max);

		glLineWidth(d_radius.getValue());
		glBegin(GL_LINES);
		if (d_drawMapBeam.getValue()) {
			RGBAColor colorL = drawColor;
			glColor4f(colorL[0], colorL[1], colorL[2], colorL[3]);
			for (unsigned int i = 0; i + 1 < positions.size(); ++i) {
				vparams->drawTool()->drawLine(positions[i], positions[i + 1], colorL);
			}
		} else {
			for (unsigned int i = 0; i + 1 < positions.size(); ++i) {
				int j = (i < this->m_frame_to_section_indices.size())
								? static_cast<int>(this->m_frame_to_section_indices[i]) - 1
								: 0;
				if (j >= 0 && j < static_cast<int>(xPos.size())) {
					RGBAColor color = _eval(xPos[j][d_deformationAxis.getValue()]);
					vparams->drawTool()->drawLine(positions[i], positions[i + 1], color);
				}
			}
		}
		glLineWidth(1);
		glEnd();
	}

} // namespace Cosserat::mapping
