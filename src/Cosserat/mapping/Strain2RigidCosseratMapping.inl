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

#include <Cosserat/mapping/Strain2RigidCosseratMapping.h>
#include <sofa/core/Multi2Mapping.inl>
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

	template<class TIn1, class TIn2, class TOut>
	Strain2RigidCosseratMapping<TIn1, TIn2, TOut>::Strain2RigidCosseratMapping() :
		Inherit(), d_deformationAxis(initData(&d_deformationAxis, (int) 1, "deformationAxis",
											  "the axis in which we want to show the deformation.\n")),
		d_max(initData(&d_max, (SReal) 1.0e-2, "max", "the maximum of the deformation.\n")),
		d_min(initData(&d_min, (SReal) 0.0, "min", "the minimum of the deformation.\n")),
		d_radius(initData(&d_radius, (SReal) 0.05, "radius", "the radius for beam visualization.\n")),
		d_drawMapBeam(initData(&d_drawMapBeam, true, "nonColored",
							   "if this parameter is false, you draw the beam with "
							   "color according to the force apply to each beam")),
		d_color(initData(&d_color, sofa::type::RGBAColor(40 / 255.0, 104 / 255.0, 137 / 255.0, 0.8), "color",
						 "The default beam color")),
		d_index(initData(&d_index, "index",
						 "if this parameter is false, you draw the beam with color "
						 "according to the force apply to each beam")),
		d_baseIndex(initData(&d_baseIndex, static_cast<unsigned int>(0), "baseIndex",
							 "This parameter defines the index of the rigid "
							 "base of Cosserat models, 0 by default this can"
							 "take another value if the rigid base is given "
							 "by another body.")) {

		// Register callback for updating frame transformations when geometry changes
		this->addUpdateCallback("updateFrames", {&d_curv_abs_section, &d_curv_abs_frames, &d_debug},
								[this](const sofa::core::DataTracker &t) {
									msg_info() << "Strain2RigidCosseratMapping updateFrames callback called";
									SOFA_UNUSED(t);
									this->updateGeometryInfo();
									msg_info_when(d_debug.getValue()) << "updateFrames callback triggered";
									const sofa::VecCoord_t<In1> &strain_state =
											m_strain_state->read(sofa::core::vec_id::read_access::position)->getValue();

									// This is also done in apply() So, no really need here !!!
									this->updateFrameTransformations(strain_state);
									return sofa::core::objectmodel::ComponentState::Valid;
								},
								{});
	}

	template<class TIn1, class TIn2, class TOut>
	void Strain2RigidCosseratMapping<TIn1, TIn2, TOut>::doBaseCosseratInit() {
		// Initialize colormap for visualization
		m_colorMap.setColorScheme("Blue to Red");
		m_colorMap.reinit();

		msg_info() << "Strain2RigidCosseratMapping initialized with liegroups SE(3) integration";
	}

	/*********************start debugging **************************/
	template<class TIn1, class TIn2, class TOut>
	void
	Strain2RigidCosseratMapping<TIn1, TIn2, TOut>::apply(const sofa::core::MechanicalParams * /* mparams */,
													  const vector<sofa::DataVecCoord_t<Out> *> &dataVecOutPos,
													  const vector<const sofa::DataVecCoord_t<In1> *> &dataVecIn1Pos,
													  const vector<const sofa::DataVecCoord_t<In2> *> &dataVecIn2Pos) {

		msg_info("Strain2RigidCosseratMapping") << "Strain2RigidCosseratMapping::apply called";

		if (dataVecOutPos.empty() || dataVecIn1Pos.empty() || dataVecIn2Pos.empty())
			return;

		// Check component state for validity
		if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
			return;

		msg_info_when(d_debug.getValue()) << "===== apply =====";

		// Get input data
		const sofa::VecCoord_t<In1> &strainState = dataVecIn1Pos[0]->getValue();
		const sofa::VecCoord_t<In2> &rigidBase = dataVecIn2Pos[0]->getValue();

		const auto frame_count = d_curv_abs_frames.getValue().size();
		sofa::VecCoord_t<Out> &output_frames = *dataVecOutPos[0]->beginEdit();
		output_frames.resize(frame_count);
		const auto baseIndex = d_baseIndex.getValue();

		// Debug output if enabled
		if (d_debug.getValue()) {
			displayMappingState("apply - start");
			displayStrainState(strainState, "apply - input");
			displayRigidState(rigidBase, "apply - input");
			displaySectionProperties("apply - before update");
		}

		// Recompute all local SE(3) exponentials from current strain state.
		// Fills m_section_properties[k].getTransformation() = expCosserat(ξ_k, L_k)
		// and m_frameProperties[i].getTransformation() = expCosserat(ξ, d_i) for each frame.
		updateFrameTransformations(strainState);

		// Get base frame transformation from rigid input
		const auto &base_rigid = rigidBase[baseIndex];
		Vector3 base_trans(base_rigid.getCenter()[0], base_rigid.getCenter()[1], base_rigid.getCenter()[2]);

		// Convert SOFA quaternion to Eigen quaternion (SOFA: x,y,z,w; Eigen: w,x,y,z)
		const auto &base_sofa_quat = base_rigid.getOrientation();
		Eigen::Quaternion<double> base_rot(base_sofa_quat[3], base_sofa_quat[0], base_sofa_quat[1], base_sofa_quat[2]);

		// Create SE3 transformation
		SE3Types base_frame(SE3Types::SO3Type(base_rot), base_trans);

		
		// Cache the printLog value out of the loop
		bool doPrintLog = this->f_printLog.getValue();

		// Pre-compute cumulative node poses in O(N) — avoids restarting from base_frame
		// for every output frame, which would cost O(F·N) overall.
		//
		//   node_poses[0] = base_frame                              (before any section)
		//   node_poses[k] = base_frame · g_0 · g_1 · … · g_{k-1}  (after k sections)
		//
		// where g_j = m_section_properties[j].getTransformation().
		// Allocate n_sections+1 so that related_beam_idx == n_sections is safe
		// (mirrors the original assert: related_beam_idx <= m_section_properties.size()).
		const size_t n_sections = m_section_properties.size();
		std::vector<SE3Types> node_poses(n_sections + 1);
		node_poses[0] = base_frame;
		for (size_t j = 1; j <= n_sections; ++j)
			node_poses[j] = node_poses[j - 1] * m_section_properties[j - 1].getTransformation();

		// Apply transformations to compute output frames — O(F) with pre-computed poses.
		for (unsigned int i = 0; i < frame_count; i++) {
			// Bounds checking
			assert(i < m_frameProperties.size() && "Frame index out of bounds");
			assert(i < output_frames.size() && "Output frames index out of bounds");

			// Each frame sits inside a section; related_beam_idx is the number of
			// full sections accumulated before reaching the frame's origin node.
			const auto related_beam_idx = m_frameProperties[i].get_related_beam_index_();
			assert(related_beam_idx <= n_sections && "Invalid beam index");

			// G_frame = node_poses[related_beam_idx] · g_frame(x)
			// where g_frame(x) = expCosserat(ξ, distance_to_nearest_node)
			const SE3Types current_frame =
				node_poses[related_beam_idx] * m_frameProperties[i].getTransformation();

			msg_info_when(d_debug.getValue()) << "Frame " << i << " = " << current_frame;
				
			// Save current rigid frame transformation into frame's properties
			//m_frameProperties[i].setTransformation(current_frame);

			// Convert SE3 to SOFA rigid coordinates
			const auto &translation = current_frame.translation();
			const auto &rotation = current_frame.rotation();

			// Convert rotation matrix to quaternion for SOFA
			Eigen::Quaternion<double> quat(rotation.matrix());
			sofa::type::Quat<SReal> sofa_quat(quat.x(), quat.y(), quat.z(), quat.w());
			sofa::type::Vec3 sofa_trans(translation[0], translation[1], translation[2]);

			output_frames[i] = sofa::Coord_t<Out>(sofa_trans, sofa_quat);

			if (doPrintLog) {
				msg_info() << "Frame " << i << " transformation applied";
			}
		}

		// Print distances between frames for debugging
		if (doPrintLog) {
			std::stringstream tmp;
			for (unsigned int i = 0; i < output_frames.size() - 1; i++) {
				sofa::type::Vec3 diff = output_frames[i + 1].getCenter() - output_frames[i].getCenter();
				tmp << "dist " << i << "  : " << diff.norm() << std::endl;
			}
			msg_info() << tmp.str();
		}

		// Debug output if enabled
		if (d_debug.getValue()) {
			displaySectionProperties("apply - after update");
			displayFrameProperties("apply - after update");
			displayOutputFrames(output_frames, "apply - output");
		}

		dataVecOutPos[0]->endEdit();
	}

	template<class TIn1, class TIn2, class TOut>
	void Strain2RigidCosseratMapping<TIn1, TIn2, TOut>::updateFrameTransformations(
			const sofa::type::vector<Coord1> &vec_of_strains) {
		auto nb_node = vec_of_strains.size();

		// Update node properties with current strain values
		// std::cout<<"[Loop]: Update node properties with current strain values"<<std::endl;
		for (size_t i = 0; i < nb_node; ++i) {
			// Extract strain components based on input type
			TangentVector strain = TangentVector::Zero();
			// No need anymore, this is already done in the constructor
			// Handle different strain input types (Vec3 or Vec6)
			if constexpr (std::is_same_v<Coord1, sofa::type::Vec3>) {
				// For Vec3 input, assume first 3 components are curvature
				strain.head<3>() = Vector3(vec_of_strains[i][0], vec_of_strains[i][1], vec_of_strains[i][2]);
			} else {
				// For Vec6 input, use all components
				for (int j = 0; j < 6 && j < vec_of_strains[i].size(); ++j) {
					strain[j] = vec_of_strains[i][j];
				}
			}
			// Update node info with strain values
			// i+1, since m_section_properties is 0-indexed
			m_section_properties[i + 1].setStrain(strain);

			// Compute SE(3) exponential for this section
			// Change input and give as input of the function m_section_properties[i]
			// SE3Types _gx = computeSE3Exponential(m_section_properties[i+1].getLength(),

			auto section_length = m_section_properties[i + 1].getLength();
			SE3Types _gx = SE3Types::expCosserat(strain, section_length);
			m_section_properties[i + 1].setTransformation(_gx);
		}

		// ── Rebuild CosseratBodyJacobian ─────────────────────────────────────────
		// After all section exponentials have been computed, cache the result in
		// m_bodyJacobian so that applyJ / applyJT can call applyForward() /
		// applyTranspose() without re-doing the section-level recurrence.
		//
		// Convention: body Jacobian section k  ↔  m_section_properties[k+1]
		//             (index 0 in section_properties is the fixed base, not a section)
		m_bodyJacobian.clear();
		m_bodyJacobian.reserve(static_cast<int>(nb_node));
		for (size_t i = 0; i < nb_node; ++i) {
			m_bodyJacobian.pushSection(
				m_section_properties[i + 1].getTransformation(),
				m_section_properties[i + 1].getTangAdjointMatrix()
			);
		}
		// ─────────────────────────────────────────────────────────────────────────

		// Update frame properties based on their position within sections
		for (size_t i = 0; i < m_frameProperties.size(); ++i) {
			if (i < m_indices_vectors.size()) {
				int sectionIndex = m_frameProperties[i].get_related_beam_index_();
				if (sectionIndex >= 0 && sectionIndex < static_cast<int>(vec_of_strains.size() + 1)) {
					// Compute frame transformation at its specific position
					SE3Types frame_gx = SE3Types::expCosserat(m_section_properties[sectionIndex].getStrainsVec(),
															  m_frameProperties[i].getDistanceToNearestBeamNode());
					m_frameProperties[i].setTransformation(frame_gx);
				}
			}
		}

	}


	template<class TIn1, class TIn2, class TOut>
	void
	Strain2RigidCosseratMapping<TIn1, TIn2, TOut>::applyJ(const sofa::core::MechanicalParams * /* mparams */,
													   const vector<sofa::DataVecDeriv_t<Out> *> &dataVecOutVel,
													   const vector<const sofa::DataVecDeriv_t<In1> *> &dataVecIn1Vel,
													   const vector<const sofa::DataVecDeriv_t<In2> *> &dataVecIn2Vel) {

		if (dataVecOutVel.empty() || dataVecIn1Vel.empty() || dataVecIn2Vel.empty())
			return;

		if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
			return;

		msg_info_when(d_debug.getValue()) << "===== applyJ =====";

		const sofa::VecDeriv_t<In1> &strain_vel = dataVecIn1Vel[0]->getValue();
		const sofa::VecDeriv_t<In2> &base_vel = dataVecIn2Vel[0]->getValue();
		sofa::VecDeriv_t<Out> &frame_vel = *dataVecOutVel[0]->beginEdit();

		const sofa::VecCoord_t<Out> &framePositions =
				this->m_frames->read(sofa::core::vec_id::read_access::position)->getValue();

		// Debug input velocities if enabled
		if (d_debug.getValue()) {
			displayVelocities(strain_vel, base_vel, frame_vel, "applyJ - input");
		}

		const auto base_index = d_baseIndex.getValue();
		const auto frame_count = d_curv_abs_frames.getValue().size();
		frame_vel.resize(frame_count);
		for (auto &vel : frame_vel){
			vel.clear();
		}

		// 1. Compute current tangent exponential SE3 matrices
		this->updateTangExpSE3();

		// 2. Compute the base velocity in SE(3) tangent space
		// 2.1 Convert base velocity to se(3) tangent vector
		TangentVector base_vel_local = TangentVector::Zero();
		for (auto u = 0; u < 6; u++)
			base_vel_local[u] = base_vel[base_index][u];

		// 2.2 Apply the local transform from SOFA's frame to Cosserat's frame

		auto frame0 = framePositions[0];		
		Vector3 trans0(frame0.getCenter()[0], frame0.getCenter()[1], frame0.getCenter()[2]);
		const auto &quat0 = frame0.getOrientation();
		Eigen::Quaternion<double> rot0(quat0[3], quat0[0], quat0[1], quat0[2]);
			
		SE3Types absoluteFrame0(SE3Types::SO3Type(rot0), trans0);
		SE3Types absoluteFrame0_inv = absoluteFrame0.inverse();
		//projection de la force globale dans le repère local
		AdjointMatrix base_projector = absoluteFrame0_inv.buildProjectionMatrix(absoluteFrame0_inv.rotation().matrix());


		// ── 3. Forward propagation via CosseratBodyJacobian ─────────────────────
		//
		// Replaces the manual section-level recurrence:
		//   η_k = Ad_{g_k⁻¹} · (η_{k-1} + J_local_k · ξ̇_k)
		//
		// m_bodyJacobian was built by updateFrameTransformations() earlier in
		// apply(), so section transforms and tangent-exp matrices are up-to-date.
		//
		// Strain DOF count is resolved at compile time (Vec3 → 3, Vec6 → 6).
		static constexpr int NStrainDOF_J = static_cast<int>(Coord1::total_size);

		// Build one TwistType per section from the input strain velocities.
		const int N_sections = m_bodyJacobian.size();
		std::vector<TwistType> strain_rates(static_cast<size_t>(N_sections), TwistType::Zero());
		for (int k = 0; k < N_sections; ++k) {
			TangentVector sv = TangentVector::Zero();
			if (k < static_cast<int>(strain_vel.size())) {
				for (int j = 0; j < NStrainDOF_J; ++j)
					sv[j] = strain_vel[k][j];
			}
			strain_rates[k] = TwistType(sv);
		}

		// Forward pass: returns N+1 twists [η_0, η_1, …, η_N]
		//   η_0 = base_projector · v_base   (root body twist)
		//   η_k = Ad_{g_k⁻¹} · (η_{k-1} + J_local_k · ξ̇_k)
		const TwistType base_twist(base_projector * base_vel_local);
		const auto body_twists = m_bodyJacobian.applyForward(strain_rates, base_twist);

		// Unpack into std::vector<TangentVector> for the frame-level step below.
		// body_twists[0] = root, body_twists[k] = after k sections.
		std::vector<TangentVector> node_velocities(body_twists.size());
		for (size_t i = 0; i < body_twists.size(); ++i)
			node_velocities[i] = body_twists[i].data();

		if (d_debug.getValue()) {
			for (size_t i = 0; i < node_velocities.size(); ++i)
				msg_info() << "Node velocity [" << i << "]: " << node_velocities[i].transpose();
		}
		// ─────────────────────────────────────────────────────────────────────────

		// 4. Compute velocity at each output frame
		for (size_t i = 0; i < frame_count; ++i) {
			const auto &frame = m_frameProperties[i];
			const auto &tang_adj = frame.getTangAdjointMatrix();

			// Get the section node index this frame belongs to (1-based → 0-based).
			// Clamp to a valid range to avoid out-of-bounds access.
			const int section_idx = [&]() {
				int idx = (i < m_indices_vectors.size()) ? static_cast<int>(m_indices_vectors[i]) - 1 : 0;
				if (idx < 0 || idx >= static_cast<int>(node_velocities.size())) idx = 0;
				return idx;
			}();

			// Extract frame strain velocity (same DOFs as the parent section).
			// NStrainDOF_J is resolved at compile time (Vec3 → 3, Vec6 → 6).
			TangentVector frame_strain_vel = TangentVector::Zero();
			if (section_idx < static_cast<int>(strain_vel.size())) {
				for (int j = 0; j < NStrainDOF_J; ++j)
					frame_strain_vel[j] = strain_vel[section_idx][j];
			}

			// Compute frame velocity: η_frame = Ad_{g_frame^{-1}} · (η_node + T_frame · ξ̇_frame)
			// See docs/mapping_explanation.md §3 — applyJ recurrence.
			const auto g_inv = frame.getInverseTransformation();
			const AdjointMatrix Ad_gm1 = g_inv.computeAdjoint();

			const TangentVector eta_frame =
					Ad_gm1 * (node_velocities[section_idx] + tang_adj * frame_strain_vel);

			// Project to output frame (convert from local to SOFA global frame)
			const auto &pos = framePositions[i]; 
			Vector3 translation(pos.getCenter()[0], pos.getCenter()[1], pos.getCenter()[2]);
			const auto &quat = pos.getOrientation();
			Eigen::Quaternion<double> rotation(quat[3], quat[0], quat[1], quat[2]);
				
			SE3Types absoluteFrame(SE3Types::SO3Type(rotation), translation);
			AdjointMatrix frame_projector = absoluteFrame.buildProjectionMatrix(absoluteFrame.rotation().matrix());

			
			TangentVector output_vel = frame_projector * eta_frame;

			// Convert to SOFA format
			for (int k = 0; k < 6; ++k) {
				frame_vel[i][k] = output_vel[k];
			}
			msg_info_when(d_debug.getValue()) << "Frame velocity [" << i << "]: " << output_vel.transpose();
		}

		// Debug output velocities if enabled
		if (d_debug.getValue()) {
			displayVelocities(strain_vel, base_vel, frame_vel, "applyJ - output");
		}

		dataVecOutVel[0]->endEdit();
	}

	template<class TIn1, class TIn2, class TOut>
	void Strain2RigidCosseratMapping<TIn1, TIn2, TOut>::applyJT(
			const sofa::core::MechanicalParams * /*mparams*/,
			const vector<sofa::DataVecDeriv_t<In1> *> &dataVecOut1Force,
			const vector<sofa::DataVecDeriv_t<In2> *> &dataVecOut2Force,
			const vector<const sofa::DataVecDeriv_t<Out> *> &dataVecInForce) {

		if (dataVecOut1Force.empty() || dataVecInForce.empty() || dataVecOut2Force.empty())
			return;

		if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
			return;

		msg_info_when(d_debug.getValue()) << "===== applyJT (VecDeriv) =====";



		const sofa::VecDeriv_t<Out> &inputForces = dataVecInForce[0]->getValue();
		sofa::VecDeriv_t<In1> &strainForces = *dataVecOut1Force[0]->beginEdit();
		sofa::VecDeriv_t<In2> &baseForces = *dataVecOut2Force[0]->beginEdit();
		const auto baseIndex = d_baseIndex.getValue();

		// Get current positions to compute transformations
		const sofa::VecCoord_t<Out> &framePositions =
				this->m_frames->read(sofa::core::vec_id::read_access::position)->getValue();
		
		// Initialize output forces to zero before accumulation.
		// Size = number of sections = m_bodyJacobian.size() (built in apply()).
		// resize() may not zero-initialise existing elements on repeated calls.
		strainForces.resize(static_cast<size_t>(m_bodyJacobian.size()));
		for (auto &f : strainForces) f.clear();

		// B^T : selects NStrainDOF active strain DOFs from the 6D body-twist.
		// Vec3Types → [I_3 | 0_3] (angular only, 3 DOF).
		// Vec6Types → I_6 (angular + linear, 6 DOF).
		// See docs/mapping_explanation.md §12.6
		static constexpr int NStrainDOF = static_cast<int>(Coord1::total_size);
		using MatB = Eigen::Matrix<double, NStrainDOF, 6>;
		MatB matB_trans = MatB::Zero();
		for (int k = 0; k < NStrainDOF; k++) matB_trans(k, k) = 1.0;

		// ── Phase 0: project each frame force from SOFA global frame to body frame ──
		const size_t tab_size = inputForces.size();
		const size_t sz       = m_indices_vectors.size();

		if (sz == 0 || tab_size == 0) {
			dataVecOut1Force[0]->endEdit();
			dataVecOut2Force[0]->endEdit();
			return;
		}

		vector<TangentVector> localForces;
		localForces.reserve(tab_size);
		for (size_t i = 0; i < tab_size; ++i) {
			TangentVector frameForce = TangentVector::Zero();
			for (unsigned j = 0; j < 6; j++) frameForce[j] = inputForces[i][j];

			const auto &pos = framePositions[i];
			Vector3 translation(pos.getCenter()[0], pos.getCenter()[1], pos.getCenter()[2]);
			const auto &quat = pos.getOrientation();
			Eigen::Quaternion<double> rotation(quat[3], quat[0], quat[1], quat[2]);
			SE3Types absoluteFrame(SE3Types::SO3Type(rotation), translation);
			AdjointMatrix P_trans = absoluteFrame.buildProjectionMatrix(absoluteFrame.rotation().matrix());
			localForces.push_back(P_trans.transpose() * frameForce);
		}

		// ── Phase 1: per-frame direct contributions ────────────────────────────────
		//
		// For each frame i in section s (1-based):
		//   w_body_i = Ad_{g_frame_i^{-T}} · w_local_i   (co-adjoint pull-back)
		//
		//   a) Direct frame strain force (J_frame^T term):
		//         strainForces[s-1] += B^T · J_frame_i^T · w_body_i
		//
		//   b) Aggregated node wrench:
		//         external_wrenches[s-1] += w_body_i
		//
		// external_wrenches[k] is then fed into CosseratBodyJacobian::applyTranspose()
		// which handles the section-to-section backward propagation.
		//
		// Duality proof (2-section example in docs/impro_mapping.md §A2).
		const int N_bj = m_bodyJacobian.size();  // number of sections
		std::vector<WrenchType> external_wrenches(static_cast<size_t>(N_bj + 1), WrenchType::Zero());

		for (size_t i = 0; i < sz; ++i) {
			const int s1  = m_indices_vectors[i];  // section index, 1-based
			const int nidx = s1 - 1;               // node index, 0-based

			const FrameInfo &frame = m_frameProperties[i];
			const TangentVector w_body = frame.getCoAdjoint() * localForces[i];

			// a) Direct strain force for this frame's section
			const Eigen::Matrix<double, NStrainDOF, 1> f =
				matB_trans * frame.getTangAdjointMatrix().transpose() * w_body;
			for (int j = 0; j < NStrainDOF; j++)
				strainForces[s1 - 1][j] += f[j];

			// b) Aggregate node wrench (bounds-checked)
			if (nidx >= 0 && nidx < N_bj + 1)
				external_wrenches[nidx] = WrenchType(external_wrenches[nidx].data() + w_body);
		}

		// ── Phase 2: section-to-section backward transport via CosseratBodyJacobian ──
		//
		// strain_bj[k] = J_section_k^T · Ad_{g_k^{-T}} · (sum of downstream wrenches)
		// base_wrench   = accumulated root wrench
		//
		// These section contributions are ADDED to strainForces — they are distinct
		// from the direct J_frame^T contributions computed in Phase 1.
		WrenchType base_wrench_out;
		const auto strain_bj = m_bodyJacobian.applyTranspose(external_wrenches, base_wrench_out);

		for (int k = 0; k < N_bj; ++k) {
			const Eigen::Matrix<double, NStrainDOF, 1> sf = matB_trans * strain_bj[k].data();
			for (int j = 0; j < NStrainDOF; j++)
				strainForces[k][j] += sf[j];
		}

		// ── Phase 3: project root wrench back to SOFA convention ──────────────────
		const auto &frame0 = framePositions[0];
		const Vector3 trans0(frame0.getCenter()[0], frame0.getCenter()[1], frame0.getCenter()[2]);
		const auto &quat0 = frame0.getOrientation();
		SE3Types absoluteFrame0(
			SE3Types::SO3Type(Eigen::Quaternion<double>(quat0[3], quat0[0], quat0[1], quat0[2])),
			trans0);
		const AdjointMatrix M = absoluteFrame0.buildProjectionMatrix(absoluteFrame0.rotation().matrix());
		const TangentVector base_force_sofa = M * base_wrench_out.data();
		for (int j = 0; j < 6; j++)
			baseForces[baseIndex][j] += base_force_sofa[j];

		if (d_debug.getValue()) {
			msg_info() << "applyJT | frames=" << tab_size
			           << " | base wrench: [" << base_wrench_out.data().transpose() << "]"
			           << " | base index: " << baseIndex;
		}

		dataVecOut1Force[0]->endEdit();
		dataVecOut2Force[0]->endEdit();

	}

	template<class TIn1, class TIn2, class TOut>
	void Strain2RigidCosseratMapping<TIn1, TIn2, TOut>::applyJT(
			const sofa::core::ConstraintParams * /*cparams*/,
			const vector<sofa::DataMatrixDeriv_t<In1> *> &dataMatOut1Const,
			const vector<sofa::DataMatrixDeriv_t<In2> *> &dataMatOut2Const,
			const vector<const sofa::DataMatrixDeriv_t<Out> *> &dataMatInConst) {


		if (dataMatOut1Const.empty() || dataMatOut2Const.empty() || dataMatInConst.empty())
			return;

		if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
			return;

		msg_info_when(d_debug.getValue()) << "===== applyJT (MatrixDeriv) =====";

		
		// Prepare input and output data matrices
		sofa::MatrixDeriv_t<In1> &out1 = *dataMatOut1Const[0]->beginEdit();
		sofa::MatrixDeriv_t<In2> &out2 = *dataMatOut2Const[0]->beginEdit();
		const sofa::MatrixDeriv_t<Out> &in = dataMatInConst[0]->getValue();

		// Get current frame positions (needed for projection matrices).
		const sofa::VecCoord_t<Out> &framePositions =
				this->m_frames->read(sofa::core::vec_id::read_access::position)->getValue();

		// B^T : selects NStrainDOF active strain DOFs from the 6D body-twist.
		// Vec3Types → [I_3 | 0_3] (3 DOF angular); Vec6Types → I_6 (6 DOF).
		// See docs/mapping_explanation.md §12.6
		static constexpr int NStrainDOF = static_cast<int>(Coord1::total_size);
		using MatB = Eigen::Matrix<double, NStrainDOF, 6>;
		MatB matB_trans = MatB::Zero();
		for (int k = 0; k < NStrainDOF; k++) matB_trans(k, k) = 1.0;

		vector<std::tuple<int, TangentVector>> NodesInvolved;
		vector<std::tuple<int, TangentVector>> NodesInvolvedCompressed;
		// Process constraints
		for (auto rowIt = in.begin(); rowIt != in.end(); ++rowIt) {
			auto colIt = rowIt.begin();
			if (colIt == rowIt.end())
				continue;

			typename sofa::MatrixDeriv_t<In1>::RowIterator o1 = out1.writeLine(rowIt.index());
			typename sofa::MatrixDeriv_t<In2>::RowIterator o2 = out2.writeLine(rowIt.index());

			NodesInvolved.clear();

			while (colIt != rowIt.end()) {
				int frameIndex = colIt.index();

				TangentVector constraintValue = TangentVector::Zero();

				// Convert constraint value to TangentVector
				const sofa::Deriv_t<Out> val = colIt.val();
				for (unsigned int j = 0; j < 6; ++j) {
					constraintValue[j] = val[j];
				}

				int sectionIndex = m_indices_vectors[frameIndex];

				auto &pos = framePositions[frameIndex];
				Vector3 translation(pos.getCenter()[0], pos.getCenter()[1], pos.getCenter()[2]);
				auto &quat = pos.getOrientation();
				Eigen::Quaternion<double> rotation(quat[3], quat[0], quat[1], quat[2]);
				SE3Types absoluteFrame(SE3Types::SO3Type(rotation), translation);
				
				AdjointMatrix P_trans = absoluteFrame.buildProjectionMatrix(absoluteFrame.rotation().matrix());
				
				const FrameInfo &frame = m_frameProperties[frameIndex];
			
				const AdjointMatrix coAdjoint = frame.getCoAdjoint();

				// Local force: ℓ_c^s = Ad_{g_s^{-T}} · P_s^T · δ_c^s
				// (constraint direction pulled back to the body frame of element b_s)
				const TangentVector localForce = coAdjoint * P_trans.transpose() * constraintValue;

				// Project into strain space: f = B^T · J_local^T · ℓ_c^s
				const Eigen::Matrix<double, 6, 6> temp = frame.getTangAdjointMatrix().transpose();
				const Eigen::Matrix<double, NStrainDOF, 1> f = matB_trans * temp * localForce;

				sofa::Deriv_t<In1> f_out;
				for (int k = 0; k < NStrainDOF; k++) f_out[k] = f[k];
				o1.addCol(sectionIndex - 1, f_out);

				NodesInvolved.emplace_back(sectionIndex, localForce);
				colIt++;
			}

			// Sort by decreasing node index so that same-node entries are adjacent.
			std::sort(
				NodesInvolved.begin(), NodesInvolved.end(),
				[](const std::tuple<int, TangentVector> &a, const std::tuple<int, TangentVector> &b) {
					return std::get<0>(a) > std::get<0>(b);
				});

			// Merge consecutive entries that share the same node index (group-by).
			// Replaces the previous fragile while+break pattern.
			NodesInvolvedCompressed.clear();
			{
				size_t n = 0;
				while (n < NodesInvolved.size()) {
					int currentNode = std::get<0>(NodesInvolved[n]);
					TangentVector cumulativeF = std::get<1>(NodesInvolved[n]);
					++n;
					while (n < NodesInvolved.size() && std::get<0>(NodesInvolved[n]) == currentNode) {
						cumulativeF += std::get<1>(NodesInvolved[n]);
						++n;
					}
					NodesInvolvedCompressed.emplace_back(currentNode, cumulativeF);
				}
			}

			// Pre-compute the base-frame projector once per constraint row.
			// (frame0 does not change during the backward loop below.)
			const auto &frame0    = framePositions[0];
			const Vector3 trans0(frame0.getCenter()[0], frame0.getCenter()[1], frame0.getCenter()[2]);
			const auto   &quat0   = frame0.getOrientation();
			const SE3Types absoluteFrame0(
				SE3Types::SO3Type(Eigen::Quaternion<double>(quat0[3], quat0[0], quat0[1], quat0[2])),
				trans0);
			const AdjointMatrix M_base =
				absoluteFrame0.buildProjectionMatrix(absoluteFrame0.rotation().matrix());
			const unsigned int baseIdx = d_baseIndex.getValue();

			// Backward co-adjoint propagation toward the root.
			// Implements §4 backward recurrence: f_{k-1} = Ad_{g_k^{-T}} · f_k
			for (const auto &[startNode, startForce] : NodesInvolvedCompressed) {
				int i = startNode;
				TangentVector CumulativeF = startForce;

				while (i > 0) {
					const SectionInfo &section = m_section_properties[i - 1];
					// Transport: f_{i-1} = Ad_{g_i^{-T}} · f_i
					CumulativeF = section.getCoAdjoint() * CumulativeF;
					// Project into strain space: q_{i-1} = B^T · J_{i-1}^T · f_{i-1}
					const Eigen::Matrix<double, 6, 6> tempSection = section.getTangAdjointMatrix().transpose();
					const Eigen::Matrix<double, NStrainDOF, 1> temp_f = matB_trans * tempSection * CumulativeF;

					if (i > 1) {
						sofa::Deriv_t<In1> temp_f_out;
						for (int k = 0; k < NStrainDOF; k++) temp_f_out[k] = temp_f[k];
						o1.addCol(i - 2, temp_f_out);
					}
					i--;
				}

				// Project accumulated wrench to base-frame SOFA convention.
				const TangentVector base_force = M_base * CumulativeF;
				sofa::type::Vec6 base_force_trans;
				for (int k = 0; k < 6; ++k) base_force_trans[k] = base_force[k];
				o2.addCol(baseIdx, base_force_trans);
			}
		
		}

		dataMatOut1Const[0]->endEdit();
		dataMatOut2Const[0]->endEdit();


	}


	template<class TIn1, class TIn2, class TOut>
	void Strain2RigidCosseratMapping<TIn1, TIn2, TOut>::draw(const sofa::core::visual::VisualParams *vparams) {
		if (!vparams->displayFlags().getShowMechanicalMappings())
			return;

		// draw cable
		typedef sofa::type::RGBAColor RGBAColor;

		const auto stateLifeCycle = vparams->drawTool()->makeStateLifeCycle();

		const sofa::DataVecCoord_t<Out> *xfromData = this->m_frames->read(sofa::core::vec_id::read_access::position);
		const sofa::VecCoord_t<Out> xData = xfromData->getValue();
		vector<sofa::type::Vec3> positions;
		vector<sofa::type::Quat<SReal>> Orientation;
		positions.clear();
		Orientation.clear();
		unsigned int sz = xData.size();
		for (unsigned int i = 0; i < sz; i++) {
			positions.push_back(xData[i].getCenter());
			Orientation.push_back(xData[i].getOrientation());
		}

		// Get access articulated
		const sofa::DataVecCoord_t<In1> *artiData = this->m_strain_state->read(sofa::core::vec_id::read_access::position);
		const sofa::VecCoord_t<In1> xPos = artiData->getValue();

		RGBAColor drawColor = d_color.getValue();
		// draw each segment of the beam as a cylinder.
		for (unsigned int i = 0; i < sz - 1; i++)
			vparams->drawTool()->drawCylinder(positions[i], positions[i + 1], d_radius.getValue(), drawColor);

		// Define color map
		SReal min = d_min.getValue();
		SReal max = d_max.getValue();
		sofa::helper::ColorMap::evaluator<SReal> _eval = m_colorMap.getEvaluator(min, max);

		glLineWidth(d_radius.getValue());
		glBegin(GL_LINES);
		if (d_drawMapBeam.getValue()) {
			sofa::type::RGBAColor _color = d_color.getValue();
			RGBAColor colorL = RGBAColor(_color[0], _color[1], _color[2], _color[3]);
			glColor4f(colorL[0], colorL[1], colorL[2], colorL[3]);
			for (unsigned int i = 0; i < sz - 1; i++) {
				vparams->drawTool()->drawLine(positions[i], positions[i + 1], colorL);
			}
		} else {
			int j = 0;
			vector<int> index = d_index.getValue();
			for (unsigned int i = 0; i < sz - 1; i++) {
				j = m_indices_vectors[i] - 1; // to get the articulation on which the frame is related to
				RGBAColor color = _eval(xPos[j][d_deformationAxis.getValue()]);
				vparams->drawTool()->drawLine(positions[i], positions[i + 1], color);
			}
		}
		glLineWidth(1);
		if (!vparams->displayFlags().getShowMappings())
			if (!d_debug.getValue())
				return;

		// Debug output if needed
		if (this->f_printLog.getValue()) {
			displayOutputFrames(xData, "draw - rendering frames");
		}

		glEnd();
	}




	template<class TIn1, class TIn2, class TOut>
	void Strain2RigidCosseratMapping<TIn1, TIn2, TOut>::computeBBox(const sofa::core::ExecParams *params,
																 bool onlyVisible) {
		// Compute bounding box for visualization
		// Implementation would calculate the extent of all frames
		Inherit::computeBBox(params, onlyVisible);
	}

	// Debug display functions implementation
	template<class TIn1, class TIn2, class TOut>
	void Strain2RigidCosseratMapping<TIn1, TIn2, TOut>::displayStrainState(const sofa::type::vector<Coord1> &strainState,
																		const std::string &context) const {

		std::cout << "\n=== STRAIN STATE DEBUG" << (context.empty() ? "" : " (" + context + ")") << " ===\n";
		std::cout << "Strain state size: " << strainState.size() << std::endl;

		for (size_t i = 0; i < strainState.size(); ++i) {
			std::cout << "  Strain[" << i << "]: ";

			if constexpr (std::is_same_v<Coord1, sofa::type::Vec3>) {
				std::cout << "[" << strainState[i][0] << ", " << strainState[i][1] << ", " << strainState[i][2] << "]";
			} else {
				std::cout << "[";
				for (int j = 0; j < strainState[i].size(); ++j) {
					std::cout << strainState[i][j];
					if (j < strainState[i].size() - 1)
						std::cout << ", ";
				}
				std::cout << "]";
			}
			std::cout << std::endl;
		}
		std::cout << "================================\n";
	}

	template<class TIn1, class TIn2, class TOut>
	void Strain2RigidCosseratMapping<TIn1, TIn2, TOut>::displayRigidState(
			const sofa::type::vector<sofa::Coord_t<In2>> &rigidState, const std::string &context) const {

		std::cout << "\n=== RIGID STATE DEBUG" << (context.empty() ? "" : " (" + context + ")") << " ===\n";
		std::cout << "Rigid state size: " << rigidState.size() << std::endl;

		for (size_t i = 0; i < rigidState.size(); ++i) {
			const auto &coord = rigidState[i];
			const auto &center = coord.getCenter();
			const auto &orientation = coord.getOrientation();

			std::cout << "  Rigid[" << i << "]:";
			std::cout << " pos=[" << center[0] << ", " << center[1] << ", " << center[2] << "]";
			std::cout << " quat=[" << orientation[0] << ", " << orientation[1] << ", " << orientation[2] << ", "
					  << orientation[3] << "]" << std::endl;
		}
		std::cout << "==============================\n";
	}

	template<class TIn1, class TIn2, class TOut>
	void
	Strain2RigidCosseratMapping<TIn1, TIn2, TOut>::displayOutputFrames(const sofa::type::vector<OutCoord> &outputFrames,
																	const std::string &context) const {

		std::cout << "\n=== OUTPUT FRAMES DEBUG" << (context.empty() ? "" : " (" + context + ")") << " ===\n";
		std::cout << "Output frames size: " << outputFrames.size() << std::endl;

		for (size_t i = 0; i < outputFrames.size(); ++i) {
			const auto &frame = outputFrames[i];
			const auto &center = frame.getCenter();
			const auto &orientation = frame.getOrientation();

			std::cout << "  Frame[" << i << "]:";
			std::cout << " pos=[" << center[0] << ", " << center[1] << ", " << center[2] << "]";
			std::cout << " quat=[" << orientation[0] << ", " << orientation[1] << ", " << orientation[2] << ", "
					  << orientation[3] << "]" << std::endl;

			// Display distance to previous frame
			if (i > 0) {
				sofa::type::Vec3 diff = center - outputFrames[i - 1].getCenter();
				std::cout << "    Distance to prev: " << diff.norm() << std::endl;
			}
		}
		std::cout << "==================================\n";
	}

	template<class TIn1, class TIn2, class TOut>
	void Strain2RigidCosseratMapping<TIn1, TIn2, TOut>::displaySectionProperties(const std::string &context) const {

		std::cout << "\n=== SECTION PROPERTIES DEBUG" << (context.empty() ? "" : " (" + context + ")") << " ===\n";
		std::cout << "Section properties size: " << m_section_properties.size() << std::endl;

		for (size_t i = 0; i < m_section_properties.size(); ++i) {
			const auto &section = m_section_properties[i];
			const auto &strain = section.getStrainsVec();
			const auto &transform = section.getTransformation();

			std::cout << "  Section[" << i << "]:";
			std::cout << " length=" << section.getLength();
			std::cout << " strain=[" << strain << "]";
			std::cout << " indices=[" << section.getIndex0() << ", " << section.getIndex1() << "]" << std::endl;

			// Display transformation matrix
			const auto &translation = transform.translation();
			const auto &rotation = transform.rotation();
			std::cout << "    Transform: trans=[" << translation[0] << ", " << translation[1] << ", " << translation[2]
					  << "]";
			std::cout << " rot_det=" << rotation.matrix().determinant() << std::endl;
		}
		std::cout << "=====================================\n";
	}

	template<class TIn1, class TIn2, class TOut>
	void Strain2RigidCosseratMapping<TIn1, TIn2, TOut>::displayFrameProperties(const std::string &context) const {

		std::cout << "\n=== FRAME PROPERTIES DEBUG" << (context.empty() ? "" : " (" + context + ")") << " ===\n";
		std::cout << "Frame properties size: " << m_frameProperties.size() << std::endl;

		for (size_t i = 0; i < m_frameProperties.size(); ++i) {
			const auto &frame = m_frameProperties[i];
			const auto &transform = frame.getTransformation();

			std::cout << "  Frame[" << i << "]:";
			std::cout << " length=" << frame.getLength();
			std::cout << " frames_sect_length_="
					  << frame.getLength(); // Same as length, but explicitly named as requested

			if (i < m_indices_vectors.size()) {
				std::cout << " section_idx=" << m_indices_vectors[i];
			}

			// Display distance to nearest beam node
			std::cout << " distance_to_nearest_beam_node=" << frame.getDistanceToNearestBeamNode();

			const auto &translation = transform.translation();
			const auto &rotation = transform.rotation();
			std::cout << " trans=[" << translation[0] << ", " << translation[1] << ", " << translation[2] << "]";
			std::cout << " rot_det=" << rotation.matrix().determinant() << std::endl;

			// Display adjoint matrix (6x6 matrix)
			// std::cout << "    adjoint_=[";
			// const auto &adjoint = frame.getAdjoint();
			// for (int row = 0; row < 6; ++row) {
			// 	if (row > 0) std::cout << "             ";
			// 	std::cout << "[";
			// 	for (int col = 0; col < 6; ++col) {
			// 		std::cout << adjoint(row, col);
			// 		if (col < 5) std::cout << ", ";
			// 	}
			// 	std::cout << "]";
			// 	if (row < 5) std::cout << ",\n";
			// }
			// std::cout << "]" << std::endl;
		}
		std::cout << "===================================\n";
	}

	template<class TIn1, class TIn2, class TOut>
	void Strain2RigidCosseratMapping<TIn1, TIn2, TOut>::displaySE3Transform(const SE3Types &transform,
																		 const std::string &name) const {

		std::cout << "\n=== SE3 TRANSFORM DEBUG: " << name << " ===\n";

		const auto &translation = transform.translation();
		const auto &rotation = transform.rotation();

		std::cout << "Translation: [" << translation[0] << ", " << translation[1] << ", " << translation[2] << "]\n";
		std::cout << "Rotation matrix:\n";
		const auto &R = rotation.matrix();
		for (int i = 0; i < 3; ++i) {
			std::cout << "  [" << R(i, 0) << ", " << R(i, 1) << ", " << R(i, 2) << "]\n";
		}
		std::cout << "Rotation determinant: " << R.determinant() << std::endl;

		// Convert to quaternion and display
		Eigen::Quaternion<double> quat(R);
		std::cout << "Quaternion: [" << quat.w() << ", " << quat.x() << ", " << quat.y() << ", " << quat.z() << "]\n";

		std::cout << "Matrix form:\n";
		const auto &M = transform.matrix();
		for (int i = 0; i < 4; ++i) {
			std::cout << "  [" << M(i, 0) << ", " << M(i, 1) << ", " << M(i, 2) << ", " << M(i, 3) << "]\n";
		}
		std::cout << "==========================================\n";
	}

	template<class TIn1, class TIn2, class TOut>
	void Strain2RigidCosseratMapping<TIn1, TIn2, TOut>::displayMappingState(const std::string &context) const {

		std::cout << "\n=== MAPPING STATE DEBUG" << (context.empty() ? "" : " (" + context + ")") << " ===\n";
		std::cout << "Base index: " << d_baseIndex.getValue() << std::endl;
		std::cout << "Debug mode: " << (d_debug.getValue() ? "ON" : "OFF") << std::endl;

		// @Todo change and use m_frameProperties instead
		const auto &curvFrames = d_curv_abs_frames.getValue();
		std::cout << "Curv abs frames size: " << curvFrames.size() << std::endl;
		if (!curvFrames.empty()) {
			std::cout << "  Values: [";
			for (size_t i = 0; i < curvFrames.size(); ++i) {
				std::cout << curvFrames[i];
				if (i < curvFrames.size() - 1)
					std::cout << ", ";
			}
			std::cout << "]\n";
		}

		std::cout << "Indices vectors size: " << m_indices_vectors.size() << std::endl;
		if (!m_indices_vectors.empty()) {
			std::cout << "  Values: [";
			for (size_t i = 0; i < m_indices_vectors.size(); ++i) {
				std::cout << m_indices_vectors[i];
				if (i < m_indices_vectors.size() - 1)
					std::cout << ", ";
			}
			std::cout << "]\n";
		}

		std::cout << "Beam length vectors size: " << m_beam_length_vectors.size() << std::endl;
		if (!m_beam_length_vectors.empty()) {
			std::cout << "  Values: [";
			for (size_t i = 0; i < m_beam_length_vectors.size(); ++i) {
				std::cout << m_beam_length_vectors[i];
				if (i < m_beam_length_vectors.size() - 1)
					std::cout << ", ";
			}
			std::cout << "]\n";
		}

		std::cout << "==============================\n";
	}

	template<class TIn1, class TIn2, class TOut>
	void Strain2RigidCosseratMapping<TIn1, TIn2, TOut>::displayVelocities(
			const sofa::type::vector<Deriv1> &strainVel, const sofa::type::vector<sofa::Deriv_t<In2>> &baseVel,
			const sofa::type::vector<OutDeriv> &outputVel, const std::string &context) const {

		std::cout << "\n=== VELOCITIES DEBUG" << (context.empty() ? "" : " (" + context + ")") << " ===\n";

		std::cout << "Strain velocities (size: " << strainVel.size() << "):" << std::endl;
		for (size_t i = 0; i < strainVel.size(); ++i) {
			std::cout << "  StrainVel[" << i << "]: ";

			if constexpr (std::is_same_v<Deriv1, sofa::type::Vec3>) {
				std::cout << "[" << strainVel[i][0] << ", " << strainVel[i][1] << ", " << strainVel[i][2] << "]";
			} else {
				std::cout << "[";
				for (int j = 0; j < strainVel[i].size(); ++j) {
					std::cout << strainVel[i][j];
					if (j < strainVel[i].size() - 1)
						std::cout << ", ";
				}
				std::cout << "]";
			}
			std::cout << std::endl;
		}

		std::cout << "\nBase velocities (size: " << baseVel.size() << "):" << std::endl;
		for (size_t i = 0; i < baseVel.size(); ++i) {
			const auto &vel = baseVel[i];
			std::cout << "  BaseVel[" << i << "]: [" << vel[0] << ", " << vel[1] << ", " << vel[2] << ", " << vel[3]
					  << ", " << vel[4] << ", " << vel[5] << "]" << std::endl;
		}

		std::cout << "\nOutput velocities (size: " << outputVel.size() << "):" << std::endl;
		for (size_t i = 0; i < outputVel.size(); ++i) {
			const auto &vel = outputVel[i];
			std::cout << "  OutputVel[" << i << "]: [" << vel[0] << ", " << vel[1] << ", " << vel[2] << ", " << vel[3]
					  << ", " << vel[4] << ", " << vel[5] << "]" << std::endl;
		}

		std::cout << "==========================\n";
	}

} // namespace Cosserat::mapping
