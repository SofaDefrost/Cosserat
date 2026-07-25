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

#include <Cosserat/mapping/Strain2FramesCosseratMapping.h>
#include <Cosserat/mapping/SofaLieGroupsUtils.h>
#include <sofa/core/MechanicalParams.h>
#include <sofa/core/Multi2Mapping.inl>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/gl/template.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/helper/logging/Message.h>
#include <sofa/helper/visual/DrawTool.h>
#include <sofa/type/Quat.h>

#include <cassert>
#include <map>
#include <string>

namespace Cosserat::mapping {

	using sofa::core::objectmodel::BaseContext;
	using sofa::helper::AdvancedTimer;
	using sofa::helper::WriteAccessor;
	using sofa::type::RGBAColor;
	using sofa::type::vector;
	using namespace sofa::component::cosserat::liegroups;

	template<class TIn1, class TIn2, class TOut>
	Strain2FramesCosseratMapping<TIn1, TIn2, TOut>::Strain2FramesCosseratMapping() :
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
							 "by another body.")),
		d_integrationMethod(initData(&d_integrationMethod, 0, "integrationMethod",
							 "Integration method for the Cosserat ODE g'(s)=g(s)·hat(ξ(s)).\n"
							 "0 = Euler (order 1, identical to expCosserat — default)\n"
							 "1 = Midpoint (order 2, RKMK2)\n"
							 "2 = RKMK4 (order 4, Magnus expansion)\n"
							 "For piecewise-constant strain all methods are equivalent.\n"
							 "Use RKMK4 with Legendre-polynomial strain fields.")),
		d_smoothPathSamples(initData(&d_smoothPathSamples, 8, "smoothPathSamples",
							 "Number of Bezier samples per section for the smoothed centerline.\n"
							 "Used by computeSmoothedPath() and the smooth-draw visualization.\n"
							 "1 = section endpoints only (equivalent to discrete simulation).\n"
							 "Higher values produce a smoother rendered rod.")),
		d_strainNoiseLevel(initData(&d_strainNoiseLevel, 1.0e-6, "strainNoiseLevel",
							 "Isotropic strain noise variance σ² for uncertainty propagation.\n"
							 "Σ_ξ = σ²·I₆ for each rod section (units: [1/m]²).\n"
							 "Set to 0 to disable (returns zero covariances).\n"
							 "Used by computeUncertainties().")) {

		// Register callback for updating frame transformations when geometry changes
		this->addUpdateCallback("updateFrames", {&d_curv_abs_section, &d_curv_abs_frames, &d_debug},
								[this](const sofa::core::DataTracker &t) {
									msg_info() << "Strain2FramesCosseratMapping updateFrames callback called";
									SOFA_UNUSED(t);
									this->updateGeometryInfo();
									std::cout << "====> Update Callback <====" << std::endl;
									const sofa::VecCoord_t<In1> &strain_state =
											m_strain_state->read(sofa::core::vec_id::read_access::position)->getValue();

									// This is also done in apply() So, no really need here !!!
									this->updateFrameTransformations(strain_state);
									return sofa::core::objectmodel::ComponentState::Valid;
								},
								{});
	}

	template<class TIn1, class TIn2, class TOut>
	void Strain2FramesCosseratMapping<TIn1, TIn2, TOut>::doBaseCosseratInit() {
		// Initialize colormap for visualization.
		// (SOFA deleted setColorScheme/reinit; build the palette via the
		//  string constructor instead.)
		m_colorMap = sofa::helper::ColorMap(256, "Blue to Red");

		msg_info() << "Strain2FramesCosseratMapping initialized with liegroups SE(3) integration";
	}

	template<class TIn1, class TIn2, class TOut>
	void Strain2FramesCosseratMapping<TIn1, TIn2, TOut>::initialization() {	
			// Get the initial configuration g(X):frames and initialize FrameInfo objects
		if (m_frames) {
			auto xfromData = m_frames->read(sofa::core::vec_id::read_access::position);
			const auto &xfrom = xfromData->getValue();

			// Initialize frame properties using the initial frame states
			const auto frame_count = xfrom.size();

			m_frame_properties.clear();
			m_frame_properties.reserve(frame_count);

			for (size_t i = 0; i < frame_count; ++i) {
				// Convert SOFA Rigid3 coord -> SE3 (via shared helper)
				const SE3Type gXi = rigidCoordToSE3(xfrom[i]);

				// Length and kappa will be set later in initializeFrameProperties
				FrameInfo frameInfo;
				frameInfo.setTransformation(gXi);
				m_frame_properties.push_back(frameInfo);
			}
		}

	}

	template<class TIn1, class TIn2, class TOut>
	void
	Strain2FramesCosseratMapping<TIn1, TIn2, TOut>::apply(const sofa::core::MechanicalParams * /* mparams */,
													  const vector<sofa::DataVecCoord_t<Out> *> &dataVecOutPos,
													  const vector<const sofa::DataVecCoord_t<In1> *> &dataVecIn1Pos,
													  const vector<const sofa::DataVecCoord_t<In2> *> &dataVecIn2Pos) {


		msg_info("Strain2FramesCosseratMapping") << "Strain2FramesCosseratMapping::apply called";

		if (dataVecOutPos.empty() || dataVecIn1Pos.empty() || dataVecIn2Pos.empty())
			return;

		// Check component state for validity
		if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
			return;

		if (d_debug.getValue())
			std::cout << " ########## Apply Function ########" << std::endl;
		
		
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

		// Update frame transformations using liegroups SE(3) exponential map
		// update the Exponential matrices according to new deformation
		// Here we update m_framesExponentialSE3Vectors & m_nodesExponentialSE3Vectors
		// Which are the homogeneous matrices of the frames: name g(X), beam's configuration
		updateFrameTransformations(strainState);

		// Get base frame transformation from rigid input (via SOFA<->Eigen helper)
		const SE3Types base_frame = rigidCoordToSE3(rigidBase[baseIndex]);

		
		// Cache the printLog value out of the loop
		bool doPrintLog = this->f_printLog.getValue();

		// Apply transformations to compute output frames
		for (unsigned int i = 0; i < frame_count; i++) {
			// Bounds checking
			assert(i < m_frame_properties.size() && "Frame index out of bounds");
			assert(i < output_frames.size() && "Output frames index out of bounds");

			// Start with the base frame
			auto current_frame = base_frame;

			// Apply section transformations up to the frame
			const auto related_beam_idx = m_frame_properties[i].getRelatedSectionIndex();
			assert(related_beam_idx <= m_section_properties.size() && "Invalid beam index");

			for (unsigned int j = 0; j < related_beam_idx; j++) {
				assert(j < m_section_properties.size() && "Section index out of bounds");
				// Compose with section transformation
				//// frame = gX(L_0)*...*gX(L_{n-1})
				current_frame = current_frame * m_section_properties[j].getTransformation();
			}

			// Apply additional frame transformation
			// frame*gX(x)
			current_frame = current_frame * m_frame_properties[i].getTransformation();
			// std::cout<<"Frame: "<<current_frame<<std::endl;

			if(d_debug.getValue())
				std::cout << "Frame  : " << i << " = " << current_frame << std::endl;
				
			// Save current rigid frame transformation into frame's properties
			//m_frame_properties[i].setTransformation(current_frame);

			// Convert SE3 -> SOFA rigid coordinates (via helper)
			se3ToRigidCoord(current_frame, output_frames[i]);

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
	void Strain2FramesCosseratMapping<TIn1, TIn2, TOut>::updateFrameTransformations(
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

		// Update frame properties based on their position within sections
		for (size_t i = 0; i < m_frame_properties.size(); ++i) {
			if (i < m_frame_to_section_indices.size()) {
				int sectionIndex = m_frame_properties[i].getRelatedSectionIndex();
				if (sectionIndex >= 0 && sectionIndex < static_cast<int>(vec_of_strains.size() + 1)) {
					// Compute frame transformation at its specific position
					SE3Types frame_gx = SE3Types::expCosserat(m_section_properties[sectionIndex].getStrainsVec(),
															  m_frame_properties[i].getDistanceToSectionStart());
					m_frame_properties[i].setTransformation(frame_gx);
				}
			}
		}

	}


	template<class TIn1, class TIn2, class TOut>
	void
	Strain2FramesCosseratMapping<TIn1, TIn2, TOut>::applyJ(const sofa::core::MechanicalParams * /* mparams */,
													   const vector<sofa::DataVecDeriv_t<Out> *> &dataVecOutVel,
													   const vector<const sofa::DataVecDeriv_t<In1> *> &dataVecIn1Vel,
													   const vector<const sofa::DataVecDeriv_t<In2> *> &dataVecIn2Vel) {


		if (dataVecOutVel.empty() || dataVecIn1Vel.empty() || dataVecIn2Vel.empty())
			return;

		if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
			return;

		if (d_debug.getValue())
			std::cout << " ########## Strain2FramesCosseratMapping ApplyJ Function ########" << std::endl;

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

		// 2.2 Apply the local transform from SOFA's frame to Cosserat's frame.
		// Project global velocity into the base body frame.
		const SE3Types absoluteFrame0     = rigidCoordToSE3(framePositions[0]);
		const SE3Types absoluteFrame0_inv = absoluteFrame0.inverse();
		const AdjointMatrix base_projector =
			absoluteFrame0_inv.buildProjectionMatrix(absoluteFrame0_inv.rotation().matrix());


		// 3. Compute velocity at each section node
		std::vector<TangentVector> node_velocities;
		node_velocities.resize(m_section_properties.size());

		// Base node velocity (transformed from SOFA frame)
		node_velocities[0] = base_projector * base_vel_local;
		if (d_debug.getValue())
    		std::cout << "Base local Velocity :" << node_velocities[0].transpose() << std::endl;

		
		for (size_t i = 1; i < m_section_properties.size(); ++i) {
			const auto &section = m_section_properties[i];
			const auto &tang_adj = section.getTangAdjointMatrix();

			// Extract strain velocity for this section
			TangentVector strain_vel_i = TangentVector::Zero();
			
			if constexpr (std::is_same_v<Deriv1, sofa::type::Vec3>){ 
				for (int j = 0; j < 3; ++j)
					strain_vel_i[j] = strain_vel[i - 1][j];
			}
			else{
				for (int j = 0; j < 6; ++j)
					strain_vel_i[j] = strain_vel[i - 1][j];
			}

			// Propagate velocity: η_i = Ad_{g_i^{-1}} * (η_{i-1} + T_i * ξ̇_i)
			// where Ad_{g_i^{-1}} is the inverse adjoint (transpose for SE(3))
			auto Ad_gim1 = section.inverse().getAdjoint();
			node_velocities[i] = Ad_gim1 * (node_velocities[i - 1] + tang_adj * strain_vel_i);
			if (d_debug.getValue()) {
				std::cout << "Node velocity [" << i << "]: " << node_velocities[i].transpose()<<"\n";
			}
		}

		// 4. Compute velocity at each output frame
		for (size_t i = 0; i < frame_count; ++i) {
			const auto &frame = m_frame_properties[i];
			const auto &tang_adj = frame.getTangAdjointMatrix();

			// Get the section index this frame belongs to
			int section_idx = (i < m_frame_to_section_indices.size()) ? m_frame_to_section_indices[i] - 1 : 0;

			// Ensure valid section index
			if (section_idx < 0 || section_idx >= static_cast<int>(node_velocities.size())) {
				section_idx = 0;
			}

			// Extract frame strain velocity (same as section strain).
			//
			// SAFETY (P1 fix): m_frame_to_section_indices[i] is the 1-based section index of
			// this frame.  For a frame located exactly on the rigid base (section 0),
			// the index can be 0 and `m_frame_to_section_indices[i] - 1` would wrap to SIZE_MAX
			// on the unsigned arithmetic of the subscript, causing a heap OOB read.
			// We clamp to [0, strain_vel.size()-1] and skip strain contribution if no
			// strain section is associated with this frame.
			TangentVector frame_strain_vel = TangentVector::Zero();
			const int strain_idx = static_cast<int>(m_frame_to_section_indices[i]) - 1;
			const bool has_strain = (strain_idx >= 0) &&
									(strain_idx < static_cast<int>(strain_vel.size()));

			if (has_strain) {
				constexpr int K = std::is_same_v<Deriv1, sofa::type::Vec3> ? 3 : 6;
				for (int j = 0; j < K; ++j)
					frame_strain_vel[j] = strain_vel[strain_idx][j];
			}

			// Compute frame velocity: η_frame = Ad_{g_frame^{-1}} * (η_node + T_frame * ξ̇_frame)
			// (we explicitly invert g_frame here to get the adjoint of g^{-1}, not Ad^T)
			const auto g_inv = frame.getInverseTransformation();
			const AdjointMatrix Ad_gm1 = g_inv.computeAdjoint();

			TangentVector eta_frame =
					Ad_gm1 * (node_velocities[m_frame_to_section_indices[i]-1] + tang_adj * frame_strain_vel);

			// Project body twist to SOFA Rigid3 convention (P(R) · η_body).
			const SE3Types absoluteFrame = rigidCoordToSE3(framePositions[i]);
			const AdjointMatrix frame_projector =
				absoluteFrame.buildProjectionMatrix(absoluteFrame.rotation().matrix());

			
			TangentVector output_vel = frame_projector * eta_frame;

			// Convert to SOFA format
			for (int k = 0; k < 6; ++k) {
				frame_vel[i][k] = output_vel[k];
			}
			if (d_debug.getValue()) {
				std::cout << "Frame velocity [" << i << "]: " << output_vel.transpose() <<"\n";
			}

		}

		// Debug output velocities if enabled
		if (d_debug.getValue()) {
			displayVelocities(strain_vel, base_vel, frame_vel, "applyJ - output");
		}

		dataVecOutVel[0]->endEdit();

	}

	template<class TIn1, class TIn2, class TOut>
	void Strain2FramesCosseratMapping<TIn1, TIn2, TOut>::applyJT(
			const sofa::core::MechanicalParams * /*mparams*/,
			const vector<sofa::DataVecDeriv_t<In1> *> &dataVecOut1Force,
			const vector<sofa::DataVecDeriv_t<In2> *> &dataVecOut2Force,
			const vector<const sofa::DataVecDeriv_t<Out> *> &dataVecInForce) {


		if (dataVecOut1Force.empty() || dataVecInForce.empty() || dataVecOut2Force.empty())
			return;

		if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
			return;

		if (d_debug.getValue())
			std::cout << " ########## Strain2FramesCosseratMapping ApplyJT Force Function ########" << std::endl;



		const sofa::VecDeriv_t<Out> &inputForces = dataVecInForce[0]->getValue();
		sofa::VecDeriv_t<In1> &strainForces = *dataVecOut1Force[0]->beginEdit();
		sofa::VecDeriv_t<In2> &baseForces = *dataVecOut2Force[0]->beginEdit();
		const auto baseIndex = d_baseIndex.getValue();

		// Get current positions to compute transformations
		const sofa::VecCoord_t<Out> &framePositions =
				this->m_frames->read(sofa::core::vec_id::read_access::position)->getValue();
		
		const sofa::DataVecCoord_t<In1> *x1fromData =
     		 	m_strain_state->read(sofa::core::vec_id::read_access::position);
		const sofa::VecCoord_t<In1> strainState = x1fromData->getValue();

		// Initialize output forces
		strainForces.resize(strainState.size());

		// Convert input forces from global frame to local frame and accumulate
		vector<TangentVector> localForces;
		auto tab_size = inputForces.size();
		localForces.clear();


		for (size_t i = 0; i < tab_size; ++i) {
			// Convert SOFA force to SE(3) tangent vector
			TangentVector frameForce = TangentVector::Zero();
			for (unsigned j = 0; j < 6; j++)
				frameForce[j] = inputForces[i][j];

			// Transform from global SOFA frame to local beam frame
			const SE3Types absoluteFrame = rigidCoordToSE3(framePositions[i]);
			const AdjointMatrix P_trans =
				absoluteFrame.buildProjectionMatrix(absoluteFrame.rotation().matrix());

			localForces.push_back(P_trans.transpose() * frameForce);
		}

		// Process forces following the beam structure (similar to DiscreteCosseratMapping)
		auto sz = m_frame_to_section_indices.size();

		if (sz == 0 || localForces.empty()) {
			dataVecOut1Force[0]->endEdit();
			dataVecOut2Force[0]->endEdit();
			// Add error or warning message 
			return;
		}
		

		auto lastSectionIndex = m_frame_to_section_indices[sz - 1];
		TangentVector totalForce = TangentVector::Zero();

		constexpr int N = std::is_same_v<Deriv1, sofa::type::Vec3> ? 3 : 6;
		Eigen::Matrix<double, N, 6> matB_trans = Eigen::Matrix<double, N, 6>::Zero();
		for (int k=0; k<N; k++)
			matB_trans(k, k) = 1.0;
		
		// Process frames in reverse order to accumulate forces
		for (auto s = sz; s--;) {

			int currentSectionIndex = m_frame_to_section_indices[s];
			const FrameInfo &frame = m_frame_properties[s];
		
			AdjointMatrix coAdjoint = frame.getCoAdjoint();
			
			
			TangentVector currentLocalForce = coAdjoint * localForces[s];
			
			AdjointMatrix temp = frame.getTangAdjointMatrix().transpose();
			
			Eigen::Matrix<double, N, 1>  f = matB_trans * temp * currentLocalForce;

			// Handle section change - propagate accumulated force
			if (lastSectionIndex != currentSectionIndex) {
				lastSectionIndex--;
				// Transform accumulated force to new section reference
				
				const SectionInfo &section = m_section_properties[lastSectionIndex];
				AdjointMatrix coAdjoint = section.getCoAdjoint();
				totalForce = coAdjoint * totalForce;

				AdjointMatrix tempSection = section.getTangAdjointMatrix().transpose();

				// apply F_tot to the new beam
				Eigen::Matrix<double, N, 1>  temp_f = matB_trans * tempSection * totalForce;
								
				// Add accumulated force to strain outpute
				for (int j=0; j<N; j++){
					strainForces[lastSectionIndex-1][j] +=temp_f[j];
				}
				
			}

			totalForce += currentLocalForce;
			for (int j=0; j<N; j++){
					strainForces[currentSectionIndex-1][j] +=f[j];
			}
		}
		
		// Propagate accumulated body wrench through base frame projector
		const SE3Types absoluteFrame0 = rigidCoordToSE3(framePositions[0]);
		const AdjointMatrix M =
			absoluteFrame0.buildProjectionMatrix(absoluteFrame0.rotation().matrix());
		const TangentVector toAdd = M * totalForce;
		
		for (int j=0; j<6; j++)
				baseForces[baseIndex][j] +=toAdd[j];	

		if (d_debug.getValue()) {
			std::cout << "Strain forces computed from " << inputForces.size() << " input forces" << std::endl;
			std::cout << "Total base force: [" << totalForce.transpose() << "]" << std::endl;
			std::cout << "Applied to base index: " << baseIndex << std::endl;
		}


		dataVecOut1Force[0]->endEdit();
		dataVecOut2Force[0]->endEdit();

	}

	template<class TIn1, class TIn2, class TOut>
	void Strain2FramesCosseratMapping<TIn1, TIn2, TOut>::applyJT(
			const sofa::core::ConstraintParams * /*cparams*/,
			const vector<sofa::DataMatrixDeriv_t<In1> *> &dataMatOut1Const,
			const vector<sofa::DataMatrixDeriv_t<In2> *> &dataMatOut2Const,
			const vector<const sofa::DataMatrixDeriv_t<Out> *> &dataMatInConst) {


		if (dataMatOut1Const.empty() || dataMatOut2Const.empty() || dataMatInConst.empty())
			return;

		if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
			return;

		if (d_debug.getValue())
			std::cout << " ########## Strain2FramesCosseratMapping ApplyJT Constraint Function ########" << std::endl;

		
		// Prepare input and output data matrices
		sofa::MatrixDeriv_t<In1> &out1 = *dataMatOut1Const[0]->beginEdit();
		sofa::MatrixDeriv_t<In2> &out2 = *dataMatOut2Const[0]->beginEdit();
		const sofa::MatrixDeriv_t<Out> &in = dataMatInConst[0]->getValue();

		// Get current positions to compute transformations
		const sofa::VecCoord_t<Out> &framePositions =
				this->m_frames->read(sofa::core::vec_id::read_access::position)->getValue();
		
		const sofa::DataVecCoord_t<In1> *x1fromData =
     		 	m_strain_state->read(sofa::core::vec_id::read_access::position);
		const sofa::VecCoord_t<In1> strainState = x1fromData->getValue();

		
		constexpr int N = std::is_same_v<Deriv1, sofa::type::Vec3> ? 3 : 6;
		Eigen::Matrix<double, N, 6> matB_trans = Eigen::Matrix<double, N, 6>::Zero();
		for (int k=0; k<N; k++)
			matB_trans(k, k) = 1.0;

		// Per-section accumulator : sectionIndex -> Σ localForce.
		// std::map is ordered by ascending key, so we'll iterate it in REVERSE
		// (rbegin/rend) to reproduce the previous "decreasing node index" sweep.
		// This replaces a vector<tuple<int,TangentVector>> + std::sort + a
		// fragile while-loop-coalesce step that had a long history of off-by-one
		// patches (cf. git blame for the previous implementation).
		std::map<int, TangentVector> accumByNode;

		// Process constraints
		for (auto rowIt = in.begin(); rowIt != in.end(); ++rowIt) {
			auto colIt = rowIt.begin();
			if (colIt == rowIt.end())
				continue;

			typename sofa::MatrixDeriv_t<In1>::RowIterator o1 = out1.writeLine(rowIt.index());
			typename sofa::MatrixDeriv_t<In2>::RowIterator o2 = out2.writeLine(rowIt.index());

			accumByNode.clear();

			while (colIt != rowIt.end()) {
				int frameIndex = colIt.index();

				TangentVector constraintValue = TangentVector::Zero();

				// Convert constraint value to TangentVector
				const sofa::Deriv_t<Out> val = colIt.val();
				for (unsigned int j = 0; j < 6; ++j) {
					constraintValue[j] = val[j];
				}

				int sectionIndex = m_frame_to_section_indices[frameIndex];

				const SE3Types absoluteFrame = rigidCoordToSE3(framePositions[frameIndex]);
				const AdjointMatrix P_trans =
					absoluteFrame.buildProjectionMatrix(absoluteFrame.rotation().matrix());
				
				const FrameInfo &frame = m_frame_properties[frameIndex];
			
				AdjointMatrix coAdjoint = frame.getCoAdjoint();

				TangentVector localForce = coAdjoint * P_trans.transpose() * constraintValue; // constraint direction in local frame of the beam.


				AdjointMatrix temp = frame.getTangAdjointMatrix().transpose();

				Eigen::Matrix<double, N, 1> f = matB_trans * temp * localForce; // constraint direction in the strain space.
				
				Deriv1 f_trans;
				for(int k=0; k<N; k++){
					f_trans[k] = f[k];
				}
				
				o1.addCol(sectionIndex-1, f_trans);

				// Accumulate per-section : if the same section is hit by multiple
				// frames in this constraint row, their localForce contributions add.
				auto [it_insert, inserted] = accumByNode.try_emplace(sectionIndex, localForce);
				if (!inserted) {
					it_insert->second += localForce;
				}
				colIt++;
			}

			// Iterate per-section accumulator in DECREASING section order
			// (rbegin/rend on a std::map sorted by ascending key).
			for (auto it = accumByNode.rbegin(); it != accumByNode.rend(); ++it) {
				int i = it->first;
				TangentVector CumulativeF = it->second;

				while(i>0){
					const SectionInfo &section = m_section_properties[i-1];
					AdjointMatrix coAdjoint = section.getCoAdjoint();

					CumulativeF = coAdjoint * CumulativeF;
					// transfer to strain space (local coordinates)
					AdjointMatrix tempSection = section.getTangAdjointMatrix().transpose();
					
					Eigen::Matrix<double, N, 1>  temp_f = matB_trans * tempSection * CumulativeF;

					Deriv1 temp_f_trans;
					for(int k=0; k<N; k++){
						temp_f_trans[k] = temp_f[k];
					}

					if(i>1){
						o1.addCol(i-2, temp_f_trans);
					}
					i--;
				}
				
				const SE3Types absoluteFrame0 = rigidCoordToSE3(framePositions[0]);
				const AdjointMatrix M =
					absoluteFrame0.buildProjectionMatrix(absoluteFrame0.rotation().matrix());
				const TangentVector base_force = M * CumulativeF;


				sofa::type::Vec6 base_force_trans;
				for (int k = 0; k < 6; ++k) {
					base_force_trans[k] = base_force[k];
				}

				o2.addCol(d_baseIndex.getValue(), base_force_trans);

			}
		
		}

		dataMatOut1Const[0]->endEdit();
		dataMatOut2Const[0]->endEdit();


	}


	template<class TIn1, class TIn2, class TOut>
	void Strain2FramesCosseratMapping<TIn1, TIn2, TOut>::draw(const sofa::core::visual::VisualParams *vparams) {
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
			for (unsigned int i = 0; i < sz - 1; i++) {
				j = m_frame_to_section_indices[i] - 1; // to get the articulation on which the frame is related to
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
	void Strain2FramesCosseratMapping<TIn1, TIn2, TOut>::computeBBox(const sofa::core::ExecParams *params,
																 bool onlyVisible) {
		// Compute bounding box for visualization
		// Implementation would calculate the extent of all frames
		Inherit::computeBBox(params, onlyVisible);
	}

	// Debug display functions (displayStrainState, displayRigidState, displayOutputFrames,
	// displaySectionProperties, displayFrameProperties, displaySE3Transform,
	// displayMappingState, displayVelocities) live in a separate _debug.inl included
	// at the bottom of this file — see Strain2FramesCosseratMapping_debug.inl.

	// applyDJT — PARTIAL geometric stiffness contribution (child wrench frozen).
	//
	// ⚠ SCOPE — read before relying on this for Newton convergence analysis.
	//
	// The exact directional derivative of Jᵀ(ξ)·f_x w.r.t. ξ has three terms:
	//   (1) ∂T/∂ξ · δξ · node_F            — tangent-map variation
	//   (2) Tᵀ · ∂coAd(F)/∂ξ · δξ · …      — coAdjoint variation
	//   (3) Tᵀ · coAd(F) · ∂P/∂ξ · δξ · …  — projector variation
	// This implementation computes ONLY term (2), faithfully ported from
	// DiscreteCosseratMapping::applyDJT which has the same limitation (its
	// finite-difference validation tests StraightBeam_FD / CurvedBeam_FD are
	// GTEST_SKIP'd for exactly this reason — see
	// tests/unit/DiscreteCosseratMappingApplyDJTTest.cpp for the analysis;
	// at ξ=0 terms (2) and (3) cancel, so term (2) alone is wrong even there).
	//
	// Consequences:
	//   - The returned matrix-vector product is NOT the exact mapping tangent;
	//     it does not match central finite differences of applyJT.
	//   - It is strictly more information than the previous zero stub, and the
	//     structural properties hold: zero for f_x = 0, linear in kFactor,
	//     mirrors applyJT's backward sweep (same coAdjoint / tangent sources).
	//   - Implicit solvers remain stable in our scene tests, but quadratic
	//     Newton convergence is NOT guaranteed by this term alone.
	//
	// Sweep structure (identical index handling to applyJT), with the geometric
	// factor ad(v)ᵀ where v = T · δξ is the local twist from the strain increment:
	//   (a) frame direct:      δf_k += kFactor · Bᵀ · T_frameᵀ · ad(v_s)ᵀ · node_F
	//   (b) section transport: δf_k += kFactor · Bᵀ · T_nodeᵀ  · ad(v_n)ᵀ · F_tot
	//
	// TODO(geometric-stiffness): implement terms (1) and (3) + the cross-section
	// couplings; then re-enable the FD tests. Until then this is an approximation
	// beyond the usual "freeze T" one (Simo & Vu-Quoc 1986).
	template<class TIn1, class TIn2, class TOut>
	void Strain2FramesCosseratMapping<TIn1, TIn2, TOut>::applyDJT(
			const sofa::core::MechanicalParams* mparams,
			sofa::core::MultiVecDerivId          inForce,
			sofa::core::ConstMultiVecDerivId     /*outForce*/) {

		if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
			return;

		const SReal kFactor = mparams->kFactor();
		if (kFactor == 0.0)
			return;

		// δξ : strain displacement (In1). f_x : child wrenches (Out), held constant.
		const sofa::VecDeriv_t<In1> &dx     = mparams->readDx(m_strain_state)->getValue();
		const sofa::VecDeriv_t<Out> &childF = mparams->readF(m_frames)->getValue();

		// Parent strain force output (In1, accumulate +=).
		sofa::Data<sofa::VecDeriv_t<In1>> &out1Data = *inForce[m_strain_state].write();
		sofa::VecDeriv_t<In1> &out1 = *out1Data.beginEdit();

		const sofa::VecCoord_t<Out> &framePositions =
				this->m_frames->read(sofa::core::vec_id::read_access::position)->getValue();

		// Refresh tangent-exp matrices (applyJ already does this each step; guard
		// against a standalone applyDJT call).
		this->updateTangExpSE3();

		const auto sz = m_frame_to_section_indices.size();
		if (sz == 0) {
			out1Data.endEdit();
			return;
		}

		// Convert child wrenches to the local beam frame (same as applyJT).
		vector<TangentVector> local_F(sz);
		for (size_t s = 0; s < sz; ++s) {
			TangentVector vec = TangentVector::Zero();
			for (unsigned j = 0; j < 6; ++j)
				vec[j] = childF[s][j];
			const SE3Types absoluteFrame = rigidCoordToSE3(framePositions[s]);
			const AdjointMatrix P_trans =
				absoluteFrame.buildProjectionMatrix(absoluteFrame.rotation().matrix());
			local_F[s] = P_trans.transpose() * vec;
		}

		constexpr int N = std::is_same_v<Deriv1, sofa::type::Vec3> ? 3 : 6;
		Eigen::Matrix<double, N, 6> matB_trans = Eigen::Matrix<double, N, 6>::Zero();
		for (int k = 0; k < N; ++k)
			matB_trans(k, k) = 1.0;

		// Little adjoint ad(v) for v = [φ; ρ] (angular head), same convention as
		// CosseratBeamGeometry::computeTangExpImplementation and DiscreteCosserat.
		auto littleAdjoint = [this](const TangentVector &v) -> AdjointMatrix {
			AdjointMatrix adv;
			this->buildAdjoint(this->getTildeMatrix(Vector3(v.template head<3>())),
							   this->getTildeMatrix(Vector3(v.template tail<3>())), adv);
			return adv;
		};

		// Embed a strain increment δξ_k (Vec3 → angular head, Vec6 → full) into a twist.
		auto embedStrainDelta = [&](int k) -> TangentVector {
			TangentVector xi_dot = TangentVector::Zero();
			for (int j = 0; j < N; ++j)
				xi_dot[j] = dx[k][j];
			return xi_dot;
		};

		auto lastSectionIndex = m_frame_to_section_indices[sz - 1];
		TangentVector totalForce = TangentVector::Zero();

		for (auto s = sz; s--;) {
			const int currentSectionIndex = static_cast<int>(m_frame_to_section_indices[s]);
			const FrameInfo &frame = m_frame_properties[s];

			// node_F = coAd(g_frame) · local_F[s]  (identical to applyJT's currentLocalForce)
			const TangentVector node_F = frame.getCoAdjoint() * local_F[s];

			// (a) frame direct geometric term
			{
				const int k_s = currentSectionIndex - 1;
				const AdjointMatrix &T_frame = frame.getTangAdjointMatrix();
				const TangentVector v_s     = T_frame * embedStrainDelta(k_s);
				const AdjointMatrix adv     = littleAdjoint(v_s);
				const Eigen::Matrix<double, N, 1> delta_f =
					kFactor * (matB_trans * (T_frame.transpose() * (adv.transpose() * node_F)));
				for (int j = 0; j < N; ++j)
					out1[k_s][j] += delta_f[j];
			}

			// Section boundary: transport accumulated wrench + (b) node term
			if (lastSectionIndex != m_frame_to_section_indices[s]) {
				lastSectionIndex--;
				const SectionInfo &section = m_section_properties[lastSectionIndex];
				totalForce = section.getCoAdjoint() * totalForce;

				const int k_node = static_cast<int>(lastSectionIndex) - 1;
				const AdjointMatrix &T_node = section.getTangAdjointMatrix();
				const TangentVector v_node  = T_node * embedStrainDelta(k_node);
				const AdjointMatrix adv_n   = littleAdjoint(v_node);
				const Eigen::Matrix<double, N, 1> delta_fn =
					kFactor * (matB_trans * (T_node.transpose() * (adv_n.transpose() * totalForce)));
				for (int j = 0; j < N; ++j)
					out1[k_node][j] += delta_fn[j];
			}

			totalForce += node_F;
		}

		out1Data.endEdit();
	}

} // namespace Cosserat::mapping

// Debug helpers (display*() definitions). Kept out-of-line so this file stays
// focused on the mapping algebra.
#include <Cosserat/mapping/Strain2FramesCosseratMapping_debug.inl>
