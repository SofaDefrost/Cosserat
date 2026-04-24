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
									std::cout << "====> Update Callback <====" << std::endl;
									const sofa::VecCoord_t<In1> &strain_state =
											m_strain_state->read(sofa::core::vec_id::read_access::position)->getValue();

									// This is also done in apply() So, no really need here !!!
									// this->updateFrameTransformations(strain_state);
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

		// Apply transformations to compute output frames
		for (unsigned int i = 0; i < frame_count; i++) {
			// Bounds checking
			assert(i < m_frameProperties.size() && "Frame index out of bounds");
			assert(i < output_frames.size() && "Output frames index out of bounds");

			// Start with the base frame
			auto current_frame = base_frame;

			// Apply section transformations up to the frame
			const auto related_beam_idx = m_frameProperties[i].get_related_beam_index_();
			assert(related_beam_idx <= m_section_properties.size() && "Invalid beam index");

			for (unsigned int j = 0; j < related_beam_idx; j++) {
				assert(j < m_section_properties.size() && "Section index out of bounds");
				// Compose with section transformation
				//// frame = gX(L_0)*...*gX(L_{n-1})
				current_frame = current_frame * m_section_properties[j].getTransformation();
			}

			// Apply additional frame transformation
			// frame*gX(x)
			current_frame = current_frame * m_frameProperties[i].getTransformation();

			if(d_debug.getValue())
				std::cout << "Frame  : " << i << " = " << current_frame << std::endl;
				
			// Save current rigid frame transformation into frame's properties
			m_frameProperties[i].setTransformation(current_frame);

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
			//	m_section_properties[i+1].getStrainsVec());

			auto section_length = m_section_properties[i + 1].getLength();
			SE3Types _gx = SE3Types::expCosserat(strain, section_length);

			m_section_properties[i + 1].setTransformation(_gx);
		}

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

		if (d_debug.getValue())
			std::cout << " ########## Strain2RigidCosseratMapping ApplyJ Function ########" << std::endl;

		const sofa::VecDeriv_t<In1> &strain_vel = dataVecIn1Vel[0]->getValue();
		const sofa::VecDeriv_t<In2> &base_vel = dataVecIn2Vel[0]->getValue();
		sofa::VecDeriv_t<Out> &frame_vel = *dataVecOutVel[0]->beginEdit();

		// Debug input velocities if enabled
		if (d_debug.getValue()) {
			displayVelocities(strain_vel, base_vel, frame_vel, "applyJ - input");
		}

		const auto base_index = d_baseIndex.getValue();
		const auto frame_count = d_curv_abs_frames.getValue().size();
		frame_vel.resize(frame_count);

		// 1. Compute current tangent exponential SE3 matrices
		this->updateTangExpSE3();

		// 2. Compute the base velocity in SE(3) tangent space
		// 2.1 Convert base velocity to se(3) tangent vector
		TangentVector base_vel_local = TangentVector::Zero();
		for (auto u = 0; u < 6; u++)
			base_vel_local[u] = base_vel[base_index][u];

		// 2.2 Apply the local transform from SOFA's frame to Cosserat's frame
		const SE3Types base_sofa_pos = m_frameProperties[0].getTransformation().inverse();
		AdjointMatrix base_projector = base_sofa_pos.buildProjectionMatrix(base_sofa_pos.rotation().matrix());

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
			if (i - 1 < strain_vel.size()) {
				// Handle both Vec3 and Vec6 input types
				if constexpr (std::is_same_v<typename sofa::Deriv_t<In1>, sofa::type::Vec3>) {
					for (int j = 0; j < 3; ++j) {
						strain_vel_i[j] = strain_vel[i - 1][j];
					}
				} else {
					for (int j = 0; j < 6 && j < strain_vel[i - 1].size(); ++j) {
						strain_vel_i[j] = strain_vel[i - 1][j];
					}
				}
			}

			// Propagate velocity: η_i = Ad_{g_i^{-1}} * (η_{i-1} + T_i * ξ̇_i)
			// where Ad_{g_i^{-1}} is the inverse adjoint (transpose for SE(3))
			node_velocities[i] = section.getAdjoint() * (node_velocities[i - 1] + tang_adj * strain_vel_i);

			if (d_debug.getValue()) {
				std::cout << "Node velocity [" << i << "]: " << node_velocities[i].transpose()<<"\n";
			}
		}

		// 4. Compute velocity at each output frame
		for (size_t i = 0; i < frame_count; ++i) {
			const auto &frame = m_frameProperties[i];
			const auto &tang_adj = frame.getTangAdjointMatrix();

			// Get the section index this frame belongs to
			int section_idx = (i < m_indices_vectors.size()) ? m_indices_vectors[i] - 1 : 0;

			// Ensure valid section index
			if (section_idx < 0 || section_idx >= static_cast<int>(node_velocities.size())) {
				section_idx = 0;
			}

			// Extract frame strain velocity (same as section strain)
			TangentVector frame_strain_vel = TangentVector::Zero();
			if (section_idx >= 0 && section_idx < static_cast<int>(strain_vel.size())) {
				if constexpr (std::is_same_v<typename sofa::Deriv_t<In1>, sofa::type::Vec3>) {
					for (int j = 0; j < 3; ++j) {
						frame_strain_vel[j] = strain_vel[section_idx][j];
					}
				} else {
					for (int j = 0; j < 6 && j < strain_vel[section_idx].size(); ++j) {
						frame_strain_vel[j] = strain_vel[section_idx][j];
					}
				}
			}

			// Compute frame velocity: η_frame = Ad_{g_frame^{-1}} * (η_node + T_frame * ξ̇_frame)
			TangentVector eta_frame =
					frame.getAdjoint() * (node_velocities[section_idx] + tang_adj * frame_strain_vel);

			// Project to output frame (convert from local to SOFA global frame)
			AdjointMatrix frame_projector =
					frame.getTransformation().buildProjectionMatrix(frame.getTransformation().rotation().matrix());
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
	void Strain2RigidCosseratMapping<TIn1, TIn2, TOut>::applyJT(
			const sofa::core::MechanicalParams * /*mparams*/,
			const vector<sofa::DataVecDeriv_t<In1> *> &dataVecOut1Force,
			const vector<sofa::DataVecDeriv_t<In2> *> &dataVecOut2Force,
			const vector<const sofa::DataVecDeriv_t<Out> *> &dataVecInForce) {

		if (dataVecOut1Force.empty() || dataVecInForce.empty() || dataVecOut2Force.empty())
			return;

		if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
			return;

		if (d_debug.getValue())
			std::cout << " ########## Strain2RigidCosseratMapping ApplyJT Force Function ########" << std::endl;

		const sofa::VecDeriv_t<Out> &inputForces = dataVecInForce[0]->getValue();
		sofa::VecDeriv_t<In1> &strainForces = *dataVecOut1Force[0]->beginEdit();
		sofa::VecDeriv_t<In2> &baseForces = *dataVecOut2Force[0]->beginEdit();
		const auto baseIndex = d_baseIndex.getValue();

		// Get current positions to compute transformations
		const sofa::VecCoord_t<Out> &framePositions =
				this->m_frames->read(sofa::core::vec_id::read_access::position)->getValue();
		const sofa::VecCoord_t<In1> &strainState =
				this->m_strain_state->read(sofa::core::vec_id::read_access::position)->getValue();

		// Initialize output forces
		strainForces.resize(strainState.size());


		updateFrameTransformations(strainState); // :) On n'avait pas fait cette mise à jour ! (raison pour laquelle frameExp == frame :{, ce qui n'est pas correct)

		// Convert input forces from global frame to local frame and accumulate
		vector<TangentVector> localForces;
		auto tab_size = inputForces.size();

		localForces.reserve(tab_size);

		for (size_t i = 0; i < tab_size; ++i) {
			// Convert SOFA force to SE(3) tangent vector
			TangentVector frameForce = TangentVector::Zero();
			//ajout2: passer de for (int j = 0; j < 6 && j < static_cast<int>(inputForces[i].size()); ++j)
			//		  à for (unsigned j = 0; j < 6; j++)
			for (unsigned j = 0; j < 6; j++) {
				frameForce[j] = inputForces[i][j];
			}

			// Transform from global SOFA frame to local beam frame 
			const auto &pos = framePositions[i]; //utiliser framePositions (position dans le repère monde au lieu de frameProperties qui est expo. ds strain)
				// std::cout<<"=> frame "<<pos<<std::endl;
				
			Vector3 translation(pos.getCenter()[0], pos.getCenter()[1], pos.getCenter()[2]);
			const auto &quat = pos.getOrientation();
			Eigen::Quaternion<double> rotation(quat[3], quat[0], quat[1], quat[2]);
				
			SE3Types absoluteFrame(SE3Types::SO3Type(rotation), translation);
			//projection de la force globale dans le repère local
			AdjointMatrix P_trans = absoluteFrame.buildProjectionMatrix(absoluteFrame.rotation().matrix());
			TangentVector localForce = P_trans.transpose() * frameForce;
			
			// std::cout<<"Projection matrix:"<<std::endl;
				// std::cout<<P_trans<<std::endl;
				// std::cout<<"Local force:"<<std::endl;
				// std::cout<<localForce.transpose()<<std::endl;				
			localForces.push_back(localForce);
			// std::cout<<"#############################"<< std::endl;
			
		}
        //ajout5
		// std::cout<<"Out of force transformation loop"<<std::endl;
		// std::cout<<"Local forces size: "<<localForces.size()<<std::endl;
		// //ajout6
		// std::cout<<"Local forces:"<<std::endl;
		// for (size_t i = 0; i < localForces.size(); ++i) {
		// 	std::cout<<"Force "<<i<<": "<<localForces[i].transpose()<<std::endl;
		// }
		// Process forces following the beam structure (similar to DiscreteCosseratMapping)
		auto sz = m_indices_vectors.size();
		//ajout7
		std::cout<<"Section count: "<<sz<<std::endl;
		if (sz == 0 || localForces.empty()) {
			dataVecOut1Force[0]->endEdit();
			dataVecOut2Force[0]->endEdit();
			return;
		}
		auto lastSectionIndex = m_indices_vectors[sz - 1];
		TangentVector totalForce = TangentVector::Zero();

		// Process frames in reverse order to accumulate forces
		//ajout9: Modification de la loop.. passer de for (int s = static_cast<int>(sz) - 1; s >= 0; s--)
		//									à for (auto s = sz; s--;)

		Eigen::Matrix<double, 3, 6> matB_trans = Eigen::Matrix<double, 3, 6>::Zero();
		for (int k=0; k<3; k++)
			matB_trans(k, k) = 1.0;

		for (auto s = sz; s--;) {

			int currentSectionIndex = m_indices_vectors[s];
			const FrameInfo &frame = m_frameProperties[s];
			
			std::cout<< "=== s: " << s << " frame: "<< frame.getTransformation() << std::endl;
			std::cout<< "=== s: " << s << " framePos: "<< framePositions[s] << std::endl;

			//@TODO Swap the computation of Adjoint and CoAdjoint
			AdjointMatrix coAdjoint = frame.getAdjoint();// Etonnant getAdjoint donne exactement computeCoAdjoint de DCM !!
			
			std::cout<<"CoAdjoint Matrix frame"<<std::endl;
			std::cout<< coAdjoint <<std::endl;
			
			
			TangentVector currentLocalForce = coAdjoint * localForces[s];
			
			AdjointMatrix temp = frame.getTangAdjointMatrix().transpose();
			
			std::cout << "Exponential Tangent Matrix (transpose) frame"<<std::endl;
			std::cout << temp <<std::endl;
			
			Vector3 f = matB_trans * temp * currentLocalForce;

			//@todo : Use the seclector matrix B which is 3x6 or 6x6

			// Handle section change - propagate accumulated force
			if (lastSectionIndex != currentSectionIndex) {
				lastSectionIndex--;
				// Transform accumulated force to new section reference
				
				const SectionInfo &section = m_section_properties[lastSectionIndex];
				AdjointMatrix coAdjoint = section.getAdjoint();
				totalForce = coAdjoint * totalForce;
				
				std::cout<<"CoAdjoint Matrix section: " <<std::endl; std::cout<< coAdjoint <<std::endl;
				std::cout<<"Total force : "<<totalForce.transpose()<<std::endl;
				
				// //ajout 
				AdjointMatrix temp = section.getTangAdjointMatrix().transpose();
				std::cout<<"Tangent Exp (transpose) section"<< std::endl; std::cout<< temp<<std::endl;

				// apply F_tot to the new beam
				Vector3 temp_f = matB_trans * temp * totalForce;
								
				// Add accumulated force to strain output
				for (int j=0; j<3; j++)
					strainForces[lastSectionIndex-1][j] +=temp_f[j];
				
				
			}
			
			totalForce += currentLocalForce;
			for (int j=0; j<f.size(); j++){
					strainForces[currentSectionIndex-1][j] +=f[j];
			}

		}

		auto frame0 = framePositions[0];		
		Vector3 trans0(frame0.getCenter()[0], frame0.getCenter()[1], frame0.getCenter()[2]);
		const auto &quat0 = frame0.getOrientation();
		Eigen::Quaternion<double> rot0(quat0[3], quat0[0], quat0[1], quat0[2]);
			
		SE3Types absoluteFrame0(SE3Types::SO3Type(rot0), trans0);
		//projection de la force globale dans le repère local
		AdjointMatrix M = absoluteFrame0.buildProjectionMatrix(absoluteFrame0.rotation().matrix());
		std::cout<<"Proj Mat 0: "<<"\n"<<M <<std::endl;

		TangentVector toAdd = TangentVector::Zero();
		toAdd = M * totalForce;
		
		for (int j=0; j<6; j++)
				baseForces[baseIndex][j] +=toAdd[j];	
		

		std::cout << "Node forces " << strainForces<< std::endl;
		std::cout << "base Force: " << baseForces[baseIndex] << std::endl;

		if (d_debug.getValue()) {
			std::cout << "Strain forces computed from " << inputForces.size() << " input forces" << std::endl;
			std::cout << "Total base force: [" << totalForce.transpose() << "]" << std::endl;
			std::cout << "Applied to base index: " << baseIndex << std::endl;
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

		if (d_debug.getValue())
			std::cout << " ########## Strain2RigidCosseratMapping ApplyJT Constraint Function ########" << std::endl;

		// Prepare input and output data matrices
		sofa::MatrixDeriv_t<In1> &out1 = *dataMatOut1Const[0]->beginEdit();
		sofa::MatrixDeriv_t<In2> &out2 = *dataMatOut2Const[0]->beginEdit();
		const sofa::MatrixDeriv_t<Out> &in = dataMatInConst[0]->getValue();

		// Process constraints
		for (auto rowIt = in.begin(); rowIt != in.end(); ++rowIt) {
			auto colIt = rowIt.begin();
			if (colIt == rowIt.end())
				continue;

			typename sofa::MatrixDeriv_t<In1>::RowIterator o1 = out1.writeLine(rowIt.index());
			typename sofa::MatrixDeriv_t<In2>::RowIterator o2 = out2.writeLine(rowIt.index());

			while (colIt != rowIt.end()) {
				int frameIndex = colIt.index();
				TangentVector constraintValue;
				// Convert constraint value to TangentVector
				const auto &val = colIt.val();
				for (unsigned int j = 0; j < 6 && j < val.size(); ++j) {
					constraintValue[j] = val[j];
				}

				const FrameInfo &frame = m_frameProperties[frameIndex];
				AdjointMatrix adjoint = frame.getAdjoint();
				TangentVector localForce = adjoint.transpose() * constraintValue;

				int sectionIndex = m_indices_vectors[frameIndex] - 1;
				if (sectionIndex >= 0) {
					const SectionInfo &section = m_section_properties[sectionIndex];
					AdjointMatrix coAdjoint = section.getCoAdjoint();
					TangentVector strainSpaceForce = coAdjoint * localForce;

					// Convert Eigen vector to SOFA Vec3
					auto strainForce3D = strainSpaceForce.head<3>();
					sofa::type::Vec3 sofaStrainForce(strainForce3D[0], strainForce3D[1], strainForce3D[2]);
					o1.addCol(sectionIndex, sofaStrainForce);

					while (sectionIndex-- > 0) {
						const SectionInfo &prevSection = m_section_properties[sectionIndex];
						coAdjoint = prevSection.getCoAdjoint();
						strainSpaceForce = coAdjoint * strainSpaceForce;

						if (sectionIndex > 0) {
							auto prevStrainForce3D = strainSpaceForce.head<3>();
							sofa::type::Vec3 sofaPrevStrainForce(prevStrainForce3D[0], prevStrainForce3D[1],
																 prevStrainForce3D[2]);
							o1.addCol(sectionIndex - 1, sofaPrevStrainForce);
						}
					}

					// Convert constraintValue to SOFA format for base
					sofa::type::Vec6 sofaConstraintValue;
					for (int k = 0; k < 6; ++k) {
						sofaConstraintValue[k] = constraintValue[k];
					}
					auto baseIndex = d_baseIndex.getValue();
					o2.addCol(baseIndex, sofaConstraintValue);
				}
				++colIt;
			}
		}

		dataMatOut1Const[0]->endEdit();
		dataMatOut2Const[0]->endEdit();
	}

	// template<class TIn1, class TIn2, class TOut>
	// void Strain2RigidCosseratMapping<TIn1, TIn2, TOut>::draw(const sofa::core::visual::VisualParams *vparams) {
	// 	if (!vparams->displayFlags().getShowMappings())
	// 		return;

	// 	// Draw implementation similar to DiscreteCosseratMapping
	// 	// This would include beam visualization with colormap
	// 	if (d_drawMapBeam.getValue()) {
	// 		// Draw colored beam based on deformation
	// 		// Implementation would depend on specific visualization requirements
	// 	}
	// }

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
