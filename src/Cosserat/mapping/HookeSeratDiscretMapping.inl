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

#include <Cosserat/mapping/HookeSeratDiscretMapping.h>
#include <sofa/core/Multi2Mapping.inl>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/gl/template.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/helper/logging/Message.h>
#include <sofa/helper/visual/DrawTool.h>
#include <sofa/type/Quat.h>

#include <string>

namespace Cosserat::mapping {

	using sofa::core::objectmodel::BaseContext;
	using sofa::helper::AdvancedTimer;
	using sofa::helper::WriteAccessor;
	using sofa::type::RGBAColor;
	using sofa::type::vector;
	using namespace sofa::component::cosserat::liegroups;

	template<class TIn1, class TIn2, class TOut>
	HookeSeratDiscretMapping<TIn1, TIn2, TOut>::HookeSeratDiscretMapping() :
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
									msg_info() << "HookeSeratDiscretMapping updateFrames callback called";
									SOFA_UNUSED(t);
									std::cout << "====> Update Callback <===="<<std::endl;
									const sofa::VecCoord_t<In1> &strain_state =
											m_strain_state->read(sofa::core::vec_id::read_access::position)->getValue();

									// This is also done in apply() So, no really need here !!!
									//this->updateFrameTransformations(strain_state);
									return sofa::core::objectmodel::ComponentState::Valid;
								},
								{});
	}

	template<class TIn1, class TIn2, class TOut>
	void HookeSeratDiscretMapping<TIn1, TIn2, TOut>::doBaseCosseratInit() {
		// Initialize colormap for visualization
		m_colorMap.setColorScheme("Blue to Red");
		m_colorMap.reinit();

		msg_info() << "HookeSeratDiscretMapping initialized with liegroups SE(3) integration";
	}

	template<class TIn1, class TIn2, class TOut>
	void
	HookeSeratDiscretMapping<TIn1, TIn2, TOut>::apply(const sofa::core::MechanicalParams * /* mparams */,
													  const vector<sofa::DataVecCoord_t<Out> *> &dataVecOutPos,
													  const vector<const sofa::DataVecCoord_t<In1> *> &dataVecIn1Pos,
													  const vector<const sofa::DataVecCoord_t<In2> *> &dataVecIn2Pos) {

		msg_info("HookeSeratDiscretMapping") << "HookeSeratDiscretMapping::apply called";

		if (dataVecOutPos.empty() || dataVecIn1Pos.empty() || dataVecIn2Pos.empty())
			return;

		// Check component state for validity
		if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
			return;

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
		const auto &baseRigid = rigidBase[baseIndex];
		Vector3 base_trans(baseRigid.getCenter()[0], baseRigid.getCenter()[1], baseRigid.getCenter()[2]);

		// Convert SOFA quaternion to Eigen quaternion (SOFA: x,y,z,w; Eigen: w,x,y,z)
		const auto &base_sofa_quat = baseRigid.getOrientation();
		Eigen::Quaternion<double> base_rot(base_sofa_quat[3], base_sofa_quat[0], base_sofa_quat[1], base_sofa_quat[2]);

		std::cout << "The base config  : "<< base_trans.transpose() << " \n" << base_rot << std::endl;
		// Create SE3 transformation
		SE3Types base_frame(SE3Types::SO3Type(base_rot), base_trans);

		// Cache the printLog value out of the loop
		bool doPrintLog = this->f_printLog.getValue();

		// Apply transformations to compute output frames
		for (unsigned int i = 0; i < frame_count; i++) {
			// Start with the base frame
			auto current_frame = base_frame;

			// Apply section transformations up to the frame
			for (unsigned int j = 0; j < m_frameProperties[i].get_related_beam_index_(); j++) {
				// Compose with section transformation
				std::cout<< " section : "<< j << std::endl;
				//// frame = gX(L_0)*...*gX(L_{n-1})
				current_frame = current_frame * m_section_properties[j].getTransformation();
			}

			std::cout << "Frame "<< i << "related to section : "<< m_frameProperties[i].get_related_beam_index_() <<std::endl;

			std::cout<< "curentFrame : "<< current_frame<< "\nFrameProTran : "
															""<<m_frameProperties[i] <<std::endl;

			// Apply additional frame transformation
			// frame*gX(x)
			current_frame = current_frame * m_frameProperties[i].getTransformation();

			// Convert SE3 to SOFA rigid coordinates
			const auto &translation = current_frame.translation();
			const auto &rotation = current_frame.rotation();

			std::cout << "Frame : "<<i<<" val :"<< translation.transpose()<< " " <<
				rotation << std::endl;

			// Convert rotation matrix to quaternion for SOFA
			Eigen::Quaternion<double> quat(rotation.matrix());
			sofa::type::Quat<SReal> sofa_quat(quat.x(), quat.y(), quat.z(), quat.w());
			sofa::type::Vec3 sofa_trans(translation[0], translation[1], translation[2]);

			output_frames[i] = sofa::Coord_t<Out>(sofa_trans, sofa_quat);

			if (doPrintLog) {
				msg_info() << "Frame " << i << " transformation applied";
			}
		}

		// const auto frame0 = Frame(In2::getCPos(in2[baseIndex]), In2::getCRot(in2[baseIndex]));
		//
		// // Cache the printLog value out of the loop, otherwise it will trigger a graph
		// // update at every iteration.
		// bool doPrintLog = this->f_printLog.getValue();
		// for (unsigned int i = 0; i < sz; i++) {
		// 	auto frame = frame0;
		// 	for (unsigned int u = 0; u < m_indices_vectors[i]; u++) {
		// 		frame *= m_nodes_exponential_se3_vectors[u]; // frame = gX(L_0)*...*gX(L_{n-1})
		// 	}
		// 	frame *= m_frames_exponential_se3_vectors[i]; // frame*gX(x)
		//
		// 	// This is a lazy printing approach, so there is no time consuming action in
		// 	// the core of the loop.
		// 	msg_info_when(doPrintLog) << "Frame  : " << i << " = " << frame;
		//
		// 	Vec3 origin = frame.getOrigin();
		// 	Quat orientation = frame.getOrientation();
		// 	out[i] = sofa::Coord_t<Out>(origin, orientation);
		// }

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
	void HookeSeratDiscretMapping<TIn1, TIn2, TOut>::updateFrameTransformations(
			const sofa::type::vector<Coord1> &vec_of_strains) {

		auto nb_node = vec_of_strains.size();

		// Update node properties with current strain values
		std::cout << "Begin updateFrameTransformations "<< std::endl;
		std::cout << "==================Begin Section ====================== "<< std::endl;
		for (size_t i = 0; i < nb_node; ++i) {
			// Extract strain components based on input type
			Vector6 strain = Vector6::Zero();
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
			std::cout << "DEBUG == HookeSeratDiscretMapping == > _strain : "<< strain.transpose() <<std::endl;
			// Update node info with strain values
			// i+1, since m_section_properties is 0-indexed
			m_section_properties[i+1].setStrain(strain);

			// Compute SE(3) exponential for this section
			// Change input and give as input of the function m_section_properties[i]
			SE3Types _gx = computeSE3Exponential(m_section_properties[i+1].getLength(),
				m_section_properties[i+1].getStrainsVec());

			std::cout << "Node " << i+1 << " transformation: " << _gx<< " : "<<
				_gx.translation().transpose() << std::endl;
			m_section_properties[i+1].setTransformation(_gx);
		}
		std::cout << "================== End Section ====================== "<< std::endl;

		std::cout << "================== Begin Frames ====================== "<< std::endl;
		// Update frame properties based on their position within sections
		for (size_t i = 0; i < m_frameProperties.size(); ++i) {
			if (i < m_indices_vectors.size()) {
				int sectionIndex = m_frameProperties[i].get_related_beam_index_();
				std::cout << "Frame " << i << " related to section index: " << sectionIndex << std::endl;
				if (sectionIndex >= 0 && sectionIndex < static_cast<int>(vec_of_strains.size()+1)) {
					// Compute frame transformation at its specific position
					SE3Types frame_gx =
							computeSE3Exponential(m_frameProperties[i].getDistanceToNearestBeamNode(),
								m_section_properties[sectionIndex].getStrainsVec());
					std::cout << "Frame " << i << " gX: " << frame_gx<< " : " << std::endl;
					m_frameProperties[i].setTransformation(frame_gx);
				}
			}
		}
		std::cout << "================== End Frames ====================== "<< std::endl;
	}

	template<class TIn1, class TIn2, class TOut>
	typename HookeSeratDiscretMapping<TIn1, TIn2, TOut>::SE3Types
	HookeSeratDiscretMapping<TIn1, TIn2, TOut>::computeSE3Exponential(const double section_length, const Vector6 &strain) {

		// std::cout<< "<------------------------Begin computeSE3Exponential function-------------------->"<<std::endl;
		// std::cout << "computeSE3Exponential function called with : len = " << section_length
		// 		  << " and strain = " << strain << std::endl;
		// // Extract strain vector components
		//
		//
		// if constexpr (std::is_same_v<Coord1, sofa::type::Vec3>) {
		// 	// For Vec3 input, assume curvature only
		// 	strain_.head<3>() = Vector3(strain[0], strain[1], strain[2]); //*section_length;
		// 	std::cout << "====> Strain_ : "<< strain_<< std::endl;
		// } else {
		// 	// For Vec6 input, use all components
		// 	for (int i = 0; i < 6 && i < strain.size(); ++i) {
		// 		strain_[i] = strain[i]; //*section_length;
		// 	}
		// }
		// std::cout << "Strain vector for SE(3) exponential: " << strain_.transpose() << std::endl;
		// // Compute SE(3) exponential using liegroups library
		// std::cout<< "<--------------------------End function------------------------->"<<std::endl;
		std::cout << "Input strain : " << strain.transpose() << "\n section_length : "<< section_length << std::endl;
		return SE3Types::expCosserat(strain, section_length);
		//return SE3Types::computeExp(strain*section_length);
	}

	template<class TIn1, class TIn2, class TOut>
	void
	HookeSeratDiscretMapping<TIn1, TIn2, TOut>::applyJ(const sofa::core::MechanicalParams * /* mparams */,
													   const vector<sofa::DataVecDeriv_t<Out> *> &dataVecOutVel,
													   const vector<const sofa::DataVecDeriv_t<In1> *> &dataVecIn1Vel,
													   const vector<const sofa::DataVecDeriv_t<In2> *> &dataVecIn2Vel) {

		if (dataVecOutVel.empty() || dataVecIn1Vel.empty() || dataVecIn2Vel.empty())
			return;

		if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
			return;

		if (d_debug.getValue())
			std::cout << " ########## HookeSeratDiscretMapping ApplyJ Function ########" << std::endl;

		const sofa::VecDeriv_t<In1> &strainVel = dataVecIn1Vel[0]->getValue();
		const sofa::VecDeriv_t<In2> &baseVel = dataVecIn2Vel[0]->getValue();
		sofa::VecDeriv_t<Out> &outputVel = *dataVecOutVel[0]->beginEdit();

		// Debug input velocities if enabled
		if (d_debug.getValue()) {
			displayVelocities(strainVel, baseVel, outputVel, "applyJ - input");
		}

		const auto baseIndex = d_baseIndex.getValue();
		const auto sz = d_curv_abs_frames.getValue().size();
		outputVel.resize(sz);

		// Convert base velocity to SE(3) tangent space
		Vector6 baseVelocitySE3 = Vector6::Zero();
		for (auto u = 0; u < 6; u++)
			baseVelocitySE3[u] = baseVel[baseIndex][u];

		// Compute velocities for each output frame
		for (unsigned int i = 0; i < sz; i++) {
			Vector6 frameVelocity = baseVelocitySE3;

			// Add contributions from strain velocities
			for (unsigned int u = 0; u < m_indices_vectors[i] && u < strainVel.size(); u++) {
				Vector6 strainVelSE3 = Vector6::Zero();

				if constexpr (std::is_same_v<typename sofa::Deriv_t<In1>, sofa::type::Vec3>) {
					// For Vec3 input
					strainVelSE3.head<3>() = Vector3(strainVel[u][0], strainVel[u][1], strainVel[u][2]);
				} else {
					// For Vec6 input
					for (int j = 0; j < 6 && j < strainVel[u].size(); ++j) {
						strainVelSE3[j] = strainVel[u][j];
					}
				}

				// Scale by section length
				if (u < m_section_properties.size()) {
					strainVelSE3 *= m_section_properties[u].getLength();
				}

				frameVelocity += strainVelSE3;
			}

			// Convert back to SOFA derivative format
			outputVel[i] = sofa::Deriv_t<Out>(frameVelocity.data());

			if (d_debug.getValue())
				std::cout << "Frame velocity : " << i << " = " << frameVelocity.transpose() << std::endl;
		}

		// Debug output velocities if enabled
		if (d_debug.getValue()) {
			displayVelocities(strainVel, baseVel, outputVel, "applyJ - output");
		}

		dataVecOutVel[0]->endEdit();
	}

	template<class TIn1, class TIn2, class TOut>
	void HookeSeratDiscretMapping<TIn1, TIn2, TOut>::applyJT(
			const sofa::core::MechanicalParams * /*mparams*/,
			const vector<sofa::DataVecDeriv_t<In1> *> &dataVecOut1Force,
			const vector<sofa::DataVecDeriv_t<In2> *> &dataVecOut2Force,
			const vector<const sofa::DataVecDeriv_t<Out> *> &dataVecInForce) {

		if (dataVecOut1Force.empty() || dataVecInForce.empty() || dataVecOut2Force.empty())
			return;

		if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
			return;

		if (d_debug.getValue())
			std::cout << " ########## HookeSeratDiscretMapping ApplyJT Force Function ########" << std::endl;

		const sofa::VecDeriv_t<Out> &inputForces = dataVecInForce[0]->getValue();
		sofa::VecDeriv_t<In1> &strainForces = *dataVecOut1Force[0]->beginEdit();
		sofa::VecDeriv_t<In2> &baseForces = *dataVecOut2Force[0]->beginEdit();
		const auto baseIndex = d_baseIndex.getValue();

		// Initialize output forces
		const sofa::VecCoord_t<In1> &strainState =
				m_strain_state->read(sofa::core::vec_id::read_access::position)->getValue();
		strainForces.resize(strainState.size());

		// Accumulate forces
		Vector6 totalBaseForce = Vector6::Zero();

		for (size_t i = 0; i < inputForces.size(); ++i) {
			// Convert input force to SE(3) format
			Vector6 frameForce = Vector6::Zero();
			for (auto j = 0; j < 6 && j < inputForces[i].size(); ++j) {
				frameForce[j] = inputForces[i][j];
			}

			// Add to base force (all forces ultimately go through the base)
			totalBaseForce += frameForce;

			// Distribute force to strain components
			if (i < m_indices_vectors.size()) {
				int sectionIndex = m_indices_vectors[i] - 1;
				if (sectionIndex >= 0 && sectionIndex < static_cast<int>(strainForces.size())) {
					// Project force to strain space and scale by section length
					Vector3 strainForce3D = frameForce.head<3>();
					if (sectionIndex < static_cast<int>(m_section_properties.size())) {
						strainForce3D *= m_section_properties[sectionIndex].getLength();
					}

					if constexpr (std::is_same_v<typename sofa::Deriv_t<In1>, sofa::type::Vec3>) {
						strainForces[sectionIndex] +=
								sofa::type::Vec3(strainForce3D[0], strainForce3D[1], strainForce3D[2]);
					} else {
						// For Vec6 output
						for (int k = 0; k < 6 && k < strainForces[sectionIndex].size(); ++k) {
							strainForces[sectionIndex][k] += frameForce[k];
						}
					}
				}
			}
		}

		// Set base force
		baseForces[baseIndex] += sofa::Deriv_t<In2>(totalBaseForce.data());

		if (d_debug.getValue()) {
			std::cout << "Strain forces computed" << std::endl;
			std::cout << "Base Force: " << totalBaseForce.transpose() << std::endl;
		}

		dataVecOut1Force[0]->endEdit();
		dataVecOut2Force[0]->endEdit();
	}

	template<class TIn1, class TIn2, class TOut>
	void HookeSeratDiscretMapping<TIn1, TIn2, TOut>::applyJT(
			const sofa::core::ConstraintParams * /*cparams*/,
			const vector<sofa::DataMatrixDeriv_t<In1> *> &dataMatOut1Const,
			const vector<sofa::DataMatrixDeriv_t<In2> *> &dataMatOut2Const,
			const vector<const sofa::DataMatrixDeriv_t<Out> *> &dataMatInConst) {

		if (dataMatOut1Const.empty() || dataMatOut2Const.empty() || dataMatInConst.empty())
			return;

		if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
			return;

		if (d_debug.getValue())
			std::cout << " ########## HookeSeratDiscretMapping ApplyJT Constraint Function ########" << std::endl;

		// Get constraint matrices
		sofa::MatrixDeriv_t<In1> &out1 = *dataMatOut1Const[0]->beginEdit();
		sofa::MatrixDeriv_t<In2> &out2 = *dataMatOut2Const[0]->beginEdit();
		const sofa::MatrixDeriv_t<Out> &in = dataMatInConst[0]->getValue();

		// Process constraints (simplified implementation)
		// This would need to be expanded based on specific constraint requirements
		// for (auto it = in.begin(); it != in.end(); ++it) {
		//     int constraintId = it.index();
		//     const auto &constraintLine = it.row();
		//
		//     // Apply constraint to base (simplified)
		//     auto baseIt = out2.writeLine(constraintId);
		//     baseIt.addCol(d_baseIndex.getValue(), constraintLine[0]); // Assuming first column is base force
		//
		//     // Apply constraint to strain space (simplified)
		//     auto strainIt = out1.writeLine(constraintId);
		//     for (auto colIt = constraintLine.begin(); colIt != constraintLine.end(); ++colIt) {
		//         int frameIndex = colIt.index();
		//         if (frameIndex < static_cast<int>(m_indices_vectors.size())) {
		//             int sectionIndex = m_indices_vectors[frameIndex] - 1;
		//             if (sectionIndex >= 0) {
		//                 strainIt.addCol(sectionIndex, colIt.val());
		//             }
		//         }
		//     }
		// }

		dataMatOut1Const[0]->endEdit();
		dataMatOut2Const[0]->endEdit();
	}

	template<class TIn1, class TIn2, class TOut>
	void HookeSeratDiscretMapping<TIn1, TIn2, TOut>::draw(const sofa::core::visual::VisualParams *vparams) {
		if (!vparams->displayFlags().getShowMappings())
			return;

		// Draw implementation similar to DiscreteCosseratMapping
		// This would include beam visualization with colormap
		if (d_drawMapBeam.getValue()) {
			// Draw colored beam based on deformation
			// Implementation would depend on specific visualization requirements
		}
	}

	template<class TIn1, class TIn2, class TOut>
	void HookeSeratDiscretMapping<TIn1, TIn2, TOut>::computeBBox(const sofa::core::ExecParams *params,
																 bool onlyVisible) {
		// Compute bounding box for visualization
		// Implementation would calculate the extent of all frames
		Inherit::computeBBox(params, onlyVisible);
	}

	// Debug display functions implementation
	template<class TIn1, class TIn2, class TOut>
	void HookeSeratDiscretMapping<TIn1, TIn2, TOut>::displayStrainState(const sofa::type::vector<Coord1> &strainState,
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
	void HookeSeratDiscretMapping<TIn1, TIn2, TOut>::displayRigidState(
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
	HookeSeratDiscretMapping<TIn1, TIn2, TOut>::displayOutputFrames(const sofa::type::vector<OutCoord> &outputFrames,
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
	void HookeSeratDiscretMapping<TIn1, TIn2, TOut>::displaySectionProperties(const std::string &context) const {

		std::cout << "\n=== SECTION PROPERTIES DEBUG" << (context.empty() ? "" : " (" + context + ")") << " ===\n";
		std::cout << "Section properties size: " << m_section_properties.size() << std::endl;

		for (size_t i = 0; i < m_section_properties.size(); ++i) {
			const auto &section = m_section_properties[i];
			const auto &strain = section.getStrainsVec();
			const auto &transform = section.getTransformation();

			std::cout << "  Section[" << i << "]:";
			std::cout << " length=" << section.getLength();
			std::cout << " strain=[" << strain <<  "]";
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
void HookeSeratDiscretMapping<TIn1, TIn2, TOut>::displayFrameProperties(const std::string &context) const {

	std::cout << "\n=== FRAME PROPERTIES DEBUG" << (context.empty() ? "" : " (" + context + ")") << " ===\n";
	std::cout << "Frame properties size: " << m_frameProperties.size() << std::endl;

	for (size_t i = 0; i < m_frameProperties.size(); ++i) {
		const auto &frame = m_frameProperties[i];
		const auto &transform = frame.getTransformation();

		std::cout << "  Frame[" << i << "]:";
		std::cout << " length=" << frame.getLength();
		std::cout << " frames_sect_length_=" << frame.getLength(); // Same as length, but explicitly named as requested

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
	void HookeSeratDiscretMapping<TIn1, TIn2, TOut>::displaySE3Transform(const SE3Types &transform,
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
	void HookeSeratDiscretMapping<TIn1, TIn2, TOut>::displayMappingState(const std::string &context) const {

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
	void HookeSeratDiscretMapping<TIn1, TIn2, TOut>::displayVelocities(
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
