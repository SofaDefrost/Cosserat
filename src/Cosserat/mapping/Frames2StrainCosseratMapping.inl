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
 *******************************************************************************
 * Authors: The SOFA Team and external contributors (see Authors.txt)          *
 *                                                                             *
 * Contact information: contact@sofa-framework.org                             *
 ******************************************************************************/
#pragma once

#include <Cosserat/mapping/Frames2StrainCosseratMapping.h>
#include <Cosserat/mapping/SofaLieGroupsUtils.h>
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
    Frames2StrainCosseratMapping<TIn1, TIn2, TOut>::Frames2StrainCosseratMapping():
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
							 "by another body.")) 
        
    {
        this->addUpdateCallback("updateFrames", {&d_curv_abs_section, &d_curv_abs_frames, &d_debug},
								[this](const sofa::core::DataTracker &t) {
									SOFA_UNUSED(t);
									msg_info() << "Frames2StrainCosseratMapping updateFrames callback called";
									this->updateGeometryInfo();
									msg_info_when(d_debug.getValue())
										<< "====> Update Callback <====";
									return sofa::core::objectmodel::ComponentState::Valid;
								},
								{});
	}        

	template<class TIn1, class TIn2, class TOut>
	void Frames2StrainCosseratMapping<TIn1, TIn2, TOut>::doBaseCosseratInit() {
		// Initialize colormap for visualization
		m_colorMap.setColorScheme("Blue to Red");
		m_colorMap.reinit();

		msg_info() << "Frames2StrainCosseratMapping initialized";
	}    

	template<class TIn1, class TIn2, class TOut>
	void Frames2StrainCosseratMapping<TIn1, TIn2, TOut>::initialization(){

		if(m_strain_state){
			auto xfromData = m_strain_state->read(sofa::core::vec_id::read_access::position);
			const auto &xfrom = xfromData->getValue();

			// Initialize frame properties using the initial frame states
			const auto frame_count = xfrom.size();

			m_frameProperties.clear();
			m_frameProperties.reserve(frame_count);

			for(size_t i=0; i<frame_count; i++){
				m_frameProperties.emplace_back();
			}
		}

	}


	template<class TIn1, class TIn2, class TOut>
	void
	Frames2StrainCosseratMapping<TIn1, TIn2, TOut>::apply(const sofa::core::MechanicalParams * /* mparams */,
													  const vector<sofa::DataVecCoord_t<Out> *> &dataVecOutPos,
													  const vector<const sofa::DataVecCoord_t<In1> *> &dataVecIn1Pos,
													  const vector<const sofa::DataVecCoord_t<In2> *> &dataVecIn2Pos) {



		msg_info("Frames2StrainCosseratMapping") << "Frames2StrainCosseratMapping::apply called";
		
        
        if (dataVecOutPos.empty() || dataVecIn1Pos.empty() || dataVecIn2Pos.empty())
			return;

		// Check component state for validity
		if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
			return;

		msg_info_when(d_debug.getValue()) << " ########## Apply Function ########";
		
        // Get input data
		const sofa::VecCoord_t<In1> &frames = dataVecIn1Pos[0]->getValue(); // frames positions
		const sofa::VecCoord_t<In2> &rigidBase = dataVecIn2Pos[0]->getValue(); // Rigid base

        //Output: strain (to evaluate)
        const auto nbSections = m_section_properties.size()-1;
		sofa::VecCoord_t<Out> &strains = *dataVecOutPos[0]->beginEdit();
        strains.resize(nbSections);
		for(auto &s : strains){
			s.clear();
		}

		const auto baseIndex = d_baseIndex.getValue();
        
        //Get base config g(0)  (via SOFA<->Eigen helper)
		const SE3Types g_base = rigidCoordToSE3(rigidBase[baseIndex]);

        // g(X) = g(L_{n-1})*exp((X-L_{n-1})*strain_n)
        // for each frame (except the base frame), compute  g(L_{n-1}).inverse * g(X) = g_rel
        // then compute the strain  => strain_n = log(g_rel)/(X-L_{n-1})
    
		const auto &curv_abs_sections = d_curv_abs_section.getValue();
		const auto &curv_abs_frames = d_curv_abs_frames.getValue();

		//compute transformation of each frame (via SOFA<->Eigen helper)
		std::vector<SE3Types> g_frames(frames.size());
		for (unsigned int i = 0; i < frames.size(); i++) {
			g_frames[i] = rigidCoordToSE3(frames[i]);
		}

		for(unsigned int i=0; i<nbSections; i++){
			double section_start = curv_abs_sections[i];
			double section_end = curv_abs_sections[i+1];

			// Trouver l'indice des frames correspondants
			int left_frame_idx = -1, right_frame_idx = -1;

			for(size_t j=0; j<curv_abs_frames.size(); j++){
				if (std::abs(curv_abs_frames[j] - section_start) <1e-12){
					left_frame_idx = j;
				}
				if (std::abs(curv_abs_frames[j] - section_end) <1e-12){
					right_frame_idx = j;
				}
			}

			//si les deux frames sont trouvés, calculer le strain
			if(left_frame_idx >=0 && right_frame_idx >=0){
				double dx = section_end - section_start; //section length

				SE3Types g_rel = g_frames[left_frame_idx].computeInverse() * g_frames[right_frame_idx];

				TangentVector xi = g_rel.computeLog()/dx;

				for(int j=0; j<6; j++){
					strains[i][j] = xi[j];
				}
			}
			else{
				msg_warning() << "Could not find frames for section " << i;
			}
		}

		dataVecOutPos[0]->endEdit();
    
    }


	template<class TIn1, class TIn2, class TOut>
	Matrix3 Frames2StrainCosseratMapping<TIn1, TIn2, TOut>::buildHatMatrix(const Vector3& k){
		Matrix3 k_hat;
    	k_hat<< 0., -k(2),  k(1),
              	k(2), 0., -k(0),
             	-k(1), k(0), 0.;

		return k_hat;

	}
	template<class TIn1, class TIn2, class TOut>
	AdjointMatrix Frames2StrainCosseratMapping<TIn1, TIn2, TOut>::compute_adjoint(const TangentVector& Omega){
		Vector3 k = Vector3::Zero(); //first 3 components of Omega
		Vector3 q = Vector3::Zero(); //last components

		for(int i=0; i<3; i++){
			k(i) = Omega(i);
			q(i) = Omega(i+3);
		}

		Matrix3 k_hat = buildHatMatrix(k);
		Matrix3 q_hat = buildHatMatrix(q);
		
		AdjointMatrix ad = AdjointMatrix::Zero();

		ad.template block<3,3>(0, 0) = k_hat;
		ad.template block<3,3>(3, 3) = k_hat;
		ad.template block<3,3>(3, 0) = q_hat;

		return ad;
	}

	template<class TIn1, class TIn2, class TOut>
	AdjointMatrix Frames2StrainCosseratMapping<TIn1, TIn2, TOut>::computeInverseTangentOperator(const TangentVector& Omega){
		//computation at order 5
		//5 Bernouilli nb. : B0=1, B1=-1/2, B2=1/6, B3=0, B4=-1/30

		AdjointMatrix res = AdjointMatrix::Zero();
		
		const Vector3 phi = Omega.template head<3>();
		double theta = phi.norm();
		AdjointMatrix Id6 = AdjointMatrix::Identity();
		AdjointMatrix adOmega = compute_adjoint(Omega); //Compute adjoint
		AdjointMatrix adOmega2 = adOmega * adOmega;
		AdjointMatrix adOmega4 = adOmega2 * adOmega2;		
		
		double B0=1., B1=-1./2., B2=1./6., B4=-1./30.;

		
		if (theta < 1e-4){ //pour de petite deformation

			res = B0*Id6 + B1*adOmega + (1./2.)*B2*adOmega2 + (1./24.)*B4*adOmega4;
		}
		else{
			double cot_half = 1.0 / std::tan(theta / 2.0);
			double c4 = ((theta/2.) * cot_half - 1. + (theta*theta)/12.)/std::pow(theta, 4);
			
			res = B0*Id6 + B1*adOmega + (1./2.)*B2*adOmega2 + c4*adOmega4; 
		}

		return res;

	}

	template<class TIn1, class TIn2, class TOut>	
	void 
	Frames2StrainCosseratMapping<TIn1, TIn2, TOut>::applyJ(const sofa::core::MechanicalParams *mparams,
														const sofa::type::vector<sofa::DataVecDeriv_t<Out> *> &dataVecOutVel,
														const sofa::type::vector<const sofa::DataVecDeriv_t<In1> *> &dataVecIn1Vel,
														const sofa::type::vector<const sofa::DataVecDeriv_t<In2> *> &dataVecIn2Vel){
    															
    	if (dataVecOutVel.empty() || dataVecIn1Vel.empty() || dataVecIn2Vel.empty())
			return;

		if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
			return;

		msg_info_when(d_debug.getValue())
			<< " ########## Frames2StrainCosseratMapping ApplyJ Function ########";
															
		const sofa::VecDeriv_t<In1> &frame_vel = dataVecIn1Vel[0]->getValue();
		const sofa::VecDeriv_t<In2> &base_vel = dataVecIn2Vel[0]->getValue();
		sofa::VecDeriv_t<Out> &strain_vel = *dataVecOutVel[0]->beginEdit();

		const sofa::VecCoord_t<Out> &strain =
				this->m_frames->read(sofa::core::vec_id::read_access::position)->getValue();
	
		const sofa::DataVecCoord_t<In1> *x1fromData =
     		 	m_strain_state->read(sofa::core::vec_id::read_access::position);
		const sofa::VecCoord_t<In1> framePositions = x1fromData->getValue();				

		const auto base_index = d_baseIndex.getValue();
		const auto section_count = d_curv_abs_section.getValue().size() - 1;

		strain_vel.resize(section_count);
		for (auto &vel : strain_vel){
			vel.clear();
		}
										


		// Compute the base velocity in SE(3) tangent space
		//    Convert base velocity to se(3) tangent vector
		TangentVector base_vel_local = TangentVector::Zero();
		for (auto u = 0; u < 6; u++)
			base_vel_local[u] = base_vel[base_index][u];

		//compute transformation of each frame (via SOFA<->Eigen helper)
		std::vector<SE3Types> g_frames(framePositions.size());
		for (unsigned int i = 0; i < framePositions.size(); i++) {
			g_frames[i] = rigidCoordToSE3(framePositions[i]);
		}

		//
		// compute the Jacobians J1 and J2
		// Omega = Log(ga^-1gb
		// J1 = 1/h dexp^-1_{-Omega)} Ad_{exp(-Omega)}
		// J2 = 1/h dexp^-1_{Omega}

		for(int i=0; i<section_count; i++){
			const auto &section = m_section_properties[i+1];
			double dx = section.getLength();

			TangentVector strain_i = TangentVector::Zero();

			for(int j=0; j<6; j++)
				strain_i[j] = strain[i][j];

			TangentVector Omega_i = dx*strain_i;
			AdjointMatrix dexp_inv = computeInverseTangentOperator(Omega_i);
			AdjointMatrix J2 = (1./dx)*dexp_inv;

			SE3Types g = SE3Types::expCosserat(strain_i, -dx); // = exp(-Omega_i)
			AdjointMatrix Adg = g.computeAdjoint();
			AdjointMatrix J1 = -(1./dx)*dexp_inv*Adg;

			TangentVector eta_a = TangentVector::Zero();
			TangentVector eta_b = TangentVector::Zero();
		

			//Projection (global -> local)
			//frame a
			SE3Types ga_inv = g_frames[i].inverse();
			AdjointMatrix a_projector = ga_inv.buildProjectionMatrix(ga_inv.rotation().matrix());
			TangentVector vela_global(frame_vel[i][0], frame_vel[i][1], frame_vel[i][2], frame_vel[i][3], frame_vel[i][4], frame_vel[i][5]);
			
			eta_a = a_projector * vela_global;

			// idem pour le frame b
			SE3Types gb_inv = g_frames[i+1].inverse();
			AdjointMatrix b_projector = gb_inv.buildProjectionMatrix(gb_inv.rotation().matrix());
			TangentVector velb_global(frame_vel[i+1][0], frame_vel[i+1][1], frame_vel[i+1][2], frame_vel[i+1][3], frame_vel[i+1][4], frame_vel[i+1][5]);
			
			eta_b = b_projector * velb_global;
			
			TangentVector output_vel = J1*eta_a + J2*eta_b;
			
			for(int k=0; k<6; k++){
				strain_vel[i][k] = output_vel[k];
			}
		}

		dataVecOutVel[0]->endEdit();

	}




	template<class TIn1, class TIn2, class TOut>
	void 
	Frames2StrainCosseratMapping<TIn1, TIn2, TOut>::applyJT(const sofa::core::MechanicalParams * /*mparams*/,
					 										const sofa::type::vector<sofa::DataVecDeriv_t<In1> *> &dataVecOut1Force,
					 										const sofa::type::vector<sofa::DataVecDeriv_t<In2> *> &dataVecOut2Force,
					 										const sofa::type::vector<const sofa::DataVecDeriv_t<Out> *> &dataVecInForce){

		if (dataVecOut1Force.empty() || dataVecInForce.empty() || dataVecOut2Force.empty())
			return;

		if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
			return;

		msg_info_when(d_debug.getValue())
			<< " ########## Frames2StrainCosseratMapping ApplyJT Force Function ########";



		const sofa::VecDeriv_t<Out> &strainForces = dataVecInForce[0]->getValue();
		sofa::VecDeriv_t<In1> &frameForces = *dataVecOut1Force[0]->beginEdit();
		sofa::VecDeriv_t<In2> &baseForces = *dataVecOut2Force[0]->beginEdit();
		const auto baseIndex = d_baseIndex.getValue();

		// Get current strain of each section
		const sofa::VecCoord_t<Out> &strainState =
				this->m_frames->read(sofa::core::vec_id::read_access::position)->getValue();

		const sofa::DataVecCoord_t<In1> *x1fromData =
     		 	m_strain_state->read(sofa::core::vec_id::read_access::position);
		const sofa::VecCoord_t<In1> framePositions = x1fromData->getValue();	

		
		// Initialize output forces
		frameForces.resize(framePositions.size());

		//compute transformation of each frame (via SOFA<->Eigen helper)
		std::vector<SE3Types> g_frames(framePositions.size());
		for (unsigned int i = 0; i < framePositions.size(); i++) {
			g_frames[i] = rigidCoordToSE3(framePositions[i]);
		}

		const auto section_count = d_curv_abs_section.getValue().size() - 1;

		for(unsigned int i=0; i<section_count; i++){
			const auto& section = m_section_properties[i+1];

			double dx = section.getLength(); //section length

			//get current strain of the section
			TangentVector strain_i = TangentVector::Zero();
			TangentVector lambda = TangentVector::Zero();

			for(int j=0; j<6; j++){
				strain_i[j] = strainState[i][j];
				lambda[j] = strainForces[i][j];
			}

			TangentVector Omega_i = dx * strain_i;

			//compute Jacobians
			AdjointMatrix dexp_inv = computeInverseTangentOperator(Omega_i);
			AdjointMatrix J2 = (1./dx)*dexp_inv;

			SE3Types g = SE3Types::expCosserat(strain_i, -dx); // = exp(-Omega_i)
			AdjointMatrix AdgT = g.computeAdjoint().transpose();
			AdjointMatrix J1_transpose = -(1./dx) * AdgT * dexp_inv.transpose();

			// NB : le facteur dx est déjà appliqué dans le ForceField, on ne le
			// remultiplie pas ici (cf. J1_transpose ci-dessus qui contient déjà -1/dx)
			TangentVector fa_local = J1_transpose * lambda; //a (b): extremite gauche (droite) de la section
			TangentVector fb_local = J2.transpose() * lambda;

			//Projection (local -> global)
			//frame a
			SE3Types ga = g_frames[i];
			AdjointMatrix a_projector = ga.buildProjectionMatrix(ga.rotation().matrix());
			TangentVector fa_global = a_projector.transpose().inverse() * fa_local;

			// idem pour le frame b
			SE3Types gb = g_frames[i+1];
			AdjointMatrix b_projector = gb.buildProjectionMatrix(gb.rotation().matrix());
			TangentVector fb_global = b_projector.transpose().inverse() * fb_local;

			for(int k=0; k<6; k++){
				frameForces[i][k] +=fa_global[k];
				frameForces[i+1][k] +=fb_global[k];
			}


		}

		dataVecOut1Force[0]->endEdit();
		dataVecOut2Force[0]->endEdit();
					
	}


	template<class TIn1, class TIn2, class TOut>
	void 
	Frames2StrainCosseratMapping<TIn1, TIn2, TOut>::applyJT(const sofa::core::ConstraintParams * /*cparams*/,
															const sofa::type::vector<sofa::DataMatrixDeriv_t<In1> *> &dataMatOut1Const,
					 										const sofa::type::vector<sofa::DataMatrixDeriv_t<In2> *> &dataMatOut2Const,
					 										const sofa::type::vector<const sofa::DataMatrixDeriv_t<Out> *> &dataMatInConst){

		
		if (dataMatOut1Const.empty() || dataMatOut2Const.empty() || dataMatInConst.empty())
			return;

		if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
			return;

		msg_info_when(d_debug.getValue())
			<< " ########## Frames2StrainCosseratMapping ApplyJT Constraint Function ########";

		
		// Prepare input and output data matrices
		sofa::MatrixDeriv_t<In1> &out1 = *dataMatOut1Const[0]->beginEdit(); // frames
		sofa::MatrixDeriv_t<In2> &out2 = *dataMatOut2Const[0]->beginEdit(); // rigid base
		const sofa::MatrixDeriv_t<Out> &in = dataMatInConst[0]->getValue(); // strain

		// Get current positions to compute transformations
		const sofa::DataVecCoord_t<In1> *x1fromData =
		m_strain_state->read(sofa::core::vec_id::read_access::position);
		const sofa::VecCoord_t<In1> framePositions = x1fromData->getValue(); 

		// Get current strain of each section
		const sofa::VecCoord_t<Out> &strainState =
		this->m_frames->read(sofa::core::vec_id::read_access::position)->getValue();


		//compute transformation of each frame (via SOFA<->Eigen helper)
		std::vector<SE3Types> g_frames(framePositions.size());
		for (unsigned int i = 0; i < framePositions.size(); i++) {
			g_frames[i] = rigidCoordToSE3(framePositions[i]);
		}

		// Process constraints
		for(auto rowIt = in.begin(); rowIt != in.end(); ++rowIt){
		
			auto colIt = rowIt.begin();
			if (colIt == rowIt.end())
			continue; 

			typename sofa::MatrixDeriv_t<In1>::RowIterator o1 = out1.writeLine(rowIt.index());
			// typename sofa::MatrixDeriv_t<In2>::RowIterator o2 = out2.writeLine(rowIt.inde));

			while(colIt != rowIt.end()){
				int strainIndex = colIt.index(); //quel strain est constrain

				const auto& section = m_section_properties[strainIndex+1];

				TangentVector strain = TangentVector::Zero(); //strain
				TangentVector constraintValue = TangentVector::Zero();


				// Convert constraint value to TangentVector
				const sofa::Deriv_t<Out> val = colIt.val();
				for (unsigned int j = 0; j < 6; ++j) {
					strain[j] = strainState[strainIndex][j];
					constraintValue[j] = val[j];
				}

				//computation of the jacobians
				double dx = section.getLength();
				TangentVector Omega = dx * strain;


				//compute Jacobians
				AdjointMatrix dexp_inv = computeInverseTangentOperator(Omega);
				AdjointMatrix J2 = (1./dx)*dexp_inv;

				SE3Types g = SE3Types::expCosserat(strain, -dx); // = exp(-Omega)
				AdjointMatrix AdgT = g.computeAdjoint().transpose();
				AdjointMatrix J1_transpose = -(1./dx) * AdgT * dexp_inv.transpose();

				TangentVector fa_local = J1_transpose * constraintValue; //a (b): extremite gauche (droite) de la section
				TangentVector fb_local = J2.transpose() * constraintValue;


				//Projection (local -> global)
				//frame a
				SE3Types ga = g_frames[strainIndex];
				AdjointMatrix a_projector = ga.buildProjectionMatrix(ga.rotation().matrix());
				TangentVector fa_global = a_projector.transpose().inverse() * fa_local;
				

				// idem pour le frame b
				SE3Types gb = g_frames[strainIndex+1];
				AdjointMatrix b_projector = gb.buildProjectionMatrix(gb.rotation().matrix());
				TangentVector fb_global = b_projector.transpose().inverse() * fb_local;

				sofa::type::Vec<6, double> fa_vec, fb_vec;

				for(int k=0; k<6; k++){ 
					fa_vec[k] = fa_global[k]; 
					fb_vec[k] = fb_global[k];
				}
				
				o1.addCol(strainIndex, fa_vec); // Impact sur Frame A
				o1.addCol(strainIndex + 1, fb_vec); // Impact sur Frame B
				
				colIt++;
			}

		}

		dataMatOut1Const[0]->endEdit();
		dataMatOut2Const[0]->endEdit();														

	}
		

	template<class TIn1, class TIn2, class TOut>
	void 
	Frames2StrainCosseratMapping<TIn1, TIn2, TOut>::computeBBox(const sofa::core::ExecParams *params,
																bool onlyVisible) {
		
		// Compute bounding box for visualization
		// Implementation would calculate the extent of all frames
		Inherit::computeBBox(params, onlyVisible);																	

	}


	template<class TIn1, class TIn2, class TOut>
	void Frames2StrainCosseratMapping<TIn1, TIn2, TOut>::draw(const sofa::core::visual::VisualParams *vparams) {
		if(!vparams->displayFlags().getShowMechanicalMappings())
			return;

		// draw cable
		typedef sofa::type::RGBAColor RGBAColor;

		const auto stateLifeCycle = vparams->drawTool()->makeStateLifeCycle();

		//Strain (Output)
		const::sofa::DataVecCoord_t<Out> *artiData = this->m_frames->read(sofa::core::vec_id::read_access::position);
		const::sofa::VecCoord_t<Out> xPos = artiData->getValue();

		//Frames (In1)
		const sofa::DataVecCoord_t<In1> *xfromData = this->m_strain_state->read(sofa::core::vec_id::read_access::position);
		const sofa::VecCoord_t<In1> xData = xfromData->getValue();
		vector<sofa::type::Vec3> positions;
		vector<sofa::type::Quat<SReal>> Orientation;
		positions.clear();
		Orientation.clear();
		unsigned int sz = xData.size();
		for (unsigned int i = 0; i < sz; i++) {
			positions.push_back(xData[i].getCenter());
			Orientation.push_back(xData[i].getOrientation());
		}

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
				j = m_indices_vectors[i] - 1; // to get the articulation on which the frame is related to
				RGBAColor color = _eval(xPos[j][d_deformationAxis.getValue()]);
				vparams->drawTool()->drawLine(positions[i], positions[i + 1], color);
			}
		}
		glLineWidth(1);
		if (!vparams->displayFlags().getShowMappings())
			if (!d_debug.getValue())
				return;

		// // Debug output if needed
		// if (this->f_printLog.getValue()) {
		// 	displayOutputFrames(xData, "draw - rendering frames");
		// }

		glEnd();
	}


	// applyDJT — geometric stiffness K_G · δξ.
	// NOT YET IMPLEMENTED. The stub satisfies the vtable so the component loads,
	// but with an implicit solver (EulerImplicitSolver, StaticSolver, …) the
	// tangent stiffness will be incomplete → Newton's method will not converge
	// quadratically on large deformations.  Mirror of the same stub already in
	// Strain2RigidCosseratMapping — see the analogous warning policy there.
	template<class TIn1, class TIn2, class TOut>
	void Frames2StrainCosseratMapping<TIn1, TIn2, TOut>::applyDJT(
			const sofa::core::MechanicalParams* /*mparams*/,
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

}