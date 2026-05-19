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
//@appa
#pragma once

#include <Cosserat/mapping/Frames2StrainCosseratMapping.h>
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
									msg_info() << "Frames2StrainCosseratMapping updateFrames callback called";
									SOFA_UNUSED(t);
									this->updateGeometryInfo();
									std::cout << "====> Update Callback <====" << std::endl;
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


		std::cout<<"========In apply function========"<<std::endl;

		msg_info("Frames2StrainCosseratMapping") << "Frames2StrainCosseratMapping::apply called";
		
        
        if (dataVecOutPos.empty() || dataVecIn1Pos.empty() || dataVecIn2Pos.empty())
			return;

		// Check component state for validity
		if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
			return;

		if (d_debug.getValue())
			std::cout << " ########## Apply Function ########" << std::endl;
		
        // Get input data
		const sofa::VecCoord_t<In1> &frames = dataVecIn1Pos[0]->getValue(); // frames positions
		const sofa::VecCoord_t<In2> &rigidBase = dataVecIn2Pos[0]->getValue(); // Rigid base

        //Output: strain (to evaluate)
        const auto nbSections = m_section_properties.size()-1;
		sofa::VecCoord_t<Out> &strains = *dataVecOutPos[0]->beginEdit();
        strains.resize(nbSections);

		const auto baseIndex = d_baseIndex.getValue();
        
        //Get base config g(0)
		const auto &base_rigid = rigidBase[baseIndex];
		Vector3 base_trans(base_rigid.getCenter()[0], base_rigid.getCenter()[1], base_rigid.getCenter()[2]);

		// Convert SOFA quaternion to Eigen quaternion (SOFA: x,y,z,w; Eigen: w,x,y,z)
		const auto &base_sofa_quat = base_rigid.getOrientation();
		Eigen::Quaternion<double> base_rot(base_sofa_quat[3], base_sofa_quat[0], base_sofa_quat[1], base_sofa_quat[2]);

		// Create SE3 transformation
		SE3Types g_base(SE3Types::SO3Type(base_rot), base_trans);

        // g(X) = g(L_{n-1})*exp((X-L_{n-1})*strain_n)
        // for each frame (except the base frame), compute  g(L_{n-1}).inverse * g(X) = g_rel
        // then compute the strain  => strain_n = log(g_rel)/(X-L_{n-1})
    
		std::cout<<"g_base: "<<g_base<<std::endl;

		SE3Types g_prev = g_base;
		SE3Types g_curr, g_prev_inverse;
		std::cout<<"g_curr: "<<g_curr<<std::endl;

		double dx = 0.;
		TangentVector xi = TangentVector::Zero();
		for(unsigned int i=0; i<nbSections; i++){
			std::cout<<"i: "<<i<<std::endl;

			const auto& frame = frames[i+1];
			const auto& section = m_section_properties[i+1];
			std::cout<<"frame: "<<frame<<std::endl;

			Vector3 curr_translation(frame.getCenter()[0], frame.getCenter()[1], frame.getCenter()[2]);

			// Convert SOFA quaternion to Eigen quaternion (SOFA: x,y,z,w; Eigen: w,x,y,z)
			const auto &quat = frame.getOrientation();
			Eigen::Quaternion<double> curr_rotation(quat[3], quat[0], quat[1], quat[2]);

			// Create SE3 transformation
			g_curr = SE3Types(SE3Types::SO3Type(curr_rotation), curr_translation);
			g_prev_inverse = g_prev.computeInverse();
			std::cout<<"g_frame: "<<g_curr<<std::endl;
			std::cout<<"g_prev: "<<g_prev<<std::endl;
			std::cout<<"g_prev^{-1}: "<<g_prev_inverse<<std::endl;

			SE3Types g_rel = g_prev_inverse*g_curr;
			std::cout<<"g_rel: "<<g_rel<<std::endl;


			// find dx = L_{n} - L_{n-1}
			dx = section.getLength();
			xi = g_rel.computeLog()/dx; 

			std::cout<<"xi: "<< xi <<std::endl;

			for(int j=0; j<3; j++){
				strains[i][j] = xi[j];
			}

			g_prev = g_curr;
			
		}
    
		dataVecOutPos[0]->endEdit();

		std::cout<<"========End apply========"<<std::endl;
    
    }


	template<class TIn1, class TIn2, class TOut>
	Matrix3 Frames2StrainCosseratMapping<TIn1, TIn2, TOut>::buildHatMatrix(const Vector3& k){
		Matrix3 k_hat = Matrix3::Zero();

		k_hat(0, 1) = -k(2);
		k_hat(0, 2) = k(1);
		k_hat(1, 2) = -k(0);

		k_hat(1, 0) = -k_hat(0, 1);
		k_hat(2, 0) = -k_hat(0, 2);
		k_hat(2, 1) = -k_hat(1, 2);

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
		ad.template block<3,3>(1, 1) = k_hat;
		ad.template block<3,3>(1, 0) = q_hat;

		return ad;
	}

	template<class TIn1, class TIn2, class TOut>
	AdjointMatrix Frames2StrainCosseratMapping<TIn1, TIn2, TOut>::compute_Adjoint(const Matrix4& g){
		Vector3 p = Vector3::Zero();
		Matrix3 rot = Matrix3::Zero();
		// Vector3 q = Vector3::Zero(); //last components

		for(int i=0; i<3; i++){
			p(i) = g(i, 3);
			for(int j=0; j<3; j++){
				rot(i, j) = g(i, j);
			}
		}

		
		Matrix3 p_hat = buildHatMatrix(p);
		
		AdjointMatrix Ad = AdjointMatrix::Zero();

		Ad.template block<3,3>(0, 0) = rot;
		Ad.template block<3,3>(1, 1) = rot;
		Ad.template block<3,3>(1, 0) = p_hat;

		return Ad;
	}



	template<class TIn1, class TIn2, class TOut>
	AdjointMatrix Frames2StrainCosseratMapping<TIn1, TIn2, TOut>::computeInverseTangentOperator(const TangentVector& Omega){
		//computation at order 4
		//4 Bernouilli nb. : B0=1, B1=-1/2, B2=1/6, B3=0

		AdjointMatrix res = AdjointMatrix::Zero();
		AdjointMatrix Id6 = AdjointMatrix::Identity();
		AdjointMatrix adOmega = compute_adjoint(Omega); //Compute adjoint
		AdjointMatrix adOmega2 = adOmega * adOmega;
		AdjointMatrix adOmega3 = adOmega2 * adOmega;

		double B0=1., B1=-1./2., B2=1./6., B3=0.;

		res = B0*Id6 + B1*adOmega + (1./2.)*B2*adOmega2 + (1./6.)*B3*adOmega3;

		return res;

	}

	template<class TIn1, class TIn2, class TOut>	
	void 
	Frames2StrainCosseratMapping<TIn1, TIn2, TOut>::applyJ(const sofa::core::MechanicalParams *mparams,
														const sofa::type::vector<sofa::DataVecDeriv_t<Out> *> &dataVecOutVel,
														const sofa::type::vector<const sofa::DataVecDeriv_t<In1> *> &dataVecIn1Vel,
														const sofa::type::vector<const sofa::DataVecDeriv_t<In2> *> &dataVecIn2Vel){
    
		std::cout<<"========In applyJ function========"<<std::endl;
															
    	if (dataVecOutVel.empty() || dataVecIn1Vel.empty() || dataVecIn2Vel.empty())
			return;

		if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
			return;

		if (d_debug.getValue())
			std::cout << " ########## Frames2StrainCosseratMapping ApplyJ Function ########" << std::endl;
															
		const sofa::VecDeriv_t<In1> &frame_vel = dataVecIn1Vel[0]->getValue();
		const sofa::VecDeriv_t<In2> &base_vel = dataVecIn2Vel[0]->getValue();
		sofa::VecDeriv_t<Out> &strain_vel = *dataVecOutVel[0]->beginEdit();

		const sofa::VecCoord_t<Out> &strain =
				this->m_frames->read(sofa::core::vec_id::read_access::position)->getValue();
	
		const sofa::DataVecCoord_t<In1> *x1fromData =
     		 	m_strain_state->read(sofa::core::vec_id::read_access::position);
		const sofa::VecCoord_t<In1> framePositions = x1fromData->getValue();				

		const auto base_index = d_baseIndex.getValue();
		const auto section_count = d_curv_abs_section.getValue().size() - 1; //Pourquoi d_curv_abs_section donne le nombre de section + 1 (@appa: enlever 1)
		std::cout<<"Nb. section: "<<section_count<<std::endl;

		strain_vel.resize(section_count);
		for (auto &vel : strain_vel){
			vel.clear();
		}


		// 2. Compute the base velocity in SE(3) tangent space
		// 2.1 Convert base velocity to se(3) tangent vector
		TangentVector base_vel_local = TangentVector::Zero();
		for (auto u = 0; u < 6; u++)
			base_vel_local[u] = base_vel[base_index][u];

		// std::cout<<"Strain : [ ";
		// for(auto v : strain)
		// 	std::cout<<"["<<v<<"] ";
		// std::cout<<"]"<<std::endl;


		// std::cout<<"Frame Positions : [ ";
		// for(auto v : framePositions)
		// 	std::cout<<"["<<v<<"] ";
		// std::cout<<"]"<<std::endl;

		// std::cout<<"Base Velocity: ["<<base_vel_local.transpose()<<"]"<<std::endl;

		// std::cout<<"Frames Velocity: [ ";
		// for(auto v : frame_vel)
		// 	std::cout<<"["<<v<<"] ";
		// std::cout<<"]"<<std::endl;


		// std::cout<<"Strain Velocity: [ ";
		// for(auto v : strain_vel)
		// 	std::cout<<"["<<v<<"] ";
		// std::cout<<"]"<<std::endl;		

		//À COMPLÉTER
		// compute the Jacobians J1 and J2
		// Omega = Log(ga^-1gb
		// J1 = 1/h dexp^-1_{-Omega)} Ad_{exp(-Omega)}
		// J2 = 1/h dexp^-1_{Omega}
		AdjointMatrix J1 = AdjointMatrix::Zero();
		AdjointMatrix J2 = AdjointMatrix::Zero();

		TangentVector strain_i = TangentVector::Zero();
		TangentVector Omega_i = TangentVector::Zero();
		// strain_i(3) = 1.;
		double dx = 0.;
		AdjointMatrix dexp_inv = AdjointMatrix::Zero();
		for(int i=0; i<section_count; i++){
			const auto &section = m_section_properties[i+1];
			dx = section.getLength();

			for(int j=0; j<3; j++)
				strain_i[j] = strain[i][j];

			Omega_i = dx*strain_i;
			// std::cout<<"Omega["<<i<<"]: "<<Omega_i.transpose()<<std::endl;
			// std::cout<<"dx["<<i<<"]: "<<dx<<std::endl;
			dexp_inv = computeInverseTangentOperator(Omega_i);
			J2 = (1./dx)*dexp_inv;

			SE3Types g = SE3Types::expCosserat(strain_i, -dx); // = exp(-Omega_i)
			AdjointMatrix Adg = g.computeAdjoint();
			J1 = -(1./dx)*dexp_inv*Adg;

			TangentVector eta_a = TangentVector::Zero();
			TangentVector eta_b = TangentVector::Zero();

			if(i == 0){
				eta_a = base_vel_local;
				for(int u=0; u<6; u++){
					eta_b[u] = frame_vel[i][u];
				}
			} else{
				for(int u=0; u<6; u++){
					eta_a[u] = frame_vel[i][u];
					eta_b[u] = frame_vel[i+1][u];
				}
			}

			TangentVector output_vel = J1*eta_a + J2*eta_b;

			for(int k=0; k<3; k++){
				strain_vel[i][k] = output_vel[k];
			}

			std::cout << "Strain velocity [" << i << "]: " << output_vel.transpose() <<"\n";
		}

		dataVecOutVel[0]->endEdit();
		std::cout<<"========End applyJ========"<<std::endl;

	}

	template<class TIn1, class TIn2, class TOut>
	void 
	Frames2StrainCosseratMapping<TIn1, TIn2, TOut>::applyJT(const sofa::core::MechanicalParams *mparams,
					 										const sofa::type::vector<sofa::DataVecDeriv_t<In1> *> &dataVecOut1Force,
					 										const sofa::type::vector<sofa::DataVecDeriv_t<In2> *> &dataVecOut2RootForce,
					 										const sofa::type::vector<const sofa::DataVecDeriv_t<Out> *> &dataVecInForce){

	}


	template<class TIn1, class TIn2, class TOut>
	void 
	Frames2StrainCosseratMapping<TIn1, TIn2, TOut>::applyJT(const sofa::core::ConstraintParams *cparams,
															const sofa::type::vector<sofa::DataMatrixDeriv_t<In1> *> &dataMatOut1Const,
					 										const sofa::type::vector<sofa::DataMatrixDeriv_t<In2> *> &dataMatOut2Const,
					 										const sofa::type::vector<const sofa::DataMatrixDeriv_t<Out> *> &dataMatInConst){

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

		// // Debug output if needed
		// if (this->f_printLog.getValue()) {
		// 	displayOutputFrames(xData, "draw - rendering frames");
		// }

		glEnd();				
	}

}