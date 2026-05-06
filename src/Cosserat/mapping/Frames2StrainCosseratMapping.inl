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

		if(m_frames){
			auto xfromData = m_frames->read(sofa::core::vec_id::read_access::position);
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
		SE3Types g_curr;
		std::cout<<"g_curr: "<<g_curr<<std::endl;

		Vector3 diff = Vector3::Zero();
		double dx = 0.;
		TangentVector xi = TangentVector::Zero();

		for(unsigned int i=0; i<nbSections; i++){
			std::cout<<"i: "<<i<<std::endl;

			const auto& frame = frames[i+1];
			std::cout<<"frame: "<<frame<<std::endl;

			Vector3 curr_translation(frame.getCenter()[0], frame.getCenter()[1], frame.getCenter()[2]);

			Vector3 prev_translation(frames[i].getCenter()[0], frames[i].getCenter()[1], frames[i].getCenter()[2]);

			// Convert SOFA quaternion to Eigen quaternion (SOFA: x,y,z,w; Eigen: w,x,y,z)
			const auto &quat = frame.getOrientation();
			Eigen::Quaternion<double> curr_rotation(quat[3], quat[0], quat[1], quat[2]);

			// Create SE3 transformation
			g_curr = SE3Types(SE3Types::SO3Type(curr_rotation), curr_translation);
			std::cout<<"g_frame: "<<g_curr<<std::endl;
			std::cout<<"g_prev: "<<g_prev<<std::endl;
			std::cout<<"g_prev^{-1}: "<<g_prev.computeInverse()<<std::endl;

			SE3Types g_rel = g_prev.computeInverse()*g_curr;
			std::cout<<"g_rel: "<<g_rel<<std::endl;


			// find dx = L_{n} - L_{n-1}
			diff = curr_translation - prev_translation;
			dx = diff.norm();
			xi = g_rel.computeLog()/dx; 

			std::cout<<"xi: "<< xi <<std::endl;

			for(int j=0; j<3; j++){
				strains[i][j] = xi[j];
			}

			g_prev = g_curr;
			
		}
    
		dataVecOutPos[0]->endEdit();
    
    }

	template<class TIn1, class TIn2, class TOut>	
	void 
	Frames2StrainCosseratMapping<TIn1, TIn2, TOut>::applyJ(const sofa::core::MechanicalParams *mparams,
														const sofa::type::vector<sofa::DataVecDeriv_t<Out> *> &dataVecOutVel,
														const sofa::type::vector<const sofa::DataVecDeriv_t<In1> *> &dataVecIn1Vel,
														const sofa::type::vector<const sofa::DataVecDeriv_t<In2> *> &dataVecIn2Vel){

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
	void Frames2StrainCosseratMapping<TIn1, TIn2, TOut>::draw(const sofa::core::visual::VisualParams *vparams) {}

}