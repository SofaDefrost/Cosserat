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
#include "DiscreteDynamicCosseratMapping.h"

#include <sofa/core/Multi2Mapping.inl>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/helper/logging/Message.h>
#include <sofa/type/Quat.h>

#include <string>


namespace sofa::component::mapping
{
using sofa::core::objectmodel::BaseContext ;
using sofa::helper::AdvancedTimer;
using sofa::helper::WriteAccessor;


template <class TIn1, class TIn2, class TOut>
DiscreteDynamicCosseratMapping<TIn1, TIn2, TOut>::DiscreteDynamicCosseratMapping()
	: m_fromModel1(NULL)
	, m_fromModel2(NULL)
	, m_toModel(NULL)
{
}


// _________________________________________________________________________________________

template <class TIn1, class TIn2, class TOut>
void DiscreteDynamicCosseratMapping<TIn1, TIn2, TOut>::init()
{
    Inherit1::init();

	if(this->getFromModels1().empty())
	{
		msg_error() << "Error while initializing ; input getFromModels1 not found" ;
		return;
	}

	if(this->getFromModels2().empty())
	{
		msg_error() << "Error while initializing ; output getFromModels2 not found" ;
		return;
	}

	if(this->getToModels().empty())
	{
		msg_error() << "Error while initializing ; output Model not found" ;
		return;
	}

	printf("=================================> Init from the DiscretDynamicCosseratMapping component \n");

	m_fromModel1 = this->getFromModels1()[0];
	m_fromModel2 = this->getFromModels2()[0];
	m_toModel = this->getToModels()[0];

	// Fill the initial vector
	const OutDataVecCoord* xfromData = m_toModel->read(core::ConstVecCoordId::position());
	const OutVecCoord xfrom = xfromData->getValue();
	//    WriteAccessor<Data < helper::vector<double>>> curv_abs_output = d_curv_abs_frames;
	//    curv_abs_output.clear();

	m_vecTransform.clear();
	for (unsigned int i = 0; i < xfrom.size(); i++) {
		m_vecTransform.push_back(xfrom[i]);
	}

	//initialize the constant matrix m_matrixBi[0][0] = 1.0;
	for(size_t i = 0 ; i < 3; i++) m_matrixBi[i][i] = 1.0;

	this->initialize();
}


template <class TIn1, class TIn2, class TOut>
void DiscreteDynamicCosseratMapping<TIn1, TIn2, TOut>::bwdInit()
{

}

template <class TIn1, class TIn2, class TOut>
void DiscreteDynamicCosseratMapping<TIn1, TIn2, TOut>::reinit()
{

}

template <class TIn1, class TIn2, class TOut>
void DiscreteDynamicCosseratMapping<TIn1, TIn2, TOut>::reset()
{
	reinit();
}



template <class TIn1, class TIn2, class TOut>
void DiscreteDynamicCosseratMapping<TIn1, TIn2, TOut>::apply(
		const core::MechanicalParams* /* mparams */, const type::vector<OutDataVecCoord*>& dataVecOutPos,
		const type::vector<const In1DataVecCoord*>& dataVecIn1Pos ,
		const type::vector<const In2DataVecCoord*>& dataVecIn2Pos)
{

	if(dataVecOutPos.empty() || dataVecIn1Pos.empty() || dataVecIn2Pos.empty())
		return;

	///Do Apply
	//We need only one input In model and input Root model (if present)
	const In1VecCoord& in1 = dataVecIn1Pos[0]->getValue();
	const In2VecCoord& in2 = dataVecIn2Pos[0]->getValue();

	size_t sz = d_curv_abs_frames.getValue().size();
	OutVecCoord& out = *dataVecOutPos[0]->beginEdit();
	out.resize(sz);

	//update the Exponential Matrices according to new deformation
	this->update_ExponentialSE3(in1); // ==> update m_framesExponentialSE3Vectors & m_nodesExponentialSE3Vectors

	Transform frame0 = Transform(In2::getCPos(in2[0]),In2::getCRot(in2[0]));
	for(unsigned int i=0; i<sz; i++){
		Transform frame = frame0;
		for (unsigned int u = 0; u < m_indicesVectors[i]; u++) {
			frame *= m_nodesExponentialSE3Vectors[u];
		}
		frame *= m_framesExponentialSE3Vectors[i];

		type::Vec3 v = frame.getOrigin();
		type::Quat q = frame.getOrientation();
		out[i] = OutCoord(v,q);
	}

	m_index_input = 0;
	dataVecOutPos[0]->endEdit();
}




template <class TIn1, class TIn2, class TOut>
void DiscreteDynamicCosseratMapping<TIn1, TIn2, TOut>:: applyJ(
		const core::MechanicalParams* /* mparams */, const type::vector< OutDataVecDeriv*>& dataVecOutVel,
		const type::vector<const In1DataVecDeriv*>& dataVecIn1Vel,
		const type::vector<const In2DataVecDeriv*>& dataVecIn2Vel) {

	if(dataVecOutVel.empty() || dataVecIn1Vel.empty() ||dataVecIn2Vel.empty() )
		return;
	const In1VecDeriv& in1 = dataVecIn1Vel[0]->getValue();
	const In2VecDeriv& in2_vecDeriv = dataVecIn2Vel[0]->getValue();
	OutVecDeriv& outVel = *dataVecOutVel[0]->beginEdit();


	helper::ReadAccessor<Data<type::vector<double>>> curv_abs_input  = d_curv_abs_section; // This is the vector of X in the paper
	helper::ReadAccessor<Data<type::vector<double>>> curv_abs_output = d_curv_abs_frames;
	helper::ReadAccessor<Data<bool>> debug = d_debug;
	m_frameJacobienVector.clear();
	m_frameJacobienDotVector.clear();


	// Compute the tangent Exponential SE3 vectors
	const In1VecCoord& inDeform = m_fromModel1->read(core::ConstVecCoordId::position())->getValue();
	this->update_TangExpSE3(inDeform,curv_abs_input.ref(),curv_abs_output.ref());

	//Get base velocity as input this is also called eta
	m_nodesVelocityVectors.clear();
	Deriv2 _baseVelocity;
	if (!in2_vecDeriv.empty())
		_baseVelocity = in2_vecDeriv[0];
	//convert to Vec6
	type::Vec6 baseVelocity;
	for (size_t u=0;u<6;u++) {baseVelocity[u] = _baseVelocity[u];}

	//Apply the local transform i.e from SOFA frame to Frederico frame
	const In2VecCoord& xfrom2Data = m_fromModel2->read(core::ConstVecCoordId::position())->getValue();
	Transform Tinverse = Transform(xfrom2Data[0].getCenter(),xfrom2Data[0].getOrientation()).inversed();
	Mat6x6 P = this->build_projector(Tinverse);
	m_nodeAdjointVectors.clear();

	type::Vec6 baseLocalVelocity = P * baseVelocity;
	m_nodesVelocityVectors.push_back(baseLocalVelocity);
	if(debug)
		std::cout << "Base local Velocity :"<< baseLocalVelocity <<std::endl;

	//Compute velocity at nodes
	for (size_t i = 1 ; i < curv_abs_input.size(); i++) {
		Transform t= m_nodesExponentialSE3Vectors[i].inversed();
		Mat6x6 Adjoint; Adjoint.clear();
		this->computeAdjoint(t,Adjoint);
                //Add this line because need for the jacobian computation
		m_nodeAdjointVectors.push_back(Adjoint);

        //Compute velocity (eta) at node i != 0 eq.(13) paper
        type::Vec6 Xi_dot = Vec6(in1[i-1],type::Vec3(0.0,0.0,0.0)) ;
        Vec6 temp = Adjoint * (m_nodesVelocityVectors[i-1] +
                               m_nodesTangExpVectors[i] * Xi_dot );
        m_nodesVelocityVectors.push_back(temp);
        if(debug)
                std::cout<< "Node velocity : "<< i << " = " << temp<< std::endl;
	}

	const OutVecCoord& out = m_toModel->read(core::ConstVecCoordId::position())->getValue();
	size_t sz =curv_abs_output.size();
	outVel.resize(sz);
	for (size_t i = 0 ; i < sz; i++) {
		Transform t= m_framesExponentialSE3Vectors[i].inversed();
		Mat6x6 Adjoint; Adjoint.clear();
		this->computeAdjoint(t,Adjoint);

        type::Vec6 Xi_dot = Vec6(in1[m_indicesVectors[i]-1],
                                         type::Vec3(0.0,0.0,0.0)) ;
        // eta
        Vec6 etaFrame = Adjoint * (m_nodesVelocityVectors[m_indicesVectors[i]-1]
                                   + m_framesTangExpVectors[i] * Xi_dot );

        //Compute here jacobien and jacobien_dot
        type::vector<Mat6x3> J_i, J_dot_i;
        computeJ_Jdot_i(Adjoint, i, J_i, etaFrame, J_dot_i);
        m_frameJacobienVector.push_back(J_i);
        m_frameJacobienVector.push_back(J_dot_i);

        //Convert from Federico node to Sofa node
        Transform _T = Transform(out[i].getCenter(),out[i].getOrientation());
        Mat6x6 _P = this->build_projector(_T);
        //std::cout<< "Eta local : "<< eta << std::endl;

        outVel[i] = _P * etaFrame;
        if(debug)
          std::cout<< "Frame velocity : "<< i << " = " << etaFrame<< std::endl;
	}
	//    std::cout << "Inside the apply J, outVel after computation  :  "<< outVel << std::endl;
	dataVecOutVel[0]->endEdit();
	m_index_input = 0;
}

template <class TIn1, class TIn2, class TOut>
void DiscreteDynamicCosseratMapping<TIn1, TIn2, TOut>::computeJ_Jdot_i(const Mat6x6 &Adjoint, size_t frameId,
                                                                      type::vector<Mat6x3> &J_i, const Vec6& etaFrame,
                                                                      type::vector<Mat6x3> &J_dot_i){
	size_t sz = d_curv_abs_section.getValue().size();
	Mat6x3 Si;
	Mat6x6 M;

	Mat6x3 Si_dot;
	Mat6x6 adj_eta; //to be computed
	//std::cout << "indice vector : "<< this->m_indicesVectors <<" vecId :"<< this->m_indicesVectors[frameId] << std::endl;
	//std::cout << "this->m_framesTangExpVectors[frameId]  :"<< this->m_framesTangExpVectors[frameId] << std::endl;
	//std::cout << "m_nodesTangExpVectors[frameId]  :"<< this->m_nodesTangExpVectors << std::endl;

	bool reachNode = false;
	for (unsigned int i = 1; i < sz; i++) {
		M = Adjoint;
		unsigned int u = m_indicesVectors[frameId];
		//std::cout << "frame : "<< frameId << " ==> section :" << i << " ==> u :"<< u << std::endl;
		if(i < u ){
			for (int j = u; j>0; j--) {
				M = M * m_nodeAdjointVectors[j] ;
			}
			//std::cout << "this->m_nodesTangExpVectors[u] : "<< this->m_nodesTangExpVectors[u]<< std::endl;
			Mat6x6 temp = M * m_nodesTangExpVectors[u];
			Si = temp * m_matrixBi;
			J_i.push_back(Si);

			Vec6 etaNode = m_nodesVelocityVectors[i];
			this->compute_adjointVec6(etaNode, adj_eta);
			Si_dot = temp * adj_eta * m_matrixBi;
			J_dot_i.push_back(Si_dot);
			//std::cout << "K1 Si : "<< Si << std::endl;
		}else{
			if(!reachNode){
				Mat6x6 temp = M * m_framesTangExpVectors[frameId] ;
				Si = temp * m_matrixBi;
				J_i.push_back(Si);

				this->compute_adjointVec6(etaFrame, adj_eta);
				Si_dot = temp * adj_eta * m_matrixBi;
				J_dot_i.push_back(Si_dot);
				reachNode = true;
				//std::cout << "K2 Si : "<< Si << std::endl;
			}else {
				Si.clear();
				J_i.push_back(Si);
				Si_dot.clear();
				J_dot_i.push_back(Si);
				//std::cout << "K3 Si : "<< Si << std::endl;
			}
		}
	}
    //    printf("J_%d:\n",frameId);
	for (unsigned int k=0; k<J_i.size(); k++) {
		std::cout << J_i[k] << std::endl;
	}
	printf("_____________________________________\n");
}

template <class TIn1, class TIn2, class TOut>
void DiscreteDynamicCosseratMapping<TIn1, TIn2, TOut>:: applyJT(
		const core::MechanicalParams* /*mparams*/, const type::vector< In1DataVecDeriv*>& dataVecOut1Force,
		const type::vector< In2DataVecDeriv*>& dataVecOut2Force,
		const type::vector<const OutDataVecDeriv*>& dataVecInForce)  {

	if(dataVecOut1Force.empty() || dataVecInForce.empty() || dataVecOut2Force.empty())
		return;

	const OutVecDeriv& in = dataVecInForce[0]->getValue();

	In1VecDeriv& out1 = *dataVecOut1Force[0]->beginEdit();
	In2VecDeriv& out2 = *dataVecOut2Force[0]->beginEdit();

	//Maybe need, in case the apply funcion is not call this must be call before
	//update_ExponentialSE3(in1);

	const OutVecCoord& frame = m_toModel->read(core::ConstVecCoordId::position())->getValue();
	const In1DataVecCoord* x1fromData = m_fromModel1->read(core::ConstVecCoordId::position());
	const In1VecCoord x1from = x1fromData->getValue();
        type::vector<Vec6> local_F_Vec ;   local_F_Vec.clear();

	out1.resize(x1from.size());

	//convert the input from Deriv type to vec6 type, for the purpose of the matrix vector multiplication
	//    std::cout<< "Size of frames :"<< in.size()<< std::endl;
	for (size_t var = 0; var < in.size(); ++var) {
        type::Vec6 vec;
		for(unsigned j = 0; j < 6; j++) vec[j] = in[var][j];

		//Convert input from global frame(SOFA) to local frame
		Transform _T = Transform(frame[var].getCenter(),frame[var].getOrientation());
		Mat6x6 P_trans =(this->build_projector(_T)); P_trans.transpose();
        type::Vec6 local_F = P_trans * vec;
		local_F_Vec.push_back(local_F);
	}

	//Compute output forces
	size_t sz = m_indicesVectors.size();

	unsigned int index =  m_indicesVectors[sz-1];
	m_totalBeamForceVectors.clear();
	m_totalBeamForceVectors.resize(sz);

	Vec6 F_tot; F_tot.clear();
	m_totalBeamForceVectors.push_back(F_tot);

	Mat3x6 matB_trans; matB_trans.clear();
	for(unsigned int k=0; k<3; k++) matB_trans[k][k] = 1.0;
	for (unsigned int s = sz ; s-- ; ) {
		Mat6x6 coAdjoint;
		//
		this->compute_coAdjoint(m_framesExponentialSE3Vectors[s], coAdjoint);  // m_framesExponentialSE3Vectors[s] computed in apply
		Vec6 node_F_Vec = coAdjoint * local_F_Vec[s];
		Mat6x6 temp = m_framesTangExpVectors[s];   // m_framesTangExpVectors[s] computed in applyJ (here we transpose)
		temp.transpose();
		type::Vec3 f = matB_trans * temp * node_F_Vec;

		if(index!=m_indicesVectors[s]){ // TODO to be replaced by while
			index--;
			//bring F_tot to the reference of the new beam
			this->compute_coAdjoint(m_nodesExponentialSE3Vectors[index],coAdjoint);  //m_nodesExponentialSE3Vectors computed in apply
			F_tot = coAdjoint * F_tot;
			Mat6x6 temp = m_nodesTangExpVectors[index];
			temp.transpose();
			//apply F_tot to the new beam
			type::Vec3 temp_f = matB_trans * temp * F_tot;
			out1[index-1] += temp_f;
		}
		if(d_debug.getValue())
			std::cout << "f at s ="<< s <<" and index"<< index <<  " is : "<< f << std::endl;

		//compte F_tot
		F_tot += node_F_Vec;
		out1[m_indicesVectors[s]-1] += f;
	}
	Transform frame0 = Transform(frame[0].getCenter(),frame[0].getOrientation());
	Mat6x6 M = this->build_projector(frame0);
	out2[0] += M * F_tot;

	if(d_debug.getValue()){
		std::cout << "Node forces "<< out1 << std::endl;
		std::cout << "base Force: "<< out2[0] << std::endl;
	}

	dataVecOut1Force[0]->endEdit();
	dataVecOut2Force[0]->endEdit();
}

//___________________________________________________________________________
template <class TIn1, class TIn2, class TOut>
void DiscreteDynamicCosseratMapping<TIn1, TIn2, TOut>::applyJT(
		const core::ConstraintParams*/*cparams*/ , const type::vector< In1DataMatrixDeriv*>&  dataMatOut1Const,
		const type::vector< In2DataMatrixDeriv*>&  dataMatOut2Const ,
		const type::vector<const OutDataMatrixDeriv*>& dataMatInConst)
{
	if(dataMatOut1Const.empty() || dataMatOut2Const.empty() || dataMatInConst.empty() )
		return;

	//We need only one input In model and input Root model (if present)
	In1MatrixDeriv& out1 = *dataMatOut1Const[0]->beginEdit(); // constraints on the strain space (reduced coordinate)
	In2MatrixDeriv& out2 = *dataMatOut2Const[0]->beginEdit(); // constraints on the reference frame (base frame)
	const OutMatrixDeriv& in = dataMatInConst[0]->getValue(); // input constraints defined on the mapped frames

	const OutVecCoord& frame = m_toModel->read(core::ConstVecCoordId::position())->getValue();
	const In1DataVecCoord* x1fromData = m_fromModel1->read(core::ConstVecCoordId::position());
	const In1VecCoord x1from = x1fromData->getValue();
	helper::ReadAccessor<Data<bool>> debug = d_debug;

	Mat3x6 matB_trans; matB_trans.clear();
	for(unsigned int k=0; k<3; k++) matB_trans[k][k] = 1.0;


    type::vector< std::tuple<int,Vec6> > NodesInvolved;
    type::vector< std::tuple<int,Vec6> > NodesInvolvedCompressed;
	//helper::vector<Vec6> NodesConstraintDirection;

	typename OutMatrixDeriv::RowConstIterator rowItEnd = in.end();

	for (typename OutMatrixDeriv::RowConstIterator rowIt = in.begin(); rowIt != rowItEnd; ++rowIt)
	{
		if (debug){
			std::cout<<"************* Apply JT (MatrixDeriv) iteration on line ";
			std::cout<<rowIt.index();
			std::cout<<"*************  "<<std::endl;
		}
		typename OutMatrixDeriv::ColConstIterator colIt = rowIt.begin();
		typename OutMatrixDeriv::ColConstIterator colItEnd = rowIt.end();

		// Creates a constraints if the input constraint is not empty.

		if (colIt == colItEnd)
		{
			if (debug){
				std::cout<<"no column for this constraint"<<std::endl;
			}
			continue;
		}
		typename In1MatrixDeriv::RowIterator o1 = out1.writeLine(rowIt.index()); // we store the constraint number
		typename In2MatrixDeriv::RowIterator o2 = out2.writeLine(rowIt.index());


		NodesInvolved.clear();
		//NodesConstraintDirection.clear();

		while (colIt != colItEnd)
		{
			int childIndex = colIt.index();


			const OutDeriv valueConst_ = colIt.val();
            type::Vec6 valueConst;
			for(unsigned j = 0; j < 6; j++) valueConst[j] = valueConst_[j];

			int indexBeam =  m_indicesVectors[childIndex];

			Transform _T = Transform(frame[childIndex].getCenter(),frame[childIndex].getOrientation());
			Mat6x6 P_trans =(this->build_projector(_T));
			P_trans.transpose();

			Mat6x6 coAdjoint;
			this->compute_coAdjoint(m_framesExponentialSE3Vectors[childIndex], coAdjoint);  // m_framesExponentialSE3Vectors[s] computed in apply
			Mat6x6 temp = m_framesTangExpVectors[childIndex];   // m_framesTangExpVectors[s] computed in applyJ (here we transpose)
			temp.transpose();

            type::Vec6 local_F =  coAdjoint * P_trans * valueConst; // constraint direction in local frame of the beam.


			type::Vec3 f = matB_trans * temp * local_F; // constraint direction in the strain space.


			o1.addCol(indexBeam-1, f);
			//std::cout<< "colIt :"<< colIt.index() <<" ; indexBeam :"<< indexBeam << " childIndex :"<< childIndex << " local_F : "<< local_F << std::endl;
			std::tuple<int,Vec6> test = std::make_tuple(indexBeam, local_F);

			NodesInvolved.push_back(test);
			colIt++;

		}
		if (debug){
			std::cout<<"==> NodesInvolved : "<<std::endl;
			for (size_t i = 0; i < NodesInvolved.size(); i++)
				std::cout << "index :" <<get<0>(NodesInvolved[i]) << " force :"
						  << get<1>(NodesInvolved[i]) << "\n ";
		}


		//std::cout<<" NodesInvolved before sort "<<NodesInvolved<<std::endl;

		// sort the Nodes Invoved by decreasing order
		std::sort(begin(NodesInvolved), end(NodesInvolved),
				  [](std::tuple<int, Vec6> const &t1, std::tuple<int, Vec6> const &t2) {
			return std::get<0>(t1) > std::get<0>(t2); // custom compare function
		} );

		//        for (size_t i = 0; i < NodesInvolved.size(); i++)
		//            std::cout << "index :" <<get<0>(NodesInvolved[i]) << " force :"
		//                      << get<1>(NodesInvolved[i]) << "\n ";

		NodesInvolvedCompressed.clear();


		for (unsigned n=0; n<NodesInvolved.size(); n++)
		{
			std::tuple<int,Vec6> test_i = NodesInvolved[n];
			int numNode_i= std::get<0>(test_i);
			Vec6 cumulativeF =std::get<1>(test_i);

			if (n<NodesInvolved.size()-1)
			{
				std::tuple<int,Vec6> test_i1 = NodesInvolved[n+1];
				int numNode_i1= std::get<0>(test_i1);

				while (numNode_i == numNode_i1)
				{
					cumulativeF += std::get<1>(test_i1);
					if ((n!=NodesInvolved.size()-2)||(n==0)){
						n++;
						break;
					}
					test_i1 = NodesInvolved[n+1];
					numNode_i1= std::get<0>(test_i1);
				}

			}
			NodesInvolvedCompressed.push_back(std::make_tuple(numNode_i, cumulativeF));
		}

		if (debug){
			std::cout<<" NodesInvolved after sort and compress"<<std::endl;
			for (size_t i = 0; i < NodesInvolvedCompressed.size(); i++)
				std::cout << "index :" <<get<0>(NodesInvolvedCompressed[i]) << " force :"
						  << get<1>(NodesInvolvedCompressed[i]) << "\n ";
		}

		for (unsigned n=0; n<NodesInvolvedCompressed.size(); n++)
		{

			std::tuple<int,Vec6> test = NodesInvolvedCompressed[n];
			int numNode= std::get<0>(test);
			int i = numNode;
			Vec6 CumulativeF = std::get<1>(test);

			while(i>0)
			{
				//cumulate on beam frame
				Mat6x6 coAdjoint;
				this->compute_coAdjoint(m_nodesExponentialSE3Vectors[i-1],coAdjoint);  //m_nodesExponentialSE3Vectors computed in apply
				CumulativeF = coAdjoint * CumulativeF;
				// transfer to strain space (local coordinates)
				Mat6x6 temp = m_nodesTangExpVectors[i-1];
				temp.transpose();
				type::Vec3 temp_f = matB_trans * temp * CumulativeF;

				if(i>1) o1.addCol(i-2, temp_f);

				i--;
			}

			Transform frame0 = Transform(frame[0].getCenter(),frame[0].getOrientation());
			Mat6x6 M = this->build_projector(frame0);

			Vec6 base_force = M * CumulativeF;
			o2.addCol(0, base_force);
		}
	}

	////// END ARTICULATION SYSTEM MAPPING
	dataMatOut1Const[0]->endEdit();
	dataMatOut2Const[0]->endEdit();
}

template <class TIn1, class TIn2, class TOut>
void DiscreteDynamicCosseratMapping<TIn1, TIn2, TOut>::draw(const core::visual::VisualParams* vparams)
{
	if (!vparams->displayFlags().getShowMappings())
		return;
}

} // namespace sofa::component::mapping
