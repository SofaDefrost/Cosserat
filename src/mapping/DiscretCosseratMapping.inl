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
#ifndef SOFA_COMPONENT_MAPPING_POEMAPING_INL
#define SOFA_COMPONENT_MAPPING_POEMAPING_INL

#include <sofa/core/Multi2Mapping.inl>
#include "DiscretCosseratMapping.h"
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <string>
#include <sofa/core/visual/VisualParams.h>

#include <sofa/helper/AdvancedTimer.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/helper/logging/Message.h>



namespace sofa
{
using sofa::core::objectmodel::BaseContext ;
using sofa::helper::AdvancedTimer;
using sofa::helper::WriteAccessor;


namespace component
{

namespace mapping
{
template <class TIn1, class TIn2, class TOut>
DiscretCosseratMapping<TIn1, TIn2, TOut>::DiscretCosseratMapping()
    : d_curv_abs_input( initData( &d_curv_abs_input, "curv_abs_input", " need to be com...."))
    , d_curv_abs_output( initData( &d_curv_abs_output, "curv_abs_output", " need to be com...."))
    , d_debug( initData( &d_debug, true,"debug", "printf for the debug"))
    , m_fromModel1(NULL)
    , m_fromModel2(NULL)
    , m_toModel(NULL)
    , m_index_input(0)
{
}


// _________________________________________________________________________________________

template <class TIn1, class TIn2, class TOut>
void DiscretCosseratMapping<TIn1, TIn2, TOut>::init()
{
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

    m_fromModel1 = this->getFromModels1()[0];
    m_fromModel2 = this->getFromModels2()[0];
    m_toModel = this->getToModels()[0];

    // Fill the initial vector
    const OutDataVecCoord* xfromData = m_toModel->read(core::ConstVecCoordId::position());
    const OutVecCoord xfrom = xfromData->getValue();
    //    WriteAccessor<Data < helper::vector<double>>> curv_abs_output = d_curv_abs_output;
    //    curv_abs_output.clear();

    m_vecTransform.clear();
    for (unsigned int i = 0; i < xfrom.size(); i++) {
        m_vecTransform.push_back(xfrom[i]);
        //        if(i==0)
        //            curv_abs_output.push_back(1.0);
        //        else
        //            curv_abs_output.push_back(curv_abs_output[i-1] + (xfrom[i]-xfrom[i-1]).norm());
    }

    //    std::cout<< "=============> curv_abs_output :" << curv_abs_output << std::endl;
    //    WriteAccessor<Data<helper::vector<double>>> curv_abs_input = d_curv_abs_input;
    //    curv_abs_input.clear();

    //    const In1VecCoord x1from = m_fromModel1->read(core::ConstVecCoordId::position())->getValue();
    //    std::cout << "=====>  x1form :"<< x1from << std::endl;
    //    for (size_t i = 0 ; i < x1from.size(); i++) {
    //        std::cout << "=============> input :"<< x1from[i] << std::endl;
    //        if(i==0) curv_abs_input.push_back(0.0);
    //        else curv_abs_input.push_back(curv_abs_input[i-1] + (x1from[i] - x1from[i-1]).norm());
    //    }

    initialize();
    //Inherit::init();

}


template <class TIn1, class TIn2, class TOut>
void DiscretCosseratMapping<TIn1, TIn2, TOut>::bwdInit()
{

}

template <class TIn1, class TIn2, class TOut>
void DiscretCosseratMapping<TIn1, TIn2, TOut>::reinit()
{

}




template <class TIn1, class TIn2, class TOut>
void DiscretCosseratMapping<TIn1, TIn2, TOut>::reset()
{
    reinit();
}


template <class TIn1, class TIn2, class TOut>
void DiscretCosseratMapping<TIn1, TIn2, TOut>::computeExponentialSE3(double & x, const defaulttype::Vector3& k, Transform & Trans){
    Eigen::MatrixXd I4 = Eigen::Matrix4d::Identity(4,4);

    double theta = k.norm();

    Eigen::Matrix4d g_x;
    Eigen::Matrix4d Xi = Eigen::MatrixXd::Zero(4,4);

    Xi(0,1) = -k(2);
    Xi(0,2) = k[1];
    Xi(1,2) = -k[0];

    Xi(1,0) = -Xi(0,1);
    Xi(2,0) = -Xi(0,2);
    Xi(2,1) = -Xi(1,2);

    Xi(0,3) = 1.0;

    if(theta <= std::numeric_limits<double>::epsilon()){
        g_x = I4 + x*Xi;
    }else {
        double scalar1= (1.0 - std::cos(x*theta))/(theta*theta);
        double scalar2 = (x*theta - std::sin(x*theta))/(theta*theta*theta);
        g_x = I4 + x*Xi + scalar1*Xi*Xi + scalar2*Xi*Xi*Xi ;
    }

    defaulttype::Mat3x3 M;
    for(size_t j = 0; j <3; j++) {
        for(size_t i = 0; i <3; i++){
            M[j][i] = g_x(j,i);
        }
    }
    defaulttype::Quat R;
    R.fromMatrix(M);
    Vector3 T = Vector3(g_x(0,3),g_x(1,3),g_x(2,3));

    Trans = Transform(T,R);
}


//Fill exponential vectors
template <class TIn1, class TIn2, class TOut>
void DiscretCosseratMapping<TIn1, TIn2, TOut>::update_ExponentialSE3(const In1VecCoord & inDeform){
    //helper::ReadAccessor<Data<helper::vector<double>>> curv_abs_input = d_curv_abs_input;
    helper::ReadAccessor<Data<helper::vector<double>>> curv_abs_output = d_curv_abs_output;
    //m_index_input = 0;

    m_ExponentialSE3Vectors.clear();
    m_nodesExponentialSE3Vectors.clear();
    size_t sz = curv_abs_output.size();

    //Compute exponential at frame points
    for (size_t i = 0; i < sz; i++) {
        Transform T ;

        Vector3 k = inDeform[m_indicesVectors[i]-1];
        double  x = m_framesLenghtVectors[i];
        computeExponentialSE3(x,k,T) ;
        m_ExponentialSE3Vectors.push_back(T);

        if (d_debug.getValue()){
            printf("__________________________________________\n");
            std::cout << "x :"<< x << "; k :"<< k << std::endl;
            std::cout<< "m_ExponentialSE3Vectors :"<< m_ExponentialSE3Vectors[i] << std::endl;
        }
    }

    //Compute the exponential at the nodes
    m_nodesExponentialSE3Vectors.push_back(Transform(Vector3(0.0,0.0,0.0),defaulttype::Quat(0,0,0,1)));
    for (size_t j = 0; j < inDeform.size(); j++) {
        Vector3 k = inDeform[j];
        double  x = m_beamLenghtVectors[j];
        Transform T;         computeExponentialSE3(x,k,T) ;
        m_nodesExponentialSE3Vectors.push_back(T);
    }
    if (d_debug.getValue()){
        printf("_________________Beam Expo___________________\n");
        std::cout << "Beam Expo : "<< m_nodesExponentialSE3Vectors << std::endl;
    }
}


template <class TIn1, class TIn2, class TOut>
void DiscretCosseratMapping<TIn1, TIn2, TOut>::apply(
        const core::MechanicalParams* /* mparams */, const helper::vector<OutDataVecCoord*>& dataVecOutPos,
        const helper::vector<const In1DataVecCoord*>& dataVecIn1Pos ,
        const helper::vector<const In2DataVecCoord*>& dataVecIn2Pos)
{

    if(dataVecOutPos.empty() || dataVecIn1Pos.empty() || dataVecIn2Pos.empty())
        return;

    ///Do Apply
    //We need only one input In model and input Root model (if present)
    const In1VecCoord& in1 = dataVecIn1Pos[0]->getValue();
    const In2VecCoord& in2 = dataVecIn2Pos[0]->getValue();

    size_t sz = d_curv_abs_output.getValue().size();
    OutVecCoord& out = *dataVecOutPos[0]->beginEdit();
    out.resize(sz);

    //update the Exponential Matrices according to new deformation
    update_ExponentialSE3(in1);

    Transform frame0 = Transform(In2::getCPos(in2[0]),In2::getCRot(in2[0]));

    for(unsigned int i=0; i<sz; i++){
        Transform frame = frame0;
        for (int u = 0; u < m_indicesVectors[i]; u++) {
            frame *= m_nodesExponentialSE3Vectors[u];
        }
        frame *= m_ExponentialSE3Vectors[i];

        Vector3 v = frame.getOrigin();
        defaulttype::Quat q = frame.getOrientation();
        out[i] = OutCoord(v,q);
    }
    m_index_input = 0;

    dataVecOutPos[0]->endEdit();
}

template <class TIn1, class TIn2, class TOut>
void DiscretCosseratMapping<TIn1, TIn2, TOut>:: computeAdjoint(const Transform & frame, Mat6x6 &Adjoint)
{

    Matrix3 R = extract_rotMatrix(frame);
    Vector3 u = frame.getOrigin();
    Matrix3 tild_u = getTildMatrix(u);
    Matrix3 tild_u_R = tild_u*R;

    buildaAdjoint(R,tild_u_R, Adjoint);
}
template <class TIn1, class TIn2, class TOut>
void DiscretCosseratMapping<TIn1, TIn2, TOut>:: compute_coAdjoint(const Transform & frame, Mat6x6 &coAdjoint)
{
    Matrix3 R = extract_rotMatrix(frame);
    Vector3 u = frame.getOrigin();
    Matrix3 tild_u = getTildMatrix(u);
    Matrix3 tild_u_R = tild_u*R;

    build_coaAdjoint(R,tild_u_R, coAdjoint);
}



template <class TIn1, class TIn2, class TOut>
void DiscretCosseratMapping<TIn1, TIn2, TOut>::compute_Tang_Exp(double & x, const defaulttype::Vector3& k, Mat6x6 & TgX){

    double theta = k.norm();
    Matrix3 tild_p = getTildMatrix(Vector3(1.0,0.0,0.0));
    Matrix3 tild_k = getTildMatrix(k);

    Mat6x6 ad_Xi ;
    buildaAdjoint(tild_k,tild_p, ad_Xi);

    Mat6x6 Id6;
    for (unsigned int i =0; i< 6;i++) Id6[i][i]=1.0; //define identity 6x6

    if(theta <= std::numeric_limits<double>::epsilon()){
        double scalar0 = x*x/2.0;
        TgX = x*Id6 + scalar0 * ad_Xi;
    }else {
        double scalar1= (4.0 -4.0*cos(x*theta)-x*theta*sin(x*theta))/(2.0*theta*theta);
        double scalar2 = (4.0*x*theta + x*theta*cos(x*theta)-5.0*sin(x*theta))/(2.0*theta*theta*theta);
        double scalar3= (2.0 -2.0*cos(x*theta)-x*theta*sin(x*theta))/(2.0*theta*theta*theta*theta);
        double scalar4 = (2.0*x*theta + x*theta*cos(x*theta)-3.0*sin(x*theta))/(2.0*theta*theta*theta*theta*theta);

        TgX = x*Id6 + scalar1*ad_Xi + scalar2*ad_Xi*ad_Xi + scalar3*ad_Xi*ad_Xi*ad_Xi + scalar4*ad_Xi*ad_Xi*ad_Xi*ad_Xi ;
    }
}


//template<class In1VecCoord, class Mat6x6>
//void computeViolation(In1VecCoord& inDeform, const helper::vector<double> m_framesLenghtVectors, const
//                       size_t sz, std::function<double(int i, int j)> f)
//{


//    for (std::size_t i = 0; i < sz; i++) {
//        Mat6x6 temp ;

//        Vector3 k = inDeform[m_indicesVectors[i]-1];
//        double  x = m_framesLenghtVectors[i];
//        compute_Tang_Exp(x,k,temp) ;
//        m_framesTangExpVectors.push_back(temp);

////        if (d_debug.getValue()){
////            printf("__________________________________________\n");
////            std::cout << "x :"<< x << "; k :"<< k << std::endl;
////            std::cout<< "m_framesTangExpVectors :"<< m_framesTangExpVectors[i] << std::endl;
////        }
//    }
//}


template <class TIn1, class TIn2, class TOut>
void DiscretCosseratMapping<TIn1, TIn2, TOut>::update_TangExpSE3(const In1VecCoord & inDeform, const helper::vector<double> &curv_abs_input , const helper::vector<double> &curv_abs_output ){

    m_framesTangExpVectors.clear();
    size_t sz = curv_abs_output.size();
    //Compute tangExpo at frame points
    for (size_t i = 0; i < sz; i++) {
        Mat6x6 temp ;

        Vector3 k = inDeform[m_indicesVectors[i]-1];
        double  x = m_framesLenghtVectors[i];
        compute_Tang_Exp(x,k,temp) ;
        m_framesTangExpVectors.push_back(temp);

        if (d_debug.getValue()){
            printf("__________________________________________\n");
            std::cout << "x :"<< x << "; k :"<< k << std::endl;
            std::cout<< "m_framesTangExpVectors :"<< m_framesTangExpVectors[i] << std::endl;
        }
    }

    //Compute the TangExpSE3 at the nodes
    m_nodesTangExpVectors.clear();
    Mat6x6 tangExpO; tangExpO.clear();
    m_nodesTangExpVectors.push_back(tangExpO);

    for (size_t j = 1; j < curv_abs_input.size(); j++) {
        Vector3 k = inDeform[j-1];
        double  x = m_beamLenghtVectors[j-1];
        Mat6x6 temp; temp.clear();
        compute_Tang_Exp(x,k,temp);
        m_nodesTangExpVectors.push_back(temp);
    }
    if (d_debug.getValue()){
        printf("_________________Node TangExpo___________________\n");
        std::cout << "Node TangExpo : "<< m_nodesTangExpVectors << std::endl;
    }
}


template <class TIn1, class TIn2, class TOut>
defaulttype::Vec6 DiscretCosseratMapping<TIn1, TIn2, TOut>::compute_eta(const defaulttype::Vec6 & baseEta, const In1VecDeriv & k_dot, const double abs_input){

    // Fill the initial vector
    const In1DataVecCoord* x1fromData = m_fromModel1->read(core::ConstVecCoordId::position());
    const In1VecCoord x1from = x1fromData->getValue();

    helper::ReadAccessor<Data<helper::vector<double>>> curv_abs_input = d_curv_abs_input;
    //helper::ReadAccessor<Data<helper::vector<double>>> curv_abs_output = d_curv_abs_output;

    Transform out_Trans;
    Mat6x6 Adjoint, Tg;


    defaulttype::Vec6 Xi_dot;
    for(unsigned int i = 0; i < 3; i++) Xi_dot[i] = k_dot[m_index_input][i];


    double diff0 ;
    double _diff0 ;

    if(m_index_input == 0){
        diff0 = abs_input; //
        _diff0 = -abs_input;
    }else {
        diff0 = abs_input - curv_abs_input[m_index_input - 1];
        _diff0 = curv_abs_input[m_index_input - 1] - abs_input;
    }

    computeExponentialSE3(_diff0,x1from[m_index_input],out_Trans);
    computeAdjoint(out_Trans,Adjoint);

    compute_Tang_Exp(diff0,x1from[m_index_input],Tg);

    return Adjoint * (baseEta + Tg * Xi_dot );
}


template <class TIn1, class TIn2, class TOut>
void DiscretCosseratMapping<TIn1, TIn2, TOut>:: applyJ(
        const core::MechanicalParams* /* mparams */, const helper::vector< OutDataVecDeriv*>& dataVecOutVel,
        const helper::vector<const In1DataVecDeriv*>& dataVecIn1Vel,
        const helper::vector<const In2DataVecDeriv*>& dataVecIn2Vel) {

    if(dataVecOutVel.empty() || dataVecIn1Vel.empty() ||dataVecIn2Vel.empty() )
        return;
    const In1VecDeriv& in1 = dataVecIn1Vel[0]->getValue();
    const In2VecDeriv& in2_vecDeriv = dataVecIn2Vel[0]->getValue();
    OutVecDeriv& outVel = *dataVecOutVel[0]->beginEdit();



    helper::ReadAccessor<Data<helper::vector<double>>> curv_abs_input = d_curv_abs_input;
    helper::ReadAccessor<Data<helper::vector<double>>> curv_abs_output = d_curv_abs_output;

    // Compute the tangent Exponential SE3 vectors
    const In1VecCoord& inDeform = m_fromModel1->read(core::ConstVecCoordId::position())->getValue();
    update_TangExpSE3(inDeform,curv_abs_input.ref(),curv_abs_output.ref());

    //Get base velocity as input this is also called eta
    m_nodesVelocityVectors.clear();
    Deriv2 _baseVelocity;
    if (!in2_vecDeriv.empty())
        _baseVelocity = in2_vecDeriv[0];
    //convert to Vec6
    defaulttype::Vec6 baseVelocity;
    for (size_t u=0;u<6;u++) {baseVelocity[u] = _baseVelocity[u];}

    //Apply the local transform i.e from SOFA frame to Frederico frame
    const In2VecCoord& xfrom2Data = m_fromModel2->read(core::ConstVecCoordId::position())->getValue();
    Transform Tinverse = Transform(xfrom2Data[0].getCenter(),xfrom2Data[0].getOrientation()).inversed();
    Mat6x6 P = build_projector(Tinverse);
    defaulttype::Vec6 baseLocalVelocity = P * baseVelocity;
    m_nodesVelocityVectors.push_back(baseLocalVelocity);
    if(d_debug.getValue())
        std::cout << "Base local Velocity :"<< baseLocalVelocity <<std::endl;

    //Compute velocity at nodes
    for (size_t i = 1 ; i < curv_abs_input.size(); i++) {
        Transform t= m_nodesExponentialSE3Vectors[i].inversed();
        Mat6x6 Adjoint; Adjoint.clear();
        computeAdjoint(t,Adjoint);

        defaulttype::Vec6 Xi_dot = Vec6(in1[i-1],Vector3(0.0,0.0,0.0)) ;
        Vec6 temp = Adjoint * (m_nodesVelocityVectors[i-1] + m_nodesTangExpVectors[i] * Xi_dot );
        m_nodesVelocityVectors.push_back(temp);
        if(d_debug.getValue())
            std::cout<< "Node velocity : "<< i << " = " << temp<< std::endl;
    }

    const OutVecCoord& out = m_toModel->read(core::ConstVecCoordId::position())->getValue();
    size_t sz =curv_abs_output.size();
    outVel.resize(sz);
    for (size_t i = 0 ; i < sz; i++) {
        Transform t= m_ExponentialSE3Vectors[i].inversed();
        Mat6x6 Adjoint; Adjoint.clear();
        computeAdjoint(t,Adjoint);

        defaulttype::Vec6 Xi_dot = Vec6(in1[m_indicesVectors[i]-1],Vector3(0.0,0.0,0.0)) ;
        Vec6 temp = Adjoint * (m_nodesVelocityVectors[m_indicesVectors[i]-1] + m_framesTangExpVectors[i] * Xi_dot );

        Transform _T = Transform(out[i].getCenter(),out[i].getOrientation());
        Mat6x6 _P = build_projector(_T);
        //std::cout<< "Eta local : "<< eta << std::endl;

        outVel[i] = _P * temp;

        if(d_debug.getValue())
            std::cout<< "Frame velocity : "<< i << " = " << temp<< std::endl;
    }
    //    std::cout << "Inside the apply J, outVel after computation  :  "<< outVel << std::endl;
    dataVecOutVel[0]->endEdit();
    m_index_input = 0;
}


template <class TIn1, class TIn2, class TOut>
void DiscretCosseratMapping<TIn1, TIn2, TOut>:: applyJT(
        const core::MechanicalParams* /*mparams*/, const helper::vector< In1DataVecDeriv*>& dataVecOut1Force,
        const helper::vector< In2DataVecDeriv*>& dataVecOut2Force,
        const helper::vector<const OutDataVecDeriv*>& dataVecInForce)  {

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
    helper::vector<Vec6> local_F_Vec ;   local_F_Vec.clear();

    out1.resize(x1from.size());

    //convert the input from Deriv type to vec6 type, for the purpose of the matrix vector multiplication
    //    std::cout<< "Size of frames :"<< in.size()<< std::endl;
    for (size_t var = 0; var < in.size(); ++var) {
        defaulttype::Vec6 vec;
        for(unsigned j = 0; j < 6; j++) vec[j] = in[var][j];

        //Convert input from global frame(SOFA) to local frame
        Transform _T = Transform(frame[var].getCenter(),frame[var].getOrientation());
        Mat6x6 P_trans =(build_projector(_T)); P_trans.transpose();
        defaulttype::Vec6 local_F = P_trans * vec;
        local_F_Vec.push_back(local_F);
    }

    //Compute output forces
    size_t sz = m_indicesVectors.size();

    int index =  m_indicesVectors[sz-1];
    m_totalBeamForceVectors.clear();
    m_totalBeamForceVectors.resize(sz);

    Vec6 F_tot; F_tot.clear();
    m_totalBeamForceVectors.push_back(F_tot);

    Mat3x6 matB_trans; matB_trans.clear();
    for(unsigned int k=0; k<3; k++) matB_trans[k][k] = 1.0;
    for (size_t s = sz ; s-- ; ) {
        Mat6x6 coAdjoint;
        //
        compute_coAdjoint(m_ExponentialSE3Vectors[s],coAdjoint);  // m_ExponentialSE3Vectors[s] computed in apply
        Vec6 node_F_Vec = coAdjoint * local_F_Vec[s];
        Mat6x6 temp = m_framesTangExpVectors[s];   // m_framesTangExpVectors[s] computed in applyJ (here we transpose)
        temp.transpose();
        Vector3 f = matB_trans * temp * node_F_Vec;

        if(index!=m_indicesVectors[s]){ // TODO to be replaced by while
            index--;
            //bring F_tot to the reference of the new beam
            compute_coAdjoint(m_nodesExponentialSE3Vectors[index],coAdjoint);  //m_nodesExponentialSE3Vectors computed in apply
            F_tot = coAdjoint * F_tot;
            Mat6x6 temp = m_nodesTangExpVectors[index];
            temp.transpose();
            //apply F_tot to the new beam
            Vector3 temp_f = matB_trans * temp * F_tot;
            out1[index-1] += temp_f;
        }
        if(d_debug.getValue())
            std::cout << "f at s ="<< s <<" and index"<< index <<  " is : "<< f << std::endl;

        //compte F_tot
        F_tot += node_F_Vec;
        out1[m_indicesVectors[s]-1] += f;
    }
    Transform frame0 = Transform(frame[0].getCenter(),frame[0].getOrientation());
    Mat6x6 M = build_projector(frame0);
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
void DiscretCosseratMapping<TIn1, TIn2, TOut>::applyJT(
        const core::ConstraintParams*/*cparams*/ , const helper::vector< In1DataMatrixDeriv*>&  dataMatOut1Const,
        const helper::vector< In2DataMatrixDeriv*>&  dataMatOut2Const ,
        const helper::vector<const OutDataMatrixDeriv*>& dataMatInConst)
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


    Mat3x6 matB_trans; matB_trans.clear();
    for(unsigned int k=0; k<3; k++) matB_trans[k][k] = 1.0;


    helper::vector< std::tuple<int,Vec6> > NodesInvolved;
    helper::vector< std::tuple<int,Vec6> > NodesInvolvedCompressed;
    //helper::vector<Vec6> NodesConstraintDirection;

    typename OutMatrixDeriv::RowConstIterator rowItEnd = in.end();

    for (typename OutMatrixDeriv::RowConstIterator rowIt = in.begin(); rowIt != rowItEnd; ++rowIt)
    {
        if (d_debug.getValue()){
            std::cout<<"************* Apply JT (MatrixDeriv) iteration on line ";
            std::cout<<rowIt.index();
            std::cout<<"*************  "<<std::endl;
        }
        typename OutMatrixDeriv::ColConstIterator colIt = rowIt.begin();
        typename OutMatrixDeriv::ColConstIterator colItEnd = rowIt.end();

        // Creates a constraints if the input constraint is not empty.

        if (colIt == colItEnd)
        {
            if (d_debug.getValue()){
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
            defaulttype::Vec6 valueConst;
            for(unsigned j = 0; j < 6; j++) valueConst[j] = valueConst_[j];

            int indexBeam =  m_indicesVectors[childIndex];

            Transform _T = Transform(frame[childIndex].getCenter(),frame[childIndex].getOrientation());
            Mat6x6 P_trans =(build_projector(_T));
            P_trans.transpose();

            Mat6x6 coAdjoint;
            compute_coAdjoint(m_ExponentialSE3Vectors[childIndex],coAdjoint);  // m_ExponentialSE3Vectors[s] computed in apply
            Mat6x6 temp = m_framesTangExpVectors[childIndex];   // m_framesTangExpVectors[s] computed in applyJ (here we transpose)
            temp.transpose();

            defaulttype::Vec6 local_F =  coAdjoint * P_trans * valueConst; // constraint direction in local frame of the beam.


            Vector3 f = matB_trans * temp * local_F; // constraint direction in the strain space.


            o1.addCol(indexBeam-1, f);
            //std::cout<< "colIt :"<< colIt.index() <<" ; indexBeam :"<< indexBeam << " childIndex :"<< childIndex << " local_F : "<< local_F << std::endl;
            std::tuple<int,Vec6> test = std::make_tuple(indexBeam, local_F);

            NodesInvolved.push_back(test);
            colIt++;

        }
        if (d_debug.getValue()){
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

        if (d_debug.getValue()){
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
                compute_coAdjoint(m_nodesExponentialSE3Vectors[i-1],coAdjoint);  //m_nodesExponentialSE3Vectors computed in apply
                CumulativeF = coAdjoint * CumulativeF;
                // transfer to strain space (local coordinates)
                Mat6x6 temp = m_nodesTangExpVectors[i-1];
                temp.transpose();
                Vector3 temp_f = matB_trans * temp * CumulativeF;

                if(i>1) o1.addCol(i-2, temp_f);

                i--;
            }

            Transform frame0 = Transform(frame[0].getCenter(),frame[0].getOrientation());
            Mat6x6 M = build_projector(frame0);

            Vec6 base_force = M * CumulativeF;
            o2.addCol(0, base_force);
        }
    }

    ////// END ARTICULATION SYSTEM MAPPING
    dataMatOut1Const[0]->endEdit();
    dataMatOut2Const[0]->endEdit();
}


//___________________________________________________________________________
template <class TIn1, class TIn2, class TOut>
void DiscretCosseratMapping<TIn1, TIn2, TOut>::initialize()
{
    //find the beam on which each output is
    helper::ReadAccessor<Data<helper::vector<double>>> curv_abs_input = d_curv_abs_input;
    helper::ReadAccessor<Data<helper::vector<double>>> curv_abs_output = d_curv_abs_output;
    size_t sz = d_curv_abs_output.getValue().size();
    m_indicesVectors.clear();
    m_framesLenghtVectors.clear();
    m_beamLenghtVectors.clear();
    size_t input_index = 1;

    for (size_t i=0; i < sz; i++) {
        if (curv_abs_input[input_index] >= curv_abs_output[i]) {
            m_indicesVectors.push_back(input_index);
        }else {
            m_indicesVectors.push_back(input_index+1);
            input_index++;
        }
        //Fill the vector of distance between output and the closest beam node toward the base
        m_framesLenghtVectors.push_back(curv_abs_output[i] - curv_abs_input[m_indicesVectors[i]-1]);
    }

    for (size_t j = 0; j < curv_abs_input.size()-1; j++) {
        m_beamLenghtVectors.push_back(curv_abs_input[j+1] - curv_abs_input[j]);
    }

    if(d_debug.getValue()){
        std::cout<< "m_indicesVectors : "<< m_indicesVectors << std::endl;
        std::cout<< "m_framesLenghtVectors : "<< m_framesLenghtVectors << std::endl;
        std::cout<< "m_beamLenghtVectors : "<< m_beamLenghtVectors << std::endl;
    }
}



//template <class TIn1, class TIn2, class TOut>
//void DiscretCosseratMapping<TIn1, TIn2, TOut>::applyDJT(const core::MechanicalParams* mparams, core::MultiVecDerivId inForce, core::ConstMultiVecDerivId outForce){}


//template <class TIn1, class TIn2, class TOut>
//void DiscretCosseratMapping<TIn1, TIn2, TOut>::do_applyJT( In1MatrixDeriv& out1, const OutMatrixDeriv& in, In2MatrixDeriv* out2 ){}

//template <class TIn1, class TIn2, class TOut>
//void DiscretCosseratMapping<TIn1, TIn2, TOut>:: applyJT(
//        const core::ConstraintParams* /* cparams */, const helper::vector< In1DataMatrixDeriv*>& dataMatOut1Const ,
//        const helper::vector< In2DataMatrixDeriv*>&  dataMatOut2Const ,
//        const helper::vector<const OutDataMatrixDeriv*>& dataMatInConst) {}


template <class TIn1, class TIn2, class TOut>
void DiscretCosseratMapping<TIn1, TIn2, TOut>::draw(const core::visual::VisualParams* vparams)
{
    if (!vparams->displayFlags().getShowMappings())

        if(!d_debug.getValue()) return;
    for (unsigned int i = 0;i < m_vecTransform.size(); i++) {

        defaulttype::Quat q = m_vecTransform[i].getOrientation();
        q.normalize();

        defaulttype::Vector3 P1, x,y,z;
        P1 = m_vecTransform[i].getCenter();

        x= q.rotate(defaulttype::Vector3(1.0,0,0));
        y= q.rotate(defaulttype::Vector3(0,1.0,0));
        z= q.rotate(defaulttype::Vector3(0,0,1.0));
        double radius_arrow = 1.0/8.0;

        vparams->drawTool()->drawArrow(P1,(P1 + x)*1.0, radius_arrow, defaulttype::Vec<4,double>(1,0,0,1));
        vparams->drawTool()->drawArrow(P1,(P1 + y)*1.0, radius_arrow, defaulttype::Vec<4,double>(0,1,0,1));
        vparams->drawTool()->drawArrow(P1,(P1 + z)*1.0, radius_arrow, defaulttype::Vec<4,double>(0,0,1,1));
    }
    //return;
}


} // namespace mapping

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_MAPPING_POEMAPING_INL
