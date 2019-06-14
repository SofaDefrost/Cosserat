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
#include "POEMapping.h"
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
using sofa::defaulttype::Matrix3;
using sofa::defaulttype::Matrix4;
using sofa::defaulttype::Vector3;

namespace component
{

namespace mapping
{
template <class TIn1, class TIn2, class TOut>
POEMapping<TIn1, TIn2, TOut>::POEMapping()
    : d_applyRestPosition( initData( &d_applyRestPosition, false, "applyRestPosition", "set to true to apply this mapping to restPosition at init"))
    , d_curv_abs_input1( initData( &d_curv_abs_input1, "curv_abs_input1", " need to be com...."))
    , d_curv_abs_input2( initData( &d_curv_abs_input2, "curv_abs_input2", " need to be com...."))
    , m_fromModel1(NULL)
    , m_fromModel2(NULL)
    , m_toModel(NULL)
    , m_index_input_1(0)
{
}

template <class TIn1, class TIn2, class TOut>
void POEMapping<TIn1, TIn2, TOut>::init()
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

    m_vecTransform.clear();
    for (unsigned int i = 0; i < xfrom.size(); i++) {
        m_vecTransform.push_back(xfrom[i]);
    }
}


template <class TIn1, class TIn2, class TOut>
void POEMapping<TIn1, TIn2, TOut>::bwdInit()
{

}

template <class TIn1, class TIn2, class TOut>
void POEMapping<TIn1, TIn2, TOut>::reinit()
{
    //init();

}


template <class TIn1, class TIn2, class TOut>
void POEMapping<TIn1, TIn2, TOut>::reset()
{
    reinit();
}



template <class TIn1, class TIn2, class TOut>
void POEMapping<TIn1, TIn2, TOut>::computeExponentialSE3(double & x, const defaulttype::Vector3& k, Transform & Trans){
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
    for(unsigned int j = 0; j <3; j++) {
        for(unsigned int i = 0; i <3; i++){
            M[j][i] = g_x(j,i);
        }
    }
    defaulttype::Quat R;
    R.fromMatrix(M);
    Vector3 T = Vector3(g_x(0,3),g_x(1,3),g_x(2,3));

    Trans = Transform(T,R);
}

template <class TIn1, class TIn2, class TOut>
void POEMapping<TIn1, TIn2, TOut>::compute_frame(Transform &latestFrame, const In1VecCoord & inDeform, const unsigned int index2)
{
    helper::ReadAccessor<Data<helper::vector<double>>> curv_abs_input1 = d_curv_abs_input1;
    helper::ReadAccessor<Data<helper::vector<double>>> curv_abs_input2 = d_curv_abs_input2;

    if (curv_abs_input1[m_index_input_1]<curv_abs_input2[index2]){

        double diff1;
        if((m_index_input_1==0)&&(index2==0))
            diff1 = curv_abs_input1[m_index_input_1] ;
        else
            diff1 = curv_abs_input1[m_index_input_1] - curv_abs_input2[index2-1];

        double diff2 = curv_abs_input2[index2] - curv_abs_input1[m_index_input_1];

        Transform T0,T1;
        computeExponentialSE3(diff1,inDeform[m_index_input_1],T0) ;
        computeExponentialSE3(diff2,inDeform[m_index_input_1+1],T1);
        latestFrame= latestFrame * T0 * T1;
        m_index_input_1++;
        return;
    }

    double diff0;
    if(index2==0)
        diff0 = curv_abs_input2[index2];
    else
        diff0 = curv_abs_input2[index2] - curv_abs_input2[index2-1];

    Transform T2;
    computeExponentialSE3(diff0,inDeform[m_index_input_1],T2);

    latestFrame = latestFrame * T2;
}

template <class TIn1, class TIn2, class TOut>
void POEMapping<TIn1, TIn2, TOut>::apply(
        const core::MechanicalParams* /* mparams */, const helper::vector<OutDataVecCoord*>& dataVecOutPos,
        const helper::vector<const In1DataVecCoord*>& dataVecIn1Pos ,
        const helper::vector<const In2DataVecCoord*>& dataVecIn2Pos)
{

    if(dataVecOutPos.empty() || dataVecIn1Pos.empty() || dataVecIn2Pos.empty())
        return;

    ///Do Apply
    //We need only one input In model and input Root model (if present)
    OutVecCoord& out = *dataVecOutPos[0]->beginEdit();
    const In1VecCoord& in1 = dataVecIn1Pos[0]->getValue();
    const In2VecCoord& in2 = dataVecIn2Pos[0]->getValue();

    Transform Frame = Transform(In2::getCPos(in2[0]),In2::getCRot(in2[0]));
    //m_vecTransform.push_back(Frame);

    out.resize(d_curv_abs_input2.getValue().size());

    for(unsigned int i=0; i<d_curv_abs_input2.getValue().size(); i++){
        compute_frame(Frame, in1, i);

        Vector3 v = Frame.getOrigin();
        defaulttype::Quat q = Frame.getOrientation();
        out[i] = outCoord(v,q);
    }
    m_index_input_1 = 0;

    dataVecOutPos[0]->endEdit();
}

template <class TIn1, class TIn2, class TOut>
void POEMapping<TIn1, TIn2, TOut>:: computeAdjoint(const Transform & frame, Mat6x6 &Adjoint)
{

    Matrix3 R = extract_rotMatrix(frame);
    Vector3 u = frame.getOrigin();
    Matrix3 tild_u = getTildMatrix(u);
    Matrix3 tild_u_R = tild_u*R;

    buildaAdjoint(R,tild_u_R, Adjoint);
}


template <class TIn1, class TIn2, class TOut>
void POEMapping<TIn1, TIn2, TOut>::compute_Tang_Exp(double & x, const defaulttype::Vector3& k, Mat6x6 & TgX){

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


template <class TIn1, class TIn2, class TOut>
defaulttype::Vec6 POEMapping<TIn1, TIn2, TOut>::compute_eta(const defaulttype::Vec6 & baseEta, const In1VecDeriv & k_dot, const double abs_input){

    // Fill the initial vector
    const In1DataVecCoord* x1fromData = m_fromModel1->read(core::ConstVecCoordId::position());
    const In1VecCoord x1from = x1fromData->getValue();


    Transform out_Trans;
    Mat6x6 Adjoint, Tg;


    defaulttype::Vec6 Xi_dot;
    for(unsigned int i = 0; i < 3; i++) Xi_dot[i] = k_dot[m_index_input_1][i];


    double diff0 ;
    double _diff0 ;

    if(m_index_input_1 == 0){
        diff0 = abs_input; //
        _diff0 = -abs_input;
    }else {
        diff0 = abs_input - d_curv_abs_input1.getValue()[m_index_input_1 - 1];
        _diff0 = d_curv_abs_input1.getValue()[m_index_input_1 - 1] - abs_input;
    }

    computeExponentialSE3(_diff0,x1from[m_index_input_1],out_Trans);
    computeAdjoint(out_Trans,Adjoint);

    compute_Tang_Exp(diff0,x1from[m_index_input_1],Tg);

    return Adjoint * (baseEta + Tg * Xi_dot );
}



template <class TIn1, class TIn2, class TOut>
void POEMapping<TIn1, TIn2, TOut>:: applyJ(
        const core::MechanicalParams* /* mparams */, const helper::vector< OutDataVecDeriv*>& dataVecOutVel,
        const helper::vector<const In1DataVecDeriv*>& dataVecIn1Vel,
        const helper::vector<const In2DataVecDeriv*>& dataVecIn2Vel) {

    printf("====> inside the applyJ \n");
    if(dataVecOutVel.empty() || dataVecIn1Vel.empty())
        return;

    //We need only one input In1 model and input In2 model (if present)
    OutVecDeriv& outVel = *dataVecOutVel[0]->beginEdit();
    const In1VecDeriv& in1 = dataVecIn1Vel[0]->getValue();
    const In2VecDeriv& in2_vecDeriv = dataVecIn2Vel[0]->getValue();
    defaulttype::Vec6 eta;

    helper::ReadAccessor<Data<helper::vector<double>>> curv_abs_input2 = d_curv_abs_input2;

    //Get base velocity
    Deriv2 _base_eta;
    if (!in2_vecDeriv.empty())
        _base_eta = in2_vecDeriv[0];
    defaulttype::Vec6 base_eta;  for (unsigned p=0;p<6;p++) {base_eta[p] = _base_eta[p];}


    //For loop for computing all etas
    //applyJ(out,in, inroot);

    outVel.resize(d_curv_abs_input2.getValue().size());
    for (unsigned int j = 0; j < d_curv_abs_input2.getValue().size(); j++) {

        if(d_curv_abs_input2.getValue()[j] <= d_curv_abs_input1.getValue()[m_index_input_1]){
            //(const defaulttype::Vec6 & baseEta, const In1VecDeriv & k_dot, const unsigned int index_input_2)
            eta = compute_eta(base_eta,in1,d_curv_abs_input2.getValue()[j]);
        }else {
            base_eta = compute_eta(base_eta,in1,d_curv_abs_input1.getValue()[m_index_input_1]);
            m_index_input_1++;
            eta = compute_eta(base_eta,in1,d_curv_abs_input2.getValue()[j]);
        }
        std::cout<< "Eta : "<< eta << std::endl;
        outVel[j] = eta;
    }

    dataVecOutVel[0]->endEdit();

    m_index_input_1 = 0;
}



template <class TIn1, class TIn2, class TOut>
void POEMapping<TIn1, TIn2, TOut>:: applyJT(
        const core::MechanicalParams* /* mparams */, const helper::vector< In1DataVecDeriv*>& dataVecOut1Force,
        const helper::vector< In2DataVecDeriv*>& dataVecOut2RootForce,
        const helper::vector<const OutDataVecDeriv*>& dataVecInForce) {}

//template <class TIn1, class TIn2, class TOut>
//void POEMapping<TIn1, TIn2, TOut>::applyDJT(const core::MechanicalParams* mparams, core::MultiVecDerivId inForce, core::ConstMultiVecDerivId outForce){}


//template <class TIn1, class TIn2, class TOut>
//void POEMapping<TIn1, TIn2, TOut>::do_applyJT( In1MatrixDeriv& out1, const OutMatrixDeriv& in, In2MatrixDeriv* out2 ){}

//template <class TIn1, class TIn2, class TOut>
//void POEMapping<TIn1, TIn2, TOut>:: applyJT(
//        const core::ConstraintParams* /* cparams */, const helper::vector< In1DataMatrixDeriv*>& dataMatOut1Const ,
//        const helper::vector< In2DataMatrixDeriv*>&  dataMatOut2Const ,
//        const helper::vector<const OutDataMatrixDeriv*>& dataMatInConst) {}


template <class TIn1, class TIn2, class TOut>
void POEMapping<TIn1, TIn2, TOut>::draw(const core::visual::VisualParams* vparams)
{
    if (!vparams->displayFlags().getShowMappings())

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
