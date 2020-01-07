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
#include "BaseCosserat.h"
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
BaseCosserat<TIn1, TIn2, TOut>::BaseCosserat()
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
void BaseCosserat<TIn1, TIn2, TOut>::init()
{
    //    if(this->getFromModels1().empty())
    //    {
    //        msg_error() << "Error while initializing ; input getFromModels1 not found" ;
    //        return;
    //    }

    //    if(this->getFromModels2().empty())
    //    {
    //        msg_error() << "Error while initializing ; output getFromModels2 not found" ;
    //        return;
    //    }

    //    if(this->getToModels().empty())
    //    {
    //        msg_error() << "Error while initializing ; output Model not found" ;
    //        return;
    //    }

    //    m_fromModel1 = this->getFromModels1()[0];
    //    m_fromModel2 = this->getFromModels2()[0];
    //    m_toModel = this->getToModels()[0];

    // Fill the initial vector
    const OutDataVecCoord* xfromData = m_toModel->read(core::ConstVecCoordId::position());
    const OutVecCoord xfrom = xfromData->getValue();
    //    WriteAccessor<Data < helper::vector<double>>> curv_abs_output = d_curv_abs_output;
    //    curv_abs_output.clear();

    //    m_vecTransform.clear();
    //    for (unsigned int i = 0; i < xfrom.size(); i++) {
    //        m_vecTransform.push_back(xfrom[i]);
    //    }
    printf("=================================> Init from the BaseCosserat component \n");
    //    initialize();

}


template <class TIn1, class TIn2, class TOut>
void BaseCosserat<TIn1, TIn2, TOut>::bwdInit()
{

}

template <class TIn1, class TIn2, class TOut>
void BaseCosserat<TIn1, TIn2, TOut>::reinit()
{

}


template <class TIn1, class TIn2, class TOut>
void BaseCosserat<TIn1, TIn2, TOut>::computeExponentialSE3(double & x, const defaulttype::Vector3& k, Transform & Trans){
    Matrix4 I4; I4.identity();

    double theta = k.norm();

    Matrix4 g_X;
    Matrix4 Xi ;

    Xi[0][1] = -k(2);
    Xi[0][2] = k[1];
    Xi[1][2] = -k[0];

    Xi[1][0] = -Xi(0,1);
    Xi[2][0] = -Xi(0,2);
    Xi[2][1] = -Xi(1,2);

    Xi[0][3] = 1.0;

    if(theta <= std::numeric_limits<double>::epsilon()){
        g_X = I4 + x*Xi;
    }else {
        double scalar1= (1.0 - std::cos(x*theta))/(theta*theta);
        double scalar2 = (x*theta - std::sin(x*theta))/(theta*theta*theta);
        g_X = I4 + x*Xi + scalar1*Xi*Xi + scalar2*Xi*Xi*Xi ;
    }

    defaulttype::Mat3x3 M; g_X.getsub(0,0,M);

    defaulttype::Quat R; R.fromMatrix(M);
    Vector3 T = Vector3(g_X(0,3),g_X(1,3),g_X(2,3));

    Trans = Transform(T,R);
}


//Fill exponential vectors
template <class TIn1, class TIn2, class TOut>
void BaseCosserat<TIn1, TIn2, TOut>::update_ExponentialSE3(const In1VecCoord & inDeform){
    //helper::ReadAccessor<Data<helper::vector<double>>> curv_abs_input = d_curv_abs_input;
    helper::ReadAccessor<Data<helper::vector<double>>> curv_abs_output = d_curv_abs_output;
    //m_index_input = 0;

    m_ExponentialSE3Vectors.clear();
    m_nodesExponentialSE3Vectors.clear();
    m_nodesLogarithmeSE3Vectors.clear();
    size_t sz = curv_abs_output.size();

    //Compute exponential at frame points
    for (size_t i = 0; i < sz; i++) {
        Transform T ;

        Vector3 k = inDeform[m_indicesVectors[i]-1];
        double  x = m_framesLenghtVectors[i];
        computeExponentialSE3(x,k,T);
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
        Transform T; computeExponentialSE3(x,k,T) ;
        m_nodesExponentialSE3Vectors.push_back(T);

        //////////////////
        //        Eigen::Matrix4d gX = convertTransformToMatrix4x4(T);
        //        Eigen::Matrix4d log_gX= (1.0/x) * computeLogarithme(x, gX);
        //        std::cout << "k : \n"<< k << std::endl;
        //        std::cout << "The logarithme : \n"<< log_gX << std::endl;
        //        m_nodesLogarithmeSE3Vectors.push_back(log_gX);
    }

    if (d_debug.getValue()){
        printf("_________________Beam Expo___________________\n");
        std::cout << "Beam Expo : "<< m_nodesExponentialSE3Vectors << std::endl;
        printf("_________________Beam Expo___________________\n");
    }

}



template <class TIn1, class TIn2, class TOut>
void BaseCosserat<TIn1, TIn2, TOut>:: computeAdjoint(const Transform & frame, Mat6x6 &Adjoint)
{
    Matrix3 R = extract_rotMatrix(frame);
    Vector3 u = frame.getOrigin();
    Matrix3 tild_u = getTildMatrix(u);
    Matrix3 tild_u_R = tild_u*R;
    buildaAdjoint(R,tild_u_R, Adjoint);
}

template <class TIn1, class TIn2, class TOut>
void BaseCosserat<TIn1, TIn2, TOut>:: compute_coAdjoint(const Transform & frame, Mat6x6 &coAdjoint)
{
    Matrix3 R = extract_rotMatrix(frame);
    Vector3 u = frame.getOrigin();
    Matrix3 tild_u = getTildMatrix(u);
    Matrix3 tild_u_R = tild_u*R;
    build_coaAdjoint(R,tild_u_R, coAdjoint);
}

template <class TIn1, class TIn2, class TOut>
void BaseCosserat<TIn1, TIn2, TOut>::compute_adjointVec6(const Vec6& eta, Mat6x6 &adjoint)
{
    Matrix3 tildMat = getTildMatrix(Vector3(eta[0],eta[1],eta[2]));
    adjoint.setsub(0,0,tildMat);
    adjoint.setsub(3,3,tildMat);
    adjoint.setsub(3,0,getTildMatrix(Vector3(eta[3],eta[4],eta[5])));
}



template <class TIn1, class TIn2, class TOut>
void BaseCosserat<TIn1, TIn2, TOut>::compute_Tang_Exp(double & x, const defaulttype::Vector3& k, Mat6x6 & TgX){
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
Matrix4 BaseCosserat<TIn1, TIn2, TOut>::computeLogarithme(const double & x, const Mat4x4 &gX){

    // Compute theta before everything
    const double theta = computeTheta(x, gX);
    Mat4x4 I4; I4.identity();
    Mat4x4 log_gX;


    double csc_theta = 1.0/(sin(x * theta/2.0));
    double sec_theta = 1.0/(cos(x * theta/2.0));
    double cst = (1.0/8) * (csc_theta*csc_theta*csc_theta) * sec_theta;
    double x_theta = x*theta;
    double cos_2Xtheta = cos(2.0 * x_theta);
    double cos_Xtheta = cos(x_theta);
    double sin_2Xtheta = sin(2.0 *x_theta);
    double sin_Xtheta = sin(x_theta);

    if(theta <= std::numeric_limits<double>::epsilon()) log_gX = I4;
    else {
        log_gX  = cst * ((x_theta*cos_2Xtheta - sin_Xtheta)*I4 -
                         (x_theta*cos_Xtheta + 2.0*x_theta*cos_2Xtheta - sin_Xtheta -sin_2Xtheta)*gX +
                         (2.0*x_theta*cos_Xtheta + x_theta*cos_2Xtheta-sin_Xtheta - sin_2Xtheta) *(gX*gX)-
                         (x_theta*cos_Xtheta - sin_Xtheta)*(gX*gX*gX));
    }

    return log_gX;
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
void BaseCosserat<TIn1, TIn2, TOut>::update_TangExpSE3(const In1VecCoord & inDeform, const helper::vector<double> &curv_abs_input , const helper::vector<double> &curv_abs_output ){

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
defaulttype::Vec6 BaseCosserat<TIn1, TIn2, TOut>::compute_eta(const defaulttype::Vec6 & baseEta, const In1VecDeriv & k_dot, const double abs_input){

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



//___________________________________________________________________________
template <class TIn1, class TIn2, class TOut>
void BaseCosserat<TIn1, TIn2, TOut>::initialize()
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
        //Fill the vector m_framesLenghtVectors with the distance
        //between frame(output) and the closest beam node toward the base
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
//void BaseCosserat<TIn1, TIn2, TOut>::applyDJT(const core::MechanicalParams* mparams, core::MultiVecDerivId inForce, core::ConstMultiVecDerivId outForce){}


//template <class TIn1, class TIn2, class TOut>
//void BaseCosserat<TIn1, TIn2, TOut>::do_applyJT( In1MatrixDeriv& out1, const OutMatrixDeriv& in, In2MatrixDeriv* out2 ){}

//template <class TIn1, class TIn2, class TOut>
//void BaseCosserat<TIn1, TIn2, TOut>:: applyJT(
//        const core::ConstraintParams* /* cparams */, const helper::vector< In1DataMatrixDeriv*>& dataMatOut1Const ,
//        const helper::vector< In2DataMatrixDeriv*>&  dataMatOut2Const ,
//        const helper::vector<const OutDataMatrixDeriv*>& dataMatInConst) {}


template <class TIn1, class TIn2, class TOut>
void BaseCosserat<TIn1, TIn2, TOut>::draw(const core::visual::VisualParams* vparams)
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

#endif // SOFA_COMPONENT_MAPPING_BASEMAPPING_INL
