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
#define SOFA_COSSERAT_CPP_BaseCosserat
#include "BaseCosserat.inl"

#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/ObjectFactory.h>

namespace sofa::component::mapping
{
using namespace sofa::defaulttype;

template<>
BaseCosserat<Vec6Types, Rigid3Types, Rigid3Types>::se3 BaseCosserat<Vec6Types, Rigid3Types, Rigid3Types>::build_Xi_hat(const Coord1& strain_i){
    se3 Xi;

    //msg_info("BaseCosserat")<<" ===========> Build Xi Hat rigid is called ";

    Xi[0][1] = -strain_i(2);
    Xi[0][2] = strain_i[1];
    Xi[1][2] = -strain_i[0];

    Xi[1][0] = -Xi(0,1);
    Xi[2][0] = -Xi(0,2);
    Xi[2][1] = -Xi(1,2);

    Xi[0][3] = 1.0;

    //std::cout <<"Before the linear part : "<< Xi <<std::endl;
    for (unsigned int i=0; i<3; i++)
        Xi[i][3] += strain_i(i+3);

    //std::cout <<"After the linear part : "<< Xi <<std::endl;

//    se3 = [
//        0               -screw(3)   screw(2)        screw(4);
//        screw(3)        0           -screw(1)       screw(5);
//        -screw(2)   screw(1)        0               screw(6);
//        0                   0                 0                 0];

    return  Xi;
}



template<>
void BaseCosserat<Vec6Types, Rigid3Types, Rigid3Types>::computeExponentialSE3(const double & curv_abs_x_n, const Coord1& strain_n, Transform & Trans){
    Matrix4 I4; I4.identity();
    //Get the angular part of the
    Vector3 k = Vector3(strain_n(0), strain_n(1), strain_n(2));
    SReal theta = k.norm(); //

    SE3 g_X_n;
    se3 Xi_hat_n = build_Xi_hat(strain_n);

    if(d_debug.getValue())
        msg_info("BaseCosserat: ")<< "matrix Xi : "<< Xi_hat_n;

    if(theta <= std::numeric_limits<double>::epsilon()){
        g_X_n = I4 + curv_abs_x_n*Xi_hat_n;
    }else {
//        se3 Xi_hat_n_2 = Xi_hat_n * Xi_hat_n;
//        se3 Xi_hat_n_3 = Xi_hat_n_2 * Xi_hat_n;
//        SReal costheta =  std::cos(theta);
//        SReal sintheta =  std::cos(theta);
//        SReal theta2 = std::pow(theta,2);
//        SReal theta3 = theta2 * theta;
//        g_X_n = I4 + curv_abs_x_n*Xi_hat_n + ((1.-costheta)/(theta2))*Xi_hat_n_2 +((theta-sintheta)/theta3)*Xi_hat_n_3;
        double scalar1= (1.0 - std::cos(curv_abs_x_n*theta))/std::pow(theta,2);
        double scalar2 = (curv_abs_x_n*theta - std::sin(curv_abs_x_n*theta))/std::pow(theta,3);
        g_X_n = I4 + curv_abs_x_n*Xi_hat_n + scalar1*Xi_hat_n*Xi_hat_n + scalar2*Xi_hat_n*Xi_hat_n*Xi_hat_n ;
    }
    if(d_debug.getValue())
        msg_info("BaseCosserat: ")<< "matrix g_X : "<< g_X_n;

    type::Mat3x3 M;
    g_X_n.getsub(0,0,M); //get the rotation matrix

    // convert the rotation 3x3 matrix to a quaternion
    sofa::type::Quat<Real> R ; R.fromMatrix(M);
    Trans = Transform(type::Vec3(g_X_n(0,3),g_X_n(1,3),g_X_n(2,3)),R);
}

template<>
void BaseCosserat<Vec6Types, Rigid3Types, Rigid3Types>::compute_Tang_Exp(double & curv_abs_n, const Coord1 & strain_i, Mat6x6 & TgX){

  SReal theta = type::Vec3(strain_i(0), strain_i(1), strain_i(2)).norm(); //Sometimes this is computed over all strain
  Matrix3 tilde_k = getTildeMatrix(type::Vec3(strain_i(0), strain_i(1), strain_i(2)));
  /* Younes @23-11-27
  old version
  @Todo ???? is p the linear deformation? If so, why didn't I just put a zero vector in place of p and the first element of p is equal to 1?
  Matrix3 tilde_p = getTildeMatrix(type::Vec3(1.0, 0.0, 0.0));
  Using the new version does not bring any difference in my three reference scenes, but need more investogation
  #TECHNICAL_DEBT
  */
  Matrix3 tilde_q = getTildeMatrix(type::Vec3(strain_i(3), strain_i(4), strain_i(5)));

  Mat6x6 ad_Xi ;
  buildAdjoint(tilde_k, tilde_q, ad_Xi);

  Mat6x6 Id6 = Mat6x6::Identity() ;
  //    for (unsigned int i =0; i< 6;i++) Id6[i][i]=1.0; //define identity 6x6

  if(theta <= std::numeric_limits<double>::epsilon()){
    double scalar0 = std::pow(curv_abs_n,2)/2.0;
    TgX = curv_abs_n*Id6 + scalar0 * ad_Xi;
  }else {
    double scalar1 = (4.0 -4.0*cos(curv_abs_n*theta)-curv_abs_n*theta*sin(curv_abs_n*theta))/(2.0*theta*theta);
    double scalar2 = (4.0*curv_abs_n*theta + curv_abs_n*theta*cos(curv_abs_n*theta)-5.0*sin(curv_abs_n*theta))/(2.0*theta*theta*theta);
    double scalar3 = (2.0 -2.0*cos(curv_abs_n*theta)-curv_abs_n*theta*sin(curv_abs_n*theta))/(2.0*theta*theta*theta*theta);
    double scalar4 = (2.0*curv_abs_n*theta + curv_abs_n*theta*cos(curv_abs_n*theta)-3.0*sin(curv_abs_n*theta))/(2.0*theta*theta*theta*theta*theta);

    TgX = curv_abs_n*Id6 + scalar1*ad_Xi + scalar2*ad_Xi*ad_Xi + scalar3*ad_Xi*ad_Xi*ad_Xi + scalar4*ad_Xi*ad_Xi*ad_Xi*ad_Xi ;
  }
}


// Register in the Factory
int BaseCosseratClass = core::RegisterObject("Set the positions and velocities of points attached to a rigid parent")
        .add< BaseCosserat< sofa::defaulttype::Vec3Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types > >()
        .add< BaseCosserat< sofa::defaulttype::Vec6Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types > >()
;


template class SOFA_COSSERAT_API BaseCosserat< sofa::defaulttype::Vec3Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types >;
template class SOFA_COSSERAT_API BaseCosserat< sofa::defaulttype::Vec6Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types >;

} // namespace sofa.
