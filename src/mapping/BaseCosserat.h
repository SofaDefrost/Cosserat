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
#ifndef SOFA_COMPONENT_MAPPING_BASEMAPPING_H
#define SOFA_COMPONENT_MAPPING_BASEMAPPING_H

#include <sofa/core/BaseMapping.h>
#include <sofa/core/core.h>
#include <sofa/core/Multi2Mapping.h>
#include "../initCosserat.h"
#include <sofa/defaulttype/SolidTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <math.h>
#include <sofa/defaulttype/Vec.h>



namespace sofa

{
using sofa::defaulttype::SolidTypes ;
using sofa::core::objectmodel::BaseContext ;
using sofa::defaulttype::Matrix3;
using sofa::defaulttype::Matrix4;
using sofa::defaulttype::Vector3;
using sofa::defaulttype::Vec6;
using std::get;
using sofa::core::objectmodel::BaseObject;

namespace component
{

namespace mapping
{

/*!
 * \class BaseCosserat
 * @brief Computes and map the length of the beams
 *
 * This is a component:
 * https://www.sofa-framework.org/community/doc/programming-with-sofa/create-your-component/
 */


template <class TIn1, class TIn2, class TOut>
class BaseCosserat : public virtual sofa::core::objectmodel::BaseObject
{
public:
    SOFA_CLASS(SOFA_TEMPLATE3(BaseCosserat, TIn1,TIn2, TOut), BaseObject );
    typedef BaseObject Inherit;

    /// Input Model Type
    typedef TIn1 In1;
    typedef TIn2 In2;

    /// Output Model Type
    typedef TOut Out;

    typedef typename In1::Coord             Coord1       ;
    typedef typename In1::Deriv             Deriv1  ;
    typedef typename In1::VecCoord In1VecCoord;
    typedef typename In1::VecDeriv In1VecDeriv;
    typedef typename In1::MatrixDeriv In1MatrixDeriv;
    typedef Data<In1VecCoord> In1DataVecCoord;
    typedef Data<In1VecDeriv> In1DataVecDeriv;
    typedef Data<In1MatrixDeriv> In1DataMatrixDeriv;
    
    typedef typename In2::Coord::value_type Real          ;
    typedef typename In2::Coord             Coord2         ;
    typedef typename In2::Deriv             Deriv2         ;
    typedef typename In2::VecCoord In2VecCoord;
    typedef typename In2::VecDeriv In2VecDeriv;
    typedef typename In2::MatrixDeriv In2MatrixDeriv;
    typedef Data<In2VecCoord> In2DataVecCoord;
    typedef Data<In2VecDeriv> In2DataVecDeriv;
    typedef Data<In2MatrixDeriv> In2DataMatrixDeriv;
    typedef defaulttype::Mat<6,6,Real> Mat6x6;
    typedef defaulttype::Mat<3,6,Real> Mat3x6;
    typedef defaulttype::Mat<6,3,Real> Mat6x3;
    typedef defaulttype::Mat<4,4,Real> Mat4x4;

    typedef typename Out::VecCoord OutVecCoord;
    typedef typename Out::Coord OutCoord;
    typedef typename Out::Deriv OutDeriv;
    typedef typename Out::VecDeriv OutVecDeriv;
    typedef typename Out::MatrixDeriv OutMatrixDeriv;
    typedef Data<OutVecCoord> OutDataVecCoord;
    typedef Data<OutVecDeriv> OutDataVecDeriv;
    typedef Data<OutMatrixDeriv> OutDataMatrixDeriv;

    typedef MultiLink<BaseCosserat<In1,In2,Out>, sofa::core::State< In1 >, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> LinkFromModels1;
    typedef MultiLink<BaseCosserat<In1,In2,Out>, sofa::core::State< In2 >, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> LinkFromModels2;
    typedef MultiLink<BaseCosserat<In1,In2,Out>, sofa::core::State< Out >, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> LinkToModels;

    typedef typename SolidTypes<Real>::Transform      Transform ;

protected:
    Data<helper::vector<double>>      d_curv_abs_input ;
    Data<helper::vector<double>>      d_curv_abs_output ;
    Data<bool>                        d_debug ;

    /// Input Models container. New inputs are added through addInputModel(In* ).
    //    LinkFromModels1 m_fromModel1;
    //    LinkFromModels2 m_fromModel2;
    //    LinkToModels m_toModel;

    core::State<In1>* m_fromModel1;
    core::State<In2>* m_fromModel2;
    core::State<Out>* m_toModel;

public:
    /*===========COSSERAT VECTORS ======================*/
    unsigned int m_index_input;
    OutVecCoord m_vecTransform ;

    helper::vector<Transform> m_ExponentialSE3Vectors;
    helper::vector<Transform> m_nodesExponentialSE3Vectors;
    helper::vector<Matrix4> m_nodesLogarithmeSE3Vectors;

    helper::vector<int> m_indicesVectors;

    helper::vector<double> m_beamLenghtVectors;
    helper::vector<double> m_framesLenghtVectors;

    helper::vector<Vec6>   m_nodesVelocityVectors;
    helper::vector<Mat6x6> m_nodesTangExpVectors;
    helper::vector<Mat6x6> m_framesTangExpVectors;
    helper::vector<Vec6>   m_totalBeamForceVectors;

    helper::vector<Mat6x6> m_nodeAjointVectors;
    helper::vector<Mat6x6> m_nodeAjointEtaVectors;
    helper::vector<Mat6x6> m_frameAjointEtaVectors;
    helper::vector<Mat6x6> m_node_coAjointEtaVectors;
    helper::vector<Mat6x6> m_frame_coAjointEtaVectors;


protected:
    /// Constructor
    //BaseCosserat();
    BaseCosserat() ;
    /// Destructor
    ~BaseCosserat()  override {}
public:


    /**********************SOFA METHODS**************************/
    void init() override;
    virtual void bwdInit() override;  // get the points
    virtual void reinit() override;
    void draw(const core::visual::VisualParams* vparams) override;

protected:
    /**********************COSSERAT METHODS**************************/
    void computeExponentialSE3(double & x, const defaulttype::Vector3& k, Transform & Trans);
    void computeAdjoint(const Transform & frame, Mat6x6 &Adjoint);
    void compute_coAdjoint(const Transform & frame, Mat6x6 &coAdjoint);
    void compute_adjointVec6(const Vec6 & frame, Mat6x6 &adjoint);

    void update_ExponentialSE3(const In1VecCoord & inDeform);
    void update_TangExpSE3(const In1VecCoord & inDeform, const helper::vector<double> &curv_abs_input , const helper::vector<double> &curv_abs_output );
    void compute_Tang_Exp(double & x, const Vector3& k, Mat6x6 & TgX);

    defaulttype::Vec6 compute_eta(const Vec6 & baseEta, const In1VecDeriv & k_dot, const double abs_input);
    defaulttype::Matrix4 computeLogarithme(const double & x, const Mat4x4 &gX);
    double computeTheta(const double &x, const Mat4x4 &gX){
        double Tr_gx = 0.0;
        for (int i = 0; i<4; i++) {
            Tr_gx += gX[i][i];
        }

        double theta;
        if( x <= std::numeric_limits<double>::epsilon()) theta = 0.0;
        else theta = (1.0/x) * std::acos((Tr_gx/2.0) -1);

        return  theta;
    }


public:
    /**********************OTHER METHODS**************************/
    void initialize();
    void print_matrix(const Mat6x6 R){
        for (unsigned int k = 0; k < 6 ; k++) {
            for (unsigned int i = 0; i < 6; i++)
                printf("  %lf",R[k][i]);
            //std::cout<< " " << R[k][i];
            printf("\n");
        }
    }

    Matrix3 extract_rotMatrix(const Transform & frame){
        defaulttype::Quat q = frame.getOrientation();
        Real R[4][4];     q.buildRotationMatrix(R);
        Matrix3 mat;

        for (unsigned int k = 0; k < 3 ; k++)
            for (unsigned int i = 0; i < 3; i++)
                mat[k][i] = R[k][i];
        return  mat;
    }
    Mat6x6 build_projector(const Transform &T){
        Mat6x6 P;
        Real R[4][4]; (T.getOrientation()).buildRotationMatrix(R);

        for (unsigned int i = 0; i < 3;  i++) {
            for (unsigned int j = 0; j < 3; j++) {
                P[i][j+3] = R[i][j];
                P[i+3][j] = R[i][j];
            }
        }
        //        print_matrix(P);
        //        printf("__________________________\n");
        return  P;
    }

    Matrix3 getTildMatrix(const Vector3 & u){
        defaulttype::Matrix3 tild;
        tild[0][1] = -u[2];
        tild[0][2] = u[1];
        tild[1][2] = -u[0];

        tild[1][0] = -tild[0][1];
        tild[2][0] = -tild[0][2];
        tild[2][1] = -tild[1][2];
        return  tild;
    }

    void buildaAdjoint(const Matrix3 &A, const Matrix3 & B, Mat6x6 & Adjoint){

        Adjoint.clear();
        for (unsigned int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; ++j) {
                Adjoint [i][j]= A[i][j];
                Adjoint [i+3][j+3]= A[i][j];
                Adjoint [i+3][j]= B[i][j];
            }
        }
    }

    void build_coaAdjoint(const Matrix3 &A, const Matrix3 & B, Mat6x6 & coAdjoint){

        coAdjoint.clear();
        for (unsigned int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; ++j) {
                coAdjoint [i][j]= A[i][j];
                coAdjoint [i+3][j+3]= A[i][j];
                coAdjoint [i][j+3]= B[i][j];
            }
        }
    }

    Matrix4 convertTransformToMatrix4x4(const Transform & T){
        Matrix4 M; M.identity();
        Matrix3 R = extract_rotMatrix(T);
        Vector3 trans = T.getOrigin();

        for (size_t i = 0; i < 3; i++) {
            for (size_t j=0; j<3; j++){
                M(i,j) = R[i][j];
                M(i,3) = trans[i];
            }
        }
        return M;
    }

};


#if  !defined(SOFA_COMPONENT_MAPPING_POE_MAPING_CPP)
extern template class SOFA_COSSERAT_MAPPING_API BaseCosserat< sofa::defaulttype::Vec3Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types >;
#endif

} // mapping

} // namespace componenet

} // namespace sofa

#endif //SOFA_COMPONENT_MAPPING_POEMAPING_H

