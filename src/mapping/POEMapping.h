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
#ifndef SOFA_COMPONENT_MAPPING_POEMAPING_H
#define SOFA_COMPONENT_MAPPING_POEMAPING_H

#include <sofa/core/BaseMapping.h>
#include <sofa/core/core.h>
#include <sofa/core/Multi2Mapping.h>
#include "../initCosserat.h"
#include <sofa/defaulttype/SolidTypes.h>
#include <sofa/defaulttype/RigidTypes.h>



namespace sofa

{
using sofa::defaulttype::SolidTypes ;
using sofa::core::objectmodel::BaseContext ;
using sofa::defaulttype::Matrix3;
using sofa::defaulttype::Matrix4;
using sofa::defaulttype::Vector3;
using sofa::defaulttype::Vec6;

namespace component
{

namespace mapping
{

/*!
 * \class POEMapping
 * @brief Computes and map the length of the beams
 *
 * This is a component:
 * https://www.sofa-framework.org/community/doc/programming-with-sofa/create-your-component/
 */


template <class TIn1, class TIn2, class TOut>
class POEMapping : public core::Multi2Mapping<TIn1, TIn2, TOut>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE3(POEMapping, TIn1,TIn2, TOut), SOFA_TEMPLATE3(core::Multi2Mapping, TIn1, TIn2, TOut) );
    typedef core::Multi2Mapping<TIn1, TIn2, TOut> Inherit;

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
    typedef typename Out::Coord outCoord;
    typedef typename Out::VecDeriv OutVecDeriv;
    typedef typename Out::MatrixDeriv OutMatrixDeriv;
    typedef Data<OutVecCoord> OutDataVecCoord;
    typedef Data<OutVecDeriv> OutDataVecDeriv;
    typedef Data<OutMatrixDeriv> OutDataMatrixDeriv;

    typedef MultiLink<POEMapping<In1,In2,Out>, sofa::core::State< In1 >, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> LinkFromModels1;
    typedef MultiLink<POEMapping<In1,In2,Out>, sofa::core::State< In2 >, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> LinkFromModels2;
    typedef MultiLink<POEMapping<In1,In2,Out>, sofa::core::State< Out >, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> LinkToModels;

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

    unsigned int m_index_input;
    OutVecCoord m_vecTransform ;

    /*===========COSSERAT VECTORS ======================*/

    helper::vector<Transform> m_ExponentialSE3Vectors;
    helper::vector<Transform> m_nodesExponentialSE3Vectors;

    helper::vector<int> m_indicesVectors;

    helper::vector<double> m_beamLenghtVectors;
    helper::vector<double> m_framesLenghtVectors;

    helper::vector<Vec6> m_nodesVelocityVectors;
    helper::vector<Mat6x6> m_nodesTangExpVectors;
    helper::vector<Mat6x6> m_framesTangExpVectors;



protected:
    /// Constructor
    //POEMapping();
    POEMapping() ;
    /// Destructor
    ~POEMapping()  override {}
public:


    /**********************SOFA METHODS**************************/
    void init() override;
    virtual void bwdInit() override;  // get the points
    virtual void reset() override;
    virtual void reinit() override;
    void draw(const core::visual::VisualParams* vparams) override;

    /**********************MAPPING METHODS**************************/
    void apply(
            const core::MechanicalParams* /* mparams */, const helper::vector<OutDataVecCoord*>& dataVecOutPos,
            const helper::vector<const In1DataVecCoord*>& dataVecIn1Pos ,
            const helper::vector<const In2DataVecCoord*>& dataVecIn2Pos) override;

    void applyJ(
            const core::MechanicalParams* /* mparams */, const helper::vector< OutDataVecDeriv*>& dataVecOutVel,
            const helper::vector<const In1DataVecDeriv*>& dataVecIn1Vel,
            const helper::vector<const In2DataVecDeriv*>& dataVecIn2Vel) override;

    //ApplyJT Force
    void applyJT(
            const core::MechanicalParams* /* mparams */, const helper::vector< In1DataVecDeriv*>& dataVecOut1Force,
            const helper::vector< In2DataVecDeriv*>& dataVecOut2RootForce,
            const helper::vector<const OutDataVecDeriv*>& dataVecInForce) override;

    void applyDJT(const core::MechanicalParams* mparams, core::MultiVecDerivId inForce, core::ConstMultiVecDerivId outForce) override{}

    /// This method must be reimplemented by all mappings if they need to support constraints.
    virtual void applyJT(
        const core::ConstraintParams* /* cparams */, const helper::vector< In1DataMatrixDeriv*>& /* dataMatOut1Const */ ,
        const helper::vector< In2DataMatrixDeriv*>&  /* dataMatOut2Const */,
        const helper::vector<const OutDataMatrixDeriv*>& /* dataMatInConst */) override { }

protected:
    /**********************COSSERAT METHODS**************************/
    void computeExponentialSE3(double & x, const defaulttype::Vector3& k, Transform & Trans);
    void computeAdjoint(const Transform & frame, Mat6x6 &adjoint);

    void update_ExponentialSE3(const In1VecCoord & inDeform);
    void update_TangExpSE3(const In1VecCoord & inDeform, const helper::vector<double> &curv_abs_input , const helper::vector<double> &curv_abs_output );

    void buildaAdjoint(const Matrix3 &R, const Matrix3 & tild_u_R, Mat6x6 & Adjoint){

        for (unsigned int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; ++j) {
                Adjoint [i][j]= R[i][j];
                Adjoint [i+3][j+3]= R[i][j];
                Adjoint [i+3][j]= tild_u_R[i][j];
            }
        }
    }

    void compute_Tang_Exp(double & x, const Vector3& k, Mat6x6 & TgX);

    defaulttype::Vec6 compute_eta(const Vec6 & baseEta, const In1VecDeriv & k_dot, const double abs_input);

    void compute_Forces_6(const unsigned int index, const Vec6& inputForce, Vec6& outForce6){
        // Fill the initial vector
        const In1DataVecCoord* x1fromData = m_fromModel1->read(core::ConstVecCoordId::position());
        const In1VecCoord x1from = x1fromData->getValue();

        helper::ReadAccessor<Data<helper::vector<double>>> curv_abs_input1 = d_curv_abs_input;
        helper::ReadAccessor<Data<helper::vector<double>>> curv_abs_input2 = d_curv_abs_output;

        Transform out_Trans;
        Mat6x6 Adjoint, temp_Adjoint;
        //        std::cout<< "index : "<< index << " ";
        //        std::cout<< "m_index_input_1 : "<< m_index_input_1 << " ";
        //        std::cout<< "curv_abs_input2[index]  : "<< curv_abs_input2[index]  << " ";
        //        std::cout<< "curv_abs_input1[m_index_input_1] : "<< curv_abs_input1[m_index_input_1] << " "<< std::endl;


        //std::cout<< " index : "<< index << " m_index_input_1 :"<< m_index_input_1 << std::endl;                                                                                                                                              ""<< std::endl;

        for (unsigned int j = 0; j < 6; j++) Adjoint[j][j] = 1.0;
        double diff0 ;
        double _diff0 ;

        while(curv_abs_input2[index] > curv_abs_input1[m_index_input]){
            if(m_index_input == 0){
                diff0 = curv_abs_input1[m_index_input]; //
                _diff0 = -diff0;
            }else {
                diff0 = curv_abs_input1[m_index_input] - curv_abs_input1[m_index_input - 1];
                _diff0 = -diff0;
            }
            computeExponentialSE3(_diff0,x1from[m_index_input],out_Trans);
            computeAdjoint(out_Trans,temp_Adjoint);
            Adjoint = temp_Adjoint * Adjoint;

            m_index_input++;
        }

        if(m_index_input == 0){
            diff0 = curv_abs_input2[index];
            _diff0 = -diff0;
        }else {
            diff0 = curv_abs_input2[index] - curv_abs_input1[m_index_input - 1];
            _diff0 = -diff0;
        }
        computeExponentialSE3(_diff0,x1from[m_index_input],out_Trans);
        computeAdjoint(out_Trans,temp_Adjoint);
        Adjoint = temp_Adjoint * Adjoint;
        //        std::cout<< "Adjoint :\n";
        //        print_matrix(Adjoint);
        //        printf("\n");
        Adjoint.transpose();
        //        std::cout<< "Trans Adjoint :\n" ;
        //        print_matrix(Adjoint);
        //        printf("\n");

        outForce6 = Adjoint * inputForce;
    }
    void compute_Forces_3(const unsigned int index1, const unsigned int index2, const Vec6& inputForce, Vector3& outForce3 ){

        // Fill the initial vector
        const In1DataVecCoord* x1fromData = m_fromModel1->read(core::ConstVecCoordId::position());
        const In1VecCoord x1from = x1fromData->getValue();

        helper::ReadAccessor<Data<helper::vector<double>>> curv_abs_input1 = d_curv_abs_input;
        helper::ReadAccessor<Data<helper::vector<double>>> curv_abs_input2 = d_curv_abs_output;

        Transform out_Trans;
        Mat6x6 Adjoint, TgX;

        for (unsigned int j = 0; j < 6; j++) Adjoint[j][j] = 1.0;

        double diff0 ;
        double _diff0 ;

        while(curv_abs_input2[index1] > curv_abs_input1[m_index_input]){
            //            printf("\n____________________While ___________________\n");
            //            std::cout<< "m_index_input_1 : "<< m_index_input_1 << " ";
            //            std::cout<< "curv_abs_input2[index]  : "<< curv_abs_input2[index1]  << " ";
            //            std::cout<< "curv_abs_input1[m_index_input_1] : "<< curv_abs_input1[m_index_input_1] << " "<< std::endl;
            if(m_index_input == 0){
                diff0 = curv_abs_input1[m_index_input]; //
                _diff0 = -diff0;
            }else {
                diff0 = curv_abs_input1[m_index_input] - curv_abs_input1[m_index_input - 1];
                _diff0 = -diff0;
            }
            Mat6x6 temp_Adjoint; temp_Adjoint.clear();
            computeExponentialSE3(_diff0,x1from[m_index_input],out_Trans);
            computeAdjoint(out_Trans,temp_Adjoint);
            Adjoint = temp_Adjoint * Adjoint;

            m_index_input++;
        }

        if(m_index_input == 0){
            diff0 = curv_abs_input2[index1];
            _diff0 = -diff0;
        }else {
            diff0 = curv_abs_input2[index1] - curv_abs_input1[m_index_input - 1];
            _diff0 = -diff0;
        }

        if(_diff0 > 0.0 ) {
            //            printf("_diff > 0 \n");
            outForce3 = Vector3(0.0, 0.0, 0.0);}
        else {
            //            printf("\n____________________Else ___________________\n");
            //            std::cout<< "m_index_input_1 : "<< m_index_input_1 << " ";
            //            std::cout<< "curv_abs_input2[index]  : "<< curv_abs_input2[index1]  << " ";
            //            std::cout<< "curv_abs_input1[m_index_input_1] : "<< curv_abs_input1[m_index_input_1] << " "<< std::endl;

            Mat6x6 temp_Adjoint; temp_Adjoint.clear();
            computeExponentialSE3(_diff0,x1from[m_index_input],out_Trans);
            computeAdjoint(out_Trans,temp_Adjoint);
            Adjoint = temp_Adjoint * Adjoint;

            //            printf("Adjoint \n");
            //            print_matrix(Adjoint);
            //            printf(" \n");

            //compute_Tang_Exp(double & x, const defaulttype::Vector3& k, Mat6x6 & TgX)
            double x = std::min(curv_abs_input2[index1], curv_abs_input1[index2 - 1]);
            if(index2>1) x-=curv_abs_input1[index2 - 2] ;

            compute_Tang_Exp(x, x1from[index2 - 1], TgX);

            //            printf("TgX \n");
            //            print_matrix(TgX);
            //            printf(" \n");

            Mat6x3 matB; matB.clear(); for(unsigned int k=0; k<3; k++) matB[k][k] = 1.0;

            Mat6x3 S = Adjoint* TgX * matB;
            Mat3x6 transS = S.transposed();

            outForce3 =  transS * inputForce;

        }
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

    Matrix3 extract_rotMatrix(const Transform frame){
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
        return  P;
        //std::cout<< "Extract mat : "<< P << std::endl;
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

};

//extern template class SOFA_POE_MAPPING_API POEMapping<defaulttype::Vec3Types>;


#if  !defined(SOFA_COMPONENT_MAPPING_POE_MAPING_CPP)
extern template class SOFA_COSSERAT_MAPPING_API POEMapping< sofa::defaulttype::Vec3Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types >;
#endif

} // mapping

} // namespace componenet

} // namespace sofa

#endif //SOFA_COMPONENT_MAPPING_POEMAPING_H

