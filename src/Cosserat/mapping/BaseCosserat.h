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
#include <Cosserat/config.h>

#include <sofa/core/BaseMapping.h>
#include <sofa/core/config.h>
#include <sofa/core/Multi2Mapping.h>
#include <sofa/defaulttype/SolidTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/type/Vec.h>

#include <iostream>
#include <cmath>

#include <Eigen/Dense>
#include <cmath>

using namespace std;
using namespace Eigen;

#include <cmath>

// TODO(dmarchal, 2024/06/07): This is polluating the namespace of sofa::components
//                             plugins should be in their own namespace like
//                             eg: cosserat::component::mapping
namespace sofa::component::mapping
{

// with a private namespace the used named are not polluating the whole sofa::component::mapping
// ones.
namespace {
using sofa::defaulttype::SolidTypes ;
using sofa::core::objectmodel::BaseContext ;
using sofa::type::Matrix3;
using sofa::type::Matrix4;
using type::Vec3;
using type::Vec6;
using std::get;
using sofa::core::objectmodel::BaseObject;
}

/*!
 * \class BaseCosserat
 * @brief Computes and map the length of the beams
 *
 * This is a component:
 * https://www.sofa-framework.org/community/doc/programming-with-sofa/create-your-component/
 */


// TODO(dmarchal: 2024/06/07) This component looks like a mapping but inherit from BaseObject *
// can you clarify why is is not inhering from BaseMapping
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


    // TODO(dmarchal:2024/06/07) This very long list of public typedefs looks questionnable
    // at least this has to be justified by comment on why these typedefs are public
    typedef typename In1::Coord             Coord1       ;
    typedef typename In1::Deriv             Deriv1  ;
    typedef typename In1::VecCoord In1VecCoord;
    typedef typename In1::VecDeriv In1VecDeriv;
    [[maybe_unused]] typedef typename In1::MatrixDeriv In1MatrixDeriv;
    typedef Data<In1VecCoord> In1DataVecCoord;
    typedef Data<In1VecDeriv> In1DataVecDeriv;

    typedef typename In2::Coord::value_type Real          ;
    typedef typename In2::Coord             Coord2         ;
    typedef typename In2::Deriv             Deriv2         ;
    typedef typename In2::VecCoord In2VecCoord;
    typedef typename In2::VecDeriv In2VecDeriv;
    [[maybe_unused]] typedef typename In2::MatrixDeriv In2MatrixDeriv;
    typedef Data<In2VecCoord> In2DataVecCoord;
    typedef Data<In2VecDeriv> In2DataVecDeriv;
    typedef type::Mat<6,6,SReal> Mat6x6;
    typedef type::Mat<4,4,SReal> Mat4x4;

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
    [[maybe_unused]] typedef MultiLink<BaseCosserat<In1,In2,Out>, sofa::core::State< Out >, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> LinkToModels;

    typedef typename SolidTypes<SReal>::Transform      Transform ;
    typedef typename type::vector<SReal> List;

    typedef typename sofa::type::Matrix4   se3;
    typedef typename sofa::type::Matrix4   SE3;
    typedef typename Eigen::Matrix4d   _SE3;
    typedef typename Eigen::Matrix4d   _se3;
    typedef typename type::Mat<6,6,SReal>   Tangent;
    typedef typename Eigen::Matrix3d RotMat;
    typedef typename Eigen::Matrix<SReal, 6, 1> Vector6d;

public:
    // TODO(dmarchal: 2024/06/07): There is a lot of public attributes is this really needed ?

    /*===========COSSERAT VECTORS ======================*/
    unsigned int m_index_input;
    OutVecCoord m_vecTransform ;

    type::vector<Transform> m_framesExponentialSE3Vectors;
    type::vector<Transform> m_nodesExponentialSE3Vectors;
    type::vector<Matrix4> m_nodesLogarithmeSE3Vectors;

    // @todo comment or explain more vectors
    type::vector<unsigned int> m_indicesVectors;
    type::vector<unsigned int> m_indicesVectorsDraw;

    type::vector<double> m_BeamLengthVectors;
    type::vector<double> m_framesLengthVectors;

    type::vector<Vec6>   m_nodesVelocityVectors;
    type::vector<Mat6x6> m_nodesTangExpVectors;
    type::vector<Mat6x6> m_framesTangExpVectors;
    type::vector<Vec6>   m_totalBeamForceVectors;

    type::vector<Mat6x6> m_nodeAdjointVectors;

    // TODO(dmarchal:2024/06/07): explain why these attributes are unused
    [[maybe_unused]] type::vector<Mat6x6> m_nodeAdjointEtaVectors;
    [[maybe_unused]] type::vector<Mat6x6> m_frameAdjointEtaVectors;
    [[maybe_unused]] type::vector<Mat6x6> m_node_coAdjointEtaVectors;
    [[maybe_unused]] type::vector<Mat6x6> m_frame_coAdjointEtaVectors;

public:
    /********************** Inhertited from BaseObject   **************/
    void init() override;
    void draw(const core::visual::VisualParams* vparams) override;

    /************************* BaseCosserat **************************/
    // TODO(dmarchal:2024/06/07), so we have "initialize" and "init" co-existances of both and their
    // roles is unclear and generates ambiguities
    void initialize();

    double computeTheta(const double &x, const Mat4x4 &gX);
    void printMatrix(const Mat6x6 R);

    Matrix3 extractRotMatrix(const Transform & frame);
    Tangent buildProjector(const Transform &T);
    se3 buildXiHat(const Coord1 & strain_i);
    Matrix3 getTildeMatrix(const type::Vec3 & u);

    void buildAdjoint(const Matrix3 &A, const Matrix3 & B, Mat6x6 & Adjoint);
    void buildCoAdjoint(const Matrix3 &A, const Matrix3 & B, Mat6x6 & coAdjoint);

    Matrix4 convertTransformToMatrix4x4(const Transform & T);
    Vec6 piecewiseLogmap(const _SE3& g_x) ;

    // TODO(dmarchal: 2024/06/07), this looks like a very common utility function... it shouldn't
    // be (re)implemented in a base classe.
    Eigen::Matrix3d rotationMatrixX(double angle) {
        Eigen::Matrix3d rotation;
        rotation << 1, 0, 0,
            0, cos(angle), -sin(angle),
            0, sin(angle), cos(angle);
        return rotation;
    }

    // TODO(dmarchal: 2024/06/07), this looks like a very common utility function... it shouldn't
    // be (re)implemented in a base classe.
    Eigen::Matrix3d rotationMatrixY(double angle) {
        Eigen::Matrix3d rotation;
        rotation << cos(angle), 0, sin(angle),
            0, 1, 0,
            -sin(angle), 0, cos(angle);
        return rotation;
    }

    // TODO(dmarchal: 2024/06/07), this looks like a very common utility function... it shouldn't
    // be (re)implemented in a base classe.
    // the type of the data return should also be unified between rotationMatrixX, Y and Z
    RotMat rotationMatrixZ(double angle) {
        RotMat rotation;
        rotation << cos(angle), -sin(angle), 0,
            sin(angle), cos(angle), 0,
            0, 0, 1;
        return rotation;
    }

protected:
    Data<type::vector<double>>      d_curv_abs_section ;
    Data<type::vector<double>>      d_curv_abs_frames ;
    Data<bool>                      d_debug ;

    /// Input Models container. New inputs are added through addInputModel(In* ).
    core::State<In1>* m_fromModel1;

    // TODO(dmarchal): why this maybe_unused on a data field ?
    [[maybe_unused]] core::State<In2>* m_fromModel2;

    core::State<Out>* m_toModel;

protected:    
    /// Constructor
    BaseCosserat() ;
    /// Destructor
    ~BaseCosserat()  override = default;

    void computeExponentialSE3(const double &x, const Coord1& k, Transform & Trans);

    // TODO(dmarchal: 2024/06/07):
    //   - two naming convention
    //   - unclear the difference between computeAdjoing and buildAdjoint ... is there room for factoring things ?
    void computeAdjoint(const Transform & frame, Mat6x6 &adjoint);
    void computeAdjoint(const Vec6 & frame, Mat6x6 &adjoint);

    void computeCoAdjoint(const Transform & frame, Mat6x6 &coAdjoint);

    void updateExponentialSE3(const In1VecCoord & inDeform);
    void updateTangExpSE3(const In1VecCoord & inDeform);
    void computeTangExp(double & x, const Coord1& k, Mat6x6 & TgX);

    [[maybe_unused]] type::Vec6 computeETA(const Vec6 & baseEta, const In1VecDeriv & k_dot, double abs_input);
    type::Matrix4 computeLogarithm(const double & x, const Mat4x4 &gX);
};


#if !defined(SOFA_COSSERAT_CPP_BaseCosserat)
extern template class SOFA_COSSERAT_API BaseCosserat< sofa::defaulttype::Vec3Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types >;
extern template class SOFA_COSSERAT_API BaseCosserat< sofa::defaulttype::Vec6Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types >;
#endif

} // namespace sofa
