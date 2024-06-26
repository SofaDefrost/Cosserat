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
#include <Cosserat/mapping/BaseCosseratMapping.h>

#include <sofa/core/BaseMapping.h>
#include <sofa/core/config.h>
#include <sofa/core/Multi2Mapping.h>
#include <sofa/defaulttype/SolidTypes.h>
#include <sofa/defaulttype/RigidTypes.h>


namespace Cosserat::mapping
{
using sofa::defaulttype::SolidTypes ;
using sofa::core::objectmodel::BaseContext ;
using sofa::type::Matrix3;
using sofa::type::Matrix4;
using sofa::type::Vec3;
using sofa::type::Vec6;
using sofa::type::Quat;
using std::get;
using sofa::type::vector;
using Cosserat::mapping::BaseCosseratMapping;

/*!
 * \class DifferenceMultiMapping
 */
template <class TIn1, class TIn2, class TOut>
class DifferenceMultiMapping : public sofa::core::Multi2Mapping<TIn1, TIn2, TOut>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE3(DifferenceMultiMapping, TIn1,TIn2, TOut),
               SOFA_TEMPLATE3(sofa::core::Multi2Mapping, TIn1, TIn2, TOut) );
    typedef sofa::core::Multi2Mapping<TIn1, TIn2, TOut> Inherit;

    /// Input Model Type
    typedef TIn1 In1;
    typedef TIn2 In2;

    /// Output Model Type
    typedef TOut Out;

    typedef typename In1::Coord Coord1;
    typedef typename In1::Deriv Deriv1;
    typedef typename In1::VecCoord In1VecCoord;
    typedef typename In1::VecDeriv In1VecDeriv;
    typedef typename In1::MatrixDeriv In1MatrixDeriv;
    typedef sofa::Data<In1VecCoord> In1DataVecCoord;
    typedef sofa::Data<In1VecDeriv> In1DataVecDeriv;
    typedef sofa::Data<In1MatrixDeriv> In1DataMatrixDeriv;

    typedef sofa::defaulttype::Rigid3dTypes::Coord Rigid;
    
    typedef typename In2::Coord::value_type Real          ;
    typedef typename In2::Coord             Coord2         ;
    typedef typename In2::Deriv             Deriv2         ;
    typedef typename In2::VecCoord In2VecCoord;
    typedef typename In2::VecDeriv In2VecDeriv;
    typedef typename In2::MatrixDeriv In2MatrixDeriv;
    typedef sofa::Data<In2VecCoord> In2DataVecCoord;
    typedef sofa::Data<In2VecDeriv> In2DataVecDeriv;
    typedef sofa::Data<In2MatrixDeriv> In2DataMatrixDeriv;
    typedef sofa::type::Mat<4,4,Real> Mat4x4;

    typedef typename Out::VecCoord OutVecCoord;
    typedef typename Out::Coord OutCoord;
    typedef typename Out::Deriv OutDeriv;
    typedef typename Out::VecDeriv OutVecDeriv;
    typedef typename Out::MatrixDeriv OutMatrixDeriv;
    typedef sofa::Data<OutVecCoord> OutDataVecCoord;
    typedef sofa::Data<OutVecDeriv> OutDataVecDeriv;
    typedef sofa::Data<OutMatrixDeriv> OutDataMatrixDeriv;

    typedef typename SolidTypes<Real>::Transform      Transform ;

public:
    /********************** The component Data **************************/
    //Input data
    sofa::Data<vector<Rigid>>                 d_direction;
    sofa::Data<vector<unsigned int>>          d_indices;
    sofa::Data<double>                        d_radius;
    sofa::Data<sofa::type::Vec4f>             d_color;
    sofa::Data<bool>                          d_drawArrows;
    sofa::Data<bool>                          d_lastPointIsFixed;

protected:   
    sofa::core::State<In1>* m_fromModel1;
    sofa::core::State<In2>* m_fromModel2;
    sofa::core::State<Out>* m_toModel;

protected:
    /// Constructor
    DifferenceMultiMapping() ;
    /// Destructor
    ~DifferenceMultiMapping()  override = default;

public:
    /**********************SOFA METHODS**************************/
    void init() override;
    void bwdInit() override;  // get the points
    void reset() override;
    void reinit() override;
    void draw(const sofa::core::visual::VisualParams* vparams) override;

    /**********************MAPPING METHODS**************************/
    void apply(
        const sofa::core::MechanicalParams* /* mparams */, const vector<OutDataVecCoord*>& dataVecOutPos,
        const vector<const In1DataVecCoord*>& dataVecIn1Pos ,
        const vector<const In2DataVecCoord*>& dataVecIn2Pos) override;

    void applyJ(
        const sofa::core::MechanicalParams* /* mparams */, const vector< OutDataVecDeriv*>& dataVecOutVel,
        const vector<const In1DataVecDeriv*>& dataVecIn1Vel,
        const vector<const In2DataVecDeriv*>& dataVecIn2Vel) override;

    //ApplyJT Force
    void applyJT(
        const sofa::core::MechanicalParams* /* mparams */, const vector< In1DataVecDeriv*>& dataVecOut1Force,
        const vector< In2DataVecDeriv*>& dataVecOut2RootForce,
        const vector<const OutDataVecDeriv*>& dataVecInForce) override;

    void applyDJT(const sofa::core::MechanicalParams* /*mparams*/, sofa::core::MultiVecDerivId /*inForce*/, sofa::core::ConstMultiVecDerivId /*outForce*/) override{}

    /// This method must be reimplemented by all mappings if they need to support constraints.
    void applyJT(
        const sofa::core::ConstraintParams*  cparams , const vector< In1DataMatrixDeriv*>& dataMatOut1Const  ,
        const vector< In2DataMatrixDeriv*>&  dataMatOut2Const ,
        const vector<const OutDataMatrixDeriv*>&  dataMatInConst) override;

    /**********************MAPPING METHODS**************************/
    void initiateTopologies();
    void computeProximity(const In1VecCoord &x1, const In2VecCoord &x2);

    void computeNeedleProximity(const In1VecCoord &x1, const In2VecCoord &x2);

    /**********************Useful METHODS**************************/
    void addPointProcess(){
        msg_warning("DifferenceMultiMapping")<< "The point you are adding is :"; //<< pointPos ;
    }

private:

    typedef struct {
        double fact;
        int p1, p2;
        int eid, cid;
        Real alpha, beta;
        Coord1 proj, Q;
        OutDeriv dirAxe, t1, t2;
        Real r, r2, Q1Q2;
        Real dist; //violation
        Real thirdConstraint;
    } Constraint;

    vector<Constraint> m_constraints;
};

} // sofa::component::mapping
