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
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/type/Transform.h>


namespace Cosserat::mapping
{
using sofa::type::vector;
using sofa::Data;


template <class TIn1, class TIn2, class TOut>
class RigidDistanceMapping : public sofa::core::Multi2Mapping<TIn1, TIn2, TOut>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE3(RigidDistanceMapping, TIn1,TIn2, TOut), SOFA_TEMPLATE3(sofa::core::Multi2Mapping, TIn1, TIn2, TOut) );
    typedef sofa::core::Multi2Mapping<TIn1, TIn2, TOut> Inherit;

    /// Input Model Type
    typedef TIn1 In1;
    typedef TIn2 In2;

    /// Output Model Type
    typedef TOut Out;

    typedef typename In1::Coord       Coord1;
    typedef typename In1::Deriv       Deriv1;
    typedef typename In1::VecCoord    In1VecCoord;
    typedef typename In1::VecDeriv    In1VecDeriv;
    typedef typename In1::MatrixDeriv In1MatrixDeriv;
    typedef Data<In1VecCoord>         In1DataVecCoord;
    typedef Data<In1VecDeriv>         In1DataVecDeriv;
    typedef Data<In1MatrixDeriv>      In1DataMatrixDeriv;
    
    typedef typename In2::Coord::value_type Real;
    typedef typename In2::Coord             Coord2;
    typedef typename In2::Deriv             Deriv2;
    typedef typename In2::VecCoord          In2VecCoord;
    typedef typename In2::VecDeriv          In2VecDeriv;
    typedef typename In2::MatrixDeriv       In2MatrixDeriv;
    typedef Data<In2VecCoord>               In2DataVecCoord;
    typedef Data<In2VecDeriv>               In2DataVecDeriv;
    typedef Data<In2MatrixDeriv>            In2DataMatrixDeriv;
    typedef sofa::type::Mat<4,4,Real>       Mat4x4;

    typedef typename Out::VecCoord OutVecCoord;
    typedef typename Out::Coord OutCoord;
    typedef typename Out::Deriv OutDeriv;
    typedef typename Out::VecDeriv OutVecDeriv;
    typedef typename Out::MatrixDeriv OutMatrixDeriv;
    typedef Data<OutVecCoord> OutDataVecCoord;
    typedef Data<OutVecDeriv> OutDataVecDeriv;
    typedef Data<OutMatrixDeriv> OutDataMatrixDeriv;

    using Transform = sofa::type::Transform<Real>;
    using SpatialVector = sofa::type::SpatialVector<Real>;

protected:
    Data<vector<unsigned int> > d_index1;
    Data<vector<unsigned int> > d_index2;
    Data<Real> d_max;
    Data<Real> d_min;
    Data<Real> d_radius;
    Data<sofa::type::RGBAColor> d_color;
    Data<vector<unsigned int> > d_index;

    sofa::core::State<Out>* m_toModel;

protected:

    RigidDistanceMapping();
    ~RigidDistanceMapping() override = default;

    sofa::Index  m_minInd{};
public:

    /**********************SOFA METHODS**************************/
    void init() override;
    void draw(const sofa::core::visual::VisualParams* /*vparams*/) override {}

    /**********************MAPPING METHODS**************************/
    void apply(
        const sofa::core::MechanicalParams* /* mparams */,
        const vector<OutDataVecCoord*>& dataVecOutPos,
        const vector<const In1DataVecCoord*>& dataVecIn1Pos ,
        const vector<const In2DataVecCoord*>& dataVecIn2Pos) override;

    void applyJ(
        const sofa::core::MechanicalParams* /* mparams */,
        const vector< OutDataVecDeriv*>& dataVecOutVel,
        const vector<const In1DataVecDeriv*>& dataVecIn1Vel,
        const vector<const In2DataVecDeriv*>& dataVecIn2Vel) override;

    //ApplyJT Force
    void applyJT(
        const sofa::core::MechanicalParams* /* mparams */,
        const vector< In1DataVecDeriv*>& dataVecOut1Force,
        const vector< In2DataVecDeriv*>& dataVecOut2RootForce,
        const vector<const OutDataVecDeriv*>& dataVecInForce) override;

    void applyDJT(const sofa::core::MechanicalParams* /*mparams*/,
                  sofa::core::MultiVecDerivId /*inForce*/,
                  sofa::core::ConstMultiVecDerivId /*outForce*/) override{}

    /// This method must be reimplemented by all mappings if they need to support constraints.
    void applyJT(
        const sofa::core::ConstraintParams*  cparams ,
        const vector< In1DataMatrixDeriv*>& dataMatOut1Const  ,
        const vector< In2DataMatrixDeriv*>&  dataMatOut2Const ,
        const vector<const OutDataMatrixDeriv*>&  dataMatInConst) override;

};

#if !defined(SOFA_COSSERAT_CPP_RigidDistanceMapping)
extern template class SOFA_COSSERAT_API RigidDistanceMapping< sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types >;
#endif

} // namespace Cosserat::mapping
