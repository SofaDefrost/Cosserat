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

    using Real = sofa::Real_t<In2>;

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
        const vector<sofa::DataVecCoord_t<Out>*>& dataVecOutPos,
        const vector<const sofa::DataVecCoord_t<In1>*>& dataVecIn1Pos ,
        const vector<const sofa::DataVecCoord_t<In2>*>& dataVecIn2Pos) override;

    void applyJ(
        const sofa::core::MechanicalParams* /* mparams */,
        const vector< sofa::DataVecDeriv_t<Out>*>& dataVecOutVel,
        const vector<const sofa::DataVecDeriv_t<In1>*>& dataVecIn1Vel,
        const vector<const sofa::DataVecDeriv_t<In2>*>& dataVecIn2Vel) override;

    //ApplyJT Force
    void applyJT(
        const sofa::core::MechanicalParams* /* mparams */,
        const vector< sofa::DataVecDeriv_t<In1>*>& dataVecOut1Force,
        const vector< sofa::DataVecDeriv_t<In2>*>& dataVecOut2RootForce,
        const vector<const sofa::DataVecDeriv_t<Out>*>& dataVecInForce) override;

    void applyDJT(const sofa::core::MechanicalParams* /*mparams*/,
                  sofa::core::MultiVecDerivId /*inForce*/,
                  sofa::core::ConstMultiVecDerivId /*outForce*/) override{}

    /// This method must be reimplemented by all mappings if they need to support constraints.
    void applyJT(
        const sofa::core::ConstraintParams*  cparams ,
        const vector< sofa::DataMatrixDeriv_t<In1>*>& dataMatOut1Const  ,
        const vector< sofa::DataMatrixDeriv_t<In2>*>&  dataMatOut2Const ,
        const vector<const sofa::DataMatrixDeriv_t<Out>*>&  dataMatInConst) override;

};

#if !defined(SOFA_COSSERAT_CPP_RigidDistanceMapping)
extern template class SOFA_COSSERAT_API RigidDistanceMapping< sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types >;
#endif

} // namespace Cosserat::mapping
