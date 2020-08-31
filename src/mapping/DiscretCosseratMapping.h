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
#include "BaseCosserat.h"


namespace sofa

{
using sofa::defaulttype::SolidTypes ;
using sofa::core::objectmodel::BaseContext ;
using sofa::defaulttype::Matrix3;
using sofa::defaulttype::Matrix4;
using sofa::defaulttype::Vector3;
using sofa::defaulttype::Vec6;
using std::get;

namespace component
{

namespace mapping
{

/*!
 * \class DiscretCosseratMapping
 * @brief Computes and map the length of the beams
 *
 * This is a component:
 * https://www.sofa-framework.org/community/doc/programming-with-sofa/create-your-component/
 */
using mapping::BaseCosserat;


template <class TIn1, class TIn2, class TOut>
class DiscretCosseratMapping : public core::Multi2Mapping<TIn1, TIn2, TOut>, public component::mapping::BaseCosserat<TIn1, TIn2, TOut>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE3(DiscretCosseratMapping, TIn1,TIn2, TOut), SOFA_TEMPLATE3(core::Multi2Mapping, TIn1, TIn2, TOut) );
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
    typedef typename Out::Coord OutCoord;
    typedef typename Out::Deriv OutDeriv;
    typedef typename Out::VecDeriv OutVecDeriv;
    typedef typename Out::MatrixDeriv OutMatrixDeriv;
    typedef Data<OutVecCoord> OutDataVecCoord;
    typedef Data<OutVecDeriv> OutDataVecDeriv;
    typedef Data<OutMatrixDeriv> OutDataMatrixDeriv;

    typedef MultiLink<DiscretCosseratMapping<In1,In2,Out>, sofa::core::State< In1 >, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> LinkFromModels1;
    typedef MultiLink<DiscretCosseratMapping<In1,In2,Out>, sofa::core::State< In2 >, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> LinkFromModels2;
    typedef MultiLink<DiscretCosseratMapping<In1,In2,Out>, sofa::core::State< Out >, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> LinkToModels;

    typedef typename SolidTypes<Real>::Transform      Transform ;

protected:
    //    Data<helper::vector<double>>      d_curv_abs_input ;
    //    Data<helper::vector<double>>      d_curv_abs_output ;
    //    Data<bool>                        d_debug ;

    /// Input Models container. New inputs are added through addInputModel(In* ).
    //    LinkFromModels1 m_fromModel1;
    //    LinkFromModels2 m_fromModel2;
    //    LinkToModels m_toModel;

    core::State<In1>* m_fromModel1;
    core::State<In2>* m_fromModel2;
    core::State<Out>* m_toModel;

    ////////////////////////// Inherited attributes ////////////////////////////
    /// https://gcc.gnu.org/onlinedocs/gcc/Name-lookup.html
    /// Bring inherited attributes and function in the current lookup context.
    /// otherwise any access to the base::attribute would require
    /// the "this->" approach.
    ///
    using BaseCosserat<TIn1, TIn2, TOut>:: m_indicesVectors ;
    using BaseCosserat<TIn1, TIn2, TOut>::d_curv_abs_input  ;
    using BaseCosserat<TIn1, TIn2, TOut>::d_curv_abs_output ;
    using BaseCosserat<TIn1, TIn2, TOut>::m_nodesTangExpVectors;
    using BaseCosserat<TIn1, TIn2, TOut>::m_nodesVelocityVectors;
    using BaseCosserat<TIn1, TIn2, TOut>::m_ExponentialSE3Vectors;
    using BaseCosserat<TIn1, TIn2, TOut>::m_framesTangExpVectors ;
    using BaseCosserat<TIn1, TIn2, TOut>::m_totalBeamForceVectors ;
    using BaseCosserat<TIn1, TIn2, TOut>::m_nodesExponentialSE3Vectors ;
    using BaseCosserat<TIn1, TIn2, TOut>::d_debug;
    using BaseCosserat<TIn1, TIn2, TOut>::m_vecTransform ;
    using BaseCosserat<TIn1, TIn2, TOut>::m_nodeAjointVectors;
    using BaseCosserat<TIn1, TIn2, TOut>::m_index_input;

    /*===========COSSERAT VECTORS ======================*/
    //    helper::vector<Matrix4> m_nodesLogarithmeSE3Vectors;

protected:
    /// Constructor    
    DiscretCosseratMapping() ;
    /// Destructor
    ~DiscretCosseratMapping()  override {}
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

    void applyDJT(const core::MechanicalParams* /*mparams*/, core::MultiVecDerivId /*inForce*/, core::ConstMultiVecDerivId /*outForce*/) override{}

    /// This method must be reimplemented by all mappings if they need to support constraints.
    virtual void applyJT(
            const core::ConstraintParams*  cparams , const helper::vector< In1DataMatrixDeriv*>& dataMatOut1Const  ,
            const helper::vector< In2DataMatrixDeriv*>&  dataMatOut2Const ,
            const helper::vector<const OutDataMatrixDeriv*>&  dataMatInConst) override;

protected:
    /**********************COSSERAT METHODS**************************/
    //    defaulttype::Matrix4 computeLogarithme(const double & x, const Matrix4 &gX);


};

//extern template class SOFA_POE_MAPPING_API DiscretCosseratMapping<defaulttype::Vec3Types>;


#if  !defined(SOFA_COMPONENT_MAPPING_POE_MAPING_CPP)
extern template class SOFA_COSSERAT_MAPPING_API DiscretCosseratMapping< sofa::defaulttype::Vec3Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types >;
#endif

} // mapping

} // namespace componenet

} // namespace sofa

#endif //SOFA_COMPONENT_MAPPING_POEMAPING_H

