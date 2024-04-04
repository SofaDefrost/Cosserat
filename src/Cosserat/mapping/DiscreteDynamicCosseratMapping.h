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

#include <Cosserat/mapping/BaseCosserat.h>

#include <sofa/core/BaseMapping.h>
#include <sofa/core/Multi2Mapping.h>
#include <sofa/defaulttype/SolidTypes.h>
#include <sofa/defaulttype/RigidTypes.h>

namespace sofa::component::mapping
{
using sofa::defaulttype::SolidTypes ;
using sofa::core::objectmodel::BaseContext ;
using sofa::type::Matrix3;
using sofa::type::Matrix4;
using type::Vec3;
using type::Vec6;
using std::get;

/*!
 * \class DiscretDynamicCosseratMapping
 * @brief Computes and map the length of the beams
 *
 * This is a component:
 * https://www.sofa-framework.org/community/doc/programming-with-sofa/create-your-component/
 */


template <class TIn1, class TIn2, class TOut>
class DiscreteDynamicCosseratMapping : public core::Multi2Mapping<TIn1, TIn2, TOut>, public component::mapping::BaseCosserat<TIn1, TIn2, TOut>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE3(DiscreteDynamicCosseratMapping, TIn1,TIn2, TOut), SOFA_TEMPLATE3(core::Multi2Mapping, TIn1, TIn2, TOut) );
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
    typedef type::Mat<6,6,Real> Mat6x6;
    typedef type::Mat<3,6,Real> Mat3x6;
    typedef type::Mat<6,3,Real> Mat6x3;
    typedef type::Mat<4,4,Real> Mat4x4;

    typedef typename Out::VecCoord OutVecCoord;
    typedef typename Out::Coord OutCoord;
    typedef typename Out::Deriv OutDeriv;
    typedef typename Out::VecDeriv OutVecDeriv;
    typedef typename Out::MatrixDeriv OutMatrixDeriv;
    typedef Data<OutVecCoord> OutDataVecCoord;
    typedef Data<OutVecDeriv> OutDataVecDeriv;
    typedef Data<OutMatrixDeriv> OutDataMatrixDeriv;

    typedef MultiLink<DiscreteDynamicCosseratMapping<In1,In2,Out>, sofa::core::State< In1 >,
            BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> LinkFromModels1;
    typedef MultiLink<DiscreteDynamicCosseratMapping<In1,In2,Out>, sofa::core::State< In2 >,
            BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> LinkFromModels2;
    typedef MultiLink<DiscreteDynamicCosseratMapping<In1,In2,Out>, sofa::core::State< Out >,
            BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> LinkToModels;

    typedef typename SolidTypes<Real>::Transform      Transform ;

protected:
    core::State<In1>* m_fromModel1;
    core::State<In2>* m_fromModel2;
    core::State<Out>* m_toModel;

protected:
    /// Constructor
    DiscreteDynamicCosseratMapping() ;
    /// Destructor
    ~DiscreteDynamicCosseratMapping()  override {}

    type::vector<type::vector<Mat6x3>> m_frameJacobienVector;
    type::vector<type::vector<Mat6x3>> m_frameJacobienDotVector;
    type::vector<Mat6x6> m_nodeJacobienVector;
    type::vector<type::vector<Mat6x6>> m_nodeJacobienDotVector;
    type::vector<Matrix3> m_MassExpressionVector;
    Mat6x3 m_matrixBi; // matrixB_i is a constant matrix due to the assumption of constant strain along the section

    ////////////////////////// Inherited attributes ////////////////////////////
    /// https://gcc.gnu.org/onlinedocs/gcc/Name-lookup.html
    /// Bring inherited attributes and function in the current lookup context.
    /// otherwise any access to the base::attribute would require
    /// the "this->" approach.
    ///
    using BaseCosserat<TIn1, TIn2, TOut>:: m_indicesVectors ;
    using BaseCosserat<TIn1, TIn2, TOut>::d_curv_abs_section  ;
    using BaseCosserat<TIn1, TIn2, TOut>::d_curv_abs_frames ;
    using BaseCosserat<TIn1, TIn2, TOut>::m_nodesTangExpVectors;
    using BaseCosserat<TIn1, TIn2, TOut>::m_nodesVelocityVectors;
    using BaseCosserat<TIn1, TIn2, TOut>::m_framesExponentialSE3Vectors;
    using BaseCosserat<TIn1, TIn2, TOut>::m_framesTangExpVectors ;
    using BaseCosserat<TIn1, TIn2, TOut>::m_totalBeamForceVectors ;
    using BaseCosserat<TIn1, TIn2, TOut>::m_nodesExponentialSE3Vectors ;
    using BaseCosserat<TIn1, TIn2, TOut>::d_debug;
    using BaseCosserat<TIn1, TIn2, TOut>::m_vecTransform ;
    using BaseCosserat<TIn1, TIn2, TOut>::m_nodeAdjointVectors;
    using BaseCosserat<TIn1, TIn2, TOut>::m_index_input;


public:


    /**********************SOFA METHODS**************************/
    void init() override;
    virtual void bwdInit() override;  // get the points
    virtual void reset() override;
    virtual void reinit() override;
    void draw(const core::visual::VisualParams* vparams) override;

    /**********************MAPPING METHODS**************************/
    void apply(
            const core::MechanicalParams* /* mparams */, const type::vector<OutDataVecCoord*>& dataVecOutPos,
            const type::vector<const In1DataVecCoord*>& dataVecIn1Pos ,
            const type::vector<const In2DataVecCoord*>& dataVecIn2Pos) override;

    void applyJ(
            const core::MechanicalParams* /* mparams */, const type::vector< OutDataVecDeriv*>& dataVecOutVel,
            const type::vector<const In1DataVecDeriv*>& dataVecIn1Vel,
            const type::vector<const In2DataVecDeriv*>& dataVecIn2Vel) override;

    //ApplyJT Force
    void applyJT(
            const core::MechanicalParams* /* mparams */, const type::vector< In1DataVecDeriv*>& dataVecOut1Force,
            const type::vector< In2DataVecDeriv*>& dataVecOut2RootForce,
            const type::vector<const OutDataVecDeriv*>& dataVecInForce) override;

    void applyDJT(const core::MechanicalParams* /*mparams*/, core::MultiVecDerivId /*inForce*/,
                  core::ConstMultiVecDerivId /*outForce*/) override{}

    /// This method must be reimplemented by all mappings if they need to support constraints.
    virtual void applyJT(
            const core::ConstraintParams*  cparams , const type::vector< In1DataMatrixDeriv*>& dataMatOut1Const  ,
            const type::vector< In2DataMatrixDeriv*>&  dataMatOut2Const ,
            const type::vector<const OutDataMatrixDeriv*>&  dataMatInConst) override;

    /**********************DISCRET DYNAMIC COSSERAT METHODS**************************/

    [[maybe_unused]] void computeMassComponent(const double sectionMass){ SOFA_UNUSED(sectionMass); }
    void computeJ_Jdot_i(const Mat6x6 &Adjoint, size_t frameId, type::vector<Mat6x3> &J_i,
                         const type::Vec6 &etaFrame, type::vector<Mat6x3> &J_dot_i);

};

#if !defined(SOFA_COSSERAT_CPP_DiscreteDynamicCosseratMapping)
extern template class SOFA_COSSERAT_API DiscreteDynamicCosseratMapping< sofa::defaulttype::Vec3Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types >;
#endif

} // namespace sofa::componenet::mapping



