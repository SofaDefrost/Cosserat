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
#include <sofa/helper/ColorMap.h>

namespace Cosserat::mapping
{
namespace
{
using Mat3x6 = sofa::type::Mat<3, 6, SReal>;
using Mat6x3 = sofa::type::Mat<6, 3, SReal>;
using sofa::type::Mat4x4;
using sofa::type::Mat6x6;
using sofa::type::Vec3;
using sofa::type::Vec6;
using sofa::type::Quat;
using sofa::Data;
}

// TODO(dmarchal: 2024/10/07) Is the description valid ? I don't think so.
/*!
 * \class DiscreteCosseratMapping
 * @brief Computes and map the length of the beams
 *
 */
template <class TIn1, class TIn2, class TOut>
class DiscreteCosseratMapping : public BaseCosseratMapping<TIn1, TIn2, TOut> {
public:
    SOFA_CLASS(SOFA_TEMPLATE3(DiscreteCosseratMapping, TIn1, TIn2, TOut),
               SOFA_TEMPLATE3(Cosserat::mapping::BaseCosseratMapping, TIn1, TIn2, TOut));

    /// Input Model Type
    typedef TIn1 In1;
    typedef TIn2 In2;
    typedef TOut Out;

    //////////////////////////////////////////////////////////////////////
    /// @name Data Fields
    /// @{
    Data<int>   d_deformationAxis;
    Data<SReal> d_max;
    Data<SReal> d_min;
    Data<SReal> d_radius;
    Data<bool>  d_drawMapBeam;
    Data<sofa::type::RGBAColor> d_color;
    Data<vector<int>>  d_index;
    Data<unsigned int> d_baseIndex;
    /// @}
    //////////////////////////////////////////////////////////////////////

public:
    //////////////////////////////////////////////////////////////////////
    /// The following methods are inherited from BaseObject
    /// @{
    void doBaseCosseratInit() final override;
    void draw(const sofa::core::visual::VisualParams *vparams) override;
    /// @}
    //////////////////////////////////////////////////////////////////////


    //////////////////////////////////////////////////////////////////////
    /// The following method are inherited from MultiMapping
    /// @{
    auto apply(const sofa::core::MechanicalParams* /* mparams */,
               const vector<sofa::DataVecCoord_t<Out>*>& dataVecOutPos,
               const vector<const sofa::DataVecCoord_t<In1>*>& dataVecIn1Pos,
               const vector<const sofa::DataVecCoord_t<In2>*>& dataVecIn2Pos) ->
        void override;

    void applyJ(const sofa::core::MechanicalParams * /* mparams */,
                const vector<sofa::DataVecDeriv_t<Out> *> &dataVecOutVel,
                const vector<const sofa::DataVecDeriv_t<In1> *> &dataVecIn1Vel,
                const vector<const sofa::DataVecDeriv_t<In2> *> &dataVecIn2Vel) override;

    void applyJT(const sofa::core::MechanicalParams * /* mparams */,
                 const vector<sofa::DataVecDeriv_t<In1> *> &dataVecOut1Force,
                 const vector<sofa::DataVecDeriv_t<In2> *> &dataVecOut2RootForce,
                 const vector<const sofa::DataVecDeriv_t<Out> *> &dataVecInForce) override;

    // TODO(dmarchal:2024/06/13): Override with an empty function is a rare code pattern
    // to make it clear this is the intented and not just an "I'm too lazy to implement it"
    // please always have a precise code comment explaning with details why it is empty.
    void applyDJT(const sofa::core::MechanicalParams * /*mparams*/,
                  sofa::core::MultiVecDerivId /*inForce*/,
                  sofa::core::ConstMultiVecDerivId /*outForce*/) override {}

    /// Support for constraints.
    void applyJT(
            const sofa::core::ConstraintParams *cparams,
            const vector<sofa::DataMatrixDeriv_t<In1> *> &dataMatOut1Const,
            const vector<sofa::DataMatrixDeriv_t<In2> *> &dataMatOut2Const,
            const vector<const sofa::DataMatrixDeriv_t<Out> *> &dataMatInConst) override;
    /// @}
    /////////////////////////////////////////////////////////////////////////////


    void computeBBox(const sofa::core::ExecParams *params, bool onlyVisible) override;
    void computeLogarithm(const double &x, const Mat4x4 &gX, Mat4x4 &log_gX);

protected:
    ////////////////////////// Inherited attributes ////////////////////////////
    /// https://gcc.gnu.org/onlinedocs/gcc/Name-lookup.html
    /// Bring inherited attributes and function in the current lookup context.
    /// otherwise any access to the base::attribute would require
    /// the "this->" approach.
    using BaseCosseratMapping<TIn1, TIn2, TOut>::m_indicesVectors;
    using BaseCosseratMapping<TIn1, TIn2, TOut>::d_curv_abs_section;
    using BaseCosseratMapping<TIn1, TIn2, TOut>::d_curv_abs_frames;
    using BaseCosseratMapping<TIn1, TIn2, TOut>::m_nodesTangExpVectors;
    using BaseCosseratMapping<TIn1, TIn2, TOut>::m_nodesVelocityVectors;
    using BaseCosseratMapping<TIn1, TIn2, TOut>::m_framesExponentialSE3Vectors;
    using BaseCosseratMapping<TIn1, TIn2, TOut>::m_framesTangExpVectors;
    using BaseCosseratMapping<TIn1, TIn2, TOut>::m_totalBeamForceVectors;
    using BaseCosseratMapping<TIn1, TIn2, TOut>::m_nodesExponentialSE3Vectors;
    using BaseCosseratMapping<TIn1, TIn2, TOut>::d_debug;
    using BaseCosseratMapping<TIn1, TIn2, TOut>::m_vecTransform;
    using BaseCosseratMapping<TIn1, TIn2, TOut>::m_nodeAdjointVectors;
    using BaseCosseratMapping<TIn1, TIn2, TOut>::m_indexInput;
    using BaseCosseratMapping<TIn1, TIn2, TOut>::m_indicesVectorsDraw;
    using BaseCosseratMapping<TIn1, TIn2, TOut>::computeTheta;

    using BaseCosseratMapping<TIn1, TIn2, TOut>::m_toModel;
    using BaseCosseratMapping<TIn1, TIn2, TOut>::m_fromModel1;
    using BaseCosseratMapping<TIn1, TIn2, TOut>::m_fromModel2;

    //////////////////////////////////////////////////////////////////////////////

    sofa::helper::ColorMap m_colorMap;
protected:
    DiscreteCosseratMapping();
    ~DiscreteCosseratMapping() override {}
};

#if !defined(SOFA_COSSERAT_CPP_DiscreteCosseratMapping)
extern template class SOFA_COSSERAT_API DiscreteCosseratMapping<
        sofa::defaulttype::Vec3Types, sofa::defaulttype::Rigid3Types,
        sofa::defaulttype::Rigid3Types>;
extern template class SOFA_COSSERAT_API DiscreteCosseratMapping<
        sofa::defaulttype::Vec6Types, sofa::defaulttype::Rigid3Types,
        sofa::defaulttype::Rigid3Types>;
#endif

} // namespace Cosserat::mapping
