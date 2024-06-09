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
#include <Cosserat/initCosserat.h>

#include <Cosserat/forcefield/BeamHookeLawForceField.h>
#include <Cosserat/mapping/BaseCosserat.h>

#include <sofa/core/BaseMapping.h>
#include <sofa/core/Multi2Mapping.h>
#include <sofa/core/config.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/SolidTypes.h>
#include <sofa/helper/ColorMap.h>

namespace Cosserat::mapping {
using sofa::component::forcefield::BeamHookeLawForceField;
using sofa::core::objectmodel::BaseContext;
using sofa::defaulttype::SolidTypes;
using sofa::type::Matrix3;
using sofa::type::Matrix4;
using sofa::type::Vec3;
using sofa::type::Vec6;
using sofa::type::Quat;
using std::get;
using sofa::Data;

/*!
 * \class DiscreteCosseratMapping
 * @brief Computes and map the length of the beams
 *
 * This is a component:
 * https://www.sofa-framework.org/community/doc/programming-with-sofa/create-your-component/
 */
// using Cosserat::mapping::BaseCosserat;

template <class TIn1, class TIn2, class TOut>
class DiscreteCosseratMapping
    : public sofa::core::Multi2Mapping<TIn1, TIn2, TOut>,
      public Cosserat::mapping::BaseCosserat<TIn1, TIn2, TOut> {
public:
    SOFA_CLASS(SOFA_TEMPLATE3(DiscreteCosseratMapping, TIn1, TIn2, TOut),
               SOFA_TEMPLATE3(sofa::core::Multi2Mapping, TIn1, TIn2, TOut));
    typedef sofa::core::Multi2Mapping<TIn1, TIn2, TOut> Inherit;

    /// Input Model Type
    typedef TIn1 In1;
    typedef TIn2 In2;

    /// Output Model Type
    typedef TOut Out;
    typedef typename In2::Coord::value_type Real1;
    typedef typename In1::Coord Coord1;
    typedef typename In1::Deriv Deriv1;
    typedef typename In1::VecCoord In1VecCoord;
    typedef typename In1::VecDeriv In1VecDeriv;
    typedef typename In1::MatrixDeriv In1MatrixDeriv;
    typedef Data<In1VecCoord> In1DataVecCoord;
    typedef Data<In1VecDeriv> In1DataVecDeriv;
    typedef Data<In1MatrixDeriv> In1DataMatrixDeriv;

    typedef typename In2::Coord::value_type Real2;
    typedef typename In2::Coord Coord2;
    typedef typename In2::Deriv Deriv2;
    typedef typename In2::VecCoord In2VecCoord;
    typedef typename In2::VecDeriv In2VecDeriv;
    typedef typename In2::MatrixDeriv In2MatrixDeriv;
    typedef Data<In2VecCoord> In2DataVecCoord;
    typedef Data<In2VecDeriv> In2DataVecDeriv;
    typedef Data<In2MatrixDeriv> In2DataMatrixDeriv;
    typedef sofa::type::Mat<6, 6, Real2> Mat6x6;
    typedef sofa::type::Mat<3, 6, Real2> Mat3x6;
    typedef sofa::type::Mat<6, 3, Real2> Mat6x3;
    typedef sofa::type::Mat<4, 4, Real2> Mat4x4;

    typedef typename Out::VecCoord OutVecCoord;
    typedef typename Out::Coord OutCoord;
    typedef typename Out::Deriv OutDeriv;
    typedef typename Out::VecDeriv OutVecDeriv;
    typedef typename Out::MatrixDeriv OutMatrixDeriv;
    typedef Data<OutVecCoord> OutDataVecCoord;
    typedef Data<OutVecDeriv> OutDataVecDeriv;
    typedef Data<OutMatrixDeriv> OutDataMatrixDeriv;
    typedef typename SolidTypes<Real2>::Transform Transform;

protected:
    sofa::core::State<In1> *m_fromModel1;
    sofa::core::State<In2> *m_fromModel2;
    sofa::core::State<Out> *m_toModel;
    Data<int> d_deformationAxis;
    Data<Real2> d_max;
    Data<Real2> d_min;
    Data<Real1> d_radius;
    Data<bool> d_drawMapBeam;
    Data<sofa::type::RGBAColor> d_color;
    Data<vector<int>> d_index;
    Data<unsigned int> d_baseIndex;

public:
    typedef sofa::MultiLink<DiscreteCosseratMapping<In1, In2, Out>,
                            sofa::core::State<In1>,
                            sofa::BaseLink::FLAG_STOREPATH | sofa::BaseLink::FLAG_STRONGLINK>
        LinkFromModels1;
    typedef sofa::MultiLink<DiscreteCosseratMapping<In1, In2, Out>,
                            sofa::core::State<In2>,
                            sofa::BaseLink::FLAG_STOREPATH | sofa::BaseLink::FLAG_STRONGLINK>
        LinkFromModels2;
    typedef sofa::MultiLink<DiscreteCosseratMapping<In1, In2, Out>,
                            sofa::core::State<Out>,
                            sofa::BaseLink::FLAG_STOREPATH | sofa::BaseLink::FLAG_STRONGLINK>
        LinkToModels;

protected:
    ////////////////////////// Inherited attributes ////////////////////////////
    /// https://gcc.gnu.org/onlinedocs/gcc/Name-lookup.html
    /// Bring inherited attributes and function in the current lookup context.
    /// otherwise any access to the base::attribute would require
    /// the "this->" approach.
    ///
    using BaseCosserat<TIn1, TIn2, TOut>::m_indicesVectors;
    using BaseCosserat<TIn1, TIn2, TOut>::d_curv_abs_section;
    using BaseCosserat<TIn1, TIn2, TOut>::d_curv_abs_frames;
    using BaseCosserat<TIn1, TIn2, TOut>::m_nodesTangExpVectors;
    using BaseCosserat<TIn1, TIn2, TOut>::m_nodesVelocityVectors;
    using BaseCosserat<TIn1, TIn2, TOut>::m_framesExponentialSE3Vectors;
    using BaseCosserat<TIn1, TIn2, TOut>::m_framesTangExpVectors;
    using BaseCosserat<TIn1, TIn2, TOut>::m_totalBeamForceVectors;
    using BaseCosserat<TIn1, TIn2, TOut>::m_nodesExponentialSE3Vectors;
    using BaseCosserat<TIn1, TIn2, TOut>::d_debug;
    using BaseCosserat<TIn1, TIn2, TOut>::m_vecTransform;
    using BaseCosserat<TIn1, TIn2, TOut>::m_nodeAdjointVectors;
    using BaseCosserat<TIn1, TIn2, TOut>::m_index_input;
    using BaseCosserat<TIn1, TIn2, TOut>::m_indicesVectorsDraw;
    using BaseCosserat<TIn1, TIn2, TOut>::computeTheta;

protected:
    /// Constructor
    DiscreteCosseratMapping();
    /// Destructor
    ~DiscreteCosseratMapping() override {}

public:
    /**********************SOFA METHODS**************************/
    void init() override;
    void reinit() override;
    void draw(const sofa::core::visual::VisualParams *vparams) override;

    /**********************MAPPING METHODS**************************/
    void
    apply(const sofa::core::MechanicalParams * /* mparams */,
          const vector<OutDataVecCoord *> &dataVecOutPos,
          const vector<const In1DataVecCoord *> &dataVecIn1Pos,
          const vector<const In2DataVecCoord *> &dataVecIn2Pos) override;

    void
    applyJ(const sofa::core::MechanicalParams * /* mparams */,
           const vector<OutDataVecDeriv *> &dataVecOutVel,
           const vector<const In1DataVecDeriv *> &dataVecIn1Vel,
           const vector<const In2DataVecDeriv *> &dataVecIn2Vel) override;

    // ApplyJT Force
    void
    applyJT(const sofa::core::MechanicalParams * /* mparams */,
            const vector<In1DataVecDeriv *> &dataVecOut1Force,
            const vector<In2DataVecDeriv *> &dataVecOut2RootForce,
            const vector<const OutDataVecDeriv *> &dataVecInForce) override;

    void applyDJT(const sofa::core::MechanicalParams * /*mparams*/,
                  sofa::core::MultiVecDerivId /*inForce*/,
                  sofa::core::ConstMultiVecDerivId /*outForce*/) override {}

    /// This method must be reimplemented by all mappings if they need to support
    /// constraints.
    virtual void applyJT(
        const sofa::core::ConstraintParams *cparams,
        const vector<In1DataMatrixDeriv *> &dataMatOut1Const,
        const vector<In2DataMatrixDeriv *> &dataMatOut2Const,
        const vector<const OutDataMatrixDeriv *> &dataMatInConst) override;
    void computeBBox(const sofa::core::ExecParams *params, bool onlyVisible) override;
    void computeLogarithm(const double &x, const Matrix4 &gX, Matrix4 &log_gX);

protected:
    sofa::helper::ColorMap m_colorMap;
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
