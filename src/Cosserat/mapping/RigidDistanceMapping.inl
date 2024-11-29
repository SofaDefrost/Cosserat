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
#include <Cosserat/mapping/RigidDistanceMapping.h>
#include <sofa/core/Multi2Mapping.inl>

namespace Cosserat::mapping
{

using sofa::Coord_t;
using sofa::Deriv_t;
using sofa::VecCoord_t;
using sofa::VecDeriv_t;
using sofa::MatrixDeriv_t;
using sofa::DataVecCoord_t;
using sofa::DataVecDeriv_t;
using sofa::DataMatrixDeriv_t;

template <class TIn1, class TIn2, class TOut>
RigidDistanceMapping<TIn1, TIn2, TOut>::RigidDistanceMapping()
    : d_index1(initData(&d_index1, "first_point", "index of the first model \n"))
    , d_index2(initData(&d_index2, "second_point", "index of the second model \n"))
    , d_max(initData(&d_max, (Real)1.0e-2, "max", "the maximum of the deformation.\n"))
    , d_min(initData(&d_min, (Real)0.0, "min", "the minimum of the deformation.\n"))
    , d_radius(initData(&d_radius, (Real)3.0, "radius", "the axis in which we want to show the deformation.\n"))
    , d_color(initData(&d_color, sofa::type::RGBAColor(1.f, 0.f, 1.f, 0.8f) ,"color", "The default beam color"))
    , d_index(initData(&d_index, "index", "if this parameter is false, you draw the beam with color "
                                          "according to the force apply to each beam"))
    , m_toModel(NULL)
{}


template <class TIn1, class TIn2, class TOut>
void RigidDistanceMapping<TIn1, TIn2, TOut>::init()
{
    if(this->getFromModels1().empty() || this->getFromModels2().empty() || this->getToModels().empty())
    {
        msg_error() << "Error while initializing ; input getFromModels1/getFromModels2/output not found" ;
        return;
    }

    const sofa::type::vector<unsigned int> &m1Indices = d_index1.getValue();
    const sofa::type::vector<unsigned int> &m2Indices = d_index2.getValue();

    m_minInd = std::min(m1Indices.size(), m2Indices.size());
    if (m_minInd == 0)
    {
        msg_error() << " The size of the indices must not be equal to zero" ;
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
        return;
    }

    if (this->getToModels().empty())
    {
        msg_error() << "Output of mapping is empty";
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
        return;
    }

    m_toModel = this->getToModels()[0];
    m_toModel->resize(m_minInd);

    Inherit1::init();

    if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Invalid)
    {
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Valid);
    }
}


template <class TIn1, class TIn2, class TOut>
void RigidDistanceMapping<TIn1, TIn2, TOut>::apply(
    const sofa::core::MechanicalParams* /* mparams */,
    const vector<DataVecCoord_t<Out>*>& dataVecOutPos,
    const vector<const DataVecCoord_t<In1>*>& dataVecIn1Pos ,
    const vector<const DataVecCoord_t<In2>*>& dataVecIn2Pos)
{

    if(dataVecOutPos.empty() || dataVecIn1Pos.empty() || dataVecIn2Pos.empty())
        return;

    ///Do Apply
    //We need only one input In model and input Root model (if present)
    const VecCoord_t<In1>& in1 = dataVecIn1Pos[0]->getValue();
    const VecCoord_t<In2>& in2 = dataVecIn2Pos[0]->getValue();

    auto out = sofa::helper::getWriteOnlyAccessor(*dataVecOutPos[0]);
    out.resize(m_minInd);

    const auto& m1Indices = d_index1.getValue();
    const auto& m2Indices = d_index2.getValue();

    for (sofa::Index pid = 0; pid < m_minInd; ++pid)
    {
        const int tm1 = m1Indices[pid];
        const int tm2 = m2Indices[pid];
        const auto outCenter = in2[tm2].getCenter() - in1[tm1].getCenter();

        sofa::type::Quat outOri = in2[tm2].getOrientation()* in1[tm1].getOrientation().inverse();
        outOri.normalize();

        out[pid] = Coord_t<Out>(outCenter, outOri); // This difference is in the word space
    }
}


template <class TIn1, class TIn2, class TOut>
void RigidDistanceMapping<TIn1, TIn2, TOut>:: applyJ(
    const sofa::core::MechanicalParams* /* mparams */,
    const vector< DataVecDeriv_t<Out>*>& dataVecOutVel,
    const vector<const DataVecDeriv_t<In1>*>& dataVecIn1Vel,
    const vector<const DataVecDeriv_t<In2>*>& dataVecIn2Vel) {

    if(dataVecOutVel.empty() || dataVecIn1Vel.empty() ||dataVecIn2Vel.empty() )
        return;

    const VecDeriv_t<In1>& in1Vel = dataVecIn1Vel[0]->getValue();
    const VecDeriv_t<In2>& in2Vel = dataVecIn2Vel[0]->getValue();

    auto outVel = sofa::helper::getWriteOnlyAccessor(*dataVecOutVel[0]);

    const auto &m1Indices = d_index1.getValue();
    const auto &m2Indices = d_index2.getValue();

    for (sofa::Index index = 0; index < m_minInd; index++)
    {
        getVCenter(outVel[index]) = getVCenter(in2Vel[m2Indices[index]]) - getVCenter(in1Vel[m1Indices[index]]);
        getVOrientation(outVel[index]) =  getVOrientation(in2Vel[m2Indices[index]]) - getVOrientation(in1Vel[m1Indices[index]]) ;
    }
}


template <class TIn1, class TIn2, class TOut>
void RigidDistanceMapping<TIn1, TIn2, TOut>:: applyJT(
    const sofa::core::MechanicalParams* /*mparams*/,
    const vector< DataVecDeriv_t<In1>*>& dataVecOut1Force,
    const vector< DataVecDeriv_t<In2>*>& dataVecOut2Force,
    const vector<const DataVecDeriv_t<Out>*>& dataVecInForce)  {

    if(dataVecOut1Force.empty() || dataVecInForce.empty() || dataVecOut2Force.empty())
        return;

    const VecDeriv_t<Out>& inForce = dataVecInForce[0]->getValue();

    auto out1Force = sofa::helper::getWriteOnlyAccessor(*dataVecOut1Force[0]);
    auto out2Force = sofa::helper::getWriteOnlyAccessor(*dataVecOut2Force[0]);

    //@todo implementation of force modification
    const auto &m1Indices = d_index1.getValue();
    const auto &m2Indices = d_index2.getValue();

    for (sofa::Index index = 0; index < m_minInd; index++)
    {
        getVCenter(     out1Force[m1Indices[index]]) -= getVCenter(     inForce[index]);
        getVOrientation(out1Force[m1Indices[index]]) -= getVOrientation(inForce[index]);

        getVCenter(     out2Force[m2Indices[index]]) += getVCenter(     inForce[index]);
        getVOrientation(out2Force[m2Indices[index]]) += getVOrientation(inForce[index]);
    }
}

//___________________________________________________________________________
template <class TIn1, class TIn2, class TOut>
void RigidDistanceMapping<TIn1, TIn2, TOut>::applyJT(
    const sofa::core::ConstraintParams*/*cparams*/ ,
    const vector< DataMatrixDeriv_t<In1>*>&  dataMatOut1Const,
    const vector< DataMatrixDeriv_t<In2>*>&  dataMatOut2Const ,
    const vector<const DataMatrixDeriv_t<Out>*>& dataMatInConst)
{
    if(dataMatOut1Const.empty() || dataMatOut2Const.empty() || dataMatInConst.empty() )
        return;

    auto out1 = sofa::helper::getWriteOnlyAccessor(*dataMatOut1Const[0]); // constraints on the reference frame 1
    auto out2 = sofa::helper::getWriteOnlyAccessor(*dataMatOut2Const[0]); // constraints on the reference frame 2

    const MatrixDeriv_t<Out>& in = dataMatInConst[0]->getValue(); // input constraints defined on the mapped frames

    const auto &m1Indices = d_index1.getValue();
    const auto &m2Indices = d_index2.getValue();
    typename MatrixDeriv_t<Out>::RowConstIterator rowItEnd = in.end();

    for (typename MatrixDeriv_t<Out>::RowConstIterator rowIt = in.begin(); rowIt != rowItEnd; ++rowIt) {
        typename MatrixDeriv_t<Out>::ColConstIterator colIt = rowIt.begin();
        //        typename MatrixDeriv_t<Out>::ColConstIterator colItEnd = rowIt.end();


        auto o1 = out1->writeLine(rowIt.index()); // we store the constraint number
        auto o2 = out2->writeLine(rowIt.index());

        int childIndex = colIt.index();

        //We compute the parents indices
        auto parentIndex1 = m1Indices[childIndex];
        auto parentIndex2 = m2Indices[childIndex];

        const Deriv_t<Out> valueConst_ = colIt.val();

        // Compute the mapped Constraint on the beam nodes
        Deriv_t<In1> direction1;
        In1::setDPos(direction1,-getVCenter(valueConst_));
        In1::setDRot(direction1,-getVOrientation(valueConst_));
        Deriv_t<In2> direction2;
        In2::setDPos(direction2,getVCenter(valueConst_));
        In2::setDRot(direction2,getVOrientation(valueConst_));

        o1.addCol(parentIndex1, direction1);
        o2.addCol(parentIndex2, direction2);
    }
}
} // namespace Cosserat::mapping
