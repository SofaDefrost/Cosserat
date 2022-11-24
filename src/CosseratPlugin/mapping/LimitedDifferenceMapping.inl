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
#include "LimitedDifferenceMapping.h"

#include <sofa/core/Multi2Mapping.inl>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/helper/logging/Message.h>
#include <sofa/type/Quat.h>

#include <string>


namespace sofa::component::mapping
{

    using sofa::core::objectmodel::BaseContext ;
    using sofa::helper::AdvancedTimer;
    using sofa::helper::WriteAccessor;
    using sofa::defaulttype::SolidTypes ;
    using sofa::type::RGBAColor;

template <class TIn1, class TIn2, class TOut>
LimitedDifferenceMapping<TIn1, TIn2, TOut>::LimitedDifferenceMapping()
    : m_toModel(NULL)
    , d_mappedIndices1(initData(&d_mappedIndices1, "mappedIndices1", "Mapped indices in the first input mechanical object \n"))
    , d_mappedIndices2(initData(&d_mappedIndices2, "mappedIndices2", "Mapped indices in the second input mechanical object \n"))
    , d_mappedEntries(initData(&d_mappedEntries, "mappedEntries", "Entries in the objects' type which are taken into account by the mapping. E.g.: [0, 1] maps the first 2 entries of a Vec3.\n"))
    , d_debug(initData(&d_debug, false, "debug", "show debug output.\n"))
{
        d_debug.setValue(false);
}


template <class TIn1, class TIn2, class TOut>
void LimitedDifferenceMapping<TIn1, TIn2, TOut>::init()
{
    if(this->getFromModels1().empty() || this->getFromModels2().empty() || this->getToModels().empty())
    {
        msg_error() << "Error while initializing ; input getFromModels1/getFromModels2/output not found" ;
        return;
    }

    // TO DO: either remove the template and implement only for Vec3, or handle non-vector cases
    const auto in1Dim = Coord1::total_size;
    const auto in2Dim = Coord2::total_size;
    if (in1Dim != in2Dim)
    {
        msg_error("") << "Both inputs should have the same dimension to use this mapping" ;
        return;
    }
    m_inputDimension = in1Dim;

    // Checking that the mapped entries are coherent with the objects dimension
    const auto &mappedEntries = d_mappedEntries.getValue();
    for (unsigned int entryIndex : mappedEntries)
    {
        if (entryIndex >= in1Dim)
        {
            msg_error("") << "Mapped entries indices shouldn't be higher than the input MechanicalObject dimension" ;
            return;
        }
    }
    const auto &m1Indices = d_mappedIndices1.getValue();
    const auto &m2Indices = d_mappedIndices2.getValue();

    if (m1Indices.size() != m2Indices.size())
    {
        msg_info("") << "The same number of indices should be provided for both inputs. Please check the mappedIndices1 and mappedIndices2 parameters" ;
        return;
    }

    m_nbMappedIndices = m1Indices.size();

}


template <class TIn1, class TIn2, class TOut>
void LimitedDifferenceMapping<TIn1, TIn2, TOut>::apply(
        const core::MechanicalParams* /* mparams */,
        const type::vector<OutDataVecCoord*>& dataVecOutPos,
        const type::vector<const In1DataVecCoord*>& dataVecIn1Pos,
        const type::vector<const In2DataVecCoord*>& dataVecIn2Pos)
{

    if(dataVecOutPos.empty() || dataVecIn1Pos.empty() || dataVecIn2Pos.empty())
        return;

    ///Do Apply
    //We need only one input In model and input Root model (if present)
    const In1VecCoord& in1 = dataVecIn1Pos[0]->getValue();
    const In2VecCoord& in2 = dataVecIn2Pos[0]->getValue();

    OutVecCoord& out = *dataVecOutPos[0]->beginEdit();
    out.resize(m_nbMappedIndices);

    const auto &m1Indices = d_mappedIndices1.getValue();
    const auto &m2Indices = d_mappedIndices2.getValue();
    const auto &mappedEntries = d_mappedEntries.getValue();

    for (unsigned int dofIndexIt=0; dofIndexIt < m_nbMappedIndices; dofIndexIt++)
    {
        auto input1DofIndex = m1Indices[dofIndexIt];
        auto input2DofIndex = m2Indices[dofIndexIt];
        for (unsigned int entryIndex : mappedEntries)
        {
            // TO DO: Are other entries set to 0 ?
            out[dofIndexIt][entryIndex] = in2[input2DofIndex][entryIndex] - in1[input1DofIndex][entryIndex]; // Difference in world space coordinates
        }

        if (d_debug.getValue())
        {
            std::cout << " in1 :" << in1[input1DofIndex] << std::endl;
            std::cout << " in2 :" << in2[input2DofIndex] << std::endl;
            std::cout << " out :" << out[dofIndexIt] << std::endl;
        }
    }

    dataVecOutPos[0]->endEdit();
}


template <class TIn1, class TIn2, class TOut>
void LimitedDifferenceMapping<TIn1, TIn2, TOut>:: applyJ(
        const core::MechanicalParams* /* mparams */,
        const type::vector< OutDataVecDeriv*>& dataVecOutVel,
        const type::vector<const In1DataVecDeriv*>& dataVecIn1Vel,
        const type::vector<const In2DataVecDeriv*>& dataVecIn2Vel)
{
    if(dataVecOutVel.empty() || dataVecIn1Vel.empty() ||dataVecIn2Vel.empty() )
        return;

    const In1VecDeriv& in1Vel = dataVecIn1Vel[0]->getValue();
    const In2VecDeriv& in2Vel = dataVecIn2Vel[0]->getValue();
    OutVecDeriv& outVel = *dataVecOutVel[0]->beginEdit();

    // TO DO: is this necessary ?
    outVel.resize(m_nbMappedIndices);

    const auto &m1Indices = d_mappedIndices1.getValue();
    const auto &m2Indices = d_mappedIndices2.getValue();
    const auto &mappedEntries = d_mappedEntries.getValue();

    for (unsigned int dofIndexIt=0; dofIndexIt < m_nbMappedIndices; dofIndexIt++)
    {
        auto input1DofIndex = m1Indices[dofIndexIt];
        auto input2DofIndex = m2Indices[dofIndexIt];
        for (unsigned int entryIndex : mappedEntries)
        {
            // TO DO: Are other entries set to 0 ?
            outVel[dofIndexIt][entryIndex] = in2Vel[input2DofIndex][entryIndex] - in1Vel[input1DofIndex][entryIndex];
        }
    }

    dataVecOutVel[0]->endEdit();
    if (d_debug.getValue()){
        std::cout << " =====> outVel[m1Indices[index]] : " << outVel << std::endl;
    }
}


template <class TIn1, class TIn2, class TOut>
void LimitedDifferenceMapping<TIn1, TIn2, TOut>:: applyJT(
        const core::MechanicalParams* /*mparams*/,
        const type::vector< In1DataVecDeriv*>& dataVecOut1Force,
        const type::vector< In2DataVecDeriv*>& dataVecOut2Force,
        const type::vector<const OutDataVecDeriv*>& dataVecInForce)
{
    if(dataVecOut1Force.empty() || dataVecInForce.empty() || dataVecOut2Force.empty())
        return;

    const OutVecDeriv& inForce = dataVecInForce[0]->getValue();

    In1VecDeriv& out1Force = *dataVecOut1Force[0]->beginEdit();
    In2VecDeriv& out2Force = *dataVecOut2Force[0]->beginEdit();

    //@todo implementation of force modification
    const auto &m1Indices = d_mappedIndices1.getValue();
    const auto &m2Indices = d_mappedIndices2.getValue();
    const auto &mappedEntries = d_mappedEntries.getValue();

    for (unsigned int dofIndexIt=0; dofIndexIt < m_nbMappedIndices; dofIndexIt++)
    {
        auto input1DofIndex = m1Indices[dofIndexIt];
        auto input2DofIndex = m2Indices[dofIndexIt];
        for (unsigned int entryIndex : mappedEntries)
        {
            // TO DO: Are other entries set to 0 ?
            out1Force[input1DofIndex][entryIndex] -= inForce[dofIndexIt][entryIndex];
            out2Force[input2DofIndex][entryIndex] += inForce[dofIndexIt][entryIndex];
        }

    }

//    for (sofa::Index index = 0; index < m_minInd; index++) {
//        getVCenter(     out1Force[m1Indices[index]]) -= getVCenter(     inForce[index]);
//        getVOrientation(out1Force[m1Indices[index]]) -= getVOrientation(inForce[index]);

//        getVCenter(     out2Force[m2Indices[index]]) += getVCenter(     inForce[index]);
//        getVOrientation(out2Force[m2Indices[index]]) += getVOrientation(inForce[index]);
//    }

    dataVecOut1Force[0]->endEdit();
    dataVecOut2Force[0]->endEdit();
}

//___________________________________________________________________________
template <class TIn1, class TIn2, class TOut>
void LimitedDifferenceMapping<TIn1, TIn2, TOut>::applyJT(
        const core::ConstraintParams*/*cparams*/,
        const type::vector< In1DataMatrixDeriv*>&  dataMatOut1Const,
        const type::vector< In2DataMatrixDeriv*>&  dataMatOut2Const,
        const type::vector<const OutDataMatrixDeriv*>& dataMatInConst)
{
    if(dataMatOut1Const.empty() || dataMatOut2Const.empty() || dataMatInConst.empty() )
        return;

    In1MatrixDeriv& out1 = *dataMatOut1Const[0]->beginEdit(); // constraints on the reference frame 1
    In2MatrixDeriv& out2 = *dataMatOut2Const[0]->beginEdit(); // constraints on the reference frame 2
    const OutMatrixDeriv& in = dataMatInConst[0]->getValue(); // input constraints defined on the mapped frames

    const auto &m1Indices = d_mappedIndices1.getValue();
    const auto &m2Indices = d_mappedIndices2.getValue();
    const auto &mappedEntries = d_mappedEntries.getValue();
    typename OutMatrixDeriv::RowConstIterator rowItEnd = in.end();

    for (typename OutMatrixDeriv::RowConstIterator rowIt = in.begin(); rowIt != rowItEnd; ++rowIt)
    {
        typename OutMatrixDeriv::ColConstIterator colIt = rowIt.begin();
//        typename OutMatrixDeriv::ColConstIterator colItEnd = rowIt.end();


        typename In1MatrixDeriv::RowIterator o1 = out1.writeLine(rowIt.index()); // we store the constraint number
        typename In2MatrixDeriv::RowIterator o2 = out2.writeLine(rowIt.index());

        int childIndex = colIt.index();

        //We compute the parents indices
        auto parentIndex1 = m1Indices[childIndex];
        auto parentIndex2 = m2Indices[childIndex];

        const OutDeriv valueConst_ = colIt.val();

//        // Compute the mapped Constraint on the beam nodes
//        Deriv1 direction1;
//        In1::setDPos(direction1,-getVCenter(valueConst_));
//        In1::setDRot(direction1,-getVOrientation(valueConst_));
//        Deriv2 direction2;
//        In2::setDPos(direction2,getVCenter(valueConst_));
//        In2::setDRot(direction2,getVOrientation(valueConst_));

//        if (d_debug.getValue()){
//            printf("1. ======================================================================================= \n");
//            std::cout << "Constraint " << rowIt.index() << " ==> childIndex: "<< childIndex << std::endl;
//            std::cout << "parentIndex1 " << parentIndex1 << " ==> parentIndex2 "<< parentIndex2 << std::endl;
//            std::cout << "valueConst_: "<< valueConst_ << std::endl;
//            std::cout << "direction1: " << direction1 << std::endl;
//            std::cout << "direction2: " << direction2 << std::endl;
//        }

//        o1.addCol(parentIndex1, direction1);
//        o2.addCol(parentIndex2, direction2);
    }
    dataMatOut1Const[0]->endEdit();
    dataMatOut2Const[0]->endEdit();
}

} // namespace sofa::components::mapping
