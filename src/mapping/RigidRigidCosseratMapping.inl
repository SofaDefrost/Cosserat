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

#include <sofa/core/Multi2Mapping.inl>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <string>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/helper/logging/Message.h>

#include "sofa/defaulttype/Quat.h"
#include "BaseCosserat.inl"
#include "RigidRigidCosseratMapping.h"


namespace sofa::component::mapping
{

    using sofa::core::objectmodel::BaseContext ;
    using sofa::helper::AdvancedTimer;
    using sofa::helper::WriteAccessor;
    using sofa::defaulttype::SolidTypes ;
    using sofa::type::RGBAColor;

template <class TIn1, class TIn2, class TOut>
RigidRigidCosseratMapping<TIn1, TIn2, TOut>::RigidRigidCosseratMapping()
    : m_fromModel1(NULL)
    , m_fromModel2(NULL)
    , m_toModel(NULL)
    , d_deformationAxis(initData(&d_deformationAxis, (int)1, "deformationAxis",
                                 "the axis in which we want to show the deformation.\n"))
    , d_indice1(initData(&d_indice1, "first_point", "index of the first model \n"))
    , d_indice2(initData(&d_indice2, "second_point", "index of the second model \n"))
    , d_max(initData(&d_max, (Real)1.0e-2, "max",
                                 "the maximum of the deformation.\n"))
    , d_min(initData(&d_min, (Real)0.0, "min",
                                 "the minimum of the deformation.\n"))

    , d_radius(initData(&d_radius, (Real)3.0, "radius",
                                 "the axis in which we want to show the deformation.\n"))
    , d_drawMapBeam(initData(&d_drawMapBeam, true,"nonColored", "if this parameter is false, you draw the beam with "
                                                                "color according to the force apply to each beam"))
    , d_color(initData(&d_color, defaulttype::Vec4f (1, 0., 1., 0.8) ,"color", "The default beam color"))
    , d_index(initData(&d_index, "index", "if this parameter is false, you draw the beam with color "
                                                          "according to the force apply to each beam"))
{
}


// _________________________________________________________________________________________

template <class TIn1, class TIn2, class TOut>
void RigidRigidCosseratMapping<TIn1, TIn2, TOut>::init()
{
    if(this->getFromModels1().empty() || this->getFromModels2().empty() || this->getToModels().empty())
    {
        msg_error() << "Error while initializing ; input getFromModels1/getFromModels2/output not found" ;
        return;
    }


}


template <class TIn1, class TIn2, class TOut>
void RigidRigidCosseratMapping<TIn1, TIn2, TOut>::apply(
        const core::MechanicalParams* /* mparams */, const type::vector<OutDataVecCoord*>& dataVecOutPos,
        const type::vector<const In1DataVecCoord*>& dataVecIn1Pos ,
        const type::vector<const In2DataVecCoord*>& dataVecIn2Pos)
{

    if(dataVecOutPos.empty() || dataVecIn1Pos.empty() || dataVecIn2Pos.empty())
        return;

    ///Do Apply
    //We need only one input In model and input Root model (if present)
    const In1VecCoord& in1 = dataVecIn1Pos[0]->getValue();
    const In2VecCoord& in2 = dataVecIn2Pos[0]->getValue();

    OutVecCoord& out = *dataVecOutPos[0]->beginEdit();

    const type::vector<int> &m1Indices = d_indice1.getValue();
    const type::vector<int> &m2Indices = d_indice2.getValue();

    auto minInd = std::min(m1Indices.size(), m2Indices.size());
    if (minInd == 0) {
        return;
    }
    out.resize(minInd);

    for (unsigned pid=0; pid<minInd; pid++) {
        int tm1 = m1Indices[pid];
        int tm2 = m2Indices[pid];
//        std::cout<< " tm1: " << tm1 << " tm2: " << tm2 << std::endl;
//        std::cout << " in1 :" << in1[tm1] << std::endl;
//        std::cout << " in2 :" << in2[tm2] << std::endl;
        Vector3 outCenter = in2[tm2].getCenter()-in1[tm1].getCenter();
        defaulttype::Quat outOri = in2[tm2].getOrientation().inverse() * in1[tm1].getOrientation();
        out[pid] = OutCoord(outCenter,outOri);
//        std::cout<< "Center: " << outCenter << std::endl;
//        std::cout<< "Ori: " << outOri << std::endl;
    }

//
//    // update the Exponential Matrices according to new deformation
//    // Here we update m_framesExponentialSE3Vectors & m_nodesExponentialSE3Vectors
//    /* Go from Cossserat to SOFA frame*/
//    Transform frame0 = Transform(In2::getCPos(in2[0]),In2::getCRot(in2[0]));
//    for(auto i=0; i<sz; i++){
//        Transform frame = frame0;
//        for (auto u = 0; u < m_indicesVectors[i]; u++) {
//            frame *= m_nodesExponentialSE3Vectors[u];
//        }
//        frame *= m_framesExponentialSE3Vectors[i];
//
//        Vector3 v = frame.getOrigin();
//        defaulttype::Quat q = frame.getOrientation();
//        out[i] = OutCoord(v,q);
//    }
//    // @todo do this another place
//    m_index_input = 0;
//    dataVecOutPos[0]->endEdit();
}


template <class TIn1, class TIn2, class TOut>
void RigidRigidCosseratMapping<TIn1, TIn2, TOut>:: applyJ(
        const core::MechanicalParams* /* mparams */, const type::vector< OutDataVecDeriv*>& dataVecOutVel,
        const type::vector<const In1DataVecDeriv*>& dataVecIn1Vel,
        const type::vector<const In2DataVecDeriv*>& dataVecIn2Vel) {

//    if(dataVecOutVel.empty() || dataVecIn1Vel.empty() ||dataVecIn2Vel.empty() )
//        return;
//    const In1VecDeriv& in1 = dataVecIn1Vel[0]->getValue();
//    const In2VecDeriv& in2_vecDeriv = dataVecIn2Vel[0]->getValue();
//    OutVecDeriv& outVel = *dataVecOutVel[0]->beginEdit();

}


template <class TIn1, class TIn2, class TOut>
void RigidRigidCosseratMapping<TIn1, TIn2, TOut>:: applyJT(
        const core::MechanicalParams* /*mparams*/, const type::vector< In1DataVecDeriv*>& dataVecOut1Force,
        const type::vector< In2DataVecDeriv*>& dataVecOut2Force,
        const type::vector<const OutDataVecDeriv*>& dataVecInForce)  {

    if(dataVecOut1Force.empty() || dataVecInForce.empty() || dataVecOut2Force.empty())
        return;

//    const OutVecDeriv& in = dataVecInForce[0]->getValue();
//
//    In1VecDeriv& out1 = *dataVecOut1Force[0]->beginEdit();
//    In2VecDeriv& out2 = *dataVecOut2Force[0]->beginEdit();
//
//
//    dataVecOut1Force[0]->endEdit();
//    dataVecOut2Force[0]->endEdit();
}

//___________________________________________________________________________
template <class TIn1, class TIn2, class TOut>
void RigidRigidCosseratMapping<TIn1, TIn2, TOut>::applyJT(
        const core::ConstraintParams*/*cparams*/ , const type::vector< In1DataMatrixDeriv*>&  dataMatOut1Const,
        const type::vector< In2DataMatrixDeriv*>&  dataMatOut2Const ,
        const type::vector<const OutDataMatrixDeriv*>& dataMatInConst)
{
    if(dataMatOut1Const.empty() || dataMatOut2Const.empty() || dataMatInConst.empty() )
        return;

    //We need only one input In model and input Root model (if present)
//    In1MatrixDeriv& out1 = *dataMatOut1Const[0]->beginEdit(); // constraints on the strain space (reduced coordinate)
//    In2MatrixDeriv& out2 = *dataMatOut2Const[0]->beginEdit(); // constraints on the reference frame (base frame)
//    const OutMatrixDeriv& in = dataMatInConst[0]->getValue(); // input constraints defined on the mapped frames
//
//    const OutVecCoord& frame = m_toModel->read(core::ConstVecCoordId::position())->getValue();
//    const In1DataVecCoord* x1fromData = m_fromModel1->read(core::ConstVecCoordId::position());
//    const In1VecCoord x1from = x1fromData->getValue();
//
//
//    ////// END ARTICULATION SYSTEM MAPPING
//    dataMatOut1Const[0]->endEdit();
//    dataMatOut2Const[0]->endEdit();
}


template <class TIn1, class TIn2, class TOut>
void RigidRigidCosseratMapping<TIn1, TIn2, TOut>::draw(const core::visual::VisualParams* vparams)
{
    //if(!d_debug.getValue()) return;


}

} // namespace sofa::components::mapping
