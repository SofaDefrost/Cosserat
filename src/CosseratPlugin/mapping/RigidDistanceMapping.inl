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
#include "RigidDistanceMapping.h"

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
RigidDistanceMapping<TIn1, TIn2, TOut>::RigidDistanceMapping()
    : m_toModel(NULL)
    , d_index1(initData(&d_index1, "first_point", "index of the first model \n"))
    , d_index2(initData(&d_index2, "second_point", "index of the second model \n"))
    , d_max(initData(&d_max, (Real)1.0e-2, "max", "the maximum of the deformation.\n"))
    , d_min(initData(&d_min, (Real)0.0, "min", "the minimum of the deformation.\n"))
    , d_radius(initData(&d_radius, (Real)3.0, "radius", "the axis in which we want to show the deformation.\n"))
    , d_color(initData(&d_color, type::Vec4f (1, 0., 1., 0.8) ,"color", "The default beam color"))
    , d_index(initData(&d_index, "index", "if this parameter is false, you draw the beam with color "
                                                          "according to the force apply to each beam"))
    , d_debug(initData(&d_debug, false, "debug", "show debug output.\n"))
{
        d_debug.setValue(false);

        this->addUpdateCallback("updateMappedIndices", {&d_index1, &d_index2}, [this](const core::DataTracker& t)
        {
            SOFA_UNUSED(t);
            this->init();
            // Resize the output MechanicalObject
            // TO DO: This callback is developped specifically to answer changes made during a navigation
            // scene of a Cosserat coaxial model. At the moment, the dynamic rediscretisation required by
            // this scenario is done entirely on the curvilinear abscissas of the Cosserat mapping component.
            // A better implementation would be to delegate the discretisation to a topology component, and
            // handle 'topological' changes (on the Mechanical Object in output of this mapping) with the
            // appropriate topology modifications API.
            core::State<Out>* toModel = this->getToModels()[0];
            toModel->resize(m_minInd);

            return sofa::core::objectmodel::ComponentState::Valid;
        }, {});
}


template <class TIn1, class TIn2, class TOut>
void RigidDistanceMapping<TIn1, TIn2, TOut>::init()
{
    if(this->getFromModels1().empty() || this->getFromModels2().empty() || this->getToModels().empty())
    {
        msg_error() << "Error while initializing ; input getFromModels1/getFromModels2/output not found" ;
        return;
    }

    const type::vector<unsigned int> &m1Indices = d_index1.getValue();
    const type::vector<unsigned int> &m2Indices = d_index2.getValue();

    m_minInd = std::min(m1Indices.size(), m2Indices.size());
    if (m_minInd == 0) {
        msg_info("") << " The size of the indices must not be equal to zero" ;
        return;
    }
}


template <class TIn1, class TIn2, class TOut>
void RigidDistanceMapping<TIn1, TIn2, TOut>::apply(
        const core::MechanicalParams* /* mparams */, const type::vector<OutDataVecCoord*>& dataVecOutPos,
        const type::vector<const In1DataVecCoord*>& dataVecIn1Pos ,
        const type::vector<const In2DataVecCoord*>& dataVecIn2Pos)
{

    if(dataVecOutPos.empty() || dataVecIn1Pos.empty() || dataVecIn2Pos.empty())
        return;

    // Checking the componentState, to trigger a callback if other data fields (specifically
    // d_index1 or d_index2) were changed
    if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
        return;

    ///Do Apply
    //We need only one input In model and input Root model (if present)
    const In1VecCoord& in1 = dataVecIn1Pos[0]->getValue();
    const In2VecCoord& in2 = dataVecIn2Pos[0]->getValue();

    OutVecCoord& out = *dataVecOutPos[0]->beginEdit();
    out.resize(m_minInd);

    auto &m1Indices = d_index1.getValue();
    auto &m2Indices = d_index2.getValue();

    size_t baseIndex = 0; // index of the first point of the beam, add this to the data
    //@TODO call the function to Compute frames or update theme here
    //@TODO : m_objects1Frames and m_objects2Frames
    Transform global_Obj1_local = Transform(In1::getCPos(in1[baseIndex]),In1::getCRot(in1[baseIndex]));
    Transform global_Obj2_local = Transform(In2::getCPos(in2[baseIndex]),In2::getCRot(in2[baseIndex]));

    for (sofa::Index pid=0; pid<m_minInd; pid++) {
        int tm1 = m1Indices[pid];
        int tm2 = m2Indices[pid];
        m_vecH.clear();
        //compute the transformation between the two points
        Transform global_H_local1 = Transform(In1::getCPos(in1[tm1]),In1::getCRot(in1[tm1]));
        Transform global_H_local2 = Transform(In2::getCPos(in2[tm2]),In2::getCRot(in2[tm2]));
        Transform Object1_H_Object2 = global_H_local1.inversed()*global_H_local2;

        m_vecH.push_back(Object1_H_Object2);

        std::cout << "The transform is :" << Object1_H_Object2 << std::endl;
        Vec3 outCenter = Object1_H_Object2.getOrigin();
        std::cout << " The center is : " <<  outCenter << std::endl;
        type::Quat outOri = Object1_H_Object2.getOrientation();

        outOri.normalize();
        out[pid] = OutCoord(outCenter,outOri);
        if (d_debug.getValue()){
            std::cout << " in1 :" << in1[tm1] << std::endl;
            std::cout << " in2 :" << in2[tm2] << std::endl;
            std::cout << " out :" << out[pid] << std::endl;
        }
    }
    dataVecOutPos[0]->endEdit();
}


template <class TIn1, class TIn2, class TOut>
void RigidDistanceMapping<TIn1, TIn2, TOut>:: applyJ(
        const core::MechanicalParams* /* mparams */, const type::vector< OutDataVecDeriv*>& dataVecOutVel,
        const type::vector<const In1DataVecDeriv*>& dataVecIn1Vel,
        const type::vector<const In2DataVecDeriv*>& dataVecIn2Vel) {

    if(dataVecOutVel.empty() || dataVecIn1Vel.empty() ||dataVecIn2Vel.empty() )
        return;

    // Checking the componentState, to trigger a callback if other data fields (specifically
    // d_index1 or d_index2) were changed
    if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
        return;

    const In1VecDeriv& in1Vel = dataVecIn1Vel[0]->getValue();
    const In2VecDeriv& in2Vel = dataVecIn2Vel[0]->getValue();
    OutVecDeriv& outVel = *dataVecOutVel[0]->beginEdit();

    const auto &m1Indices = d_index1.getValue();
    const auto &m2Indices = d_index2.getValue();

    SpatialVector vDOF1, vDOF2;
    std::cout << "the size of the Ind is : " << m_minInd << std::endl;
    std::cout << "the size of the outVel is : " << outVel.size() << std::endl;

    for (sofa::Index index = 0; index < m_minInd; index++) {
        getVCenter(outVel[index]) = getVCenter(in2Vel[m2Indices[index]]) - getVCenter(in1Vel[m1Indices[index]]);
        getVOrientation(outVel[index]) =  getVOrientation(in2Vel[m2Indices[index]]) - getVOrientation(in1Vel[m1Indices[index]]) ;
    }
    dataVecOutVel[0]->endEdit();

    // old version
    /*for (sofa::Index index = 0; index < m_minInd; index++) {
        getVCenter(outVel[index]) = getVCenter(in2Vel[m2Indices[index]]) - getVCenter(in1Vel[m1Indices[index]]);
        getVOrientation(outVel[index]) =  getVOrientation(in2Vel[m2Indices[index]]) - getVOrientation(in1Vel[m1Indices[index]]) ;
    }
    dataVecOutVel[0]->endEdit();*/

    if (d_debug.getValue()){
        std::cout << " =====> outVel[m1Indices[index]] : " << outVel << std::endl;
    }
}


template <class TIn1, class TIn2, class TOut>
void RigidDistanceMapping<TIn1, TIn2, TOut>:: applyJT(
        const core::MechanicalParams* /*mparams*/, const type::vector< In1DataVecDeriv*>& dataVecOut1Force,
        const type::vector< In2DataVecDeriv*>& dataVecOut2Force,
        const type::vector<const OutDataVecDeriv*>& dataVecInForce)  {

    if(dataVecOut1Force.empty() || dataVecInForce.empty() || dataVecOut2Force.empty())
        return;

    // Checking the componentState, to trigger a callback if other data fields (specifically
    // d_index1 or d_index2) were changed
    if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
        return;

    const OutVecDeriv& inForce = dataVecInForce[0]->getValue();

    In1VecDeriv& out1Force = *dataVecOut1Force[0]->beginEdit();
    In2VecDeriv& out2Force = *dataVecOut2Force[0]->beginEdit();

    //@todo implementation of force modification
    const auto &m1Indices = d_index1.getValue();
    const auto &m2Indices = d_index2.getValue();

    // Safety check
    // TO DO: is it necessary to raise a warning or an error?
    if (inForce.size() != m_minInd)
        return;

    for (sofa::Index index = 0; index < m_minInd; index++) {
        getVCenter(     out1Force[m1Indices[index]]) -= getVCenter(     inForce[index]);
        getVOrientation(out1Force[m1Indices[index]]) -= getVOrientation(inForce[index]);

        getVCenter(     out2Force[m2Indices[index]]) += getVCenter(     inForce[index]);
        getVOrientation(out2Force[m2Indices[index]]) += getVOrientation(inForce[index]);
    }
    dataVecOut1Force[0]->endEdit();
    dataVecOut2Force[0]->endEdit();
}

//___________________________________________________________________________
template <class TIn1, class TIn2, class TOut>
void RigidDistanceMapping<TIn1, TIn2, TOut>::applyJT(
        const core::ConstraintParams*/*cparams*/ , const type::vector< In1DataMatrixDeriv*>&  dataMatOut1Const,
        const type::vector< In2DataMatrixDeriv*>&  dataMatOut2Const ,
        const type::vector<const OutDataMatrixDeriv*>& dataMatInConst)
{
    if(dataMatOut1Const.empty() || dataMatOut2Const.empty() || dataMatInConst.empty() )
        return;

    // Checking the componentState, to trigger a callback if other data fields (specifically
    // d_index1 or d_index2) were changed
    if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
        return;

    In1MatrixDeriv& out1 = *dataMatOut1Const[0]->beginEdit(); // constraints on the reference frame 1
    In2MatrixDeriv& out2 = *dataMatOut2Const[0]->beginEdit(); // constraints on the reference frame 2
    const OutMatrixDeriv& in = dataMatInConst[0]->getValue(); // input constraints defined on the mapped frames

    const auto &m1Indices = d_index1.getValue();
    const auto &m2Indices = d_index2.getValue();
    typename OutMatrixDeriv::RowConstIterator rowItEnd = in.end();

    for (typename OutMatrixDeriv::RowConstIterator rowIt = in.begin(); rowIt != rowItEnd; ++rowIt) {
        typename OutMatrixDeriv::ColConstIterator colIt = rowIt.begin();
//        typename OutMatrixDeriv::ColConstIterator colItEnd = rowIt.end();


        typename In1MatrixDeriv::RowIterator o1 = out1.writeLine(rowIt.index()); // we store the constraint number
        typename In2MatrixDeriv::RowIterator o2 = out2.writeLine(rowIt.index());

        int childIndex = colIt.index();

        //We compute the parents indices
        auto parentIndex1 = m1Indices[childIndex];
        auto parentIndex2 = m2Indices[childIndex];

        const OutDeriv valueConst_ = colIt.val();

        // Compute the mapped Constraint on the beam nodes
        Deriv1 direction1;
        In1::setDPos(direction1,-getVCenter(valueConst_));
        In1::setDRot(direction1,-getVOrientation(valueConst_));
        Deriv2 direction2;
        In2::setDPos(direction2,getVCenter(valueConst_));
        In2::setDRot(direction2,getVOrientation(valueConst_));

        if (d_debug.getValue()){
            printf("1. ======================================================================================= \n");
            std::cout << "Constraint " << rowIt.index() << " ==> childIndex: "<< childIndex << std::endl;
            std::cout << "parentIndex1 " << parentIndex1 << " ==> parentIndex2 "<< parentIndex2 << std::endl;
            std::cout << "valueConst_: "<< valueConst_ << std::endl;
            std::cout << "direction1: " << direction1 << std::endl;
            std::cout << "direction2: " << direction2 << std::endl;
        }

        o1.addCol(parentIndex1, direction1);
        o2.addCol(parentIndex2, direction2);
    }
    dataMatOut1Const[0]->endEdit();
    dataMatOut2Const[0]->endEdit();
}


template <class TIn1, class TIn2, class TOut>
int RigidDistanceMapping<TIn1, TIn2, TOut>::computeTransform(Transform &global_H0_local,
                                                     Transform &global_H1_local,
                                                     Transform &local0_H_local1,
                                                     Quat<Real> &local_R_local0,
                                                     const Coord1 &x1, const Coord2 &x2)
{
        /// 1. Get the indices of element and nodes
        // @todo: check if we need to do this every time !

        /// 2. Computes the optional rigid transformation of DOF0_Transform_node0 and DOF1_Transform_node1
        // @todo: This part depend on the previous step, so it should be done in the same loop
        Transform OBJ0_H_local0 = Transform(type::Vec3(0,0,0), Rot::identity());
        Transform OBJ1_H_local1 = Transform(type::Vec3(0,0,0), Rot::identity());
        //getDOFtoLocalTransform(edgeInList, DOF0_H_local0,  DOF1_H_local1);

        /// 3. Computes the transformation global To local for both nodes
        Transform global_H_OBJ0(x1.getCenter(), x1.getOrientation());
        Transform global_H_OBJ1(x2.getCenter(), x2.getOrientation());

        /// - add a optional transformation
        Transform global_H_local0 = global_H_OBJ0*OBJ0_H_local0;
        Transform global_H_local1 = global_H_OBJ1*OBJ1_H_local1;


        /// 4. Compute the local frame
        /// SIMPLIFICATION: local = local0:
        local_R_local0.clear();

        global_H_OBJ0.set(type::Vec3(0,0,0), x1.getOrientation());
        global_H_OBJ1.set(type::Vec3(0,0,0), x2.getOrientation());

        /// - rotation due to the optional transformation
        global_H_local0 = global_H_OBJ0*OBJ0_H_local0;
        global_H_local1 = global_H_OBJ1*OBJ1_H_local1;

        global_H0_local = global_H_local0;
        sofa::type::Quat local0_R_local1 = local0_H_local1.getOrientation();
        Transform local0_HR_local1(type::Vec3(0,0,0), local0_R_local1);

        global_H1_local = global_H_local1 * local0_HR_local1.inversed();

        return 1;
    }


template <class TIn1, class TIn2, class TOut>
int RigidDistanceMapping<TIn1, TIn2, TOut>::computeTransform2(unsigned int edgeInList,
                                                      Transform &global_H_local0,
                                                      Transform &global_H_local1,
                                                      const OutVecCoord &x)
{
        /// 1. Get the indices of element and nodes
        unsigned int node0Idx=0, node1Idx=0;
        /*if ( getNodeIndices( edgeInList,  node0Idx, node1Idx ) == -1)
        {
            dmsg_error() << "[computeTransform2] Error in getNodeIndices(). (Aborting)" ;
            return -1;
        }*/

        /// 2. Computes the optional rigid transformation of DOF0_Transform_node0 and DOF1_Transform_node1
        Transform DOF0_H_local0, DOF1_H_local1;
        getDOFtoLocalTransform(edgeInList, DOF0_H_local0,  DOF1_H_local1);

        /// 3. Computes the transformation global To local for both nodes
        Transform global_H_DOF0(x[node0Idx].getCenter(),x[node0Idx].getOrientation());
        Transform global_H_DOF1(x[node1Idx].getCenter(),x[node1Idx].getOrientation());
        /// - add a optional transformation
        global_H_local0 = global_H_DOF0*DOF0_H_local0;
        global_H_local1 = global_H_DOF1*DOF1_H_local1;

        return 1; /// no error
    }
} // namespace sofa::components::mapping
