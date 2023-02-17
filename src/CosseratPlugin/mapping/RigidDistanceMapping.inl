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
    ,d_newVersionOfFrameComputation(initData(&d_newVersionOfFrameComputation, false, "newVersionOfFrameComputation", "if true, the frame is computed with the new version of the code"))
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
            this->Inherit::init();

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
void RigidDistanceMapping<TIn1, TIn2, TOut>::reinit()
{
    // TO DO: This method is implemented as a quickfix, to only trigger an artificial call
    // to the component's callback from a Python code. This is only a quickfix as the correct
    // way to do this would be to implement a Python binding for the callback method.

    // Checking the componentState, to trigger a callback if other data fields (specifically
    // d_index1 or d_index2) were changed
    if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
        return;
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
    const In1VecCoord& input1MOPos = dataVecIn1Pos[0]->getValue();
    const In2VecCoord& input2MOPos = dataVecIn2Pos[0]->getValue();

    OutVecCoord& outputMOPos = *dataVecOutPos[0]->beginEdit();
    outputMOPos.resize(m_minInd);
    pointsR0.resize(m_minInd);

    auto &m1Indices = d_index1.getValue();
    auto &m2Indices = d_index2.getValue();

    size_t baseIndex = 0; // index of the first point of the beam, add this to the data
    mVect_global_H_input2Frame.clear();
    mVect_input1Frame_H_input2Frame.clear();

    m_vecRotation.clear();

    for (sofa::Index pid=0; pid<m_minInd; pid++)
    {
        int tm1 = m1Indices[pid];
        int tm2 = m2Indices[pid];

        Vec3 outCenter ;
        Rot outOri;

        if (d_newVersionOfFrameComputation.getValue())
        {
            // The goal here is to compute the distance in local coordinates
            // We did the abritary choice to compute the distance in the local coordinates of input2Frame

            /* The trasnform global_H_input2Frame converts coordinates from the input2 local frame to
             * the global frame.
             * Here, by using this Transform constructor, with the position and orientation of input2
             * expressed in global coordinates, we imply that the Transform parent frame is the global
             * frame, and that the child frame is input2. See explanation below.
             *
             * NB: we use the naming convention global_H_input2Frame with input2Frame being the child frame,
             * and global the parent frame, as, in this configuration,
             * global_H_input2Frame.getOrientation().rotate() allows to pass from input2Frame axis
             * orientation to the global frame axis orientation. Additionnaly, the multiplication of two
             * Transforms is more readable this way.
             *
             * Explanation: A Transform object is defined based on two fields (cf SolidTypes.h):
             *   orientation_
             *   origin_
             * These two fields correspond to a transformation between a parent frame, and a child frame.
             * /!\ The convention in SOFA is that:
             *   > the orientation_ field stores the orientation of the *child* with respect to the *parent*
             *   > the origin_      field stores the position of the *parent* with respect to the *child*
             * This corresponds to the *Featherstone's convention*.
             * /!\ However, when using the standard constructor signature (providing a position
             * as first parameter, and a rotation as second parameter), the constructor expects
             * arguments which respect the *standard convention*, that is:
             *   > an orientation parameter describing the orientation of the *child* with respect to the *parent*
             *   > a position parameter, also describing the position of the *child* with respect to the *parent*
             * This constructor then sets the internal fields doing :
             *   > orientation_ = orientationParameter
             *   > origin_      = - orientation_.inverseRotate(positionParameter).
             * The only thing left is to access correctly the Transform parameters:
             *   > Transform::getOrientation() returns orientation_ as it is, which means
             *     it describes the orientation of the *child* with respect to the *parent*
             *   > Transform::getOrigin() returns -orientation_.rotate(origin_), which means
             *     it describes the position of the *child* with respect to the *parent* (standard convention)
             *   > Transform::getOriginOfParentInChild() returns origin_ as it is, which means
             *     it describes the position of the *parent* with respect to the *child*
             *
             * Concrete example : here, if input2Frame is simply rotated by +90 degrees around the y axis,
             * the corresponding Rigid, expressed in global frame coordinates, is (rounded):
             * input2MOPos = [0, 0, 0, 0, 0.707107, 0, 0.707107]
             * Considering a unit Vec3 along the X axis, expressed in input2Frame local coordinates :
             * xAxisVectInInput2Frame = [1, 0, 0]
             * Due to the rotation, the expression of this vector is [0, 0, -1] in global frame coordinates.
             * We can obtain these global frame coordinates by doing:
             *   global_H_input2Frame.getOrientation.rotate(xAxisVectInInput2Frame)
             * which will return [0, 0, -1].
             * If we now consider a unit Vec3 along the X axis, expressed in global frame coordinates:
             * xAxisVectInGlobalFrame = [1, 0, 0]
             * Calling:
             *   global_H_input2Frame.getOrientation.inverseRotate(xAxisVectInGlobalFrame)
             * will return the coordinates of the vector in the input2Frame local frame, that is [0, 0, 1]
             *
             * If instead of applying a rotation to input2Frame, we simply translate it so that,
             * in global coordinates : input2MOPos = [1, 0, 0, 0, 0, 0, 1]
             * Then calling :
             *   global_H_input2Frame.getOrientation()
             * will simply return the position of input2Frame (child frame) with respect to the global
             * frame (parent frame), that is : [1, 0, 0, 0, 0, 0, 1]
             */
            Transform global_H_input2Frame = Transform(In2::getCPos(input2MOPos[tm2]), In2::getCRot(input2MOPos[tm2]));
            // We compute a similar Transform for input1Frame
            Transform global_H_input1Frame = Transform(In1::getCPos(input1MOPos[tm1]), In1::getCRot(input1MOPos[tm1]));

            /* Now, we compute a new Transform from the local coordinates of input1Frame to the local
             * coordinates of input2Frame. This Transform can then be used directly to set the output
             * MechanicalObject position and orientation (which are the difference in position and
             * orientation between input1Frame and input2Frame, expressed in input2Frame)
             *
             * NB: in Transform, the * operator is defined so that the result of transform1 * transform2 is
             * a Transform in which the parent frame is transform1's parent frame, and the child frame is
             * transform2's child frame. To make sense, transform1's child frame and transform2's parent frame
             * should be the same.
             */
            Transform input2Frame_H_input1Frame = global_H_input2Frame.inversed() * global_H_input1Frame;

            // Updating the output
            outCenter = input2Frame_H_input1Frame.getOrigin();
            outOri = input2Frame_H_input1Frame.getOrientation();
            outOri.normalize();
            outputMOPos[pid] = OutCoord(outCenter,outOri);

            // Saving the computed Transforms to be used in the rest of the mapping methods
            mVect_global_H_input2Frame.push_back(global_H_input2Frame);
            mVect_input1Frame_H_input2Frame.push_back(input2Frame_H_input1Frame);
        }
        else
        {
            outCenter = (In2::getCPos(input2MOPos[tm2]) - In1::getCPos(input1MOPos[tm1]));
            outOri = (In1::getCRot(input1MOPos[tm1]).inverse()*In2::getCRot(input2MOPos[tm2]));
            outOri.normalize();
            outputMOPos[pid] = OutCoord(outCenter,outOri);
        }

        if (d_debug.getValue())
        {
            std::cout << " in1 :" << input1MOPos[tm1] << std::endl;
            std::cout << " in2 :" << input2MOPos[tm2] << std::endl;
            std::cout << " out :" << outputMOPos[pid] << std::endl;
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

    Vector omega;

    const In1VecDeriv& input1MOVelocities = dataVecIn1Vel[0]->getValue();
    const In2VecDeriv& input2MOVelocities = dataVecIn2Vel[0]->getValue();
    OutVecDeriv& outputVelocities = *dataVecOutVel[0]->beginEdit();

    const auto &m1Indices = d_index1.getValue();
    const auto &m2Indices = d_index2.getValue();

    SpatialVector vDOF1, vDOF2;

    for (sofa::Index index = 0; index < m_minInd; index++)
    {
        if (d_newVersionOfFrameComputation.getValue())
        {
            // Let's compute the velocity of the input1Frame in the local coordinates of the input2Frame

            // Getting the velocities of input1Frame and input2Frame in global coordinates
            auto input1MOVel = getVCenter(input1MOVelocities[m1Indices[index]]);
            auto input2MOVel = getVCenter(input2MOVelocities[m2Indices[index]]);

            // The angular velocity `omega`, of the input1Frame with respect to the input2Frame in global frame is :
            auto omega = getVOrientation(input1MOVelocities[m1Indices[index]]) - getVOrientation(input2MOVelocities[m2Indices[index]]);
            // auto omega = mVect_input1Frame_H_input2Frame[index].getOrientation().rotate(getVOrientation(input1MOVelocities[m1Indices[index]])) ;

            // V1_2 = 1_R_2 * V1_1 + omega x p1
            // auto input1MOVelin2 = mVect_global_H_input2Frame[index].getOrientation().rotate(input2MOVel - input1MOVel)  + cross(omega, getVCenter(input1MOVelocities[m1Indices[index]]));
            auto diffVelIn2 = mVect_global_H_input2Frame[index].getOrientation().inverseRotate(input1MOVel - input2MOVel); //  + cross(omega, mVect_input1Frame_H_input2Frame[index].getOrigin());

            // Due to the fact that the velocity of input2Frame in its frame is zero, we can write :
            getVCenter(outputVelocities[index]) = diffVelIn2;
            getVOrientation(outputVelocities[index]) = mVect_global_H_input2Frame[index].getOrientation().inverseRotate(omega) ;
        }
        else
        {
            getVCenter(outputVelocities[index]) = getVCenter(input2MOVelocities[m2Indices[index]]) - getVCenter(input1MOVelocities[m1Indices[index]]);
            getVOrientation(outputVelocities[index]) = getVOrientation(input2MOVelocities[m2Indices[index]]) - getVOrientation(input1MOVelocities[m1Indices[index]]) ;
        }
    }

    dataVecOutVel[0]->endEdit();

    if (d_debug.getValue())
        std::cout << " =====> outVel[m1Indices[index]] : " << outputVelocities << std::endl;
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

    const OutVecDeriv& outputForces = dataVecInForce[0]->getValue();

    In1VecDeriv& input1MOForces = *dataVecOut1Force[0]->beginEdit();
    In2VecDeriv& input2MOForces = *dataVecOut2Force[0]->beginEdit();

    //@todo implementation of force modification
    const auto &m1Indices = d_index1.getValue();
    const auto &m2Indices = d_index2.getValue();

    // Safety check
    // TO DO: is it necessary to raise a warning or an error?
    if (outputForces.size() != m_minInd)
        return;

    for (sofa::Index outputMOIndex = 0; outputMOIndex < m_minInd; outputMOIndex++)
    {
        if (d_newVersionOfFrameComputation.getValue())
        {
            // The outputForces are computed in the local frames of the input2Frame

            // Computing the forces of the input1Frame and input2Frame in the global coordinates
            Transform global_H_input2Frame = mVect_global_H_input2Frame[outputMOIndex];
            auto force = global_H_input2Frame.getOrientation().rotate(getVCenter(outputForces[outputMOIndex]));
            auto angForce = getVOrientation(outputForces[outputMOIndex]); // + cross(force,mVect_input1Frame_H_input2Frame[outputMOIndex].getOrigin());

            getVCenter(     input1MOForces[m1Indices[outputMOIndex]]) += force;
            getVOrientation(input1MOForces[m1Indices[outputMOIndex]]) += angForce;
            getVCenter(     input2MOForces[m2Indices[outputMOIndex]]) -= force;
            getVOrientation(input2MOForces[m2Indices[outputMOIndex]]) -= angForce;
        }
        else
        {
            getVCenter(input1MOForces[m1Indices[outputMOIndex]]) -= getVCenter(outputForces[outputMOIndex]);
            getVOrientation(input1MOForces[m1Indices[outputMOIndex]]) -= getVOrientation(outputForces[outputMOIndex]);
            getVCenter(input2MOForces[m2Indices[outputMOIndex]]) += getVCenter(outputForces[outputMOIndex]);
            getVOrientation(input2MOForces[m2Indices[outputMOIndex]]) += getVOrientation(outputForces[outputMOIndex]);
        }
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

    for (typename OutMatrixDeriv::RowConstIterator rowIt = in.begin(); rowIt != rowItEnd; ++rowIt)
    {
        typename OutMatrixDeriv::ColConstIterator colIt = rowIt.begin();
        typename In1MatrixDeriv::RowIterator o1 = out1.writeLine(rowIt.index()); // we store the constraint number
        typename In2MatrixDeriv::RowIterator o2 = out2.writeLine(rowIt.index());

        int outputMOIndex = colIt.index();

        //We compute the inputMO indices
        auto input1MOIndex = m1Indices[outputMOIndex];
        auto input2MOIndex = m2Indices[outputMOIndex];

        Deriv1 direction1;
        Deriv2 direction2;

        const OutDeriv valueConst_ = colIt.val();

        if (d_newVersionOfFrameComputation.getValue())
        {
            Transform global_H_input2Frame = mVect_global_H_input2Frame[outputMOIndex];
            auto force = global_H_input2Frame.getOrientation().rotate(getVCenter(valueConst_));
            auto angForce = getVOrientation(valueConst_); //+ cross(force,mVect_input1Frame_H_input2Frame[outputMOIndex].getOrigin());
            In1::setDPos(direction1, force);
            In1::setDRot(direction1, angForce);
            In2::setDPos(direction2, -force);
            In2::setDRot(direction2, -angForce);
        }
        else
        {
            // Compute the mapped Constraint on the beam nodes
            In1::setDPos(direction1,-getVCenter(valueConst_));
            In1::setDRot(direction1,-getVOrientation(valueConst_));
            In2::setDPos(direction2,getVCenter(valueConst_));
            In2::setDRot(direction2,getVOrientation(valueConst_));
        }

        if (d_debug.getValue())
        {
            printf("1. ======================================================================================= \n");
            std::cout << "Constraint " << rowIt.index() << " ==> outputMOIndex: "<< outputMOIndex << std::endl;
            std::cout << "input1MOIndex " << input1MOIndex << " ==> input2MOIndex "<< input2MOIndex << std::endl;
            std::cout << "valueConst_: "<< valueConst_ << std::endl;
            std::cout << "direction1: " << direction1 << std::endl;
            std::cout << "direction2: " << direction2 << std::endl;
        }

        o1.addCol(input1MOIndex, direction1);
        o2.addCol(input2MOIndex, direction2);
    }

    dataMatOut1Const[0]->endEdit();
    dataMatOut2Const[0]->endEdit();
}


template <class TIn1, class TIn2, class TOut>
int RigidDistanceMapping<TIn1, TIn2, TOut>::computeTransform(Transform &global_H_local1,
                                                     Transform &global_H_local2,
                                                     Transform &local1_H_local2,
                                                     Quat<Real> &local1_R_local2,
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
        global_H_local1 = global_H_OBJ0*OBJ0_H_local0;
        global_H_local2 = global_H_OBJ1*OBJ1_H_local1;


        /// 4. Compute the local frame
        /// SIMPLIFICATION: local = local0:
        local1_R_local2.clear();

        global_H_OBJ0.set(type::Vec3(0,0,0), x1.getOrientation());
        global_H_OBJ1.set(type::Vec3(0,0,0), x2.getOrientation());

        /// - rotation due to the optional transformation
        global_H_local1 = global_H_OBJ0*OBJ0_H_local0;
        global_H_local2 = global_H_OBJ1*OBJ1_H_local1;

//        global_H_local1 = global_H_local1;
        sofa::type::Quat local0_R_local1 = local1_H_local2.getOrientation();
        Transform local0_HR_local1(type::Vec3(0,0,0), local0_R_local1);

        global_H_local2 = global_H_local1 * local0_HR_local1.inversed();

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
