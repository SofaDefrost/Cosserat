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
#include "DifferenceMultiMapping.h"

#include <sofa/core/Multi2Mapping.inl>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/helper/logging/Message.h>
#include <sofa/type/RGBAColor.h>

#include <string>

namespace sofa::component::mapping
{
    using sofa::core::objectmodel::BaseContext;
    using sofa::helper::AdvancedTimer;
    using sofa::helper::WriteAccessor;
    using sofa::type::RGBAColor;

    template <class TIn1, class TIn2, class TOut>
    DifferenceMultiMapping<TIn1, TIn2, TOut>::DifferenceMultiMapping()
        : d_direction(initData(&d_direction, "direction", "The list of directions of fix points .\n")),
          d_indices(initData(&d_indices, "indices", "Indices of fixe points of the cable")),
          d_radius(initData(&d_radius, 2.0, "radius", "The size of the cable")),
          d_color(initData(&d_color, type::Vec4f(1, 0, 0, 1), "color", "The color of the cable")),
          d_drawArrows(initData(&d_drawArrows, false, "drawArrows", "The color of the cable")),
          d_lastPointIsFixed(initData(&d_lastPointIsFixed, true, "lastPointIsFixed", "This select the last point as fixed of not,"
                                                                                     "one.")),
          m_fromModel1(NULL), m_fromModel2(NULL), m_toModel(NULL)
    {
    }

    template <class TIn1, class TIn2, class TOut>
    void DifferenceMultiMapping<TIn1, TIn2, TOut>::initiateTopologies()
    {
        m_toModel = this->getToModels()[0];
        if (!m_toModel)
        {
            std::cout << " No output mechanical state found. Consider setting the "
                      << this->toModels.getName() << " attribute." << std::endl;
            return;
        }

        if (!d_direction.isSet())
            msg_warning() << "No direction nor indices is given.";

    }

    // _________________________________________________________________________________________

    template <class TIn1, class TIn2, class TOut>
    void DifferenceMultiMapping<TIn1, TIn2, TOut>::init()
    {
        Inherit1::init();

        if (this->getFromModels1().empty())
        {
            msg_error() << "Error while initializing ; input getFromModels1 not found";
            return;
        }

        if (this->getFromModels2().empty())
        {
            msg_error() << "Error while initializing ; output getFromModels2 not found";
            return;
        }

        if (this->getToModels().empty())
        {
            msg_error() << "Error while initializing ; output Model not found";
            // return;
        }
        m_fromModel1 = this->getFromModels1()[0];
        m_fromModel2 = this->getFromModels2()[0];
        m_toModel = this->getToModels()[0];

        m_toModel = m_fromModel1;

        initiateTopologies();

        printf("init\n");
    }

    template <class TIn1, class TIn2, class TOut>
    void DifferenceMultiMapping<TIn1, TIn2, TOut>::bwdInit()
    {
    }

    template <class TIn1, class TIn2, class TOut>
    void DifferenceMultiMapping<TIn1, TIn2, TOut>::reinit()
    {
    }

    template <class TIn1, class TIn2, class TOut>
    void DifferenceMultiMapping<TIn1, TIn2, TOut>::reset()
    {
        reinit();
    }

    template <class TIn1, class TIn2, class TOut>
    void DifferenceMultiMapping<TIn1, TIn2, TOut>::computeProximity(const In1VecCoord &x1, const In2VecCoord &x2)
    {

        In1VecCoord from = x1;
        In2VecCoord dst = x2;
        m_constraints.clear();

        size_t szFrom = from.size();
        size_t szDst = dst.size();
        type::vector<Rigid> direction = d_direction.getValue();

        /// get the last rigid direction, the main goal is to use it for the
        ///  3D bilateral constraint i.e the fix point of the cable in the robot structure
        // Rigid direction = d_direction.getValue()[szDst-1];

        // For each point in the FEM find the closest edge of the cable
        for (size_t i = 0; i < szFrom; i++)
        {
            Coord2 P = from[i];
            Constraint constraint;

            // find the min distance between a from mstate point, and it's projection on each edge of the cable (destination mstate)
            Real min_dist = std::numeric_limits<Real>::max();
            for (size_t j = 0; j < szDst - 1; j++)
            {
                Coord1 Q1 = dst[j];
                Coord1 Q2 = dst[j + 1];
                // the axis
                Coord1 dirAxe = Q2 - Q1;
                Real length = dirAxe.norm();
                Real fact_v = dot(P - Q1, dirAxe) / dot(dirAxe, dirAxe);

                if (std::abs(fact_v) < min_dist)
                {
                    // if(fact_v < min_dist){
                    min_dist = std::abs(fact_v);

                    // define the constraint variables
                    Deriv1 proj; // distVec;
                    Real alpha;  // dist;

                    /// To solve the case that the closest node is
                    ///  not the node 0 but the node 1 of the beam
                    if (fact_v < 0.0 && j != 0 && std::abs(fact_v) > 1e-8)
                    {
                        // if fact_v < 0.0 that means the last beam is the good beam
                        // printf("if fact_v < 0.0 that means the last beam is the good beam \n");
                        Q1 = dst[j - 1];
                        dirAxe = dst[j] - Q1;
                        length = dirAxe.norm();
                        fact_v = dot(P - Q1, dirAxe) / dot(dirAxe, dirAxe);
                        dirAxe.normalize();
                        alpha = (P - Q1) * dirAxe;

                        proj = Q1 + dirAxe * alpha;
                        // distVec = P - proj; // violation vector
                        // dist = (P - proj).norm(); // constraint violation
                        constraint.eid = j - 1;
                        // The direction of the axe or the beam
                        constraint.dirAxe = dirAxe;
                        // the node contribution to the constraint which is 1-coeff
                        alpha = alpha / length; // normalize, ensure that <1.0
                        if (alpha < 1e-8)
                            constraint.alpha = 1.0;
                        else
                            constraint.alpha = 1.0 - alpha;

                        // The projection on the axe
                        constraint.proj = proj;
                        constraint.Q = from[i];

                        /////
                        length = (dst[j] - Q1).norm();
                        constraint.Q1Q2 = length;
                        constraint.r2 = fact_v;

                        // We move the constraint point onto the projection
                        Deriv1 t1 = P - proj;        // violation vector
                        constraint.dist = t1.norm(); // constraint violation
                        t1.normalize();              // direction of the constraint

                        //// First method compute normals using projections
                        //                    if(t1.norm()<1.0e-1 && dirAxe[2] < 0.99){
                        //                        type::Vec3 temp = type::Vec3(dirAxe[0],dirAxe[1],dirAxe[2]+50.0);
                        //                        t1 = cross(dirAxe,temp);
                        //                        t1.normalize();
                        //                        constraint.t1 = t1;
                        //                    }
                        //                    if(t1.norm()<1.0e-1){
                        //                        type::Vec3 temp = type::Vec3(dirAxe[0],dirAxe[1]+50.0,dirAxe[2]);
                        //                        t1 = cross(dirAxe,temp);
                        //                        t1.normalize();
                        //                        constraint.t1 = t1;
                        //                    }

                        //                    if(t1.norm()<1.0e-1)
                        //                    {

                        //// Second method compute normals using frames directions
                        Rigid dir = direction[constraint.eid];
                        type::Vec3 vY = type::Vec3(0., 1., 0.);
                        type::Quat ori = dir.getOrientation();
                        vY = ori.rotate(vY);
                        vY.normalize();
                        t1 = vY;
                        //                    }

                        constraint.t1 = t1;
                        // tangential 2
                        Deriv1 t2 = cross(t1, dirAxe); t2.normalize();
                        constraint.t2 = t2;

                        if (i == szFrom - 1)
                        {
                            /// This handle the fix point constraint the last point of
                            ///  of cstr points indeed here we have
                            ///  3D bilateral constraint and alpha=1.0
                            // We use the given direction of fill H

                            if (!direction.empty())
                            {
                                type::Quat _ori = direction[szDst - 1].getOrientation();
                                type::Vec3 _vY = _ori.rotate(type::Vec3(0., 1., 0.)); _vY.normalize();
                                type::Vec3 _vZ = _ori.rotate(type::Vec3(0., 0., 1.)); _vZ.normalize();

                                constraint.t1 = _vY;
                                constraint.t2 = _vZ;
                            }
                            constraint.proj = dst[szDst - 1];
                            constraint.eid = szDst - 2;
                            constraint.alpha = 1.0;
                            constraint.dist = (dst[szDst - 1] - from[szFrom - 1]).norm();
                        }
                    }
                    else
                    {
                        // compute needs for constraint
                        dirAxe.normalize();
                        alpha = (P - Q1) * dirAxe;

                        proj = Q1 + dirAxe * alpha;
                        // distVec = P - proj; // violation vector
                        // dist = (P - proj).norm(); // constraint violation
                        constraint.eid = j;
                        // The direction of the axe or the beam
                        constraint.dirAxe = dirAxe;
                        // the node contribution to the constraint which is 1-coeff
                        alpha = alpha / length; // normalize, ensure that <1.0
                        if (alpha < 1e-8)
                            constraint.alpha = 1.0;
                        else
                            constraint.alpha = 1.0 - alpha;

                        // The projection on the axe
                        constraint.proj = proj;
                        constraint.Q = from[i];

                        /////
                        constraint.Q1Q2 = length;
                        constraint.r2 = fact_v;

                        // We move the constraint point onto the projection
                        Deriv1 t1 = P - proj;        // violation vector
                        constraint.dist = t1.norm(); // constraint violation
                        t1.normalize();              // direction of the constraint

                        /// If the violation is very small t1 is close to zero
                        ///
                        //// First method compute normals using projections
                        //                    if(t1.norm()<1.0e-1 && dirAxe[2] < 0.99){
                        //                        type::Vec3 temp = type::Vec3(dirAxe[0],dirAxe[1],dirAxe[2]+50.0);
                        //                        t1 = cross(dirAxe,temp);
                        //                        t1.normalize();
                        //                        constraint.t1 = t1;
                        //                    }
                        //                    if(t1.norm()<1.0e-1){
                        //                        type::Vec3 temp = type::Vec3(dirAxe[0],dirAxe[1]+50.0,dirAxe[2]);
                        //                        t1 = cross(dirAxe,temp);
                        //                        t1.normalize();
                        //                        constraint.t1 = t1;
                        //                    }

                        //// Second method compute normals using frames directions
                        type::Vec3 vY = type::Vec3(0., 1., 0.);
                        type::Quat ori = (direction[szDst - 1]).getOrientation();
                        vY = ori.rotate(vY); vY.normalize();
                        t1 = vY;
                        //                    }
                        constraint.t1 = t1;
                        // tangential 2
                        Deriv1 t2 = cross(t1, dirAxe);
                        t2.normalize();
                        constraint.t2 = t2;

                        /// This is need because we are applying the a
                        ///  bilateral constraint on the last node of the mstate
                        if (i == szFrom - 1)
                        {
                            /// This handle the fix point constraint the last point of
                            ///  of cstr points indeed here we have
                            ///  3D bilateral constraint and alpha=1.0
                            // We use the given direction of fill H
                            if (!d_direction.getValue().empty())
                            {
                                type::Quat _ori = (direction[szDst - 1]).getOrientation();
                                type::Vec3 _vY = _ori.rotate(type::Vec3(0., 1., 0.)); _vY.normalize();
                                type::Vec3 _vZ = _ori.rotate(type::Vec3(0., 0., 1.)); _vZ.normalize();

                                constraint.t1 = _vY;
                                constraint.t2 = _vZ;
                            }
                            constraint.proj = dst[szDst - 1];
                            constraint.eid = szDst - 2;
                            constraint.alpha = 1.0;
                            constraint.dist = (dst[szDst - 1] - from[szFrom - 1]).norm();
                        }
                    }
                }
            }            
            m_constraints.push_back(constraint);
        }
    }

    template <class TIn1, class TIn2, class TOut>
    void DifferenceMultiMapping<TIn1, TIn2, TOut>::computeNeedleProximity(const In1VecCoord &x1, const In2VecCoord &x2)
    {

        In1VecCoord from = x1;
        In2VecCoord dst = x2;
        m_constraints.clear();

        size_t szFrom = from.size();
        size_t szDst = dst.size();
        type::vector<Rigid> direction = d_direction.getValue();

        /// get the last rigid direction, the main goal is to use it for the
        ///  3D bilateral constraint i.e the fix point of the cable in the robot structure
        // Rigid direction = d_direction.getValue()[szDst-1];

        // For each point in the FEM find the closest edge of the cable
        for (size_t i = 0; i < szFrom; i++)
        {
            Coord2 P = from[i];
            Constraint constraint;

            // find the min distance between a from mstate point, and it's projection on each edge of the cable (destination mstate)
            Real min_dist = std::numeric_limits<Real>::max();
            for (size_t j = 0; j < szDst - 1; j++)
            {
                Coord1 Q1 = dst[j];
                Coord1 Q2 = dst[j + 1];
                // the axis
                Coord1 dirAxe = Q2 - Q1;
                Real length = dirAxe.norm();
                Real fact_v = dot(P - Q1, dirAxe) / dot(dirAxe, dirAxe);

                if (std::abs(fact_v) < min_dist)
                {
                    // if(fact_v < min_dist){
                    min_dist = std::abs(fact_v);

                    // define the constraint variables
                    Deriv1 proj;
                    Real alpha;

                    /// To solve the case that the closest node is
                    ///  not the node 0 but the node 1 of the beam
                    if (fact_v < 0.0 && j != 0 && std::abs(fact_v) > 1e-8)
                    {
                        // if fact_v < 0.0 that means the last beam is the good beam
                        // printf("if fact_v < 0.0 that means the last beam is the good beam \n");
                        Q1 = dst[j - 1];
                        dirAxe = dst[j] - Q1;
                        length = dirAxe.norm();
                        fact_v = dot(P - Q1, dirAxe) / dot(dirAxe, dirAxe);
                        dirAxe.normalize();
                        alpha = (P - Q1) * dirAxe;

                        proj = Q1 + dirAxe * alpha;
                        // distVec = P - proj; // violation vector
                        // dist = (P - proj).norm(); // constraint violation
                        constraint.eid = j - 1;
                        // The direction of the axe or the beam
                        constraint.dirAxe = dirAxe;
                        // the node contribution to the constraint which is 1-coeff
                        alpha = alpha / length; // normalize, ensure that <1.0
                        if (alpha < 1e-8)
                            constraint.alpha = 1.0;
                        else
                            constraint.alpha = 1.0 - alpha;

                        // The projection on the axe
                        constraint.proj = proj;
                        constraint.Q = from[i];

                        /////
                        length = (dst[j] - Q1).norm();
                        constraint.Q1Q2 = length;
                        constraint.r2 = fact_v;

                        // We move the constraint point onto the projection
                        Deriv1 t1 = P - proj;        // violation vector
                        constraint.dist = t1.norm(); // constraint violation
                        t1.normalize();              // direction of the constraint

                        //// Second method compute normals using frames directions
                        Rigid dir = direction[constraint.eid];
                        type::Vec3 vY = type::Vec3(0., 1., 0.);
                        type::Quat ori = dir.getOrientation();
                        vY = ori.rotate(vY); vY.normalize();
                        t1 = vY;
                        //                    }

                        constraint.t1 = t1;
                        // tangential 2
                        Deriv1 t2 = cross(t1, dirAxe);
                        t2.normalize();
                        constraint.t2 = t2;
                    }
                    else
                    {
                        // compute needs for constraint
                        dirAxe.normalize();
                        alpha = (P - Q1) * dirAxe;

                        proj = Q1 + dirAxe * alpha;
                        // distVec = P - proj; // violation vector
                        // dist = (P - proj).norm(); // constraint violation
                        constraint.eid = j;
                        // The direction of the axe or the beam
                        constraint.dirAxe = dirAxe;
                        // the node contribution to the constraint which is 1-coeff
                        alpha = alpha / length; // normalize, ensure that <1.0
                        if (alpha < 1e-8)
                            constraint.alpha = 1.0;
                        else
                            constraint.alpha = 1.0 - alpha;

                        // The projection on the axe
                        constraint.proj = proj;
                        constraint.Q = from[i];

                        /////
                        constraint.Q1Q2 = length;
                        constraint.r2 = fact_v;

                        // We move the constraint point onto the projection
                        Deriv1 t1 = P - proj;        // violation vector
                        constraint.dist = t1.norm(); // constraint violation
                        t1.normalize();              // direction of the constraint

                        //// Second method compute normals using frames directions
                        Rigid dir = direction[constraint.eid];
                        type::Vec3 vY = type::Vec3(0., 1., 0.);
                        type::Quat ori = dir.getOrientation();
                        vY = ori.rotate(vY); vY.normalize();
                        t1 = vY;
                        //                    }
                        constraint.t1 = t1;
                        // tangential 2
                        Deriv1 t2 = cross(t1, dirAxe); t2.normalize();
                        constraint.t2 = t2;
                    }
                }
            }
            m_constraints.push_back(constraint);
        }
    }

    template <class TIn1, class TIn2, class TOut>
    void DifferenceMultiMapping<TIn1, TIn2, TOut>::apply(
        const core::MechanicalParams * /* mparams */, const type::vector<OutDataVecCoord *> &dataVecOutPos,
        const type::vector<const In1DataVecCoord *> &dataVecIn1Pos,
        const type::vector<const In2DataVecCoord *> &dataVecIn2Pos)
    {

        if (dataVecOutPos.empty() || dataVecIn1Pos.empty() || dataVecIn2Pos.empty())
            return;

        // printf("///Do Apply//We need only one input In model and input Root model (if present) \n");
        const In1VecCoord &in1 = dataVecIn1Pos[0]->getValue();
        const In2VecCoord &in2 = dataVecIn2Pos[0]->getValue();

        OutVecCoord &out = *dataVecOutPos[0]->beginEdit();

        if (d_lastPointIsFixed.getValue())
        {
            computeProximity(in1, in2);

            // auto out = sofa::helper::writeOnly(*dataVecOutPos[0]);
            size_t sz = m_constraints.size();
            out.resize(sz);

            for (unsigned int i = 0; i < sz; i++)
            {
                Constraint &c = m_constraints[i];
                if (i < sz - 1)
                {
                    out[i][0] = 0.0;
                    out[i][1] = c.t1 * (in1[i] - c.proj); // c.dist;
                    out[i][2] = c.t2 * (in1[i] - c.proj); // 0.0
                }
                else
                {
                    out[sz - 1][0] = c.dirAxe * (in1[in1.size() - 1] - in2[in2.size() - 1]);
                    out[sz - 1][1] = c.t1 * (in1[in1.size() - 1] - in2[in2.size() - 1]); // std::abs(in2[in2.size()-1][1] - in1[in1.size()-1][1]);
                    out[sz - 1][2] = c.t2 * (in1[in1.size() - 1] - in2[in2.size() - 1]); // std::abs(in2[in2.size()-1][2] - in1[in1.size()-1][2]);
                }
            }
        }
        else
        {
            computeNeedleProximity(in1, in2);

            // auto out = sofa::helper::writeOnly(*dataVecOutPos[0]);
            size_t sz = m_constraints.size();
            out.resize(sz);

            for (unsigned int i = 0; i < sz; i++)
            {
                Constraint &c = m_constraints[i];
                out[i][0] = 0.0;
                out[i][1] = c.t1 * (in1[i] - c.proj); // c.dist;
                out[i][2] = c.t2 * (in1[i] - c.proj); // 0.0
            }
        }
        dataVecOutPos[0]->endEdit();
    }

    template <class TIn1, class TIn2, class TOut>
    void DifferenceMultiMapping<TIn1, TIn2, TOut>::applyJ(
        const core::MechanicalParams * /* mparams */, const type::vector<OutDataVecDeriv *> &dataVecOutVel,
        const type::vector<const In1DataVecDeriv *> &dataVecIn1Vel,
        const type::vector<const In2DataVecDeriv *> &dataVecIn2Vel)
    {
        if (dataVecOutVel.empty() || dataVecIn1Vel.empty() || dataVecIn2Vel.empty())
            return;
        const In1VecDeriv &in1 = dataVecIn1Vel[0]->getValue();
        const In2VecDeriv &in2 = dataVecIn2Vel[0]->getValue();
        OutVecDeriv &outVel = *dataVecOutVel[0]->beginEdit();

        // const OutVecCoord& out = m_toModel->read(core::ConstVecCoordId::position())->getValue();
        size_t sz = m_constraints.size();
        outVel.resize(sz);
        
        if (d_lastPointIsFixed.getValue())
        {
            for (size_t i = 0; i < sz; i++)
            {
                Constraint &c = m_constraints[i];
                int ei1 = c.eid;
                int ei2 = c.eid + 1;
                if (i < sz - 1)
                {
                    Real v0 = c.dirAxe * (in1[i] - c.alpha * in2[ei1] - (1 - c.alpha) * in2[ei2]);
                    Real v1 = c.t1 * (in1[i] - c.alpha * in2[ei1] - (1 - c.alpha) * in2[ei2]);
                    Real v2 = c.t2 * (in1[i] - c.alpha * in2[ei1] - (1 - c.alpha) * in2[ei2]);
                    outVel[i] = OutDeriv(v0, v1, v2);
                }
                else
                {
                    Real v0 = c.dirAxe * (in1[i] - in2[ei2]);
                    Real v1 = c.t1 * (in1[i] - in2[ei2]);
                    Real v2 = c.t2 * (in1[i] - in2[ei2]);
                    outVel[i] = OutDeriv(v0, v1, v2);
                }
            }
        }
        else
        {
            for (size_t i = 0; i < sz; i++)
            {
                Constraint &c = m_constraints[i];

                int ei1 = c.eid;
                int ei2 = c.eid + 1;
                Real v0 = c.dirAxe * (in1[i] - c.alpha * in2[ei1] - (1 - c.alpha) * in2[ei2]);
                Real v1 = c.t1 * (in1[i] - c.alpha * in2[ei1] - (1 - c.alpha) * in2[ei2]);
                Real v2 = c.t2 * (in1[i] - c.alpha * in2[ei1] - (1 - c.alpha) * in2[ei2]);

                outVel[i] = OutDeriv(v0, v1, v2);
            }
        }
        dataVecOutVel[0]->endEdit();
    }

    template <class TIn1, class TIn2, class TOut>
    void DifferenceMultiMapping<TIn1, TIn2, TOut>::applyJT(
        const core::MechanicalParams * /*mparams*/, const type::vector<In1DataVecDeriv *> &dataVecOut1Force,
        const type::vector<In2DataVecDeriv *> &dataVecOut2Force,
        const type::vector<const OutDataVecDeriv *> &dataVecInForce)
    {
        if (dataVecOut1Force.empty() || dataVecInForce.empty() || dataVecOut2Force.empty())
            return;
        
        const OutVecDeriv &in = dataVecInForce[0]->getValue();
        if (in.empty())
            return;

        In1VecDeriv &out1 = *dataVecOut1Force[0]->beginEdit();
        In2VecDeriv &out2 = *dataVecOut2Force[0]->beginEdit();
        // Compute output forces
        size_t sz = m_constraints.size();
        if (d_lastPointIsFixed.getValue())
        {
            for (size_t i = 0; i < sz; i++)
            {
                Constraint &c = m_constraints[i];
                int ei1 = c.eid;
                int ei2 = c.eid + 1;
                OutDeriv f = in[i];
                if (i < sz - 1)
                {
                    Deriv2 f1 = (f[0] * c.dirAxe) + (f[1] * c.t1) + (f[2] * c.t2);
                    Deriv1 f2_1 = (c.alpha * f[0] * c.dirAxe) + (c.alpha * f[1] * c.t1) + (c.alpha * f[2] * c.t2);
                    Deriv1 f2_2 = ((1 - c.alpha) * f[0] * c.dirAxe) + ((1 - c.alpha) * f[1] * c.t1) + ((1 - c.alpha) * f[2] * c.t2);
                    out1[i] += f1; out2[ei1] -= f2_1; out2[ei2] -= f2_2;
                }
                else
                {
                    Deriv2 f1 = (f[0] * c.dirAxe) + (f[1] * c.t1) + (f[2] * c.t2);
                    out1[i] += f1; out2[ei2] -= f1;
                }
            }
        }
        else
        {
            for (size_t i = 0; i < sz; i++)
            {
                Constraint &c = m_constraints[i];
                int ei1 = c.eid;
                int ei2 = c.eid + 1;
                OutDeriv f = in[i];
                Deriv2 f1 = (f[0] * c.dirAxe) + (f[1] * c.t1) + (f[2] * c.t2);
                Deriv1 f2_1 = (c.alpha * f[0] * c.dirAxe) + (c.alpha * f[1] * c.t1) + (c.alpha * f[2] * c.t2);
                Deriv1 f2_2 = ((1 - c.alpha) * f[0] * c.dirAxe) + ((1 - c.alpha) * f[1] * c.t1) + ((1 - c.alpha) * f[2] * c.t2);
                out1[i] += f1;
                out2[ei1] -= f2_1;
                out2[ei2] -= f2_2;
            }
        }
        dataVecOut1Force[0]->endEdit();
        dataVecOut2Force[0]->endEdit();
    }

    //___________________________________________________________________________
    template <class TIn1, class TIn2, class TOut>
    void DifferenceMultiMapping<TIn1, TIn2, TOut>::applyJT(
        const core::ConstraintParams * /*cparams*/, const type::vector<In1DataMatrixDeriv *> &dataMatOut1Const,
        const type::vector<In2DataMatrixDeriv *> &dataMatOut2Const,
        const type::vector<const OutDataMatrixDeriv *> &dataMatInConst)
    {
        if (dataMatOut1Const.empty() || dataMatOut2Const.empty() || dataMatInConst.empty())
            return;

        // We need only one input In model and input Root model (if present)
        In1MatrixDeriv &out1 = *dataMatOut1Const[0]->beginEdit(); // constraints on the FEM cable points
        In2MatrixDeriv &out2 = *dataMatOut2Const[0]->beginEdit(); // constraints on the frames cable points
        const OutMatrixDeriv &in = dataMatInConst[0]->getValue(); // input constraints defined on the mapped point
        const In1DataVecCoord *x1fromData = m_fromModel1->read(core::ConstVecCoordId::position());
        const In1VecCoord x1from = x1fromData->getValue();

        typename OutMatrixDeriv::RowConstIterator rowIt = in.begin();
        typename OutMatrixDeriv::RowConstIterator rowItEnd = in.end();

        for (rowIt = in.begin(); rowIt != rowItEnd; ++rowIt)
        {
            typename OutMatrixDeriv::ColConstIterator colIt = rowIt.begin();
            typename OutMatrixDeriv::ColConstIterator colItEnd = rowIt.end();

            // Creates a constraints if the input constraint is not empty.
            if (colIt == colItEnd)
            {
                continue;
            }
            typename In1MatrixDeriv::RowIterator o1 = out1.writeLine(rowIt.index()); // we store the constraint number
            typename In2MatrixDeriv::RowIterator o2 = out2.writeLine(rowIt.index());

            if (d_lastPointIsFixed.getValue())
            {
                if ((rowIt.index() / 2) < (int)(x1from.size() - 1))
                {
                    while (colIt != colItEnd)
                    {
                        int childIndex = colIt.index();
                        Constraint c = m_constraints[childIndex];
                        const OutDeriv h = colIt.val();
                        int indexBeam = c.eid;

                        Deriv2 h1 = (h[0] * c.dirAxe) + (h[1] * c.t1) + (h[2] * c.t2);
                        Deriv1 h2_1 = (c.alpha * h[0] * c.dirAxe) + (c.alpha * h[1] * c.t1) + (c.alpha * h[2] * c.t2);
                        Deriv1 h2_2 = ((1.0 - c.alpha) * h[0] * c.dirAxe) + ((1.0 - c.alpha) * h[1] * c.t1) + ((1.0 - c.alpha) * h[2] * c.t2);

                        o1.addCol(childIndex, h1);
                        o2.addCol(indexBeam, -h2_1);
                        o2.addCol(indexBeam + 1, -h2_2);

                        colIt++;
                    }
                }
                else
                {
                    while (colIt != colItEnd)
                    {
                        int childIndex = colIt.index();
                        Constraint c = m_constraints[childIndex];
                        const OutDeriv h = colIt.val();
                        int indexBeam = c.eid;

                        Deriv2 h1 = (h[0] * c.dirAxe) + (h[1] * c.t1) + (h[2] * c.t2);
                        Deriv1 h2 = (h[0] * c.dirAxe) + (h[1] * c.t1) + (h[2] * c.t2);

                        o1.addCol(childIndex, h1);
                        o2.addCol(indexBeam + 1, -h2);
                        colIt++;
                    }
                }
            }
            else
            {
                while (colIt != colItEnd)
                {
                    int childIndex = colIt.index();
                    Constraint c = m_constraints[childIndex];
                    const OutDeriv h = colIt.val();
                    int indexBeam = c.eid;

                    Deriv2 h1 = (h[0] * c.dirAxe) + (h[1] * c.t1) + (h[2] * c.t2);
                    Deriv1 h2_1 = (c.alpha * h[0] * c.dirAxe) + (c.alpha * h[1] * c.t1) + (c.alpha * h[2] * c.t2);
                    Deriv1 h2_2 = ((1.0 - c.alpha) * h[0] * c.dirAxe) + ((1.0 - c.alpha) * h[1] * c.t1) + ((1.0 - c.alpha) * h[2] * c.t2);

                    o1.addCol(childIndex, h1);
                    o2.addCol(indexBeam, -h2_1);
                    o2.addCol(indexBeam + 1, -h2_2);

                    colIt++;
                }
            }
        }
        dataMatOut1Const[0]->endEdit();
        dataMatOut2Const[0]->endEdit();
    }

    template <class TIn1, class TIn2, class TOut>
    void DifferenceMultiMapping<TIn1, TIn2, TOut>::draw(const core::visual::VisualParams *vparams)
    {
        /// draw cable
        if (!vparams->displayFlags().getShowInteractionForceFields())
            return;

        typedef sofa::type::RGBAColor RGBAColor;
        vparams->drawTool()->saveLastState();
        vparams->drawTool()->disableLighting();

        std::vector<type::Vec3> vertices;
        RGBAColor color = RGBAColor::magenta();

        if (d_drawArrows.getValue() && d_lastPointIsFixed.getValue())
        {
            for (size_t i = 0; i < m_constraints.size(); i++)
            {
                color = RGBAColor::green();
                vertices.push_back(m_constraints[i].proj);
                vertices.push_back(m_constraints[i].Q);
                vparams->drawTool()->drawLines(vertices, 4.0, color);
                if (i == (m_constraints.size() - 1))
                {
                    Coord2 P1 = m_constraints[i].Q;
                    Real radius_arrow = 0.30;
                    Coord2 x = m_constraints[i].dirAxe * 5.0;
                    Coord2 y = m_constraints[i].t1 * 5.0;
                    Coord2 z = m_constraints[i].t2 * 5.0;

                    vparams->drawTool()->drawArrow(P1, P1 + x, radius_arrow, RGBAColor::red());
                    vparams->drawTool()->drawArrow(P1, P1 + y, radius_arrow, RGBAColor::green());
                    vparams->drawTool()->drawArrow(P1, P1 + z, radius_arrow, RGBAColor::blue());
                }
                else
                {
                    Coord2 P1 = m_constraints[i].Q;
                    Real radius_arrow = 0.30;
                    // Coord2 x = m_constraints[i].dirAxe * 5.0;
                    Coord2 y = m_constraints[i].t1 * 5.0;
                    Coord2 z = m_constraints[i].t2 * 5.0;
                    vparams->drawTool()->drawArrow(P1, P1 + y, radius_arrow, RGBAColor::blue());
                    vparams->drawTool()->drawArrow(P1, P1 + z, radius_arrow, RGBAColor::blue());
                }
            }
            const In1DataVecDeriv *xDestData = m_fromModel1->read(core::ConstVecCoordId::position());
            const In1VecCoord &fromPos = xDestData[0].getValue();
            vparams->drawTool()->draw3DText_Indices(fromPos, 6, RGBAColor(0.0, 1.0, 0.0, 1.0));
        }
        vparams->drawTool()->restoreLastState();
    }

} // namespace sofa
