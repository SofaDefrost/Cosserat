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
#include "CosseratSlidingConstraint.h"

#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/behavior/BaseConstraint.h>
#include <sofa/type/RGBAColor.h>
#include <sofa/type/Vec.h>

namespace sofa::component::constraintset
{

template<class DataTypes>
CosseratSlidingConstraint<DataTypes>::CosseratSlidingConstraint()
    : CosseratSlidingConstraint(nullptr, nullptr)
{
}

template<class DataTypes>
CosseratSlidingConstraint<DataTypes>::CosseratSlidingConstraint(MechanicalState* m_from)
    : CosseratSlidingConstraint(m_from, m_from)
{
}

template<class DataTypes>
CosseratSlidingConstraint<DataTypes>::CosseratSlidingConstraint(MechanicalState* m_from, MechanicalState* m_dst)
    : Inherit(m_from, m_dst)
    , d_m1(initData(&d_m1, 0, "sliding_point","index of the spliding point on the first model"))
    , d_m2a(initData(&d_m2a, 0, "axis_1","index of one end of the sliding axis"))
    , d_m2b(initData(&d_m2b, 0, "axis_2","index of the other end of the sliding axis"))
    , d_force(initData(&d_force,"force","force (impulse) used to solve the constraint"))
    //, m_yetIntegrated(false)
    //, m_step(0)
{
}

template<class DataTypes>
void CosseratSlidingConstraint<DataTypes>::init()
{
    assert(this->mstate1);
    assert(this->mstate2);
}


template<class DataTypes>
void CosseratSlidingConstraint<DataTypes>::computeProximity(const DataVecCoord &x1, const DataVecCoord &x2){

    VecCoord from = x1.getValue();
    VecCoord dst  = x2.getValue();
    m_constraints.clear();

    size_t szFrom = from.size();
    size_t szDst = dst.size();
    //For each point in the FEM find the closest edge of the cable
    for (unsigned int i = 0 ; i < szFrom-1; i++) {
        Coord P = from[i];
        Constraint constraint;

        //std::cout << "P "<< P << std::endl;

        // min dist between the projection and the projected point
        Real min_dist = std::numeric_limits<Real>::max();
        for (unsigned int j = 0; j < szDst-1; j++) {
            Coord Q1 = dst[j];
            Coord Q2 = dst[j+1];

            // the axis
            Coord v = Q2 -Q1;
            Real fact_v = dot(P-Q1,v) / dot(v,v);

            //std::cout<< " i: " << i<< " \t fact_v :"<< fact_v << std::endl;
            if(fact_v <= 0.0) continue;

            if(fact_v < min_dist){

                Deriv dirAxe = v;
                dirAxe.normalize();

                // projection of the point on the axis
                Real r = (P-Q1) * dirAxe;
                Deriv proj = Q1 + dirAxe * r;

                min_dist = fact_v;
                constraint.P = proj;
                constraint.Q = from[i];
                constraint.eid = j;
                constraint.r = r;
                constraint.dirAxe = dirAxe;
                /////
                constraint.Q1Q2 = v.norm();
                constraint.r2 = fact_v;
                //std::cout << "fact_v :"<< fact_v << " ; r2 :"<< constraint.r2 << std::endl;

                // We move the constraint point onto the projection
                Deriv dirProj = P - proj; // violation vector
                constraint.dist = dirProj.norm(); // constraint violation
                dirProj.normalize(); // direction of the constraint
                constraint.dirProj = dirProj;

                Deriv dirOrtho = cross(dirProj, dirAxe);
                dirOrtho.normalize();
                constraint.dirOrtho = dirOrtho;

            }
        }
        //std::cout << i << " Closets edge is "<< constraint.eid << std::endl;
        m_constraints.push_back(constraint);
    }
    //printf("CosseratSlidingConstraint<DataTypes>::computeProximity After\n");
}

template<class DataTypes>
void CosseratSlidingConstraint<DataTypes>::buildConstraintMatrix(const core::ConstraintParams*, DataMatrixDeriv &c1_d, DataMatrixDeriv &c2_d, unsigned int &cIndex
                                                                 , const DataVecCoord &x1, const DataVecCoord &x2)
{
    computeProximity(x1,x2);
    //printf("=================================\n");

    MatrixDeriv &c1 = *c1_d.beginEdit();
    MatrixDeriv &c2 = *c2_d.beginEdit();
    VecCoord dst  = x2.getValue();

    for (size_t i = 0 ; i <  m_constraints.size(); i++) {
        Constraint& c = m_constraints[i];

        int ei1 = c.eid;
        int ei2 = (ei1==0) ? c.eid+1 : c.eid+1;

        unsigned int cid = cIndex;
        cIndex += 2;

        MatrixDerivRowIterator c1_it = c1.writeLine(cid);
        c1_it.addCol(i, c.dirProj);

        c1_it = c1.writeLine(cid + 1);
        c1_it.setCol(i, c.dirOrtho);

        MatrixDerivRowIterator c2_it = c2.writeLine(cid);
        c2_it.addCol(ei1, -c.dirProj * (1-c.r2));
        c2_it.addCol(ei2, -c.dirProj * c.r2);


        c2_it = c2.writeLine(cid + 1);
        c2_it.addCol(ei1, -c.dirOrtho * (1-c.r2));
        c2_it.addCol(ei2, -c.dirOrtho * c.r2);

        c.thirdConstraint = 0;

        if (c.r < 0)
        {
            c.thirdConstraint = c.r;
            cIndex++;

            c1_it = c1.writeLine(cid + 2);
            c1_it.setCol(i, c.dirAxe);

            c2_it = c2.writeLine(cid + 2);
            c2_it.addCol(ei1, -c.dirAxe);
        }
        else if (c.r > c.Q1Q2)
        {
            c.thirdConstraint = c.r - c.Q1Q2;
            cIndex++;

            c1_it = c1.writeLine(cid + 2);
            c1_it.setCol(i, -c.dirAxe);

            c2_it = c2.writeLine(cid + 2);
            c2_it.addCol(ei2, c.dirAxe);
        }
        c.cid = cid;
    }

    c1_d.endEdit();
    c2_d.endEdit();

}


template<class DataTypes>
void CosseratSlidingConstraint<DataTypes>::getConstraintViolation(const core::ConstraintParams *, sofa::linearalgebra::BaseVector *v, const DataVecCoord &, const DataVecCoord &
                                                                  , const DataVecDeriv &, const DataVecDeriv &)
{
    for (size_t i = 0; i < m_constraints.size(); i++) {
        Constraint& c = m_constraints[i];
        //std::cout << " c.cid :"<< c.cid << " c.dist :"<< c.dist << std::endl;

        v->set(c.cid, c.dist);
        v->set(c.cid+1, 0.0);

        if(c.thirdConstraint)
        {
            if(c.thirdConstraint>0)
                v->set(c.cid+2, - c.thirdConstraint);
            else
                v->set(c.cid+2, c.thirdConstraint);
        }
    }
}

using sofa::linearalgebra::BaseVector ;
template<class DataTypes>
void CosseratSlidingConstraint<DataTypes>::getConstraintResolution(const ConstraintParams*,
                                                                   std::vector<core::behavior::ConstraintResolution*>& resTab,
                                                                   unsigned int& offset)
{
    for (size_t i = 0; i < m_constraints.size(); i++) {
        Constraint& c = m_constraints[i];
        resTab[offset++] = new BilateralConstraintResolution();
        resTab[offset++] = new BilateralConstraintResolution();

        if(c.thirdConstraint)
            resTab[offset++] = new UnilateralConstraintResolution();
    }
}


template<class DataTypes>
void CosseratSlidingConstraint<DataTypes>::storeLambda(const ConstraintParams* /*cParams*/, sofa::core::MultiVecDerivId /*res*/, const BaseVector* lambda)
{
    Real lamb1,lamb2, lamb3;
    for (size_t i = 0; i < m_constraints.size(); i++) {
        Constraint& c = m_constraints[i];

        lamb1 = lambda->element(c.cid);
        lamb2 = lambda->element(c.cid+1);

        if(c.thirdConstraint)
        {
            lamb3 = lambda->element(c.cid+2);
            d_force.setValue( c.dirProj* lamb1 + c.dirOrtho * lamb2 + c.dirProj * lamb3);
        }
        else
        {
            d_force.setValue( c.dirProj* lamb1 + c.dirOrtho * lamb2 );
        }
    }
}

template<class DataTypes>
void CosseratSlidingConstraint<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
    //printf("CosseratSlidingConstraint<DataTypes>::draw(const core::visual::VisualParams* vparams) before \n");
    if (!vparams->displayFlags().getShowInteractionForceFields())
        return;

    vparams->drawTool()->saveLastState();

    vparams->drawTool()->disableLighting();

    sofa::type::RGBAColor color;
    //    Constraint& c = m_constraints[0];

    //    if(c.thirdConstraint<0)
    //        color = sofa::defaulttype::RGBAColor::yellow();
    //    else if(c.thirdConstraint>0)
    //        color = sofa::defaulttype::RGBAColor::green();
    //    else
    color = sofa::type::RGBAColor::magenta();

    std::vector<sofa::type::Vec3> vertices;
    //    vertices.push_back(DataTypes::getCPos((this->mstate1->read(core::ConstVecCoordId::position())->getValue())[d_m1.getValue()]));

    //    vparams->drawTool()->drawPoints(vertices, 10, color);
    //    vertices.clear();

    //    color = sofa::defaulttype::RGBAColor::blue();
    //    vertices.push_back(DataTypes::getCPos((this->mstate2->read(core::ConstVecCoordId::position())->getValue())[d_m2a.getValue()]));
    //    vertices.push_back(DataTypes::getCPos((this->mstate2->read(core::ConstVecCoordId::position())->getValue())[d_m2b.getValue()]));
    vparams->drawTool()->drawLines(vertices, 1, color);

    for (size_t i =0 ; i < m_constraints.size(); i++) {
        color = sofa::type::RGBAColor::green();
        vertices.push_back(m_constraints[i].P);
        vertices.push_back(m_constraints[i].Q);
        vparams->drawTool()->drawLines(vertices, 1, color);
    }
    //printf("CosseratSlidingConstraint<DataTypes>::draw(const core::visual::VisualParams* vparams) After \n");

    drawLinesBetweenPoints(vparams);
    vparams->drawTool()->restoreLastState();
}

template<class DataTypes>
void CosseratSlidingConstraint<DataTypes>::drawLinesBetweenPoints(const core::visual::VisualParams* vparams)
{
    const VecCoord & positions  = this->mstate2->read(core::ConstVecCoordId::position())->getValue();
    sofa::type::RGBAColor color;
    color = sofa::type::RGBAColor::magenta();
    std::vector<sofa::type::Vec3> vertices;
    for (unsigned int i=0; i<positions.size()-1; i++)
    {
        vertices.push_back(positions[i]);
        vertices.push_back(positions[i+1]);
    }

    vparams->drawTool()->drawLines(vertices, 1.5, color);
}

} // namespace sofa
