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
#ifndef COSSERAT_ProjectionEngine_INL
#define COSSERAT_ProjectionEngine_INL

#include "ProjectionEngine.h"

#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/behavior/BaseConstraint.h>
#include <sofa/type/RGBAColor.h>
#include <sofa/type/Vec.h>


namespace sofa::component::constraintset
{


template<class DataTypes>
ProjectionEngine<DataTypes>::ProjectionEngine()
    : d_from ( initData (&d_from, "fromPos", "The position of the mstate we are deforming, here the points mapped inside the FEM ") )
    , d_dest ( initData (&d_dest, "destination", "The position of the cable points") )
    , d_output( initData (&d_output, "output", "output information ") )
{
    f_listening.setValue(true);
}

template<class DataTypes>
void ProjectionEngine<DataTypes>::init()
{
    addInput(&d_from);
    addInput(&d_dest);
    addOutput(&d_output);

    Inherit1::init();
}

template<class DataTypes>
void ProjectionEngine<DataTypes>::reinit()
{

}


template<class DataTypes>
void ProjectionEngine<DataTypes>::computeProximity(){

    VecCoord from = d_from.getValue();
    VecCoord dst  = d_dest.getValue();
    m_constraints.clear();

    size_t szFrom = from.size();
    size_t szDst = dst.size();
    //For each point in the FEM find the closest edge of the cable
    for (size_t i = 0 ; i < szFrom-1; i++) {
        Coord P = from[i];
        Constraint constraint;

        //std::cout << "P "<< P << std::endl;

        // find the min distance between a from mstate point and it's projection on each edge of the cable (destination mstate)
        Real min_dist = std::numeric_limits<Real>::max();
        for (size_t j = 0; j < szDst-1; j++) {
            Coord Q1 = dst[j];
            Coord Q2 = dst[j+1];

            // the axis
            Coord v = Q2 -Q1;
            Real fact_v = dot(P-Q1,v) / dot(v,v);

            // std::cout<< " i: " << i<< " ; j" << j <<  " \t fact_v :"<< fact_v << std::endl;
            if(fact_v <= 0.0) continue;

            if(fact_v < min_dist){

                Deriv dirAxe = v;
                dirAxe.normalize();

                // projection of the point on the axis
                // std::cout << " =========>>> (P-Q1) * dirAxe : " << (P-Q1) * dirAxe << std::endl;
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
                // std::cout << "fact_v :"<< fact_v << " ; r2 :"<< constraint.r2 << std::endl;

                // We move the constraint point onto the projection
                Deriv t1 = P - proj; // violation vector
                constraint.dist = t1.norm(); // constraint violation
                t1.normalize(); // direction of the constraint
                constraint.t1 = t1;

                //tangential 2
                Deriv t2 = cross(t1, dirAxe);  t2.normalize();
                constraint.t2 = t2;

            }
        }
        //std::cout << i << " Closets edge is "<< constraint.eid << std::endl;
        m_constraints.push_back(constraint);
    }
    //printf("ProjectionEngine<DataTypes>::computeProximity After\n");
}

template <class DataTypes>
void ProjectionEngine<DataTypes>::handleEvent(core::objectmodel::Event* event)
{
    // *****************************
    // Update bending at benginEvent
    if (dynamic_cast<sofa::simulation::AnimateBeginEvent *>(event))
        computeProximity();
    // printf("==================================> \n");
}

//template<class DataTypes>
//void ProjectionEngine<DataTypes>::buildConstraintMatrix(const core::ConstraintParams*, DataMatrixDeriv &c1_d, DataMatrixDeriv &c2_d, unsigned int &cIndex
//                                                                 , const DataVecCoord &x1, const DataVecCoord &x2)
//{
//    computeProximity(x1,x2);
//    //printf("=================================\n");

//    MatrixDeriv &c1 = *c1_d.beginEdit();
//    MatrixDeriv &c2 = *c2_d.beginEdit();
//    VecCoord dst  = x2.getValue();

//    for (size_t i = 0 ; i <  m_constraints.size(); i++) {
//        Constraint& c = m_constraints[i];

//        int ei1 = c.eid;
//        int ei2 = (ei1==0) ? c.eid+1 : c.eid+1;

//        unsigned int cid = cIndex;
//        cIndex += 2;

//        MatrixDerivRowIterator c1_it = c1.writeLine(cid);
//        c1_it.addCol(i, c.dirProj);

//        c1_it = c1.writeLine(cid + 1);
//        c1_it.setCol(i, c.dirOrtho);

//        MatrixDerivRowIterator c2_it = c2.writeLine(cid);
//        c2_it.addCol(ei1, -c.dirProj * (1-c.r2));
//        c2_it.addCol(ei2, -c.dirProj * c.r2);


//        c2_it = c2.writeLine(cid + 1);
//        c2_it.addCol(ei1, -c.dirOrtho * (1-c.r2));
//        c2_it.addCol(ei2, -c.dirOrtho * c.r2);

//        c.thirdConstraint = 0;

//        if (c.r < 0)
//        {
//            c.thirdConstraint = c.r;
//            cIndex++;

//            c1_it = c1.writeLine(cid + 2);
//            c1_it.setCol(i, c.dirAxe);

//            c2_it = c2.writeLine(cid + 2);
//            c2_it.addCol(ei1, -c.dirAxe);
//        }
//        else if (c.r > c.Q1Q2)
//        {
//            c.thirdConstraint = c.r - c.Q1Q2;
//            cIndex++;

//            c1_it = c1.writeLine(cid + 2);
//            c1_it.setCol(i, -c.dirAxe);

//            c2_it = c2.writeLine(cid + 2);
//            c2_it.addCol(ei2, c.dirAxe);
//        }
//        c.cid = cid;
//    }

//    c1_d.endEdit();
//    c2_d.endEdit();

//}


template<class DataTypes>
void ProjectionEngine<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
    //printf("ProjectionEngine<DataTypes>::draw(const core::visual::VisualParams* vparams) before \n");
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
    //printf("ProjectionEngine<DataTypes>::draw(const core::visual::VisualParams* vparams) After \n");

    //drawLinesBetweenPoints(vparams);
    vparams->drawTool()->restoreLastState();
}

template<class DataTypes>
void ProjectionEngine<DataTypes>::drawLinesBetweenPoints(const core::visual::VisualParams* vparams)
{
    const VecCoord & positions  = d_dest.getValue(); // this->mstate2->read(core::ConstVecCoordId::position())->getValue();
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


} // namespace sofa::component::constraintset


#endif
