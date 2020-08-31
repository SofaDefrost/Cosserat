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
#ifndef SOFA_COMPONENT_MAPPING_POEMAPING_INL
#define SOFA_COMPONENT_MAPPING_POEMAPING_INL

#include <sofa/core/Multi2Mapping.inl>
#include "DifferenceMultiMapping.h"
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <string>
#include <sofa/core/visual/VisualParams.h>

#include <sofa/helper/AdvancedTimer.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/helper/logging/Message.h>



namespace sofa
{
using sofa::core::objectmodel::BaseContext ;
using sofa::helper::AdvancedTimer;
using sofa::helper::WriteAccessor;


namespace component
{

namespace mapping
{
template <class TIn1, class TIn2, class TOut>
DifferenceMultiMapping<TIn1, TIn2, TOut>::DifferenceMultiMapping()
    : m_fromModel1(NULL)
    , m_fromModel2(NULL)
    , m_toModel(NULL)
{}


template <class TIn1, class TIn2, class TOut>
void DifferenceMultiMapping<TIn1, TIn2, TOut>::initiatTopologies()
{
    m_toModel = this->getToModels()[0];
    if (! m_toModel) {
        std::cout << " No output mechanical state found. Consider setting the " << this->toModels.getName() << " attribute."<< std::endl;
    }else {
        msg_warning() << " Need to initialize everything " ;
        //        BaseMeshTopology * topology = nullptr;
        //        this->toModel->getContext()->get(topology);
        //        if (topology) {
        //            m_toModel.set(topology);
        //        }
    }
}


// _________________________________________________________________________________________

template <class TIn1, class TIn2, class TOut>
void DifferenceMultiMapping<TIn1, TIn2, TOut>::init()
{
    if(this->getFromModels1().empty())
    {
        msg_error() << "Error while initializing ; input getFromModels1 not found" ;
        return;
    }

    if(this->getFromModels2().empty())
    {
        msg_error() << "Error while initializing ; output getFromModels2 not found" ;
        return;
    }

    if(this->getToModels().empty())
    {
        msg_error() << "Error while initializing ; output Model not found" ;
        //return;
    }
    m_fromModel1 = this->getFromModels1()[0];
    m_fromModel2 = this->getFromModels2()[0];
    m_toModel = this->getToModels()[0];

    m_toModel = m_fromModel1;

    //initiatTopologies();
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
void DifferenceMultiMapping<TIn1, TIn2, TOut>::computeProximity(const In1VecCoord &x1, const In2VecCoord &x2){

    In1VecCoord from = x1;
    In2VecCoord dst  = x2;
    m_constraints.clear();

    size_t szFrom = from.size();
    size_t szDst = dst.size();

    printf ("###################### The size from : %d \n", szFrom);
    printf ("###################### The size szDst : %d \n", szDst);
    //For each point in the FEM find the closest edge of the cable
    for (size_t i = 0 ; i < szFrom; i++) {
        Coord2 P = from[i];
        Constraint constraint;

        //std::cout << "P "<< P << std::endl;

        // find the min distance between a from mstate point and it's projection on each edge of the cable (destination mstate)
        Real min_dist = std::numeric_limits<Real>::max();
        Real min_dist2 = std::numeric_limits<Real>::max();
        for (size_t j = 0; j < szDst-1; j++) {
            Coord1 Q1 = dst[j];
            Coord1 Q2 = dst[j+1];

            // the axis
            Coord1 v = Q2 -Q1;
            Real fact_v = dot(P-Q1,v) / dot(v,v);
            //if(fact_v <= 0.0) continue;

            // projection of the point on the axis
            Deriv1 dirAxe = v;
            dirAxe.normalize();
            Real alpha = (P-Q1) * dirAxe;
            Deriv1 proj = Q1 + dirAxe * alpha;
            Deriv1 distVec = P - proj; // violation vector
            Real dist = distVec.norm(); // constraint violation

            //std::cout << i << "====>The distance between is : "<< dist << "   j: "<< j << " fact_v : "<< fact_v<<  std::endl;


            if(std::abs(fact_v) < min_dist){
                //printf ("#################Begin \n");
                min_dist = std::abs(fact_v) ; //fact_v ;

                //                // projection of the point on the axis
                //                Deriv1 dirAxe = v;
                //                dirAxe.normalize();
                //                Real alpha = (P-Q1) * dirAxe;
                //                Deriv1 proj = Q1 + dirAxe * alpha;

                constraint.proj = proj;
                constraint.Q = from[i];
                constraint.eid = j;
                constraint.alpha = alpha;
                constraint.dirAxe = dirAxe;
                /////
                constraint.Q1Q2 = v.norm();
                constraint.r2 = fact_v;
                //std::cout << "alpha :"<< alpha << " ; r2 :"<< constraint.alpha << std::endl;

                // We move the constraint point onto the projection
                Deriv1 t1 = P - proj; // violation vector
                constraint.dist = t1.norm(); // constraint violation
                t1.normalize(); // direction of the constraint
                if(t1.norm()<1.0e-1 && dirAxe[2] < 0.99){
                    Vector3 temp = Vector3(dirAxe[0],dirAxe[1],dirAxe[2]+50.0);
                    t1 = cross(dirAxe,temp);
                    t1.normalize();
                    constraint.t1 = t1;
                }
                if(t1.norm()<1.0e-1){
                    Vector3 temp = Vector3(dirAxe[0],dirAxe[1]+50.0,dirAxe[2]);
                    t1 = cross(dirAxe,temp);
                    t1.normalize();
                    constraint.t1 = t1;
                }

                //printf("1 ########################################################### \n");
                //std::cout << " i : "<< i << " ==> j :"<< j << " ==> eid : "<< constraint.eid << std::endl;
                constraint.t1 = t1;
                //                std::cout << "t1 : "<< t1 << std::endl;
                //                std::cout << "=============> t1 : "<< constraint.t1 << std::endl;
                //tangential 2
                Deriv1 t2 = cross(t1, dirAxe);  t2.normalize();
                constraint.t2 = t2;

            }
        }
        //        std::cout << i << " Closets edge is "<< constraint.eid << std::endl;
        //std::cout << " i : "<< i << " ==> dist: :"<< constraint.dist << " ==> eid : "<< constraint.eid << std::endl;
        m_constraints.push_back(constraint);
    }
    //printf("ProjectionEngine<DataTypes>::computeProximity After\n");
}


template <class TIn1, class TIn2, class TOut>
void DifferenceMultiMapping<TIn1, TIn2, TOut>::apply(
        const core::MechanicalParams* /* mparams */, const helper::vector<OutDataVecCoord*>& dataVecOutPos,
        const helper::vector<const In1DataVecCoord*>& dataVecIn1Pos ,
        const helper::vector<const In2DataVecCoord*>& dataVecIn2Pos)
{

    if(dataVecOutPos.empty() || dataVecIn1Pos.empty() || dataVecIn2Pos.empty())
        return;

    printf("///Do Apply//We need only one input In model and input Root model (if present) \n");
    const In1VecCoord& in1 = dataVecIn1Pos[0]->getValue();
    const In2VecCoord& in2 = dataVecIn2Pos[0]->getValue();

    computeProximity(in1,in2);

    OutVecCoord& out = *dataVecOutPos[0]->beginEdit();
    size_t sz = m_constraints.size();
    out.resize(sz);

    for(unsigned int i=0; i<sz; i++){
        out[i][0] = 0.0;
        out[i][1] = m_constraints[i].dist;
        out[i][2] = 0.0;
    }
    dataVecOutPos[0]->endEdit();
}



template <class TIn1, class TIn2, class TOut>
void DifferenceMultiMapping<TIn1, TIn2, TOut>:: applyJ(
        const core::MechanicalParams* /* mparams */, const helper::vector< OutDataVecDeriv*>& dataVecOutVel,
        const helper::vector<const In1DataVecDeriv*>& dataVecIn1Vel,
        const helper::vector<const In2DataVecDeriv*>& dataVecIn2Vel) {
    if(dataVecOutVel.empty() || dataVecIn1Vel.empty() ||dataVecIn2Vel.empty() )
        return;
    const In1VecDeriv& in1 = dataVecIn1Vel[0]->getValue();
    const In2VecDeriv& in2 = dataVecIn2Vel[0]->getValue();
    OutVecDeriv& outVel = *dataVecOutVel[0]->beginEdit();

    //const OutVecCoord& out = m_toModel->read(core::ConstVecCoordId::position())->getValue();
    size_t sz = m_constraints.size();
    outVel.resize(sz);
    for (size_t i = 0 ; i < sz; i++) {
        Constraint& c = m_constraints[i];
        int ei1 = c.eid;
        int ei2 = c.eid+1;
        //        std::cout << " ei1 : " << ei1 << " ei2 : "<< ei2 << std::endl;
        Real v0 = c.dirAxe * (c.alpha * in2[ei1] + (1-c.alpha) * in2[ei2] - in1[i]);
        Real v1 = c.t1     * (c.alpha * in2[ei1] + (1-c.alpha) * in2[ei2] - in1[i]);
        Real v2 = c.t2     * (c.alpha * in2[ei1] + (1-c.alpha) * in2[ei2] - in1[i]);

        outVel[i] = OutDeriv(v0,v1,v2);
    }
    dataVecOutVel[0]->endEdit();
}


template <class TIn1, class TIn2, class TOut>
void DifferenceMultiMapping<TIn1, TIn2, TOut>::applyJT(
        const core::MechanicalParams* /*mparams*/, const helper::vector< In1DataVecDeriv*>& dataVecOut1Force,
        const helper::vector< In2DataVecDeriv*>& dataVecOut2Force,
        const helper::vector<const OutDataVecDeriv*>& dataVecInForce)  {

    if(dataVecOut1Force.empty() || dataVecInForce.empty() || dataVecOut2Force.empty())
        return;

    const OutVecDeriv& in = dataVecInForce[0]->getValue();

    In1VecDeriv& out1 = *dataVecOut1Force[0]->beginEdit();
    In2VecDeriv& out2 = *dataVecOut2Force[0]->beginEdit();

    //Compute output forces
    size_t sz = m_constraints.size();

    for (size_t i = 0 ; i < sz; i++) {
        Constraint& c = m_constraints[i];
        int ei1 = c.eid;
        int ei2 = c.eid+1;
        OutDeriv f = in[i];
        std::cout << " ================+++++++++>>>>> The force : " << f << std::endl;
        Deriv2 f1   = (f[0] * c.dirAxe) + (f[1] * c.t1) + (f[2] * c.t2) ;
        Deriv1 f2   = (c.alpha * f[0]*c.dirAxe) + (c.alpha* f[1] * c.t1) + (c.alpha * f[2] * c.t2);
        Deriv1 f2_1 = ((1-c.alpha) * f[0]*c.dirAxe )+ ((1-c.alpha) * f[1]*c.t1) + ((1-c.alpha) * f[2]*c.t2);

        std::cout << " f1 : " << f1 << "   f&_1: " << f2_1 << " ; f2 : "<< f2 << std::endl;

        out1[i] -= f1;
        out2[ei1] += f2;
        out2[ei2] += f2_1;
    }

    dataVecOut1Force[0]->endEdit();
    dataVecOut2Force[0]->endEdit();
}

//___________________________________________________________________________
template <class TIn1, class TIn2, class TOut>
void DifferenceMultiMapping<TIn1, TIn2, TOut>::applyJT(
        const core::ConstraintParams*/*cparams*/ , const helper::vector< In1DataMatrixDeriv*>&  dataMatOut1Const,
        const helper::vector< In2DataMatrixDeriv*>&  dataMatOut2Const ,
        const helper::vector<const OutDataMatrixDeriv*>& dataMatInConst)
{
    printf("######################## Inside the applyJ constrainte level 1 \n");
    if(dataMatOut1Const.empty() || dataMatOut2Const.empty() || dataMatInConst.empty() )
        return;


    //We need only one input In model and input Root model (if present)
    In1MatrixDeriv& out1 = *dataMatOut1Const[0]->beginEdit(); // constraints on the FEM cable points
    In2MatrixDeriv& out2 = *dataMatOut2Const[0]->beginEdit(); // constraints on the frames cable points
    const OutMatrixDeriv& in = dataMatInConst[0]->getValue(); // input constraints defined on the mapped point
    const OutVecCoord& mappedPoints = m_toModel->read(core::ConstVecCoordId::position())->getValue();
    const In1DataVecCoord* x1fromData = m_fromModel1->read(core::ConstVecCoordId::position());
    const In1VecCoord x1from = x1fromData->getValue();

    typename OutMatrixDeriv::RowConstIterator rowIt    = in.begin()  ;
    typename OutMatrixDeriv::RowConstIterator rowItEnd = in.end();


    if(!rowIt.index()) printf( "===++++==========>  The rowIt is empty \n");



    for (rowIt = in.begin(); rowIt != rowItEnd; ++rowIt)
    {
        std::cout<<"************* Apply JT (MatrixDeriv) iteration on line ";
        std::cout<<rowIt.index();
        std::cout<<"*************  "<<std::endl;

        typename OutMatrixDeriv::ColConstIterator colIt = rowIt.begin();
        typename OutMatrixDeriv::ColConstIterator colItEnd = rowIt.end();

        // Creates a constraints if the input constraint is not empty.
        if (colIt == colItEnd)
        {
            std::cout<<"no column for this constraint"<<std::endl;
            continue;
        }
        typename In1MatrixDeriv::RowIterator o1 = out1.writeLine(rowIt.index()); // we store the constraint number
        typename In2MatrixDeriv::RowIterator o2 = out2.writeLine(rowIt.index());


        while (colIt != colItEnd)
        {
            int childIndex = colIt.index();
            Constraint c = m_constraints[childIndex];
            const OutDeriv h = colIt.val();
            int indexBeam =  c.eid;
            std::cout << " ===> childIndex : " << childIndex<< " ===> indexBeam : " << indexBeam << std::endl;
            std::cout << " ===> h : " << h << std::endl;

            Deriv2 h1 = (h[0] * c.dirAxe) + (h[1] * c.t1) + (h[2] * c.t2) ;
            Deriv1 h2 = (c.alpha * h[0]*c.dirAxe) + (c.alpha* h[1] * c.t1) + (c.alpha * h[2] * c.t2);
            Deriv1 h2_1 = ((1-c.alpha) * h[0]*c.dirAxe )+ ((1-c.alpha) * h[1]*c.t1) + ((1-c.alpha) * h[2]*c.t2);

            std::cout << " ==> t1 : "<< c.t1 << " ==> t2 : "<< c.t2 << std::endl;

            std::cout << " ===> h1 : " << h1<< " ===> h2 : " << h2 << " ===> h2_1 : " << h2_1 << std::endl;

            o1.addCol(childIndex, -h1);
            o2.addCol(indexBeam, h2);
            o2.addCol(indexBeam+1, h2_1);

            colIt++;
        }
    }

    ////// END ARTICULATION SYSTEM MAPPING
    dataMatOut1Const[0]->endEdit();
    dataMatOut2Const[0]->endEdit();
    printf("######################## Inside the applyJ constrainte level 2 \n");
}


template <class TIn1, class TIn2, class TOut>
void DifferenceMultiMapping<TIn1, TIn2, TOut>::draw(const core::visual::VisualParams* vparams)
{
    //printf("CosseratSlidingConstraint<DataTypes>::draw(const core::visual::VisualParams* vparams) before \n");
    if (!vparams->displayFlags().getShowInteractionForceFields())
        return;

    vparams->drawTool()->saveLastState();

    vparams->drawTool()->disableLighting();

    sofa::defaulttype::RGBAColor color;
    //    Constraint& c = m_constraints[0];

    //    if(c.thirdConstraint<0)
    //        color = sofa::defaulttype::RGBAColor::yellow();
    //    else if(c.thirdConstraint>0)
    //        color = sofa::defaulttype::RGBAColor::green();
    //    else
    color = sofa::defaulttype::RGBAColor::magenta();

    std::vector<sofa::defaulttype::Vector3> vertices;
    //    vertices.push_back(DataTypes::getCPos((this->mstate1->read(core::ConstVecCoordId::position())->getValue())[d_m1.getValue()]));

    //    vparams->drawTool()->drawPoints(vertices, 10, color);
    //    vertices.clear();

    //    color = sofa::defaulttype::RGBAColor::blue();
    //    vertices.push_back(DataTypes::getCPos((this->mstate2->read(core::ConstVecCoordId::position())->getValue())[d_m2a.getValue()]));
    //    vertices.push_back(DataTypes::getCPos((this->mstate2->read(core::ConstVecCoordId::position())->getValue())[d_m2b.getValue()]));
    vparams->drawTool()->drawLines(vertices, 1, color);

    for (size_t i =0 ; i < m_constraints.size(); i++) {
        color = sofa::defaulttype::RGBAColor::green();
        //        std::cout << " Projection : "<< m_constraints[i].proj << " ;  Q :" << m_constraints[i].Q << std::endl;
        vertices.push_back(m_constraints[i].proj);
        vertices.push_back(m_constraints[i].Q);
        vparams->drawTool()->drawLines(vertices, 1, color);
    }
    //printf("CosseratSlidingConstraint<DataTypes>::draw(const core::visual::VisualParams* vparams) After \n");

    //    drawLinesBetweenPoints(vparams);
    vparams->drawTool()->restoreLastState();
}

template <class TIn1, class TIn2, class TOut>
void DifferenceMultiMapping<TIn1, TIn2, TOut>::drawLinesBetweenPoints(const core::visual::VisualParams* vparams)
{
    //    const In2VecCoord & positions  = this->mstate2->read(core::ConstVecCoordId::position())->getValue();
    //    sofa::defaulttype::RGBAColor color;
    //    color = sofa::defaulttype::RGBAColor::magenta();
    //    std::vector<sofa::defaulttype::Vector3> vertices;
    //    for (unsigned int i=0; i<positions.size()-1; i++)
    //    {
    //        vertices.push_back(positions[i]);
    //        vertices.push_back(positions[i+1]);
    //    }

    //    vparams->drawTool()->drawLines(vertices, 1.5, color);
}


} // namespace mapping

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_MAPPING_POEMAPING_INL
