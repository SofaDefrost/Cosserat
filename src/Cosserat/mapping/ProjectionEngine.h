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
#ifndef COSSERAT_ProjectionEngine_H
#define COSSERAT_ProjectionEngine_H

#include <Cosserat/config.h>

#include <sofa/core/behavior/PairInteractionConstraint.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/DataEngine.h>
#include <sofa/type/Vec.h>
#include <sofa/simulation/AnimateBeginEvent.h>

#include <iostream>


namespace sofa::component::constraintset
{


using sofa::core::ConstraintParams;

template<class DataTypes>
class ProjectionEngine : public core::DataEngine
{
public:
    //    SOFA_CLASS(SOFA_TEMPLATE(ProjectionEngine,DataTypes), SOFA_TEMPLATE(core::behavior::PairInteractionConstraint,DataTypes));
    SOFA_CLASS(SOFA_TEMPLATE(ProjectionEngine,DataTypes),core::DataEngine);

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::MatrixDeriv MatrixDeriv;
    typedef typename DataTypes::MatrixDeriv::RowIterator MatrixDerivRowIterator;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename Coord::value_type Real;
    typedef typename core::behavior::MechanicalState<DataTypes> MechanicalState;
    //typedef typename core::behavior::PairInteractionConstraint<DataTypes> Inherit;
    typedef sofa::type::vector<DataTypes> VecData;

    typedef core::objectmodel::Data<VecCoord>		DataVecCoord;
    typedef core::objectmodel::Data<VecDeriv>		DataVecDeriv;
    typedef core::objectmodel::Data<MatrixDeriv>    DataMatrixDeriv;

protected:
    Data<VecCoord> d_from; ///< input vector
    Data<VecCoord> d_dest; ///< vector to substract to input
    Data<VecCoord> d_output;


    //Data<Deriv> d_force; ///< interaction force

    //    Real m_dist;	// constraint violation
    //    Real m_thirdConstraint; // 0 if A<proj<B, -1 if proj<A, 1 if B<proj
    //    bool m_yetIntegrated;
    //    unsigned int m_cid;


    ProjectionEngine();

    virtual ~ProjectionEngine(){}

public:
    void init() override;
    void reinit() override;

    void handleEvent(core::objectmodel::Event *ev) override;

    void computeProximity();

    void draw(const core::visual::VisualParams* vparams) override;
    void drawLinesBetweenPoints(const core::visual::VisualParams* vparams);



private:
    // storage of force
    Deriv  m_dirAxe, m_dirProj, m_dirOrtho;

    typedef struct {
        double fact;
        int p1, p2;
        int eid, cid;
        Real alpha, beta;
        Coord P, Q;
        Deriv dirAxe, t1, t2;
        Real r,r2, Q1Q2;
        Real dist; //violation
        Real thirdConstraint;
        //Deriv  m_dirAxe, m_dirProj, m_dirOrtho;
    } Constraint;

    type::vector<Constraint> m_constraints;
    unsigned int m_step;

    void doUpdate() override
    {
        computeProximity();
    }

};

//#if !defined(SOFA_COSSERAT_CPP_ProjectionEngine)
//extern template class SOFA_CONSTRAINT_API ProjectionEngine< sofa::defaulttype::Vec3Types >;

//#endif

} // namespace sofa::component::constraintset


#endif // COSSERAT_ProjectionEngine_H
