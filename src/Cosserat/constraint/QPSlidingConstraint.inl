/******************************************************************************
*               SOFA, Simulation Open-Framework Architecture                  *
*                (c) 2006-2018 INRIA, USTL, UJF, CNRS, MGH                    *
*                                                                             *
* This library is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This library is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this library; if not, write to the Free Software Foundation,     *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
*******************************************************************************
*                           Plugin Cosserat  v1.0                             *
*                                                                             *
* This plugin is also distributed under the GNU LGPL (Lesser General          *
* Public License) license with the same conditions than SOFA.                 *
*                                                                             *
* Contributors: Defrost team  (INRIA, University of Lille, CNRS,              *
*               Ecole Centrale de Lille)                                      *
*                                                                             *
* Contact information: https://project.inria.fr/softrobot/contact/            *
*                     adagolodjo@protonmail.com                               *
******************************************************************************/

#pragma once

#include <sofa/core/visual/VisualParams.h>
#include <sofa/type/Vec.h>
#include <sofa/component/constraint/lagrangian/model/BilateralLagrangianConstraint.h>

#include "QPSlidingConstraint.h"

namespace sofa::component::constraintset
{

using sofa::core::objectmodel::ComponentState;
using sofa::helper::WriteAccessor;
using sofa::core::objectmodel::ComponentState;
using sofa::core::visual::VisualParams;
using sofa::linearalgebra::BaseVector;
using sofa::helper::ReadAccessor;
using sofa::type::Vec4f;
using sofa::type::Vec3;
using sofa::type::vector;
using sofa::helper::OptionsGroup;
using sofa::component::constraint::lagrangian::model::BilateralConstraintResolution;

template<class DataTypes>
QPSlidingConstraint<DataTypes>::QPSlidingConstraint(MechanicalState* object)
    : Inherit1(object)

    , d_value(initData(&d_value, "value",
                       "Displacement or force to impose.\n"))

    , d_valueIndex(initData(&d_valueIndex, (unsigned int) 0, "valueIndex",
                            "Index of the value (in InputValue vector) that we want to impose \n"
                            "If unspecified the default value is {0}"))

    , d_valueType(initData(&d_valueType, {"displacement","force"}, "valueType",
                           "displacement = the contstraint will impose the displacement provided in data value[valueIndex] \n"
                           "force = the contstraint will impose the force provided in data value[valueIndex] \n"
                           "If unspecified, the default value is displacement"))
{

}

template<class DataTypes>
QPSlidingConstraint<DataTypes>::~QPSlidingConstraint()= default;


template<class DataTypes>
void QPSlidingConstraint<DataTypes>::init()
{
    Inherit1::init();

    internalInit();
}


template<class DataTypes>
void QPSlidingConstraint<DataTypes>::reinit()
{
    internalInit();
}

template<class DataTypes>
void QPSlidingConstraint<DataTypes>::internalInit()
{
    if(d_value.getValue().size()==0)
    {
        WriteAccessor<Data<vector<Real>>> value = d_value;
        value.resize(1);
        value[0] = 0.;
    }

    // check for errors in the initialization
    if(d_value.getValue().size()<d_valueIndex.getValue())
    {
        msg_warning() << "Bad size for data value (size="<< d_value.getValue().size()<<"), or wrong value for data valueIndex (valueIndex="<<d_valueIndex.getValue()<<"). Set default valueIndex=0.";
        d_valueIndex.setValue(0);
    }
}



template<class DataTypes>
void QPSlidingConstraint<DataTypes>::buildConstraintMatrix(const ConstraintParams* cParams, DataMatrixDeriv &cMatrix, unsigned int &cIndex, const DataVecCoord &x)
{
    if(d_componentState.getValue() != ComponentState::Valid)
        return ;

    SOFA_UNUSED(cParams);
    MatrixDeriv& matrix = *cMatrix.beginEdit();
    VecCoord positions = x.getValue();
    m_constraintIndex.setValue(cIndex);
    const auto& constraintIndex = sofa::helper::getReadAccessor(m_constraintIndex);

    for (unsigned int i=0; i<positions.size(); i++)
    {
        if(i < positions.size()-1){
            MatrixDerivRowIterator c_it = matrix.writeLine(cIndex);
            c_it.addCol(i, Coord(0,1,0));
            MatrixDerivRowIterator c_it_1 = matrix.writeLine(cIndex+1);
            c_it_1.addCol(i, Coord(0,0,1));
            cIndex +=2;
        }else{
            MatrixDerivRowIterator c_it = matrix.writeLine(cIndex);
            c_it.addCol(i, Coord(1,0,0));

            MatrixDerivRowIterator c_it_1 = matrix.writeLine(cIndex+1);
            c_it_1.addCol(i, Coord(0,1,0));
            MatrixDerivRowIterator c_it_2 = matrix.writeLine(cIndex+2);
            c_it_2.addCol(i, Coord(0,0,1));
            cIndex +=3;
        }

    }
    cMatrix.endEdit();
    m_nbLines = cIndex - constraintIndex;
}


template<class DataTypes>
void QPSlidingConstraint<DataTypes>::getConstraintViolation(const ConstraintParams* cParams,
                                                            BaseVector *resV,
                                                            const BaseVector *Jdx)
{
    if(d_componentState.getValue() != ComponentState::Valid)
        return ;

    SOFA_UNUSED(cParams);
    ReadAccessor<Data<VecCoord>> positions = m_state->readPositions();
    const auto& constraintIndex = sofa::helper::getReadAccessor(m_constraintIndex);

    if(Jdx->size()==0){
        //std::cout << "Size JDX = "<< Jdx->size() << std::endl;
        for (size_t i = 0; i < positions.size(); i++){
            if( i < positions.size()-1){
                Real dfree1 =  positions[i][1];
                Real dfree2 =  positions[i][2];

                resV->set(constraintIndex + 2*i   , dfree1);
                resV->set(constraintIndex + 2*i +1, dfree2);
            }else{
                Real dfree0 =  positions[i][0];
                Real dfree1 =  positions[i][1];
                Real dfree2 =  positions[i][2];

                resV->set(constraintIndex + 2*i   , dfree0);
                resV->set(constraintIndex + 2*i +1, dfree1);
                resV->set(constraintIndex + 2*i +2, dfree2);
            }
        }

    }else{
        //std::cout << "Size JDX = "<< Jdx->size() << std::endl;
        for (size_t i = 0; i < positions.size(); i++){
            if( i < positions.size()-1){
                Real dfree1 = Jdx->element(2*i)   + positions[i][1];
                Real dfree2 = Jdx->element(2*i+1) + positions[i][2];

                resV->set(constraintIndex + 2*i   , dfree1);
                resV->set(constraintIndex + 2*i +1, dfree2);
            }else{
                //std::cout << " The laste position : "<< positions[i] << "; Jdx->element() "<< Jdx->element(2*i) <<" "<< Jdx->element(2*i+1) <<" "<< Jdx->element(2*i+2) << std::endl;
                Real dfree0 = Jdx->element(2*i)   + positions[i][0];
                Real dfree1 = Jdx->element(2*i+1) + positions[i][1];
                Real dfree2 = Jdx->element(2*i+2) + positions[i][2];

                //std::cout << " m_constraintId + 2*i : "<< m_constraintId + 2*i << std::endl;
                resV->set(constraintIndex + 2*i   , dfree0);
                resV->set(constraintIndex + 2*i +1, dfree1);
                resV->set(constraintIndex + 2*i +2, dfree2);
            }
        }
    }
}

template<class DataTypes>
void QPSlidingConstraint<DataTypes>::getConstraintResolution(const ConstraintParams*,
                                                             std::vector<core::behavior::ConstraintResolution*>& resTab,
                                                             unsigned int& offset)
{
    ReadAccessor<Data<VecCoord>> positions = m_state->readPositions();

    for (size_t i = 0; i < positions.size(); i++){

        resTab[offset++] = new BilateralConstraintResolution();
        resTab[offset++] = new BilateralConstraintResolution();
        if(i == positions.size()-1){
            resTab[offset++] = new BilateralConstraintResolution();
        }
    }
}


template<class DataTypes>
void QPSlidingConstraint<DataTypes>::draw(const VisualParams* vparams)
{
    SOFA_UNUSED(vparams);
    if(d_componentState.getValue() != ComponentState::Valid)
        return ;

}

} // namespace sofa

