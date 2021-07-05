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
#include <sofa/defaulttype/Vec.h>
#include <SofaConstraint/BilateralInteractionConstraint.h>

#include "CosseratNeedleSlidingConstraint.h"

namespace sofa::component::constraintset
{

using sofa::core::objectmodel::ComponentState;
using sofa::helper::WriteAccessor;

using sofa::core::objectmodel::ComponentState;
using sofa::core::visual::VisualParams;
using sofa::defaulttype::BaseVector;
using sofa::helper::ReadAccessor;
using sofa::defaulttype::Vec4f;
using sofa::defaulttype::Vector3;
using sofa::type::vector;
using sofa::helper::OptionsGroup;

template<class DataTypes>
CosseratNeedleSlidingConstraint<DataTypes>::CosseratNeedleSlidingConstraint(MechanicalState* object)
    : Inherit1(object)

    , d_value(initData(&d_value, "value",
                       "Displacement or force to impose.\n"))

    , d_valueIndex(initData(&d_valueIndex, (unsigned int) 0, "valueIndex",
                            "Index of the value (in InputValue vector) that we want to impose \n"
                            "If unspecified the default value is {0}"))

    , d_valueType(initData(&d_valueType, OptionsGroup(2,"displacement","force"), "valueType",
                           "displacement = the contstraint will impose the displacement provided in data value[valueIndex] \n"
                           "force = the contstraint will impose the force provided in data value[valueIndex] \n"
                           "If unspecified, the default value is displacement"))
{

}

template<class DataTypes>
CosseratNeedleSlidingConstraint<DataTypes>::~CosseratNeedleSlidingConstraint()
{
}


template<class DataTypes>
void CosseratNeedleSlidingConstraint<DataTypes>::init()
{
    Inherit1::init();

    internalInit();
}


template<class DataTypes>
void CosseratNeedleSlidingConstraint<DataTypes>::reinit()
{
    internalInit();
}

template<class DataTypes>
void CosseratNeedleSlidingConstraint<DataTypes>::internalInit()
{
    if(d_value.getValue().size()==0)
    {
        WriteAccessor<Data<vector<Real>>> value = d_value;
        value.resize(1,0.);
    }

    // check for errors in the initialization
    if(d_value.getValue().size()<d_valueIndex.getValue())
    {
        msg_warning() << "Bad size for data value (size="<< d_value.getValue().size()<<"), or wrong value for data valueIndex (valueIndex="<<d_valueIndex.getValue()<<"). Set default valueIndex=0.";
        d_valueIndex.setValue(0);
    }
}



template<class DataTypes>
void CosseratNeedleSlidingConstraint<DataTypes>::buildConstraintMatrix(const ConstraintParams* cParams, DataMatrixDeriv &cMatrix, unsigned int &cIndex, const DataVecCoord &x)
{
    if(d_componentState.getValue() != ComponentState::Valid)
        return ;

    SOFA_UNUSED(cParams);
    MatrixDeriv& matrix = *cMatrix.beginEdit();
    VecCoord positions = x.getValue();
    m_constraintId= cIndex;

    for (unsigned int i=0; i<positions.size(); i++)
    {
        MatrixDerivRowIterator c_it = matrix.writeLine(cIndex);
        c_it.addCol(i, Coord(0,1,0)); // instead of vector3(0,1,0) use the directtion of the projection
        MatrixDerivRowIterator c_it_1 = matrix.writeLine(cIndex+1);
        c_it_1.addCol(i, Coord(0,0,1)); // instead of vector3(0,1,0) use the directtion of the projection
        cIndex +=2;
    }
    cMatrix.endEdit();
    m_nbLines = cIndex - m_constraintId;
}


template<class DataTypes>
void CosseratNeedleSlidingConstraint<DataTypes>::getConstraintViolation(const ConstraintParams* cParams,
                                                                        BaseVector *resV,
                                                                        const BaseVector *Jdx)
{
    if(d_componentState.getValue() != ComponentState::Valid)
        return ;

    SOFA_UNUSED(cParams);
    ReadAccessor<Data<VecCoord>> positions = m_state->readPositions();

    msg_info("CosseratNeedleSlidingConstraint") << "The size of constraint is "<< positions.size();

    if(Jdx->size()==0){
        for (unsigned int i = 0; i < positions.size(); i++){
            Real dfree1 =  positions[i][1];
            Real dfree2 =  positions[i][2];

            resV->set(m_constraintId + 2*i   , dfree1);
            resV->set(m_constraintId + 2*i +1, dfree2);
        }

    }else{
        for (unsigned int i = 0; i < positions.size(); i++){
            Real dfree1 = Jdx->element(2*i)   + positions[i][1];
            Real dfree2 = Jdx->element(2*i+1) + positions[i][2];
            resV->set(m_constraintId + 2*i   , dfree1);
            resV->set(m_constraintId + 2*i +1, dfree2);

        }
    }
}

template<class DataTypes>
void CosseratNeedleSlidingConstraint<DataTypes>::getConstraintResolution(const ConstraintParams*,
                                                                         std::vector<core::behavior::ConstraintResolution*>& resTab,
                                                                         unsigned int& offset)
{
    ReadAccessor<Data<VecCoord>> positions = m_state->readPositions();
    double imposedValue= 1.0;
    for (size_t i = 0; i < positions.size(); i++){

        resTab[offset++] = new BilateralConstraintResolution();
        resTab[offset++] = new BilateralConstraintResolution();
    }
}


template<class DataTypes>
void CosseratNeedleSlidingConstraint<DataTypes>::draw(const VisualParams* vparams)
{
    if(d_componentState.getValue() != ComponentState::Valid)
        return ;

}

} // namespace sofa

