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
*                           Plugin Cosserat v1.0                              *
*                                                                             *
* This plugin is also distributed under the GNU LGPL (Lesser General          *
* Public License) license with the same conditions than SOFA.                 *
*                                                                             *
* Contributors: Defrost team  (INRIA, University of Lille, CNRS,              *
*               Ecole Centrale de Lille)                                      *
*                                                                             *
* Contact information: https://project.inria.fr/softrobot/contact/            *
*                                                                             *
******************************************************************************/

#pragma once

#include "CosseratActuatorConstraint.h"

using sofa::helper::OptionsGroup;

namespace sofa::component::constraintset
{

using sofa::core::objectmodel::ComponentState;
using sofa::helper::WriteAccessor;

template<class DataTypes>
CosseratActuatorConstraint<DataTypes>::CosseratActuatorConstraint(MechanicalState* object)
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
    , d_integral(initData(&d_integral,"integral","helper vector of H_i ()"))
{
}

template<class DataTypes>
CosseratActuatorConstraint<DataTypes>::~CosseratActuatorConstraint()
{
}

template<class DataTypes>
void CosseratActuatorConstraint<DataTypes>::init()
{
    Inherit1::init();
    internalInit();
}

template<class DataTypes>
void CosseratActuatorConstraint<DataTypes>::reinit()
{
    internalInit();
}

template<class DataTypes>
void CosseratActuatorConstraint<DataTypes>::internalInit()
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
void CosseratActuatorConstraint<DataTypes>::buildConstraintMatrix(const ConstraintParams* cParams, DataMatrixDeriv &cMatrix, unsigned int &cIndex, const DataVecCoord &x)
{
    if(d_componentState.getValue() != ComponentState::Valid)
        return ;

    SOFA_UNUSED(cParams);

    m_constraintIndex.setValue(cIndex);
    const auto& constraintIndex = sofa::helper::getReadAccessor(m_constraintIndex);

    MatrixDeriv& matrix = *cMatrix.beginEdit();

    MatrixDerivRowIterator rowIterator = matrix.writeLine(constraintIndex);

    VecCoord positions = x.getValue();

    const SetIndexArray &indices = d_indices.getValue();
    helper::ReadAccessor<Data<type::vector<Coord>>> integral = d_integral;

    for (unsigned int i=0; i<indices.size(); i++)
    {
        /*
        Coord previousPosition;
        Coord currentPosition;
        Coord nextPosition;

        int previousIndex = indices[i-1];
        int currentIndex  = indices[i];
        int nextIndex     = indices[i+1];

        previousPosition = positions[previousIndex];
        currentPosition  = positions[currentIndex];
        nextPosition     = positions[nextIndex];

        Deriv directionBeyond = previousPosition - currentPosition;
        directionBeyond.normalize();

        Deriv directionAhead  = currentPosition - nextPosition;
        directionAhead.normalize();

        Deriv slidingDirection = directionBeyond - directionAhead;
        rowIterator.setCol(currentIndex, slidingDirection);*/

        //        std::cout << "Integral["<<i<<"] : "<< integral[i] << std::endl;
        rowIterator.setCol(i, integral[i]);
    }
    cIndex++;
    cMatrix.endEdit();
    m_nbLines = cIndex - constraintIndex;
}


template<class DataTypes>
void CosseratActuatorConstraint<DataTypes>::getConstraintViolation(const ConstraintParams* cParams,
                                                                   BaseVector *resV,
                                                                   const BaseVector *Jdx)
{
    if(d_componentState.getValue() != ComponentState::Valid)
        return ;

    SOFA_UNUSED(cParams);

    d_cableLength.setValue(this->getCableLength(m_state->readPositions().ref()));
    //    std::cout << "===> d_cableInitialLength : "<< d_cableInitialLength.getValue() << std::endl;
    //    std::cout << "===> d_cableLength : "<< d_cableLength.getValue() << std::endl;
    //    std::cout << "===> Jdx->element(0) : "<< Jdx->element(0) << std::endl;
    //    std::cout << "&&&&&&&===> d_indices.getValue().size() : "<< d_indices.getValue().size() << std::endl;

    Real dfree = Jdx->element(0) + d_cableInitialLength.getValue() - d_cableLength.getValue();


    for (unsigned i=0;i<d_indices.getValue().size();i++) {
        resV->set(m_constraintIndex.getValue(), dfree);
    }
}


template<class DataTypes>
void CosseratActuatorConstraint<DataTypes>::getConstraintResolution(const ConstraintParams* cParam,
                                                                    std::vector<ConstraintResolution*>& resTab,
                                                                    unsigned int& offset)
{
    if(d_componentState.getValue() != ComponentState::Valid)
        return ;

    SOFA_UNUSED(cParam);

    double imposedValue=d_value.getValue()[d_valueIndex.getValue()];

    double maxDisplacement = std::numeric_limits<double>::max();
    double minDisplacement = -maxDisplacement;

    setUpForceLimits(imposedValue,minDisplacement,maxDisplacement);

    MyCableForceConstraintResolution *cr=  new MyCableForceConstraintResolution(imposedValue, minDisplacement, maxDisplacement);
    resTab[offset++] =cr;
}


template<class DataTypes>
void CosseratActuatorConstraint<DataTypes>::setUpForceLimits(double& imposedValue, double& minDisplacement, double& maxDisplacement)
{
    if(d_maxForce.isSet() && imposedValue>d_maxForce.getValue())
        imposedValue = d_maxForce.getValue();

    if(d_minForce.isSet() && imposedValue<d_minForce.getValue())
        imposedValue = d_minForce.getValue();

    if(d_maxNegativeDisplacement.isSet())
        minDisplacement=-d_maxNegativeDisplacement.getValue();
    if(d_maxPositiveDisplacement.isSet())
        maxDisplacement=d_maxPositiveDisplacement.getValue();
}

} // namespace sofa

