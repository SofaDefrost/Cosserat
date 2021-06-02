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

#include "CosseratUnilateralInteractionConstraint.h"

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
    using sofa::helper::vector;
    using sofa::helper::OptionsGroup;

    template<class DataTypes>
    CosseratUnilateralInteractionConstraint<DataTypes>::CosseratUnilateralInteractionConstraint(MechanicalState* object)
            : Inherit1(object)

            , d_value(initData(&d_value, "value",
                               "Displacement or force to impose.\n"))
            , d_force_dumping(initData(&d_force_dumping, "force_dumping",
                                 "The dumping coefficient to slow down the insertion speed.\n"))
            , d_valueIndex(initData(&d_valueIndex, (unsigned int) 0, "valueIndex",
                                    "Index of the value (in InputValue vector) that we want to impose \n"
                                    "If unspecified the default value is {0}"))
            , d_vectorOfIndices(initData(&d_vectorOfIndices, "vectorOfIndices",
                                         "vector of indices on witch we have to apply the constraint.\n"))
            , d_entryPoint(initData(&d_entryPoint, "entryPoint",
                                    "The predefined entry point, this point can also be determined automatically"
                                    "but not implemented here.\n"))
    {

    }

    template<class DataTypes>
    CosseratUnilateralInteractionConstraint<DataTypes>::~CosseratUnilateralInteractionConstraint()
    {}


    template<class DataTypes>
    void CosseratUnilateralInteractionConstraint<DataTypes>::init()
    {
        Inherit1::init();
        UpdateList();
    }


    template<class DataTypes>
    void CosseratUnilateralInteractionConstraint<DataTypes>::reinit()
    {
        UpdateList();
    }

    template<class DataTypes>
    defaulttype::Vector3 CosseratUnilateralInteractionConstraint<DataTypes>::findEntryPoint()
    {
        /// @todo:1- find the entry automatically, this is need in the case of needle insertion
        /// can also be necessary when the volume is deforming
        /// @todo:2- Build unitest function

        return defaulttype::Vector3(0.0f,0.0f,0.0f);
    }

    template<class DataTypes>
    void CosseratUnilateralInteractionConstraint<DataTypes>::UpdateList()
    {
        /// @todo:1- Update the list of points beyond the entry point
        /// @todo:2- Build unitest function

    }


    template<class DataTypes>
    void CosseratUnilateralInteractionConstraint<DataTypes>::buildConstraintMatrix(const ConstraintParams* cParams, DataMatrixDeriv &cMatrix, unsigned int &cIndex, const DataVecCoord &x)
    {
        if(d_componentState.getValue() != ComponentState::Valid)
            return ;

        SOFA_UNUSED(cParams);
        MatrixDeriv& matrix = *cMatrix.beginEdit();
        VecCoord positions = x.getValue();
        m_constraintId= cIndex;
        size_t sz = d_vectorOfIndices.getValue().size();

        for (auto i=0; i<sz; i++)
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
    void CosseratUnilateralInteractionConstraint<DataTypes>::getConstraintViolation(const ConstraintParams* cParams,
                                                                BaseVector *resV,
                                                                const BaseVector *Jdx)
    {
        if(d_componentState.getValue() != ComponentState::Valid)
            return ;

        SOFA_UNUSED(cParams);
        ReadAccessor<Data<VecCoord>> positions = m_state->readPositions();
        auto indices = d_vectorOfIndices.getValue();

        if(Jdx->size()==0){
            //std::cout << "Size JDX = "<< Jdx->size() << std::endl;

            for (auto index : indices){
                Real dfree1 =  positions[index][1];
                Real dfree2 =  positions[index][2];

                resV->set(m_constraintId + 2*index   , dfree1);
                resV->set(m_constraintId + 2*index +1, dfree2);
            }

        }else{
            for (auto index : indices){
                Real dfree1 = Jdx->element(2*index)   + positions[index][1];
                Real dfree2 = Jdx->element(2*index+1) + positions[index][2];

                resV->set(m_constraintId + 2*index   , dfree1);
                resV->set(m_constraintId + 2*index +1, dfree2);
            }
        }
    }

    template<class DataTypes>
    void CosseratUnilateralInteractionConstraint<DataTypes>::getConstraintResolution(const ConstraintParams*,
                                                                 std::vector<core::behavior::ConstraintResolution*>& resTab,
                                                                 unsigned int& offset)
    {
        ReadAccessor<Data<VecCoord>> positions = m_state->readPositions();
        //    std::cout << "The position size is : " << positions.size()<< std::endl;
        double imposedValue= 1.0;
        for (size_t i = 0; i < positions.size(); i++){

            resTab[offset++] = new BilateralConstraintResolution();
            resTab[offset++] = new BilateralConstraintResolution();
            if(i == positions.size()-1){
                resTab[offset++] = new BilateralConstraintResolution();
            }
        }
        //    std::cout << "The position size is END " << std::endl;
    }


    template<class DataTypes>
    void CosseratUnilateralInteractionConstraint<DataTypes>::draw(const VisualParams* vparams)
    {
        if(d_componentState.getValue() != ComponentState::Valid)
            return ;

    }

} // namespace sofa

