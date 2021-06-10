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

            , d_value(initData(&d_value, "value","Displacement or force to impose.\n"))
            , d_dampingCoefficient(initData(&d_dampingCoefficient, "dampingCoefficient",
                                 "The dumping coefficient to slow down the insertion speed.\n"))
            , d_valueIndex(initData(&d_valueIndex, (unsigned int) 0, "valueIndex",
                                    "Index of the value (in InputValue vector) that we want to impose \n"
                                    "If unspecified the default value is {0}"))
            , d_vectorOfIndices(initData(&d_vectorOfIndices, "vectorOfIndices",
                                         "vector of indices on witch we have to apply the constraint.\n"))
            , d_entryPoint(initData(&d_entryPoint, defaulttype::Vector3(24.95, 0., 0), "entryPoint",
                                    "The predefined entry point, this point can also be determined automatically"
                                    "but not implemented here.\n"))
            , d_direction(initData(&d_direction, defaulttype::Quat(0,0,0,1), "direction",
                                   "direction of insertion, if this is not given "
                                   "the insertion is direct along X.\n"))
    {}

    template<class DataTypes>
    CosseratUnilateralInteractionConstraint<DataTypes>::~CosseratUnilateralInteractionConstraint()
    = default;


    template<class DataTypes>
    void CosseratUnilateralInteractionConstraint<DataTypes>::init()
    {
        Inherit1::init();
        internalInit();
        // update the list of points inside the cochlea
        UpdateList();
    }


    template<class DataTypes>
    void CosseratUnilateralInteractionConstraint<DataTypes>::reinit()
    {
        internalInit();
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
    bool CosseratUnilateralInteractionConstraint<DataTypes>::UpdateList()
    {
        /// @todo:1- Update the list of points beyond the entry point
        /// @todo:2- Build unitest function
        /// @todo:3- use the direction

        printf("1- Inside the UpdateList \n");
        ReadAccessor<Data<VecCoord>> positions = m_state->readPositions();
        auto entryPoint = d_entryPoint.getValue();
//        WriteAccessor<Data<vector<size_t>>> innerPoints = d_vectorOfIndices;
        WriteAccessor<Data<vector<size_t>>> indices = d_vectorOfIndices;
        indices.clear();

        unsigned int index=0;
        for(auto position : positions){
//            std::cout << "1. position "<< position[0]<< " ==> entryPoint : "<< entryPoint[0] << std::endl;
            if (position[0] >= entryPoint[0] ){
                std::cout << index << " ==> position "<< position[0]<< " ==> entryPoint : "<< entryPoint[0] << std::endl;
                indices.push_back(index);
            }
            index++;
        }
        return true;
    }


    template<class DataTypes>
    void CosseratUnilateralInteractionConstraint<DataTypes>::buildConstraintMatrix(const ConstraintParams* cParams, DataMatrixDeriv &cMatrix, unsigned int &cIndex, const DataVecCoord &x)
    {
        if(d_componentState.getValue() != ComponentState::Valid)
            return ;
        printf("1- Inside the buildConstraintMatrix \n");
        UpdateList();

        SOFA_UNUSED(cParams);
        MatrixDeriv& matrix = *cMatrix.beginEdit();
        VecCoord positions = x.getValue();
        auto indices = d_vectorOfIndices.getValue();

        std::cout << " the size of d_vectorOfIndices: "<< indices <<std::endl;

        m_constraintId = cIndex;
        size_t i=0 ;
        for (auto index : indices)
        {
            std::cout << " ===> 1.1 inside the build constraint;  " << positions[index] << std::endl;
            MatrixDerivRowIterator c_it_x = matrix.writeLine(cIndex++);
            c_it_x.addCol(index, Coord(1.,0,0)); // @todo instead of vector3(0,1,0) use the direction of the projection
            MatrixDerivRowIterator c_it_y = matrix.writeLine(cIndex++);
            c_it_y.addCol(index, Coord(0,1.,0)); // @todo instead of vector3(0,1,0) use the direction of the projection
            MatrixDerivRowIterator c_it_z = matrix.writeLine(cIndex++);
            c_it_z.addCol(index, Coord(0,0,1.)); // @todo instead of vector3(0,1,0) use the direction of the projection
        }
        cMatrix.endEdit();
        m_nbLines = cIndex - m_constraintId;
        std::cout << "Number of lines :"<< m_nbLines << std::endl;
    }


    template<class DataTypes>
    void CosseratUnilateralInteractionConstraint<DataTypes>::getConstraintViolation(const ConstraintParams* cParams,
                                                                BaseVector *resV,
                                                                const BaseVector *Jdx)
    {
        if(d_componentState.getValue() != ComponentState::Valid)
            return ;
        printf("1- Inside the getConstraintViolation \n");
        SOFA_UNUSED(cParams);
        ReadAccessor<Data<VecCoord>> positions = m_state->readPositions();
        auto indices = d_vectorOfIndices.getValue();

        unsigned int iter = 0;
        if(Jdx->size()==0){
            for (auto index : indices ){
                //@todo this need to be recompute, use dfree = newPos-oldPos
                Real dfree1 = m_dx; // positions[indices[index]][1];
                Real dfree2 = m_dx; // positions[indices[index]][2];
                Real dfree3 = m_dx; // positions[indices[index]][2];

                resV->set(m_constraintId + 3*iter+0, dfree1);
                resV->set(m_constraintId + 3*iter+1, dfree2);
                resV->set(m_constraintId + 3*iter+2, dfree3);
                iter++;
            }
        }else{
            for (auto index : indices){
                Real dfree1 = Jdx->element(3*iter+0) + m_dx;
                Real dfree2 = Jdx->element(3*iter+1) + m_dx;
                Real dfree3 = Jdx->element(3*iter+2) + m_dx;

                resV->set(m_constraintId + 3*iter+0, dfree1);
                resV->set(m_constraintId + 3*iter+1, dfree2);
                resV->set(m_constraintId + 3*iter+2, dfree2);
                iter++;
            }
        }
        printf("2- Inside the getConstraintViolation \n");
    }

    template<class DataTypes>
    void CosseratUnilateralInteractionConstraint<DataTypes>::getConstraintResolution(const ConstraintParams*,
                                                                 std::vector<core::behavior::ConstraintResolution*>& resTab,
                                                                 unsigned int& offset)
    {
        //        ReadAccessor<Data<VecCoord>> positions = m_state->readPositions();
        std::cout << "1- Inside the getConstraintResolution, offset: "<< offset << std::endl;
        double dampingCoefficient = d_dampingCoefficient.getValue();
        std::cout << "the size resTab : "<< resTab.size() << std::endl;
        for (auto index : d_vectorOfIndices.getValue()){
            resTab[offset+0] = new MyUnilateralConstraintResolutionWithFriction(dampingCoefficient);
            resTab[offset+1] = new MyUnilateralConstraintResolutionWithFriction(dampingCoefficient);
            resTab[offset+2] = new MyUnilateralConstraintResolutionWithFriction(dampingCoefficient);
            offset += 3;
        }
        std::cout << " 2- Inside the getConstraintResolution, offset: "<< offset << std::endl;
    }


    template<class DataTypes>
    void CosseratUnilateralInteractionConstraint<DataTypes>::draw(const VisualParams* /*vparams*/)
    {
        if(d_componentState.getValue() != ComponentState::Valid)
            return ;

    }

} // namespace sofa

