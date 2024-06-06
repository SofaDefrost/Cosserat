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

#include <sofa/type/Vec.h>
#include <sofa/helper/OptionsGroup.h>
#include<sofa/defaulttype/VecTypes.h>
#include <sofa/core/behavior/BaseConstraintSet.h>
#include <sofa/core/behavior/PairInteractionConstraint.h>
#include <SoftRobots/component/constraint/model/CableModel.h>
#include <sofa/core/behavior/ConstraintResolution.h>
namespace sofa::component::constraintset
{

using sofa::type::vector;
using sofa::type::Vec;
using sofa::type::Vec3d;
using sofa::helper::WriteAccessor;
using sofa::core::ConstraintParams;
using sofa::linearalgebra::BaseVector;
using sofa::core::visual::VisualParams ;
using sofa::core::behavior::ConstraintResolution;


class MyCableForceConstraintResolution : public ConstraintResolution
{
public:
    //CableForceConstraintResolution(const double& imposedForce, const double& min, const double& max);

    //////////////////// Inherited from ConstraintResolution ////////////////////
    //    void init(int line, double** w, double *force) override;
    //    void resolution(int line, double** w, double* d, double* force, double* dfree) override;
    /////////////////////////////////////////////////////////////////////////////

protected:
    double      m_wActuatorActuator;
    double      m_imposedForce;
    double      m_minDisplacement;
    double      m_maxDisplacement;

public:
    //--------------- Force constraint -------------
    MyCableForceConstraintResolution(const double &imposedForce, const double& min, const double& max)
        : ConstraintResolution(1)
        , m_imposedForce(imposedForce)
        , m_minDisplacement(min)
        , m_maxDisplacement(max)
    {
        //        printf("The constructor is called \n");
    }


    void init(int line, double** w, double * lambda) override
    {
        SOFA_UNUSED(lambda);
        m_wActuatorActuator = w[line][line];
    }

    void resolution(int line, double** w, double* d, double* lambda, double* dfree) override
    {
        //        std::cout << "&&&&&&&&&&&&&&&&&&&&&===========================> W1 : "<< w[line][line] << std::endl;
        SOFA_UNUSED(dfree);
        SOFA_UNUSED(w);

        double displacement = m_wActuatorActuator*m_imposedForce + d[line];

        if (displacement<m_minDisplacement)
        {
            displacement=m_minDisplacement;
            lambda[line] -= (d[line]-displacement) / m_wActuatorActuator;
        }
        else if (displacement>m_maxDisplacement)
        {
            displacement=m_maxDisplacement;
            lambda[line] -= (d[line]-displacement) / m_wActuatorActuator;
        }
        else
            lambda[line] = m_imposedForce;
    }
};



/**
 * This component simulates a force exerted by a cable to solve an effector constraint.
 * Description can be found at:
 *
*/
template< class DataTypes >
class CosseratActuatorConstraint : public CableModel<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(CosseratActuatorConstraint,DataTypes), SOFA_TEMPLATE(CableModel,DataTypes));

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;

    typedef typename DataTypes::MatrixDeriv MatrixDeriv;
    typedef typename core::behavior::MechanicalState<DataTypes> MechanicalState;

    typedef Data<VecCoord>		DataVecCoord;
    typedef Data<VecDeriv>		DataVecDeriv;
    typedef Data<MatrixDeriv>    DataMatrixDeriv;
    typedef typename DataTypes::MatrixDeriv::RowIterator MatrixDerivRowIterator;
    typedef type::vector<unsigned int> SetIndexArray;



public:
    CosseratActuatorConstraint(MechanicalState* object = nullptr);
    ~CosseratActuatorConstraint() override;

    /////////////// Inherited from BaseObject //////////////////////
    void init() override;
    void reinit() override;
    ///////////////////////////////////////////////////////////////


    /////////////////// Inherited from BaseConstraint ///////////////
    void getConstraintResolution(const core::ConstraintParams *cParam,
                                 std::vector<ConstraintResolution*>& resTab,
                                 unsigned int& offset) override;
    ////////////////////////// Inherited from CableModel //////////////////////
    void buildConstraintMatrix(const ConstraintParams* cParams,
                               DataMatrixDeriv &cMatrix,
                               unsigned int &cIndex,
                               const DataVecCoord &x) override;

    void getConstraintViolation(const ConstraintParams* cParams,
                                BaseVector *resV,
                                const BaseVector *Jdx) override;
    ////////////////////////////////////////////////////////////////

    ////////////////////////// Inherited attributes ////////////////////////////
    /// https://gcc.gnu.org/onlinedocs/gcc/Name-lookup.html
    /// Bring inherited attributes and function in the current lookup context.
    /// otherwise any access to the base::attribute would require
    /// the "this->" approach.
    using CableModel<DataTypes>::d_maxDispVariation ;
    using CableModel<DataTypes>::d_maxPositiveDisplacement ;
    using CableModel<DataTypes>::d_maxNegativeDisplacement ;
    using CableModel<DataTypes>::d_maxForce ;
    using CableModel<DataTypes>::d_minForce ;
    using CableModel<DataTypes>::d_displacement ;
    using CableModel<DataTypes>::d_componentState ;
    using CableModel<DataTypes>::m_nbLines ;
    using CableModel<DataTypes>::m_constraintIndex ;
    using CableModel<DataTypes>::m_state ;
    using CableModel<DataTypes>::d_indices ;
    using CableModel<DataTypes>::d_pullPoint;
    using CableModel<DataTypes>::d_hasPullPoint;
    using CableModel<DataTypes>::d_cableInitialLength;
    using CableModel<DataTypes>::d_cableLength;
    using CableModel<DataTypes>::d_force;

protected:
    //Input data
    Data<type::vector< Real > >       d_value;
    Data<unsigned int>                  d_valueIndex;
    Data<helper::OptionsGroup>          d_valueType;
    //    Data<SetIndexArray>                 d_indices;
    Data<type::vector<Coord>>         d_integral;

    void internalInit();
private:
    void setUpForceLimits(double& imposedValue, double& minDisplacement, double& maxDisplacement);
};

// Declares template as extern to avoid the code generation of the template for
// each compilation unit. see: http://www.stroustrup.com/C++11FAQ.html#extern-templates
//extern template class CosseratActuatorConstraint<sofa::defaulttype::Vec3Types>;

} // namespace sofa
