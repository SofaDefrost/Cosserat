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
*                           Plugin SoftRobots v1.0                            *
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

#include <sofa/defaulttype/Vec.h>
#include <sofa/helper/OptionsGroup.h>
#include<sofa/defaulttype/VecTypes.h>
#include <sofa/core/behavior/BaseConstraintSet.h>
#include <sofa/core/behavior/PairInteractionConstraint.h>
#include <SoftRobots/component/constraint/model/CableModel.h>
#include <sofa/core/behavior/ConstraintResolution.h>
namespace sofa::component::constraintset
{

using sofa::helper::vector;
using sofa::defaulttype::Vec;
using sofa::defaulttype::Vec3d;
using sofa::helper::WriteAccessor;
using sofa::core::ConstraintParams;
using sofa::defaulttype::BaseVector;
using sofa::core::visual::VisualParams ;
using sofa::core::behavior::ConstraintResolution;


class MyUnilateralConstraintResolutionWithFriction : public ConstraintResolution
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
    MyUnilateralConstraintResolutionWithFriction(double mu, PreviousForcesContainer* prev = nullptr, bool* active = nullptr)
        : ConstraintResolution(3)
        , _mu(mu)
        , _prev(prev)
        , _active(active)
    {
        printf("Cosserat UnilateralConstraintResolutionWithFriction \n");
    }


    void init(int line, double** w, double * lambda) override
    {
        _W[0]=w[line  ][line  ];
        _W[1]=w[line  ][line+1];
        _W[2]=w[line  ][line+2];
        _W[3]=w[line+1][line+1];
        _W[4]=w[line+1][line+2];
        _W[5]=w[line+2][line+2];

        ////////////////// christian : the following does not work ! /////////
        if(_prev)
        {
            force[line] = _prev->popForce();
            force[line+1] = _prev->popForce();
            force[line+2] = _prev->popForce();
        }
    }

    void resolution(int line, double** w, double* d, double* lambda, double* dfree) override
    {
        double f[2];
        double normFt;

        f[0] = force[line]; f[1] = force[line+1];
        force[line] -= d[line] / _W[0];

        if(force[line] < 0)
        {
            force[line]=0; force[line+1]=0; force[line+2]=0;
            return;
        }

        d[line+1] += _W[1] * (force[line]-f[0]);
        d[line+2] += _W[2] * (force[line]-f[0]);
        force[line+1] -= 2*d[line+1] / (_W[3] +_W[5]) ;
        force[line+2] -= 2*d[line+2] / (_W[3] +_W[5]) ;

        normFt = sqrt(force[line+1]*force[line+1] + force[line+2]*force[line+2]);

        double fN = _mu*force[line];
        if(normFt > fN)
        {
            double factor = fN / normFt;
            force[line+1] *= factor;
            force[line+2] *= factor;
        }
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
    typedef helper::vector<unsigned int> SetIndexArray;



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
    using CableModel<DataTypes>::m_constraintId ;
    using CableModel<DataTypes>::m_state ;
    using CableModel<DataTypes>::d_indices ;
    using CableModel<DataTypes>::d_pullPoint;
    using CableModel<DataTypes>::d_hasPullPoint;
    using CableModel<DataTypes>::d_cableInitialLength;
    using CableModel<DataTypes>::d_cableLength;
    using CableModel<DataTypes>::d_force;

protected:
    //Input data
    Data<helper::vector< Real > >       d_value;
    Data<unsigned int>                  d_valueIndex;
    Data<helper::OptionsGroup>          d_valueType;
    //    Data<SetIndexArray>                 d_indices;
    Data<helper::vector<Coord>>         d_integral;

    void internalInit();
private:
    void setUpForceLimits(double& imposedValue, double& minDisplacement, double& maxDisplacement);
};

// Declares template as extern to avoid the code generation of the template for
// each compilation unit. see: http://www.stroustrup.com/C++11FAQ.html#extern-templates
//extern template class CosseratActuatorConstraint<sofa::defaulttype::Vec3Types>;

} // namespace sofa
