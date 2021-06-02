/******************************************************************************
*               SOFA, Simulation Open-Framework Architecture                  *
*                (c) 2020 INRIA, USTL, UJF, CNRS, MGH                    *
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
*                      Plugin Cosserat v1.0                                   *
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

#include<sofa/defaulttype/VecTypes.h>
#include <sofa/helper/OptionsGroup.h>
#include <sofa/defaulttype/defaulttype.h>
#include <SofaConstraint/UnilateralInteractionConstraint.h>

#include "../../../SoftRobots/src/SoftRobots/component/constraint/model/CableModel.h"
#include "../../../SoftRobots/src/SoftRobots/component/behavior/SoftRobotsConstraint.h"


template class SOFA_SOFACONSTRAINT_API sofa::component::constraintset::UnilateralInteractionConstraint<sofa::defaulttype::Vec3Types>;

//class SOFA_SOFACONSTRAINT_API MyUnilateralConstraintResolutionWithFriction : public UnilateralConstraintResolutionWithFriction
//{
//public:
//    MyUnilateralConstraintResolutionWithFriction(double mu, PreviousForcesContainer* prev = nullptr, bool* active = nullptr)
//            :core::behavior::ConstraintResolution(3)
//            , _mu(mu)
//            , _prev(prev)
//            , _active(active)
//    {
//    }
//
//    void init(int line, double** w, double* force) override
//    {
//        _W[0]=w[line  ][line  ];
//        _W[1]=w[line  ][line+1];
//        _W[2]=w[line  ][line+2];
//        _W[3]=w[line+1][line+1];
//        _W[4]=w[line+1][line+2];
//        _W[5]=w[line+2][line+2];
//
//        ////////////////// christian : the following does not work ! /////////
//        if(_prev)
//        {
//            force[line] = _prev->popForce();
//            force[line+1] = _prev->popForce();
//            force[line+2] = _prev->popForce();
//        }
//
//    }
//    void resolution(int line, double** w, double* d, double* force, double *dFree) override
//    {
//        double f[2];
//        double normFt;
//
//        f[0] = force[line]; f[1] = force[line+1];
//        force[line] -= d[line] / _W[0];
//
//        if(force[line] < 0)
//        {
//            force[line]=0; force[line+1]=0; force[line+2]=0;
//            return;
//        }
//
//        d[line+1] += _W[1] * (force[line]-f[0]);
//        d[line+2] += _W[2] * (force[line]-f[0]);
//        force[line+1] -= 2*d[line+1] / (_W[3] +_W[5]) ;
//        force[line+2] -= 2*d[line+2] / (_W[3] +_W[5]) ;
//
//        normFt = sqrt(force[line+1]*force[line+1] + force[line+2]*force[line+2]);
//
//        double fN = _mu*force[line];
//        if(normFt > fN)
//        {
//            double factor = fN / normFt;
//            force[line+1] *= factor;
//            force[line+2] *= factor;
//        }
//    }
//    void store(int line, double* force, bool /*convergence*/) override
//    {
//        if(_prev)
//        {
//            _prev->pushForce(force[line]);
//            _prev->pushForce(force[line+1]);
//            _prev->pushForce(force[line+2]);
//        }
//
//        if(_active)
//        {
//            *_active = (force[line] != 0);
//            _active = nullptr; // Won't be used in the haptic thread
//        }
//    }
//
//    protected:
//    double _mu;
//    double _W[6];
//    PreviousForcesContainer* _prev;
//    bool* _active; // Will set this after the resolution
//};



namespace sofa::component::constraintset
{

    using sofa::core::behavior::SoftRobotsConstraint ;
    using sofa::core::visual::VisualParams ;
    using sofa::core::objectmodel::Data ;
    using sofa::defaulttype::Vec3dTypes ;
    using sofa::defaulttype::Vec3fTypes ;
    using sofa::defaulttype::BaseVector ;
    using sofa::core::ConstraintParams ;
    using sofa::helper::ReadAccessor ;
    using sofa::core::VecCoordId ;

    using sofa::core::behavior::ConstraintResolution ;


    class SlidingForceConstraintResolution : public ConstraintResolution
    {
    public:
        SlidingForceConstraintResolution(const double& imposedForce, const double& min, const double& max);

        //////////////////// Inherited from ConstraintResolution ////////////////////
        virtual void init(int line, double** w, double *force) override;
        virtual void resolution(int line, double** w, double* d, double* force, double* dfree) override;
        /////////////////////////////////////////////////////////////////////////////

    protected:

        double      m_wActuatorActuator;
        double      m_imposedForce;
        double      m_minDisplacement;
        double      m_maxDisplacement;

    };


/**
 * This class contains common implementation of cable constraints
*/
    template< class DataTypes >
    class CosseratUnilateralInteractionConstraint : public CableModel<DataTypes>
    {
    public:

        SOFA_CLASS(SOFA_TEMPLATE(CosseratUnilateralInteractionConstraint,DataTypes), SOFA_TEMPLATE(CableModel,DataTypes));
        typedef typename DataTypes::VecCoord VecCoord;
        typedef typename DataTypes::VecDeriv VecDeriv;
        typedef typename DataTypes::Coord Coord;
        typedef typename DataTypes::Deriv Deriv;
        typedef typename DataTypes::MatrixDeriv MatrixDeriv;
        typedef typename Coord::value_type Real;
        typedef typename core::behavior::MechanicalState<DataTypes> MechanicalState;

        typedef typename DataTypes::MatrixDeriv::RowIterator MatrixDerivRowIterator;
        typedef Data<VecCoord>		DataVecCoord;
        typedef Data<VecDeriv>		DataVecDeriv;
        typedef Data<MatrixDeriv>    DataMatrixDeriv;
        typedef helper::vector<unsigned int> SetIndexArray;


    public:
        CosseratUnilateralInteractionConstraint(MechanicalState* object = nullptr);

        ~CosseratUnilateralInteractionConstraint() override;

        /*********** Inherited from BaseObject ************/
        void init() override;
        void reinit() override;
        void draw(const VisualParams* vparams) override;

        /*********** Inherited from Actuator ************/
        void buildConstraintMatrix(const ConstraintParams* cParams,
                                   DataMatrixDeriv &cMatrix,
                                   unsigned int &cIndex,
                                   const DataVecCoord &x) override;

        void getConstraintViolation(const ConstraintParams* cParams,
                                    BaseVector *resV,
                                    const BaseVector *Jdx) override;
        void getConstraintResolution(const ConstraintParams*,
                                     std::vector<core::behavior::ConstraintResolution*>& resTab,
                                     unsigned int& offset) override;
    protected:
        //Input data
        Data<helper::vector< Real > >   d_value;
        Data<Real>                      d_force_dumping;
        Data<unsigned int>              d_valueIndex;
        Data<helper::vector<size_t>>    d_vectorOfIndices;
        Data<defaulttype::Vector3>      d_entryPoint;

    protected:
        using SoftRobotsConstraint<DataTypes>::m_state ;
        using CableModel<DataTypes>::d_componentState ;
        using SoftRobotsConstraint<DataTypes>::m_nbLines ;
        using SoftRobotsConstraint<DataTypes>::m_constraintId ;

    public:
        void UpdateList();
        defaulttype::Vector3 findEntryPoint();
    };

// Declares template as extern to avoid the code generation of the template for
// each compilation unit. see: http://www.stroustrup.com/C++11FAQ.html#extern-templates
//extern template class CosseratUnilateralInteractionConstraint<defaulttype::Vec3Types>;

} // namespace sofa


