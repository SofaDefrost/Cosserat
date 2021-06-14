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
#include <sofa/defaulttype/config.h>
#include <SofaConstraint/UnilateralInteractionConstraint.h>

#include "../../../SoftRobots/src/SoftRobots/component/constraint/model/CableModel.h"
#include "../../../SoftRobots/src/SoftRobots/component/behavior/SoftRobotsConstraint.h"

template class SOFA_SOFACONSTRAINT_API sofa::component::constraintset::UnilateralInteractionConstraint<sofa::defaulttype::Vec3Types>;





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
    using sofa::helper::WriteAccessor;
    using sofa::core::VecCoordId ;

    using sofa::core::behavior::ConstraintResolution ;

    class MyPreviousForcesContainer
    {
    public:
        MyPreviousForcesContainer() : resetFlag(true) {}
        double popForce()
        {
            resetFlag = true;
            if(forces.empty()) return 0;
            double f = forces.front();
            forces.pop_front();
            return f;
        }

        void pushForce(double f)
        {
            if(resetFlag)
            {
                forces.clear();
                resetFlag = false;
            }

            forces.push_back(f);
        }

    protected:
        std::deque<double> forces;
        bool resetFlag; // We delete all forces that were not read
    };

    class SOFA_SOFACONSTRAINT_API MyUnilateralConstraintResolutionWithFriction : public ConstraintResolution
    {
    public:
        MyUnilateralConstraintResolutionWithFriction(double dampingFactor, MyPreviousForcesContainer* prev = nullptr, bool* active = nullptr)
                :sofa::core::behavior::ConstraintResolution(3)
                , _dampingFactor(dampingFactor)
                , _prev(prev)
                , _active(active)
        {
        }

        void init(int line, double** w, double* force) override
        {
            // for methode 1
            sofa::defaulttype::Mat<3,3,double> temp;
            temp[0][0] = w[line][line];
            temp[0][1] = w[line][line+1];
            temp[0][2] = w[line][line+2];
            temp[1][0] = w[line+1][line];
            temp[1][1] = w[line+1][line+1];
            temp[1][2] = w[line+1][line+2];
            temp[2][0] = w[line+2][line];
            temp[2][1] = w[line+2][line+1];
            temp[2][2] = w[line+2][line+2];

            sofa::type::invertMatrix(invW, temp);
            // for method 2
            _W[0]=w[line  ][line  ];
            _W[1]=w[line  ][line+1];
            _W[2]=w[line  ][line+2];
            _W[3]=w[line+1][line+1];
            _W[4]=w[line+1][line+2];
            _W[5]=w[line+2][line+2];
        }
        void resolution(int line, double** /*w*/, double* d, double* force, double *dFree) override
        {
            double f[2];
            double normFt;
            defaulttype::Vector3 _f ;
            for(int i=0; i<3; i++)
            {
                for(int j=0; j<3; j++)
                    _f[i] -= d[line+j] * invW[i][j];
            }

            //printf("d1 %f; d2 %f  d3 %f \n", d[line], d[line+1], d[line+2]);
            //std::cout<< " 1==> force :"<< _f << std::endl;
            //std::cout<< " 1==> Force :"<< force[line+0] << " "<< force[line+1] << " " <<force[line+2] << std::endl;

            force[line] -= _f[0]*_dampingFactor;
            force[line+1] -= _f[1]*_dampingFactor;
            force[line+2] -= _f[2]* _dampingFactor;
            //std::cout<< " 2==> Force :"<< force[line+0] << " "<< force[line+1] << " " <<force[line+2] << std::endl;
            //printf("===================================\n");
        }
        void store(int line, double* force, bool /*convergence*/) override
        {
            if(_prev)
            {
                _prev->pushForce(force[line+0]);
                _prev->pushForce(force[line+1]);
                _prev->pushForce(force[line+2]);
            }

            if(_active)
            {
                *_active = (force[line] != 0);
                _active = nullptr; // Won't be used in the haptic thread
            }
        }

    protected:
        double _dampingFactor;
        double _W[6];
        sofa::defaulttype::Mat<3,3,double> invW;
        MyPreviousForcesContainer* _prev;
        bool* _active; // Will set this after the resolution
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
        bool UpdateList();



        void storeLambda(const ConstraintParams* cParams,
                                                core::MultiVecDerivId res,
                                                const sofa::defaulttype::BaseVector* /*lambda*/) override
        {
            SOFA_UNUSED(res);
            SOFA_UNUSED(cParams);
        }

    protected:
        //Input data
        Data<helper::vector< Real > >   d_value;
        Data<Real>                      d_dampingCoefficient;
        Data<unsigned int>              d_valueIndex;
        Data<helper::vector<size_t>>    d_vectorOfIndices;
        Data<defaulttype::Vector3>      d_entryPoint;
        Data<helper::vector<helper::vector<double>> >         d_direction;
    protected:
        using SoftRobotsConstraint<DataTypes>::m_state ;
        using CableModel<DataTypes>::d_componentState ;
        using SoftRobotsConstraint<DataTypes>::m_nbLines ;

        using CableModel<DataTypes>::d_maxDispVariation ;

        using CableModel<DataTypes>::d_maxPositiveDisplacement ;
        using CableModel<DataTypes>::d_maxNegativeDisplacement ;
        using CableModel<DataTypes>::d_maxForce ;
        using CableModel<DataTypes>::d_minForce ;
        using CableModel<DataTypes>::d_displacement ;

        using SoftRobotsConstraint<DataTypes>::m_constraintId ;

        void internalInit()
        {
            if(d_value.getValue().size()==0)
            {
                WriteAccessor<Data<helper::vector<Real>>> value = d_value;
                value.resize(1,0.);
            }
            // check for errors in the initialization
            if(d_value.getValue().size()<d_valueIndex.getValue())
            {
                msg_warning() << "Bad size for data value (size="<< d_value.getValue().size()<<"), or wrong value for data valueIndex (valueIndex="<<d_valueIndex.getValue()<<"). Set default valueIndex=0.";
                d_valueIndex.setValue(0);
            }
        }
    public:
        defaulttype::Vector3 findEntryPoint();

    protected:
         Real m_dx = 0.025;



    };

// Declares template as extern to avoid the code generation of the template for
// each compilation unit. see: http://www.stroustrup.com/C++11FAQ.html#extern-templates
//extern template class CosseratUnilateralInteractionConstraint<defaulttype::Vec3Types>;

} // namespace sofa

