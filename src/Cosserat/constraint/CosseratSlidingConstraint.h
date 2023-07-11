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
#pragma once
#include <Cosserat/config.h>

#include <sofa/core/behavior/PairInteractionConstraint.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/behavior/ConstraintResolution.h>

#include <iostream>

namespace sofa::component::constraintset
{

using sofa::core::behavior::ConstraintResolution ;
class UnilateralConstraintResolution : public ConstraintResolution {
public:

    /*!
     * \brief UnilateralConstraintResolution Constructor
     * \param m : double, maxForce value
     */
    explicit UnilateralConstraintResolution(double m = std::numeric_limits<double>::max()) :
        sofa::core::behavior::ConstraintResolution(1),
        m_maxForce(m)
    {}

    /*!
     * \brief resolution : updates and keeps the force applied on 'line'
     * in the interval [0, maxForce]
     * \param line
     * \param w
     * \param d
     * \param force
     */
    void resolution (int line, double** w, double* d, double* force, double * /*dFree*/) override {
        force[line] -= d[line] / w[line][line];

        if (force[line]>m_maxForce)
            force[line] = m_maxForce;
        else if (force[line]<0)
            force[line] = 0.0;
    }

    double m_maxForce;
};

class BilateralConstraintResolution : public ConstraintResolution
{
public:
    explicit BilateralConstraintResolution(double* initF=nullptr)
        : ConstraintResolution(1)
        , _f(initF) {}
    void resolution(int line, double** w, double* d, double* force, double *dfree) override
    {
        SOFA_UNUSED(dfree);
        force[line] -= d[line] / w[line][line];
    }

    void init(int line, double** /*w*/, double* force) override
    {
        if(_f) { force[line] = *_f; }
    }

    void initForce(int line, double* force) override
    {
        if(_f) { force[line] = *_f; }
    }

    void store(int line, double* force, bool /*convergence*/) override
    {
        if(_f) *_f = force[line];
    }

protected:
    double* _f;
};


using sofa::core::ConstraintParams;

template<class DataTypes>
class CosseratSlidingConstraint : public core::behavior::PairInteractionConstraint<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(CosseratSlidingConstraint,DataTypes), SOFA_TEMPLATE(core::behavior::PairInteractionConstraint,DataTypes));

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::MatrixDeriv MatrixDeriv;
    typedef typename DataTypes::MatrixDeriv::RowIterator MatrixDerivRowIterator;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename Coord::value_type Real;
    typedef typename core::behavior::MechanicalState<DataTypes> MechanicalState;
    typedef typename core::behavior::PairInteractionConstraint<DataTypes> Inherit;

    typedef core::objectmodel::Data<VecCoord>		DataVecCoord;
    typedef core::objectmodel::Data<VecDeriv>		DataVecDeriv;
    typedef core::objectmodel::Data<MatrixDeriv>    DataMatrixDeriv;

protected:

    Data<int> d_m1; ///< index of the 3D point on the first model
    Data<int> d_m2a; ///< index of one end of the sliding axis
    Data<int> d_m2b; ///< index of the other end of the sliding axis
    Data<Deriv> d_force; ///< interaction force

    CosseratSlidingConstraint();
    explicit CosseratSlidingConstraint(MechanicalState* object);
    CosseratSlidingConstraint(MechanicalState* object1, MechanicalState* object2);

    virtual ~CosseratSlidingConstraint()=default;

public:
    void init() override;

    void computeProximity(const DataVecCoord &x1, const DataVecCoord &x2);

    void buildConstraintMatrix(const core::ConstraintParams* cParams, DataMatrixDeriv &c1, DataMatrixDeriv &c2, unsigned int &cIndex
                               , const DataVecCoord &x1, const DataVecCoord &x2) override;

    void getConstraintViolation(const core::ConstraintParams* cParams, sofa::linearalgebra::BaseVector *v, const DataVecCoord &x1, const DataVecCoord &x2
                                , const DataVecDeriv &v1, const DataVecDeriv &v2) override;

    void getConstraintResolution(const core::ConstraintParams*,
                                 std::vector<core::behavior::ConstraintResolution*>& resTab,
                                 unsigned int& offset) override;
    void storeLambda(const ConstraintParams* cParams, sofa::core::MultiVecDerivId res, const sofa::linearalgebra::BaseVector* lambda) override;

    void draw(const core::visual::VisualParams* vparams) override;

    void drawLinesBetweenPoints(const core::visual::VisualParams* vparams);

private:
    // storage of force

    typedef struct {
        double fact;
        int p1, p2;
        int eid, cid;
        Real alpha, beta;
        Coord P, Q;
        Deriv dirProj, dirAxe, dirOrtho;
        Real r,r2, Q1Q2;
        Real dist; //violation
        Real thirdConstraint;
    } Constraint;

    unsigned int m_step;
    type::vector<Constraint> m_constraints;

};

}
