//
// Created by younes on 17/11/2021.
//

#pragma once

#include "LegendrePolynomialsMapping.h"
#include <SofaRigid/RigidMapping.h>
#include <sofa/core/visual/VisualParams.h>

#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/State.h>

#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>

#include <sofa/helper/io/XspLoader.h>
#include <sofa/helper/io/SphereLoader.h>
#include <sofa/helper/io/Mesh.h>
#include <sofa/helper/decompose.h>
#include <sofa/core/MechanicalParams.h>

namespace sofa::component::mapping {

    template <class TIn, class TOut>
    LegendrePolynomialsMapping<TIn, TOut>::LegendrePolynomialsMapping()
        : Inherit()
        , index(initData(&index, (unsigned)0, "index", "input DOF index"))
    {}

    template <class TIn, class TOut>
    void LegendrePolynomialsMapping<TIn, TOut>::reinit()
    {
        const VecCoord& xTo =this->toModel->read(core::ConstVecCoordId::position())->getValue();
        unsigned int i = 0;
    }

    template <class TIn, class TOut>
    void LegendrePolynomialsMapping<TIn, TOut>::init()
    {
        this->reinit();
        this->Inherit::init();
    }

    template <class TIn, class TOut>
    void LegendrePolynomialsMapping<TIn, TOut>::apply(const core::MechanicalParams * /*mparams*/, Data<VecCoord>& dOut, const Data<InVecCoord>& dIn)
    {
        helper::WriteOnlyAccessor< Data<VecCoord> > out = dOut;
        helper::ReadAccessor< Data<InVecCoord> > in = dIn;
    }

    template <class TIn, class TOut>
    void LegendrePolynomialsMapping<TIn, TOut>::applyJ(const core::MechanicalParams * /*mparams*/, Data<VecDeriv>& dOut, const Data<InVecDeriv>& dIn)
    {
        helper::WriteOnlyAccessor< Data<VecDeriv> > out = dOut;
        helper::ReadAccessor< Data<InVecDeriv> > in = dIn;

        for(sofa::Index i=0 ; i<out.size() ; ++i)
        {
        }
    }

    template <class TIn, class TOut>
    void LegendrePolynomialsMapping<TIn, TOut>::applyJT(const core::MechanicalParams * /*mparams*/, Data<InVecDeriv>& dOut, const Data<VecDeriv>& dIn)
    {
        helper::WriteAccessor< Data<InVecDeriv> > out = dOut;
        helper::ReadAccessor< Data<VecDeriv> > in = dIn;

        for(sofa::Index i=0 ; i<in.size() ; ++i)
        {

        }
    }

// RigidMapping::applyJT( InMatrixDeriv& out, const OutMatrixDeriv& in ) //
// this function propagate the constraint through the rigid mapping :
// if one constraint along (vector n) with a value (v) is applied on the childModel (like collision model)
// then this constraint is transformed by (Jt.n) with value (v) for the rigid model
// There is a specificity of this propagateConstraint: we have to find the application point on the childModel
// in order to compute the right constaint on the rigidModel.
template <class TIn, class TOut>
void LegendrePolynomialsMapping<TIn, TOut>::applyJT(const core::ConstraintParams * /*cparams*/, Data<InMatrixDeriv>& dOut, const Data<OutMatrixDeriv>& dIn)
{
        InMatrixDeriv& out = *dOut.beginEdit();
        const OutMatrixDeriv& in = dIn.getValue();

        dmsg_info() << "J on mapped DOFs == " << in << msgendl
                    << "J on input  DOFs == " << out ;

        const unsigned int numDofs = this->getFromModel()->getSize();

        // TODO the implementation on the new data structure could maybe be optimized
        typename Out::MatrixDeriv::RowConstIterator rowItEnd = in.end();

        for (typename Out::MatrixDeriv::RowConstIterator rowIt = in.begin(); rowIt != rowItEnd; ++rowIt)
        {
            for (unsigned int ito = 0; ito < numDofs; ito++)
            {

            }
        }

        dmsg_info() << "new J on input  DOFs = " << out ;

        dOut.endEdit();
}


template <class TIn, class TOut>
void LegendrePolynomialsMapping<TIn, TOut>::draw(const core::visual::VisualParams* /*vparams*/)
{
    // draw cable

}

}