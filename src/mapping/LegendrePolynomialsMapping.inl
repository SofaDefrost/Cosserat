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
        , d_order(initData(&d_order, (unsigned)3, "order", "The order of Legendre polynomials"))
        , d_vectorOfCurvilinearAbscissa(initData(&d_vectorOfCurvilinearAbscissa, "curvAbscissa", "Vector of curvilinear Abscissa element of [0, 1]"))
        , d_vectorOfContrePointsAbs(initData(&d_vectorOfContrePointsAbs, "controlPointsAbs", "Vector of curvilinear Abscissa of control points element of [0, 1]"))
    {}

    template <class TIn, class TOut>
    double LegendrePolynomialsMapping<TIn, TOut>::legendrePoly(unsigned int n, const double x) {
        if (n == 0)
            return 1. ;
        else if (n == 1 )
            return x ;
        else
            return (((2 * x ) - 1) * x * legendrePoly(n - 1, x) - (n - 1) * legendrePoly(n - 2, x)) / double(n);
    }

    template <class TIn, class TOut>
    void LegendrePolynomialsMapping<TIn, TOut>::reinit() {
        m_matOfCoeffs.clear();
        auto curvAbs = d_vectorOfCurvilinearAbscissa.getValue();
        auto  sz = curvAbs.size();
        std::cout << " curvAbs :" << curvAbs << std::endl;
        for (unsigned int i = 1; i < sz; i++){
            type::vector<double> coeffsOf_i;
            coeffsOf_i.clear();
            for (unsigned int order = 0; order < d_order.getValue(); order++)
                coeffsOf_i.push_back(legendrePoly(order, curvAbs[i]));

            // std::cout << " = = = >coeffsOf_i: " << coeffsOf_i << std::endl;
            m_matOfCoeffs.push_back(coeffsOf_i);
        }
    }


    template <class TIn, class TOut>
    void LegendrePolynomialsMapping<TIn, TOut>::init()
    {
//        this->Inherit::init();
        //Compute the coefficients for each curv_abs at all orders of the polynomials
        reinit();
    }


    template <class TIn, class TOut>
    void LegendrePolynomialsMapping<TIn, TOut>::apply(const core::MechanicalParams * /*mparams*/, Data<VecCoord>& dOut, const Data<InVecCoord>& dIn)
    {
        helper::ReadAccessor< Data<InVecCoord> > in = dIn;
        helper::WriteOnlyAccessor< Data<VecCoord> > out = dOut;
        const auto sz = d_vectorOfCurvilinearAbscissa.getValue().size();
        out.resize(sz-1);

        for (unsigned int i = 0; i < sz-1; i++){
            type::Vector3 Xi ;
            for (unsigned int j = 0; j < in.size(); j++)
                Xi += m_matOfCoeffs[i][j] * in[j];

            // std::cout << "= = = > Xi: " << Xi << std::endl;
            out[i] = Xi;
        }
    }

    template <class TIn, class TOut>
    void LegendrePolynomialsMapping<TIn, TOut>::applyJ(const core::MechanicalParams * /*mparams*/, Data<VecDeriv>& dOut, const Data<InVecDeriv>& dIn)
    {
        helper::WriteOnlyAccessor< Data<VecDeriv> > velOut = dOut;
        helper::ReadAccessor< Data<InVecDeriv> > velIn = dIn;

        helper::WriteOnlyAccessor< Data<VecDeriv> > out = dOut;
        helper::ReadAccessor< Data<InVecDeriv> > in = dIn;

        const auto sz = d_vectorOfCurvilinearAbscissa.getValue().size();
        out.resize(sz-1);
        for(sofa::Index i=0 ; i<sz-1 ; ++i)
        {
            Vector3 vel ;
            for (unsigned int j = 0; j < velIn.size(); j++)
                vel += m_matOfCoeffs[i][j] * velIn[j];

            // std::cout << " vel :" << vel << std::endl;
            velOut[i] = vel;
        }
    }

    template <class TIn, class TOut>
    void LegendrePolynomialsMapping<TIn, TOut>::applyJT(const core::MechanicalParams * /*mparams*/, Data<InVecDeriv>& dOut, const Data<VecDeriv>& dIn)
    {
        helper::WriteAccessor< Data<InVecDeriv> > out = dOut;
        helper::ReadAccessor< Data<VecDeriv> > in = dIn;

//        std::cout << "J on mapped DOFs == " << in[0] << "; size :"<< in.size()
//                  << "\nJ on input  DOFs == " << out[0] << "; size :"<< out.size()  << std::endl;

        for (unsigned int cI = 0; cI < out.size(); cI++){
            for(sofa::Index i=0 ; i<in.size() ; ++i){
                //std::cout << " cI:" << cI << " i:"<< i <<" m_matOfCoeffs[i][cI] : "<< (m_matOfCoeffs[i][cI]) * in[i]<< std::endl;
                out[cI] += (1./m_matOfCoeffs[i][cI]) * in[i];
            }
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

//         std::cout << "J on mapped DOFs == " << in << "; size :"<< in.size()
//            << "\nJ on input  DOFs == " << out << "; size :"<< out.size()  << std::endl;

        const unsigned int numDofs = this->getFromModel()->getSize();

        std::cout << " - - > numDofs : "<< numDofs << std::endl;

        // TODO the implementation on the new data structure could maybe be optimized
        typename Out::MatrixDeriv::RowConstIterator rowItEnd = in.end();

        for (typename Out::MatrixDeriv::RowConstIterator rowIt = in.begin(); rowIt != rowItEnd; ++rowIt)
        {
            typename OutMatrixDeriv::ColConstIterator colIt = rowIt.begin();
            typename OutMatrixDeriv::ColConstIterator colItEnd = rowIt.end();

            std::cout << "rowIndexe : "<<rowIt.index() << std::endl;

            if (colIt == colItEnd)
                continue;

            typename InMatrixDeriv::RowIterator o = out.writeLine(rowIt.index()); // we store the constraint number
            int indexBeam =  0;
            while (colIt != colItEnd) {
                int childIndex = colIt.index();
                const OutDeriv f_It = colIt.val();
                /* @todo get parentIndex: here we are contributing as parent for everything right ?
                    todo get parentIndex: If not, one have to find a way to get this (help @Christian) */
                //int indexBeam =  0; //m_indicesVectors[childIndex];

                o.addCol(indexBeam, f_It);
                indexBeam++;
                indexBeam = (indexBeam == 3) ? 0 : indexBeam;
                colIt++;
            }
//            for (unsigned int ito = 0; ito < numDofs; ito++)
//            {
//                OutDeriv f;
//                bool needToInsert = false;
//
//                for (typename Out::MatrixDeriv::ColConstIterator colIt = rowIt.begin(); colIt != rowIt.end(); ++colIt) {
//                    // this used to get the parentIndex, the node on which we are computing the force
//                    auto childIndex = colIt.index();
//                    std::cout << "childIndex :"<< childIndex<<std::endl;
//                    /* @todo get parentIndex: here we are contributing as parent for everything right ?
//                    todo get parentIndex: If not, one have to find a way to get this (help @Christian) */
//
//                    const OutDeriv fIt = colIt.val();
//                    f += fIt;
//                }
//                // std::cout << "The out force :"<< f << std::endl;
//                o.addCol(ito,f);
//            }
        }

        std::cout << "new J on input  DOFs = \n" << out.size() << std::endl;

        dOut.endEdit();
}


template <class TIn, class TOut>
void LegendrePolynomialsMapping<TIn, TOut>::draw(const core::visual::VisualParams* /*vparams*/)
{
    // draw cable

}

}