/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture                          *
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
*                           Plugin Cosserat    v1.0                           *
*				                                                              *
* This plugin is also distributed under the GNU LGPL (Lesser General          *
* Public License) license with the same conditions than SOFA.                 *
*                                                                             *
* Contributors: Defrost team  (INRIA, University of Lille, CNRS,              *
*               Ecole Centrale de Lille)                                      *
*                                                                             *
* Contact information: https://www.inria.fr/fr/defrost                        *
*                                                                             *
******************************************************************************/
#pragma once

#include "LegendrePolynomialsMapping.h"
#include <boost/math/special_functions/legendre.hpp>

namespace sofa::component::mapping
{
    template <class TIn, class TOut>
    LegendrePolynomialsMapping<TIn, TOut>::LegendrePolynomialsMapping()
        : Inherit1()
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
        // std::cout << " curvAbs :" << curvAbs << std::endl;
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
        Inherit1::init();

        //Compute the coefficients for each curv_abs at all orders of the polynomials
        reinit();
    }


    template <class TIn, class TOut>
    void LegendrePolynomialsMapping<TIn, TOut>::apply(const core::MechanicalParams * /*mparams*/,
                                                      DataVecCoord_t<TOut>& dOut, const DataVecCoord_t<TIn>& dIn)
    {
        auto in = sofa::helper::getReadAccessor(dIn);
        auto out = sofa::helper::getWriteAccessor(dOut);
        const auto sz = d_vectorOfCurvilinearAbscissa.getValue().size();
        out.resize(sz-1);

        for (unsigned int i = 0; i < sz-1; i++){
            type::Vec3 Xi ;
            for (unsigned int j = 0; j < in.size(); j++)
                Xi += m_matOfCoeffs[i][j] * in[j];
            out[i] = Xi;
        }
    }

    template <class TIn, class TOut>
    void LegendrePolynomialsMapping<TIn, TOut>::applyJ(const core::MechanicalParams * /*mparams*/,
                                                       DataVecDeriv_t<TOut>& dOut, const DataVecDeriv_t<TIn>& dIn)
    {
        auto velOut = sofa::helper::getWriteAccessor(dOut);
        auto velIn = sofa::helper::getReadAccessor(dIn);
        auto out = sofa::helper::getWriteOnlyAccessor(dOut);

        const auto sz = d_vectorOfCurvilinearAbscissa.getValue().size();
        out.resize(sz-1);
        for(sofa::Index i=0 ; i<sz-1 ; ++i)
        {
            type::Vec3 vel ;
            for (unsigned int j = 0; j < velIn.size(); j++)
                vel += m_matOfCoeffs[i][j] * velIn[j];
            velOut[i] = vel;
        }
    }

    template <class TIn, class TOut>
    void LegendrePolynomialsMapping<TIn, TOut>::applyJT(const core::MechanicalParams * /*mparams*/,
                                                        DataVecDeriv_t<TOut>& dOut, const DataVecDeriv_t<TIn>& dIn)
    {
        auto out = sofa::helper::getWriteAccessor(dOut);
        auto in = sofa::helper::getReadAccessor(dIn);

        const unsigned int numDofs = this->getFromModel()->getSize();
        out.resize(numDofs);
        for (unsigned int cI = 0; cI < out.size(); cI++){
            for(sofa::Index i=0 ; i<in.size() ; ++i){
                out[cI] += m_matOfCoeffs[i][cI] * in[i];
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
void LegendrePolynomialsMapping<TIn, TOut>::applyJT(const core::ConstraintParams * /*cparams*/, DataMatrixDeriv_t<TOut>& dOut, const DataMatrixDeriv_t<TIn>& dIn)
{
        auto out = sofa::helper::getWriteAccessor(dOut);
        auto in = sofa::helper::getReadAccessor(dIn);

        const unsigned int numDofs = this->getFromModel()->getSize();

        type::vector<Deriv_t<TIn>> tabF;
        tabF.resize(numDofs);

        auto rowIt = in->begin();
        auto rowEnd = in->end();

        for (;rowIt!=rowEnd;++rowIt)
        {
            auto colIt = rowIt.begin();
            auto colItEnd = rowIt.end();

            if (colIt == colItEnd)
                continue;

            auto o = out->writeLine(rowIt.index()); // we store the constraint number

            while (colIt != colItEnd)
            {
                int childIndex = colIt.index();
                auto f_It = colIt.val();
                for (unsigned int order = 0; order < numDofs; order++){
                    Deriv_t<TIn> f;
                    f = m_matOfCoeffs[childIndex][order] * f_It;
                    tabF[order] += f;
                    o.addCol(order, f);
                }
                colIt++;
            }
        }
}

}
