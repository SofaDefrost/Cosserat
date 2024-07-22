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

#include <Cosserat/config.h>
#include <sofa/core/Mapping.h>
#include <sofa/core/trait/DataTypes.h>

namespace sofa::component::mapping
{
    namespace{
        using namespace sofa::core::trait;
    }

/*!
 * \class LegendrePolynomialsMapping
 * @brief Computes and map the length of the beams
 *
 */
template<class TIn, class TOut>
class LegendrePolynomialsMapping : public core::Mapping<TIn, TOut>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE2(LegendrePolynomialsMapping,TIn,TOut),
               SOFA_TEMPLATE2(core::Mapping,TIn,TOut));

    Data<sofa::Index> index; ///< input DOF index
    Data<unsigned int> d_order; ///< input DOF index
    Data<type::vector<double>> d_vectorOfCurvilinearAbscissa;
    Data<type::vector<double>> d_vectorOfContrePointsAbs;

protected:
    LegendrePolynomialsMapping();
    virtual ~LegendrePolynomialsMapping() = default;
    type::vector<type::vector<double>> m_matOfCoeffs;

public:

    /// Compute the local coordinates based on the current output coordinates.

    void init() override;
    void reinit() override;
    double legendrePoly(unsigned int n, const double x);

    void apply(const core::MechanicalParams *mparams,
               DataVecCoord_t<TOut>& out, const DataVecCoord_t<TIn>& in) override;

    void applyJ(const core::MechanicalParams *mparams,
                DataVecDeriv_t<TOut>& out, const DataVecDeriv_t<TIn>& in) override;

    void applyJT(const core::MechanicalParams *mparams,
                 DataVecDeriv_t<TOut>& out, const DataVecDeriv_t<TIn>& in) override;

    void applyJT(const core::ConstraintParams *cparams,
                 DataMatrixDeriv_t<TOut>& out, const DataMatrixDeriv_t<TIn>& in) override;

};

#if !defined(SOFA_COSSERAT_CPP_LegendrePolynomialsMapping)
extern template class SOFA_COSSERAT_API LegendrePolynomialsMapping< sofa::defaulttype::Vec3Types, sofa::defaulttype::Vec3Types >;
#endif

}
