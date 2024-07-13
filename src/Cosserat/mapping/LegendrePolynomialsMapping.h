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
#include <Cosserat/mapping/BaseCosseratMapping.h>

#include <sofa/core/BaseMapping.h>
#include <sofa/core/config.h>
#include <sofa/core/Mapping.h>
#include <sofa/defaulttype/SolidTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/Mapping.h>
#include <sofa/helper/ColorMap.h>

#include <boost/math/special_functions/legendre.hpp>


namespace sofa::component::mapping {
    using sofa::defaulttype::SolidTypes;
    using sofa::core::objectmodel::BaseContext;
    using sofa::type::Matrix3;
    using sofa::type::Matrix4;
    using type::Vec3;
    using type::Vec6;
    using std::get;

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

    typedef core::Mapping<TIn, TOut> Inherit;
    typedef TIn In;
    typedef TOut Out;
    typedef Out OutDataTypes;
    typedef typename Out::VecCoord VecCoord;
    typedef typename Out::VecDeriv VecDeriv;
    typedef typename Out::Coord Coord;
    typedef typename Out::Deriv OutDeriv;
    typedef typename Out::MatrixDeriv OutMatrixDeriv;
    typedef typename In::Real InReal;
    typedef typename In::Deriv InDeriv;
    typedef typename In::VecCoord InVecCoord;
    typedef typename In::VecDeriv InVecDeriv;
    typedef typename In::MatrixDeriv InMatrixDeriv;
    typedef typename Coord::value_type Real;

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

    void apply(const core::MechanicalParams *mparams, Data<VecCoord>& out, const Data<InVecCoord>& in) override;

    void applyJ(const core::MechanicalParams *mparams, Data<VecDeriv>& out, const Data<InVecDeriv>& in) override;

    void applyJT(const core::MechanicalParams *mparams, Data<InVecDeriv>& out, const Data<VecDeriv>& in) override;

    void applyJT(const core::ConstraintParams *cparams, Data<InMatrixDeriv>& out, const Data<OutMatrixDeriv>& in) override;

    void draw(const core::visual::VisualParams* vparams) override;

};

#if !defined(SOFA_COSSERAT_CPP_LegendrePolynomialsMapping)
extern template class SOFA_COSSERAT_API LegendrePolynomialsMapping< sofa::defaulttype::Vec3Types, sofa::defaulttype::Vec3Types >;
#endif

}
