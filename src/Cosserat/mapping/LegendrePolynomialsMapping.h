//
// Created by younes on 17/11/2021.
//
#pragma once
#include <Cosserat/config.h>

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
 * This is a component:
 * https://www.sofa-framework.org/community/doc/programming-with-sofa/create-your-component/
 */

template<class TIn, class TOut>
class LegendrePolynomialsMapping : public core::Mapping<TIn, TOut>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE2(LegendrePolynomialsMapping,TIn,TOut), SOFA_TEMPLATE2(core::Mapping,TIn,TOut));
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
extern template class SOFA_SOFARIGID_API LegendrePolynomialsMapping< sofa::defaulttype::Vec3Types, sofa::defaulttype::Vec3Types >;
#endif

}
