//
// Created by younes on 17/11/2021.
//

#pragma once
#include <sofa/core/BaseMapping.h>
#include <sofa/core/config.h>
#include <sofa/core/Mapping.h>
#include "../initCosserat.h"
#include <sofa/defaulttype/SolidTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <SofaOpenglVisual/OglColorMap.h>


namespace sofa::component::mapping {
    using sofa::defaulttype::SolidTypes;
    using sofa::core::objectmodel::BaseContext;
    using sofa::type::Matrix3;
    using sofa::type::Matrix4;
    using sofa::defaulttype::Vector3;
    using sofa::defaulttype::Vec6;
    using std::get;

/*!
 * \class LegendrePolynomialsMapping
 * @brief Computes and map the length of the beams
 *
 * This is a component:
 * https://www.sofa-framework.org/community/doc/programming-with-sofa/create-your-component/
 */

template<class TIn, class TOut>
class LegendrePolynomialsMapping : public core::Mapping<TIn, TOut>, public component::mapping<TIn, TOut>
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
    typedef typename Out::Deriv Deriv;
    typedef typename Out::MatrixDeriv OutMatrixDeriv;
    typedef typename In::Real InReal;
    typedef typename In::Deriv InDeriv;
    typedef typename InDeriv::Pos DPos;
    typedef typename InDeriv::Rot DRot;
    typedef typename In::VecCoord InVecCoord;
    typedef typename In::VecDeriv InVecDeriv;
    typedef typename In::MatrixDeriv InMatrixDeriv;
    typedef typename Coord::value_type Real;

protected:
    LegendrePolynomialsMapping();
    virtual ~LegendrePolynomialsMapping() {}

public:
    sofa::Size addPoint(const Coord& c);
    sofa::Size addPoint(const Coord& c, sofa::Index indexFrom);

    void init() override;

    /// Compute the local coordinates based on the current output coordinates.
    void reinit() override;

    void apply(const core::MechanicalParams *mparams, Data<VecCoord>& out, const Data<InVecCoord>& in) override;

    void applyJ(const core::MechanicalParams *mparams, Data<VecDeriv>& out, const Data<InVecDeriv>& in) override;

    void applyJT(const core::MechanicalParams *mparams, Data<InVecDeriv>& out, const Data<VecDeriv>& in) override;

    void applyJT(const core::ConstraintParams *cparams, Data<InMatrixDeriv>& out, const Data<OutMatrixDeriv>& in) override;

    void applyDJT(const core::MechanicalParams* mparams, core::MultiVecDerivId parentForce, core::ConstMultiVecDerivId  childForce ) override;

    void draw(const core::visual::VisualParams* vparams) override;

};
#if  !defined(SOFA_COMPONENT_MAPPING_SOFA_LEGENDREPOLYNOMIALSMAPPING_CPP)
    extern template class SOFA_SOFARIGID_API LegendrePolynomialsMapping< sofa::defaulttype::Vec3Types, sofa::defaulttype::Vec3Types >;
#endif

}
#endif //SOFA_LEGENDREPOLYNOMIALSMAPPING_H
