//
// Created by younes on 17/11/2021.
//
#define SOFA_COSSERAT_CPP_LegendrePolynomialsMapping
#include <Cosserat/mapping/LegendrePolynomialsMapping.inl>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/core/ObjectFactory.h>


using namespace sofa::defaulttype;
namespace Cosserat::mapping
{
    template class SOFA_COSSERAT_API LegendrePolynomialsMapping< sofa::defaulttype::Vec3Types, sofa::defaulttype::Vec3Types >;
}

namespace Cosserat
{
    // Register in the Factory
    void registerLegendrePolynomialsMapping(sofa::core::ObjectFactory* factory) {
        factory->registerObjects(sofa::core::ObjectRegistrationData(
            "In contrast to the DiscreteCosseratMapping component, which employs a constant curvature "
            "formulation, this component utilizes Legendre polynomials to enable a variable curvature formulation.")
                .add<mapping::LegendrePolynomialsMapping<Vec3Types, Vec3Types>>());
    }

} // namespace cosserat