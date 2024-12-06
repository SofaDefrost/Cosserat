//
// Created by younes on 17/11/2021.
//
#define SOFA_COSSERAT_CPP_LegendrePolynomialsMapping
#include "LegendrePolynomialsMapping.inl"

#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/ObjectFactory.h>

namespace Cosserat
{
using namespace sofa::defaulttype;

// Register in the Factory
void registerLegendrePolynomialsMapping(sofa::core::ObjectFactory* factory) {
    factory->registerObjects(sofa::core::ObjectRegistrationData("Set the positions and velocities of points attached to a rigid parent")
            .add<sofa::component::mapping::LegendrePolynomialsMapping<sofa::defaulttype::Vec3Types, sofa::defaulttype::Vec3Types>>());
}

} // namespace sofa::component::mapping

namespace sofa::component::mapping
{
    template class SOFA_COSSERAT_API LegendrePolynomialsMapping< sofa::defaulttype::Vec3Types, sofa::defaulttype::Vec3Types >;
}