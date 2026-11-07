#define SOFA_COSSERAT_CPP_KinematicChainMapping
#include <Cosserat/mapping/KinematicChainMapping.inl>

#include <sofa/core/ObjectFactory.h>

namespace Cosserat::mapping {
template class SOFA_COSSERAT_API KinematicChainMapping<sofa::defaulttype::Vec1Types, sofa::defaulttype::Rigid3Types,
                                                       sofa::defaulttype::Rigid3Types>;
} // namespace Cosserat::mapping

namespace Cosserat {
void registerKinematicChainMapping(sofa::core::ObjectFactory *factory) {
    factory->registerObjects(
        sofa::core::ObjectRegistrationData(
            "Maps joint coordinates (Vec1d: revolute / prismatic joints) to the world poses of a kinematic tree of rigid "
            "segments (Rigid3d) by product of exponentials on SE(3) (Lie group code of the plugin). Optional input2: pose "
            "of a moving base. Counterpart of ArticulatedSystemMapping, with trees and full rest poses.")
            .add<mapping::KinematicChainMapping<sofa::defaulttype::Vec1Types, sofa::defaulttype::Rigid3Types,
                                                sofa::defaulttype::Rigid3Types>>());
}
} // namespace Cosserat
