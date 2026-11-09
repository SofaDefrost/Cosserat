#define SOFA_COSSERAT_CPP_KinematicChainMapping
#include <Cosserat/mapping/KinematicChainMapping.inl>
#include <Cosserat/mapping/KinematicChainBiasForceField.h>
#include <Cosserat/mapping/RigidBodyMass.h>

#include <sofa/core/ObjectFactory.h>

namespace Cosserat::mapping {
template class SOFA_COSSERAT_API KinematicChainMapping<sofa::defaulttype::Vec1Types, sofa::defaulttype::Rigid3Types,
                                                       sofa::defaulttype::Rigid3Types>;
} // namespace Cosserat::mapping

namespace Cosserat::mapping {
KinematicChainBiasForceField::KinematicChainBiasForceField()
    : l_mapping(initLink("mapping", "The KinematicChainMapping whose Jdot qdot is used (default: the one of this node)")),
      l_mass(initLink("mass", "The mass of the frames (default: the Rigid3d mass of this node)")),
      d_scale(initData(&d_scale, 1.0, "scale", "Factor on the bias force (1 = exact, 0 = off)")) {}

void KinematicChainBiasForceField::init() {
    Inherit::init();
    if (!l_mapping) l_mapping.set(this->getContext()->get<Mapping>());
    if (!l_mass) l_mass.set(this->getContext()->get<Mass>());
    if (!l_mapping || !l_mass) {
        msg_error() << "needs a KinematicChainMapping and a mass on the frames (links 'mapping' and 'mass').";
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
        return;
    }
    this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Valid);
}

void KinematicChainBiasForceField::addForce(const sofa::core::MechanicalParams *mparams, DataVecDeriv &f, const DataVecCoord &,
                                            const DataVecDeriv &) {
    if (!this->isComponentStateValid() || d_scale.getValue() == 0.0) return;
    DataVecDeriv bias;
    l_mapping->computeBiasAcceleration(*bias.beginEdit());
    bias.endEdit();
    if (bias.getValue().size() != f.getValue().size()) return;
    l_mass->addMDx(mparams, f, bias, -d_scale.getValue()); // f -= scale * M (Jdot qdot)
}
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
    factory->registerObjects(
        sofa::core::ObjectRegistrationData(
            "Inertial force -M (Jdot qdot) on the frames of a KinematicChainMapping: with the mapped mass it completes the "
            "exact joint-space dynamics (J^T M J) qddot = J^T (f - M Jdot qdot).")
            .add<mapping::KinematicChainBiasForceField>());
    factory->registerObjects(
        sofa::core::ObjectRegistrationData(
            "Rigid3d mass with centre of mass anywhere in the body frame and a general inertia tensor: mass matrix R I R^T "
            "that turns with the body, weight and gyroscopic / centrifugal force, exact tangent. Unlike UniformMass it is "
            "exact for anisotropic inertia (use with KinematicChainBiasForceField).")
            .add<mapping::RigidBodyMass>());
}
} // namespace Cosserat
