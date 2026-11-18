/******************************************************************************
 * Cosserat plugin — KinematicChainBiasForceField
 *
 * Inertial force  -M (Jdot qdot)  on the frames of a KinematicChainMapping.
 *
 * With a mass defined on the mapped frames (Rigid3d) SOFA builds J^T M J and J^T f, but nothing provides the
 * velocity-dependent part of the frame accelerations a_i = J qddot + (Jdot qdot)_i. The exact joint-space equation is
 *
 *       (J^T M J) qddot = J^T ( f - M (Jdot qdot) - gyroscopic terms ),
 *
 * so this force field adds  -M (Jdot qdot)  to the force of the frames; the mapping then carries it to the joints through
 * J^T like any other force. The mass is applied through its own addMDx, hence any mass component works (UniformMass,
 * DiagonalMass, RigidBodyMass...). SOFA's own rigid masses use a constant inertia matrix and no gyroscopic term: they are
 * exact for isotropic inertia (I = i*Id) with the centre of mass on the frame origin. For a general inertia and an offset
 * centre of mass use RigidBodyMass (same plugin): with it route A equals the spatial-inertia dynamics of liegroups.
 *
 * Put it in the node of the frames (next to the mapping and the mass). The explicit-in-velocity evaluation means that no
 * damping-like derivative d(bias)/d(qdot) is given to implicit solvers (addDForce is empty).
 ******************************************************************************/
#pragma once

#include <Cosserat/config.h>
#include <Cosserat/mapping/KinematicChainMapping.h>

#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/behavior/Mass.h>
#include <sofa/defaulttype/RigidTypes.h>

namespace Cosserat::mapping {

class SOFA_COSSERAT_API KinematicChainBiasForceField : public sofa::core::behavior::ForceField<sofa::defaulttype::Rigid3Types> {
public:
    SOFA_CLASS(KinematicChainBiasForceField, SOFA_TEMPLATE(sofa::core::behavior::ForceField, sofa::defaulttype::Rigid3Types));
    using DataTypes = sofa::defaulttype::Rigid3Types;
    using Inherit = sofa::core::behavior::ForceField<DataTypes>;
    using VecCoord = DataTypes::VecCoord;
    using VecDeriv = DataTypes::VecDeriv;
    using DataVecCoord = sofa::Data<VecCoord>;
    using DataVecDeriv = sofa::Data<VecDeriv>;
    using Mapping = KinematicChainMapping<sofa::defaulttype::Vec1Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types>;
    using Mass = sofa::core::behavior::Mass<DataTypes>;

    sofa::SingleLink<KinematicChainBiasForceField, Mapping, sofa::BaseLink::FLAG_STOREPATH | sofa::BaseLink::FLAG_STRONGLINK> l_mapping;
    sofa::SingleLink<KinematicChainBiasForceField, Mass, sofa::BaseLink::FLAG_STOREPATH | sofa::BaseLink::FLAG_STRONGLINK> l_mass;
    sofa::Data<double> d_scale; ///< factor on the bias force (1 = exact, 0 = off)

protected:
    KinematicChainBiasForceField();

public:
    void init() override;

    void addForce(const sofa::core::MechanicalParams *mparams, DataVecDeriv &f, const DataVecCoord &x, const DataVecDeriv &v) override;
    void addDForce(const sofa::core::MechanicalParams *, DataVecDeriv &, const DataVecDeriv &) override {}
    SReal getPotentialEnergy(const sofa::core::MechanicalParams *, const DataVecCoord &) const override { return 0; }
    void buildStiffnessMatrix(sofa::core::behavior::StiffnessMatrix *) override {}
    void buildDampingMatrix(sofa::core::behavior::DampingMatrix *) override {}
};

} // namespace Cosserat::mapping
