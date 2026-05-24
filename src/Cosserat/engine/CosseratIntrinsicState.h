#pragma once

#include <liegroups/SO3.h>
#include <sofa/core/behavior/BaseMechanicalState.h>
#include <sofa/core/objectmodel/Data.h>
#include <sofa/type/Vec.h>
#include <sofa/type/Mat.h>
#include <sofa/type/vector.h>
#include <sofa/core/objectmodel/SingleLink.h>
#include <sofa/component/topology/container/dynamic/EdgeSetTopologyContainer.h>

namespace sofa::component::cosserat::engine {

/**
 * @brief Intrinsic state for a Cosserat beam.
 * Manages positions (N+1 nodes) and orientations (N segments).
 */
class CosseratIntrinsicState : public sofa::core::behavior::BaseMechanicalState {
   public:
    SOFA_CLASS(CosseratIntrinsicState, sofa::core::behavior::BaseMechanicalState);

    using Vec3d = sofa::type::Vec3d;
    using SO3 = sofa::component::cosserat::liegroups::SO3<double>;

    CosseratIntrinsicState();
    virtual ~CosseratIntrinsicState() = default;

    // Data containers
    sofa::core::objectmodel::Data<sofa::type::vector<Vec3d>> d_positions;
    sofa::core::objectmodel::Data<sofa::type::vector<SO3>> d_orientations;

    // Velocity and acceleration in tangent space
    sofa::core::objectmodel::Data<sofa::type::vector<Vec3d>> d_velocities;  // R3 for nodes
    sofa::core::objectmodel::Data<sofa::type::vector<Vec3d>>
        d_angularVelocities;  // R3 (so3) for segments
    sofa::core::objectmodel::Data<sofa::type::vector<Vec3d>> d_accelerations;  // R3 for nodes
    sofa::core::objectmodel::Data<sofa::type::vector<Vec3d>>
        d_angularAccelerations;  // R3 (so3) for segments

    // Secure accessors
    const sofa::type::vector<Vec3d> &getPositions() const { return d_positions.getValue(); }
    const sofa::type::vector<SO3> &getOrientations() const { return d_orientations.getValue(); }
    size_t getNbSegments() const { return d_orientations.getValue().size(); }

    /**
     * @brief Returns the rest lengths of the N segments, computed at init().
     *
     * rest_lengths[i] = ||x_{i+1} - x_i|| in the reference configuration.
     * Used by PainlessBeamForceField to normalise strains.
     */
    const std::vector<double> &getRestLengths() const { return rest_lengths; }

    // Virtual methods from BaseMechanicalState
    unsigned int getSpaceDimensions() const override { return 3; }
    unsigned int getDoFCount() const override { return getPositions().size(); }

    /**
     * @brief Computes the linear strain on segment i.
     * @param i Segment index (0 to N-1).
     * @return Elongation and shear strain vector.
     */
    sofa::type::Vec3d getLinearStrain(size_t i);

    /**
     * @brief Computes the angular strain (curvature + torsion) at node i.
     *
     * Staggered convention (Romanyà-Serrasolsas et al., SIGGRAPH 2025):
     * @code
     *   Omega_i = log(R_{i-1}^T * R_i) / h_dual
     *   h_dual  = (rest_lengths[i-1] + rest_lengths[i]) / 2
     * @endcode
     * where @c h_dual is the length of the dual edge connecting the midpoints of
     * segments i-1 and i.  For uniform spacing h_dual = h.
     *
     * @param i Node index (1 to N-1).  Returns (0,0,0) for i=0 and i=N
     *          (boundary strains are enforced at the ForceField level).
     * @return Angular strain vector in the body frame [torsion, bending_y, bending_z].
     */
    sofa::type::Vec3d getAngularStrain(size_t i);

    /**
     * @brief Implements the inverse of the right Lie Jacobian of SO(3), Jr^{-1}(omega).
     * @param omega The angular velocity vector.
     * @return The 3x3 inverse Jacobian matrix.
     */
    static sofa::type::Mat3x3d getInverseLieJacobian(const sofa::type::Vec3d& omega);

    void updateStrainsCache();

protected:
    void init() override;

private:
    std::vector<sofa::type::Vec3d> cached_linear_strains;
    std::vector<sofa::type::Vec3d> cached_angular_strains;
    std::vector<double> rest_lengths;
    unsigned int last_positions_counter = 0;
    unsigned int last_orientations_counter = 0;

    // As a custom state, it is important to implement resize/clear methods if needed by solvers.
    void resize(unsigned int count) override {
        // Unclear how this translates to N and N+1, just a stub for now.
    }
};

}  // namespace sofa::component::cosserat::engine
