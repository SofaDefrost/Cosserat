#pragma once

#include <liegroups/SO3.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/objectmodel/Data.h>
#include <sofa/type/Vec.h>
#include <sofa/type/Mat.h>
#include <sofa/type/vector.h>
#include <sofa/core/objectmodel/Link.h>

namespace sofa::component::cosserat::engine {

/**
 * @brief Intrinsic state for a staggered Cosserat beam.
 *
 * Stores N+1 node positions (Vec3) and N segment orientations (SO3) as
 * SOFA Data fields.  Time integration is performed externally (Python
 * explicit-Euler controller on the `painless/python-explicit` branch).
 *
 * ## DOF layout
 *   d_positions    : vector<Vec3d>  size N+1
 *   d_orientations : vector<SO3>    size N
 *
 * Velocities and accelerations (tangent-space, R³) are stored for the
 * Python integrator but are not connected to SOFA's solver pipeline.
 */
class CosseratIntrinsicState : public sofa::core::objectmodel::BaseObject {
   public:
    SOFA_CLASS(CosseratIntrinsicState, sofa::core::objectmodel::BaseObject);

    using Vec3d = sofa::type::Vec3d;
    using SO3 = sofa::component::cosserat::liegroups::SO3<double>;

    /**
     * @brief SOFA template name for factory lookup.
     *
     * Registers this component under the "Vec3d" template so that Python scenes
     * can use `template="Vec3d"` when calling `node.addObject("CosseratIntrinsicState", ...)`.
     */
    static std::string templateName(const CosseratIntrinsicState* = nullptr) {
        return "Vec3d";
    }

    CosseratIntrinsicState();
    virtual ~CosseratIntrinsicState() = default;

    // Data containers
    sofa::core::objectmodel::Data<sofa::type::vector<Vec3d>> d_positions;
    /// Orientations stored as rotation vectors ω = log(R) ∈ ℝ³ (SOFA-accessible from Python).
    /// Identity rotation = [0, 0, 0].  Convert to SO3 via SO3::exp(ω).
    sofa::core::objectmodel::Data<sofa::type::vector<Vec3d>> d_orientations;

    // Velocity and acceleration in tangent space
    sofa::core::objectmodel::Data<sofa::type::vector<Vec3d>> d_velocities;  // R3 for nodes
    sofa::core::objectmodel::Data<sofa::type::vector<Vec3d>>
        d_angularVelocities;  // R3 (so3) for segments
    sofa::core::objectmodel::Data<sofa::type::vector<Vec3d>> d_accelerations;  // R3 for nodes
    sofa::core::objectmodel::Data<sofa::type::vector<Vec3d>>
        d_angularAccelerations;  // R3 (so3) for segments

    // Secure accessors
    const sofa::type::vector<Vec3d>& getPositions() const { return d_positions.getValue(); }

    /// Returns segment orientations as SO3 objects (converted from stored rotation vectors).
    sofa::type::vector<SO3> getOrientations() const;

    size_t getNbSegments() const { return d_orientations.getValue().size(); }

    /**
     * @brief Returns the rest lengths of the N segments, computed at init().
     *
     * rest_lengths[i] = ||x_{i+1} - x_i|| in the reference configuration.
     * Used by PainlessBeamForceField to normalise strains.
     */
    const std::vector<double> &getRestLengths() const { return rest_lengths; }

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

};

}  // namespace sofa::component::cosserat::engine
