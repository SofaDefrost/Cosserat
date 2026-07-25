#pragma once

#include <Cosserat/engine/CosseratIntrinsicState.h>
#include <liegroups/SO3.h>

#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/objectmodel/Data.h>
#include <sofa/core/objectmodel/Link.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/type/RGBAColor.h>
#include <sofa/type/Vec.h>
#include <sofa/type/vector.h>

namespace sofa::component::cosserat::mapping {

/**
 * @brief Converts a staggered Cosserat state (Vec3 + SO3) into Rigid3d frames.
 *
 * The staggered Cosserat discretisation stores:
 *   - N+1 node positions  x_0 … x_N   ∈ R³
 *   - N   segment frames  R_0 … R_{N-1} ∈ SO(3)  (at edge midpoints)
 *
 * This component computes two families of output frames:
 *
 * ### Node frames (d_nodeFrames, size N+1)
 * Placed at each node position with an orientation interpolated from the
 * adjacent segment frames:
 * @code
 *   node 0       → orientation = R_0
 *   node i ∈ [1,N-1] → orientation = R_{i-1} · exp(0.5 · log(R_{i-1}^T · R_i))
 *                                     (geodesic midpoint between R_{i-1} and R_i)
 *   node N       → orientation = R_{N-1}
 * @endcode
 *
 * ### Segment frames (d_segmentFrames, size N)
 * Placed at the midpoint of each edge with orientation equal to the segment frame:
 * @code
 *   segment i → position    = (x_i + x_{i+1}) / 2
 *             → orientation = R_i
 * @endcode
 *
 * Both outputs are updated at the end of each animation step
 * (AnimateEndEvent) and can be connected to visual or collision models.
 *
 * A built-in draw() function renders the beam centerline as cylinders and
 * optionally shows the segment trihedra (X→red, Y→green, Z→blue).
 */
class StaggeredCosseratMapping : public sofa::core::objectmodel::BaseObject {
   public:
    SOFA_CLASS(StaggeredCosseratMapping, sofa::core::objectmodel::BaseObject);

    using Vec3d    = sofa::type::Vec3d;
    using Quatd    = sofa::type::Quat<double>;
    using Rigid3d  = sofa::defaulttype::Rigid3dTypes::Coord;   // position + quat
    using VecRigid = sofa::type::vector<Rigid3d>;
    using SO3      = sofa::component::cosserat::liegroups::SO3<double>;

    // ── Link to the staggered mechanical state ─────────────────────────────────
    sofa::core::objectmodel::SingleLink<
        StaggeredCosseratMapping,
        sofa::component::cosserat::engine::CosseratIntrinsicState,
        sofa::core::objectmodel::BaseLink::FLAG_STOREPATH |
        sofa::core::objectmodel::BaseLink::FLAG_STRONGLINK
    > l_state;

    // ── Computed output frames ──────────────────────────────────────────────────
    /// Rigid3d frames at the N+1 node positions (orientation interpolated).
    sofa::core::objectmodel::Data<VecRigid> d_nodeFrames;
    /// Rigid3d frames at the N segment midpoints (orientation = R_i).
    sofa::core::objectmodel::Data<VecRigid> d_segmentFrames;

    // ── Visualisation parameters ────────────────────────────────────────────────
    /// Draw the beam centerline and frames.
    sofa::core::objectmodel::Data<bool> d_drawBeam;
    /// Radius of the drawn cylinders (beam tube).
    sofa::core::objectmodel::Data<double> d_drawRadius;
    /// Length of the frame axes arrows (0 = auto = radius × 3).
    sofa::core::objectmodel::Data<double> d_drawAxisLength;
    /// Color of the beam tube.
    sofa::core::objectmodel::Data<sofa::type::RGBAColor> d_color;
    /// Draw the segment trihedra (X/Y/Z axes in R/G/B).
    sofa::core::objectmodel::Data<bool> d_drawFrames;

    // ─────────────────────────────────────────────────────────────────────────
    StaggeredCosseratMapping();
    ~StaggeredCosseratMapping() override = default;

    void init() override;

    /**
     * @brief Recompute node and segment frames from the current state.
     *
     * Called automatically via AnimateEndEvent.  Can also be called explicitly
     * after modifying the state.
     */
    void update();

    void handleEvent(sofa::core::objectmodel::Event* event) override;

    /**
     * @brief Draw the beam: centerline cylinders + optional trihedra.
     */
    void draw(const sofa::core::visual::VisualParams* vparams) override;

   private:
    /// Convert an SO3 rotation to a SOFA quaternion (x, y, z, w).
    static Quatd so3ToQuat(const SO3& R);

    /**
     * @brief Geodesic interpolation between two SO3 elements.
     * @param Ra  Left rotation.
     * @param Rb  Right rotation.
     * @param t   Parameter in [0, 1].
     * @return    Ra · exp(t · log(Ra^T · Rb))
     */
    static SO3 slerp(const SO3& Ra, const SO3& Rb, double t);
};

}  // namespace sofa::component::cosserat::mapping
