#include <Cosserat/mapping/StaggeredCosseratMapping.h>

#include <sofa/core/ObjectFactory.h>
#include <sofa/core/objectmodel/Event.h>
#include <sofa/simulation/AnimateEndEvent.h>
#include <sofa/helper/logging/Messaging.h>

// Registration
namespace Cosserat {
void registerStaggeredCosseratMapping(sofa::core::ObjectFactory* factory) {
    factory->registerObjects(
        sofa::core::ObjectRegistrationData("StaggeredCosseratMapping")
            .add<sofa::component::cosserat::mapping::StaggeredCosseratMapping>());
}
}  // namespace Cosserat

namespace sofa::component::cosserat::mapping {

// ── Constructor ───────────────────────────────────────────────────────────────

StaggeredCosseratMapping::StaggeredCosseratMapping()
    : l_state(initLink("state",
                       "CosseratIntrinsicState that holds the staggered beam DOFs")),
      d_nodeFrames(initData(&d_nodeFrames, "nodeFrames",
                            "OUTPUT — Rigid3d frames at the N+1 node positions")),
      d_segmentFrames(initData(&d_segmentFrames, "segmentFrames",
                               "OUTPUT — Rigid3d frames at the N segment midpoints")),
      d_drawBeam(initData(&d_drawBeam, true, "drawBeam",
                          "Draw the beam centerline and frames")),
      d_drawRadius(initData(&d_drawRadius, 0.05, "drawRadius",
                            "Radius of the beam tube [m]")),
      d_drawAxisLength(initData(&d_drawAxisLength, 0.0, "drawAxisLength",
                                "Length of frame axes (0 = 3 × drawRadius)")),
      d_color(initData(&d_color,
                       sofa::type::RGBAColor(0.2f, 0.5f, 0.9f, 0.85f),
                       "color", "Color of the beam tube (RGBA)")),
      d_drawFrames(initData(&d_drawFrames, true, "drawFrames",
                            "Draw the segment trihedra (X/Y/Z axes in R/G/B)")) {
    // Subscribe to animation events so update() runs each step.
    f_listening.setValue(true);
}

// ── Lifecycle ─────────────────────────────────────────────────────────────────

void StaggeredCosseratMapping::init() {
    BaseObject::init();

    if (!l_state.get()) {
        msg_error() << "No CosseratIntrinsicState linked. "
                       "StaggeredCosseratMapping will produce no output.";
        return;
    }

    // Run an initial update so d_nodeFrames / d_segmentFrames are populated
    // even before the first animation step.
    update();
}

// ── Event handling ────────────────────────────────────────────────────────────

void StaggeredCosseratMapping::handleEvent(sofa::core::objectmodel::Event* event) {
    if (sofa::simulation::AnimateEndEvent::checkEventType(event)) {
        update();
    }
}

// ── Internal helpers ──────────────────────────────────────────────────────────

StaggeredCosseratMapping::Quatd
StaggeredCosseratMapping::so3ToQuat(const SO3& R) {
    // SO3 stores an Eigen unit quaternion; SOFA Quat is (x, y, z, w).
    const auto& eq = R.quaternion();   // Eigen::Quaternion<double>
    return Quatd(eq.x(), eq.y(), eq.z(), eq.w());
}

StaggeredCosseratMapping::SO3
StaggeredCosseratMapping::slerp(const SO3& Ra, const SO3& Rb, double t) {
    // Geodesic interpolation: Ra · exp(t · log(Ra^T · Rb))
    const SO3 rel = Ra.inverse() * Rb;
    const SO3::TangentVector log_rel = rel.log();
    return Ra * SO3::exp(log_rel * t);
}

// ── Core update ───────────────────────────────────────────────────────────────

void StaggeredCosseratMapping::update() {
    const auto* state = l_state.get();
    if (!state) return;

    const auto& pos = state->getPositions();
    const auto& R   = state->getOrientations();
    const size_t N   = R.size();    // number of segments
    const size_t Np1 = pos.size();  // N+1 nodes

    if (Np1 != N + 1 || N == 0) {
        msg_warning() << "update() skipped — state sizes inconsistent "
                      << "(pos=" << Np1 << ", R=" << N << ").";
        return;
    }

    // ── Node frames (N+1) ─────────────────────────────────────────────────────
    //
    //  Orientation is interpolated from adjacent segment frames via geodesic
    //  mid-point in SO(3):
    //    node 0         → R_0
    //    node i ∈[1,N-1]→ slerp(R_{i-1}, R_i, 0.5)
    //    node N         → R_{N-1}
    //
    VecRigid& node_frames = *d_nodeFrames.beginEdit();
    node_frames.resize(Np1);

    // Boundary nodes
    node_frames[0].getCenter()      = pos[0];
    node_frames[0].getOrientation() = so3ToQuat(R[0]);

    node_frames[N].getCenter()      = pos[N];
    node_frames[N].getOrientation() = so3ToQuat(R[N - 1]);

    // Interior nodes: geodesic mid-point between the two adjacent segment frames
    for (size_t i = 1; i < N; ++i) {
        node_frames[i].getCenter()      = pos[i];
        node_frames[i].getOrientation() = so3ToQuat(slerp(R[i - 1], R[i], 0.5));
    }
    d_nodeFrames.endEdit();

    // ── Segment frames (N) ────────────────────────────────────────────────────
    //
    //  Each frame sits at the midpoint of edge i with the segment's own R_i.
    //
    VecRigid& seg_frames = *d_segmentFrames.beginEdit();
    seg_frames.resize(N);

    for (size_t i = 0; i < N; ++i) {
        const Vec3d mid = (pos[i] + pos[i + 1]) * 0.5;
        seg_frames[i].getCenter()      = mid;
        seg_frames[i].getOrientation() = so3ToQuat(R[i]);
    }
    d_segmentFrames.endEdit();
}

// ── Visualisation ─────────────────────────────────────────────────────────────

void StaggeredCosseratMapping::draw(
    const sofa::core::visual::VisualParams* vparams) {
    if (!vparams->displayFlags().getShowVisualModels()) return;
    if (!d_drawBeam.getValue()) return;

    const auto* state = l_state.get();
    if (!state) return;

    const auto& pos = state->getPositions();
    const auto& R   = state->getOrientations();
    const size_t N   = R.size();
    if (pos.size() != N + 1 || N == 0) return;

    const double radius  = d_drawRadius.getValue();
    double axisLen       = d_drawAxisLength.getValue();
    if (axisLen <= 0.0) axisLen = radius * 3.0;

    const sofa::type::RGBAColor col = d_color.getValue();

    auto* dt = vparams->drawTool();

    // ── Beam centerline: cylinders between adjacent nodes ─────────────────────
    for (size_t i = 0; i < N; ++i) {
        dt->drawCylinder(
            sofa::type::Vec3f(
                static_cast<float>(pos[i].x()),
                static_cast<float>(pos[i].y()),
                static_cast<float>(pos[i].z())),
            sofa::type::Vec3f(
                static_cast<float>(pos[i + 1].x()),
                static_cast<float>(pos[i + 1].y()),
                static_cast<float>(pos[i + 1].z())),
            static_cast<float>(radius),
            col);
    }

    // ── Node spheres ──────────────────────────────────────────────────────────
    for (size_t i = 0; i <= N; ++i) {
        dt->drawSphere(
            sofa::type::Vec3f(
                static_cast<float>(pos[i].x()),
                static_cast<float>(pos[i].y()),
                static_cast<float>(pos[i].z())),
            static_cast<float>(radius * 1.2f));
    }

    // ── Segment trihedra ──────────────────────────────────────────────────────
    if (!d_drawFrames.getValue()) return;

    const sofa::type::RGBAColor red   (1.0f, 0.1f, 0.1f, 1.0f);
    const sofa::type::RGBAColor green (0.1f, 0.8f, 0.1f, 1.0f);
    const sofa::type::RGBAColor blue  (0.1f, 0.1f, 1.0f, 1.0f);

    for (size_t i = 0; i < N; ++i) {
        const Vec3d mid = (pos[i] + pos[i + 1]) * 0.5;
        const sofa::type::Vec3f origin(
            static_cast<float>(mid.x()),
            static_cast<float>(mid.y()),
            static_cast<float>(mid.z()));

        // Build the three axis directions from R_i
        // R_i.act(e_k) gives the k-th column of the rotation matrix
        const auto xdir = R[i].act(SO3::Vector(1, 0, 0));
        const auto ydir = R[i].act(SO3::Vector(0, 1, 0));
        const auto zdir = R[i].act(SO3::Vector(0, 0, 1));

        const auto tip = [&](const SO3::Vector& d) {
            return sofa::type::Vec3f(
                static_cast<float>(mid.x() + axisLen * d.x()),
                static_cast<float>(mid.y() + axisLen * d.y()),
                static_cast<float>(mid.z() + axisLen * d.z()));
        };

        dt->drawArrow(origin, tip(xdir), static_cast<float>(radius * 0.3), red);
        dt->drawArrow(origin, tip(ydir), static_cast<float>(radius * 0.3), green);
        dt->drawArrow(origin, tip(zdir), static_cast<float>(radius * 0.3), blue);
    }
}

}  // namespace sofa::component::cosserat::mapping
