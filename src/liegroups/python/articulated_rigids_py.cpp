// Standalone pybind11 module `articulated_rigids` (no SOFA dependency): joints between two rigid bodies
// built on the Cosserat plugin Lie group code (SE3 product of exponentials).
//
// Pose convention: (x, y, z, qx, qy, qz, qw) -- same layout as a SOFA Rigid3d (quaternion xyzw).
// Twist convention (se3_exp / se3_log): (phi, rho), angular part first, as in SE3.h.
//
// Build with build_python_binding.sh.

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <array>
#include <cmath>
#include <stdexcept>

#include <liegroups/KinematicChain.h>
#include <liegroups/KinematicChainDynamics.h>
#include <liegroups/MultiDofJoints.h>
#include <liegroups/RevoluteJoint.h>

namespace py = pybind11;
using namespace sofa::component::cosserat::liegroups;

namespace {

using Vec3 = Eigen::Vector3d;
using Vec2 = Eigen::Vector2d;
using Pose = std::array<double, 7>;
using Triple = std::array<double, 3>;

Vec3 toVec3(const py::object &o, const char *what) {
    py::sequence s = py::reinterpret_borrow<py::sequence>(o);
    if (!PySequence_Check(o.ptr()) || s.size() != 3)
        throw py::value_error(std::string(what) + " must be a sequence of 3 numbers");
    return Vec3(s[0].cast<double>(), s[1].cast<double>(), s[2].cast<double>());
}

Vec2 toVec2(const py::object &o, const char *what) {
    py::sequence s = py::reinterpret_borrow<py::sequence>(o);
    if (!PySequence_Check(o.ptr()) || s.size() != 2)
        throw py::value_error(std::string(what) + " must be a sequence of 2 numbers");
    return Vec2(s[0].cast<double>(), s[1].cast<double>());
}

SE3d toSE3(const py::object &o, const char *what) {
    if (!PySequence_Check(o.ptr()))
        throw py::value_error(std::string(what) + " must be a sequence (x, y, z, qx, qy, qz, qw)");
    py::sequence s = py::reinterpret_borrow<py::sequence>(o);
    if (s.size() != 7)
        throw py::value_error(std::string(what) + " must have 7 components (x, y, z, qx, qy, qz, qw)");
    double v[7];
    for (int i = 0; i < 7; ++i) v[i] = s[i].cast<double>();
    Eigen::Quaterniond q(v[6], v[3], v[4], v[5]); // Eigen ctor order: w, x, y, z
    if (q.norm() < 1e-12) throw py::value_error(std::string(what) + ": null quaternion");
    q.normalize();
    return SE3d(SO3d(q), Vec3(v[0], v[1], v[2]));
}

Pose fromSE3(const SE3d &T) {
    Eigen::Quaterniond q = T.rotation().quaternion();
    if (q.w() < 0) q.coeffs() = -q.coeffs(); // canonical hemisphere
    const Vec3 &t = T.translation();
    return {t.x(), t.y(), t.z(), q.x(), q.y(), q.z(), q.w()};
}

py::tuple toTuple(const Pose &p) {
    return py::make_tuple(p[0], p[1], p[2], p[3], p[4], p[5], p[6]);
}

py::tuple toTuple3(const Vec3 &v) { return py::make_tuple(v.x(), v.y(), v.z()); }

} // namespace

PYBIND11_MODULE(articulated_rigids, m) {
    m.doc() = "Articulated rigid bodies (revolute, spherical, two-axis joints) on SE(3). "
              "Pose = (x, y, z, qx, qy, qz, qw); twist = (phi, rho).";

    py::class_<RevoluteJoint<double>>(m, "RevoluteJoint",
        "Hinge, 1 DoF, q = angle (float). Child pose in parent = Exp(q xi) * rest_pose. "
        "Also models the trochoid (pivot) joint with the axis along the long axis of the bones.")
        .def(py::init([](const py::object &axis, const py::object &pivot, const py::object &rest) {
                 return RevoluteJoint<double>(toVec3(axis, "axis"), toVec3(pivot, "pivot"), toSE3(rest, "rest_pose"));
             }),
             py::arg("axis"), py::arg("pivot"), py::arg("rest_pose") = py::make_tuple(0., 0., 0., 0., 0., 0., 1.))
        .def_property_readonly("dof", [](const RevoluteJoint<double> &) { return 1; })
        .def("relative_pose", [](const RevoluteJoint<double> &j, double q) { return toTuple(fromSE3(j.relativePose(q))); }, py::arg("q"))
        .def("child_pose", [](const RevoluteJoint<double> &j, const py::object &parent, double q) {
                 return toTuple(fromSE3(j.childPose(toSE3(parent, "parent_pose"), q)));
             }, py::arg("parent_pose"), py::arg("q"))
        .def("joint_angle", [](const RevoluteJoint<double> &j, const py::object &relative) {
                 return j.jointAngle(toSE3(relative, "relative_pose"));
             }, py::arg("relative_pose"));

    py::class_<SphericalJoint<double>>(m, "SphericalJoint",
        "Ball-and-socket, 3 DoF, q = rotation vector (rx, ry, rz) about pivot.")
        .def(py::init([](const py::object &pivot, const py::object &rest) {
                 return SphericalJoint<double>(toVec3(pivot, "pivot"), toSE3(rest, "rest_pose"));
             }),
             py::arg("pivot"), py::arg("rest_pose") = py::make_tuple(0., 0., 0., 0., 0., 0., 1.))
        .def_property_readonly("dof", [](const SphericalJoint<double> &) { return 3; })
        .def("relative_pose", [](const SphericalJoint<double> &j, const py::object &q) {
                 return toTuple(fromSE3(j.relativePose(toVec3(q, "q"))));
             }, py::arg("q"))
        .def("child_pose", [](const SphericalJoint<double> &j, const py::object &parent, const py::object &q) {
                 return toTuple(fromSE3(j.childPose(toSE3(parent, "parent_pose"), toVec3(q, "q"))));
             }, py::arg("parent_pose"), py::arg("q"))
        .def("joint_angles", [](const SphericalJoint<double> &j, const py::object &relative) {
                 return toTuple3(j.jointAngles(toSE3(relative, "relative_pose")));
             }, py::arg("relative_pose"));

    py::class_<TwoAxisJoint<double>>(m, "TwoAxisJoint",
        "Two successive hinges, 2 DoF, q = (q1, q2): Exp(q1 xi1) Exp(q2 xi2) rest_pose, screws in the parent "
        "frame at q = 0. point1 == point2: condylar joint; perpendicular offset axes: saddle joint.")
        .def(py::init([](const py::object &a1, const py::object &p1, const py::object &a2, const py::object &p2,
                         const py::object &rest) {
                 return TwoAxisJoint<double>(toVec3(a1, "axis1"), toVec3(p1, "point1"), toVec3(a2, "axis2"),
                                             toVec3(p2, "point2"), toSE3(rest, "rest_pose"));
             }),
             py::arg("axis1"), py::arg("point1"), py::arg("axis2"), py::arg("point2"),
             py::arg("rest_pose") = py::make_tuple(0., 0., 0., 0., 0., 0., 1.))
        .def_property_readonly("dof", [](const TwoAxisJoint<double> &) { return 2; })
        .def("relative_pose", [](const TwoAxisJoint<double> &j, const py::object &q) {
                 return toTuple(fromSE3(j.relativePose(toVec2(q, "q"))));
             }, py::arg("q"))
        .def("child_pose", [](const TwoAxisJoint<double> &j, const py::object &parent, const py::object &q) {
                 return toTuple(fromSE3(j.childPose(toSE3(parent, "parent_pose"), toVec2(q, "q"))));
             }, py::arg("parent_pose"), py::arg("q"));

    m.def("se3_exp", [](const py::object &phi, const py::object &rho) {
              SE3d::TangentVector xi;
              xi.head<3>() = toVec3(phi, "phi");
              xi.tail<3>() = toVec3(rho, "rho");
              return toTuple(fromSE3(SE3d::computeExp(xi)));
          }, py::arg("phi"), py::arg("rho"),
          "Exponential map se(3) -> SE(3), twist ordered (phi angular, rho linear). Returns a pose (x,y,z,qx,qy,qz,qw).");
    m.def("se3_log", [](const py::object &pose) {
              const SE3d::TangentVector xi = toSE3(pose, "pose").computeLog();
              return py::make_tuple(toTuple3(xi.head<3>()), toTuple3(xi.tail<3>()));
          }, py::arg("pose"), "Logarithm SE(3) -> se(3). Returns (phi, rho).");

    // ---- Generic SE(3) operations (poses as 7-tuples, matrices as numpy arrays) ----
    m.def("se3_compose", [](const py::object &a, const py::object &b) {
              return toTuple(fromSE3(toSE3(a, "a") * toSE3(b, "b")));
          }, py::arg("a"), py::arg("b"), "Group product a * b.");
    m.def("se3_inverse", [](const py::object &a) { return toTuple(fromSE3(toSE3(a, "pose").inverse())); },
          py::arg("pose"));
    m.def("se3_act", [](const py::object &a, const py::object &p) {
              const Vec3 r = toSE3(a, "pose").act(toVec3(p, "point"));
              return toTuple3(r);
          }, py::arg("pose"), py::arg("point"), "Apply the rigid transform to a point.");
    m.def("se3_matrix", [](const py::object &a) { return Eigen::Matrix4d(toSE3(a, "pose").matrix()); },
          py::arg("pose"), "4x4 homogeneous matrix.");
    m.def("se3_adjoint", [](const py::object &a) { return Eigen::Matrix<double, 6, 6>(toSE3(a, "pose").computeAdjoint()); },
          py::arg("pose"), "6x6 adjoint Ad_g in the [phi; rho] layout.");
    m.def("se3_ad", [](const Eigen::Matrix<double, 6, 1> &v) { return Eigen::Matrix<double, 6, 6>(SE3d::computeAd(v)); },
          py::arg("twist"), "6x6 Lie-algebra adjoint ad_v, [phi; rho] layout.");
    m.def("se3_bracket", [](const Eigen::Matrix<double, 6, 1> &v, const Eigen::Matrix<double, 6, 1> &w) {
              return Eigen::Matrix<double, 6, 1>(SE3d::lieBracket(v, w));
          }, py::arg("v"), py::arg("w"), "Lie bracket [v, w].");
    m.def("se3_left_jacobian_inverse", [](const Eigen::Matrix<double, 6, 1> &xi) {
              return Eigen::Matrix<double, 6, 6>(SE3d::computeLeftJacobianInverse(xi));
          }, py::arg("twist"));

    // ---- Serial chain (product of exponentials) ----
    using Chain = KinematicChain<double>;
    py::class_<Chain>(m, "KinematicChain",
                      "Kinematic tree of revolute / prismatic / fixed segments (product of exponentials on SE(3)). "
                      "parent = -1 attaches a segment to the base, -2 (default) to the previous segment.")
        .def(py::init<>())
        .def("add_revolute", [](Chain &c, const py::object &axis, const py::object &pivot, const py::object &rest,
                                int parent) {
                 return c.addRevolute(toVec3(axis, "axis"), toVec3(pivot, "pivot"),
                                      rest.is_none() ? SE3d() : toSE3(rest, "rest"), parent);
             }, py::arg("axis"), py::arg("pivot") = py::make_tuple(0., 0., 0.), py::arg("rest") = py::none(),
             py::arg("parent") = Chain::Previous)
        .def("add_prismatic", [](Chain &c, const py::object &axis, const py::object &rest, int parent) {
                 return c.addPrismatic(toVec3(axis, "axis"), rest.is_none() ? SE3d() : toSE3(rest, "rest"), parent);
             }, py::arg("axis"), py::arg("rest") = py::none(), py::arg("parent") = Chain::Previous)
        .def("add_fixed", [](Chain &c, const py::object &rest, int parent) { return c.addFixed(toSE3(rest, "rest"), parent); },
             py::arg("rest"), py::arg("parent") = Chain::Previous)
        .def_property_readonly("n_dof", &Chain::nDof)
        .def_property_readonly("n_segments", &Chain::nSegments)
        .def("ancestor_dofs", [](const Chain &c, int i) { return c.ancestorDofs(i); }, py::arg("i"))
        .def("forward", [](const Chain &c, const Eigen::VectorXd &q) { return toTuple(fromSE3(c.forward(q))); },
             py::arg("q"))
        .def("segment_pose", [](const Chain &c, const Eigen::VectorXd &q, int i) {
                 return toTuple(fromSE3(c.segmentPose(q, i)));
             }, py::arg("q"), py::arg("i"))
        .def("all_poses", [](const Chain &c, const Eigen::VectorXd &q) {
                 py::list out;
                 for (const auto &T : c.allPoses(q)) out.append(toTuple(fromSE3(T)));
                 return out;
             }, py::arg("q"))
        .def("joint_screws", [](const Chain &c, const Eigen::VectorXd &q) { return c.jointScrews(c.allPoses(q)); },
             py::arg("q"), "6 x nDof unit screws of every joint seen from the base frame.")
        .def("spatial_jacobian", [](const Chain &c, const Eigen::VectorXd &q, int i) { return c.spatialJacobian(q, i); },
             py::arg("q"), py::arg("i") = -1)
        .def("body_jacobian", [](const Chain &c, const Eigen::VectorXd &q, int i) { return c.bodyJacobian(q, i); },
             py::arg("q"), py::arg("i") = -1)
        .def("hybrid_jacobian", [](const Chain &c, const Eigen::VectorXd &q, int i) { return c.hybridJacobian(q, i); },
             py::arg("q"), py::arg("i") = -1)
        .def("world_screws", [](const Chain &c, const Eigen::VectorXd &q, const py::object &base) {
                 return c.worldFrame(q, base.is_none() ? SE3d() : toSE3(base, "base_pose")).screws;
             }, py::arg("q"), py::arg("base_pose") = py::none(),
             "6 x nDof unit joint screws in the world frame (base pose given), [w ; r], reference point = world origin.")
        .def("bias_acceleration", [](const Chain &c, const Eigen::VectorXd &q, const Eigen::VectorXd &dq,
                                     const py::object &base, const Eigen::Matrix<double, 6, 1> &baseTwist,
                                     const Eigen::Matrix<double, 6, 1> &baseBias) {
                 const auto wf = c.worldFrame(q, base.is_none() ? SE3d() : toSE3(base, "base_pose"));
                 const auto a = c.biasAcceleration(wf.screws, wf.positions, dq, baseTwist, baseBias);
                 Eigen::MatrixXd out(a.size(), 6);
                 for (std::size_t i = 0; i < a.size(); ++i) out.row(i) = a[i].transpose();
                 return out;
             }, py::arg("q"), py::arg("dq"), py::arg("base_pose") = py::none(),
             py::arg("base_twist") = Eigen::Matrix<double, 6, 1>::Zero().eval(),
             py::arg("base_bias") = Eigen::Matrix<double, 6, 1>::Zero().eval(),
             "n_segments x 6 velocity-dependent accelerations [alpha ; a_origin] (world, hybrid), i.e. d/dt(J dx) without J xddot.")
        .def("geometric_stiffness", [](const Chain &c, const Eigen::VectorXd &q, const Eigen::MatrixXd &wrenches,
                                       const Eigen::VectorXd &dq, const py::object &base,
                                       const Eigen::Matrix<double, 6, 1> &dBaseTwist) {
                 const auto wf = c.worldFrame(q, base.is_none() ? SE3d() : toSE3(base, "base_pose"));
                 std::vector<Chain::TangentVector> w(wrenches.rows());
                 for (Eigen::Index i = 0; i < wrenches.rows(); ++i) w[i] = wrenches.row(i).transpose();
                 const auto r = c.geometricStiffness(wf.screws, wf.positions, w, dq, dBaseTwist);
                 return py::make_tuple(r.dQ, r.dMoment);
             }, py::arg("q"), py::arg("wrenches"), py::arg("dq"), py::arg("base_pose") = py::none(),
             py::arg("d_base_twist") = Eigen::Matrix<double, 6, 1>::Zero().eval(),
             "(dQ, dMoment): change of J^T W (and of the total moment about the world origin) for fixed world wrenches "
             "[torque about the segment origin ; force] per segment and a displacement (dq, d_base_twist).");

    // ---- Rigid-body dynamics of a chain (RNEA, mass matrix, forward dynamics) ----
    using Dyn = KinematicChainDynamics<double>;
    using Wrenches = std::vector<Dyn::TangentVector>;
    auto toWrenches = [](const Eigen::MatrixXd &w) {
        Wrenches out(w.rows());
        for (Eigen::Index i = 0; i < w.rows(); ++i) out[i] = w.row(i).transpose();
        return out;
    };
    py::class_<Dyn>(m, "KinematicChainDynamics",
                    "Rigid-body dynamics of a KinematicChain in the world frame with a fixed base: RNEA, M(q), h(q, qdot), "
                    "forward dynamics, energies. Inertias are given per segment in the segment frame.")
        .def(py::init<const Chain &>(), py::arg("chain"))
        .def("set_inertia", [](Dyn &d, int segment, double mass, const py::object &com, const Eigen::Matrix3d &inertia_com) {
                 Dyn::BodyInertia b;
                 b.mass = mass;
                 b.com = toVec3(com, "com");
                 b.inertiaCom = inertia_com;
                 d.setInertia(segment, b);
             }, py::arg("segment"), py::arg("mass"), py::arg("com") = py::make_tuple(0., 0., 0.),
             py::arg("inertia_com") = Eigen::Matrix3d::Zero().eval(),
             "mass, centre of mass and rotational inertia about the centre of mass (segment axes) of one segment")
        .def("spatial_inertia", [](const Dyn &d, int segment) { return Eigen::Matrix<double, 6, 6>(d.inertia(segment).matrix()); },
             py::arg("segment"), "6x6 inertia about the segment origin, layout [angular ; linear]")
        .def("inverse_dynamics", [&](const Dyn &d, const Eigen::VectorXd &q, const Eigen::VectorXd &dq, const Eigen::VectorXd &ddq,
                                     const py::object &gravity, const py::object &external) {
                 const Eigen::Vector3d g = gravity.is_none() ? Eigen::Vector3d(Eigen::Vector3d::Zero()) : toVec3(gravity, "gravity");
                 if (external.is_none()) return d.inverseDynamics(q, dq, ddq, g);
                 const Wrenches w = toWrenches(external.cast<Eigen::MatrixXd>());
                 return d.inverseDynamics(q, dq, ddq, g, &w);
             }, py::arg("q"), py::arg("dq"), py::arg("ddq"), py::arg("gravity") = py::none(), py::arg("external") = py::none(),
             "RNEA: joint forces for the given acceleration. external: n_segments x 6 [torque about the origin ; force], world frame.")
        .def("mass_matrix", &Dyn::massMatrix, py::arg("q"))
        .def("bias_forces", [&](const Dyn &d, const Eigen::VectorXd &q, const Eigen::VectorXd &dq, const py::object &gravity,
                                const py::object &external) {
                 const Eigen::Vector3d g = gravity.is_none() ? Eigen::Vector3d(Eigen::Vector3d::Zero()) : toVec3(gravity, "gravity");
                 if (external.is_none()) return d.biasForces(q, dq, g);
                 const Wrenches w = toWrenches(external.cast<Eigen::MatrixXd>());
                 return d.biasForces(q, dq, g, &w);
             }, py::arg("q"), py::arg("dq"), py::arg("gravity") = py::none(), py::arg("external") = py::none())
        .def("forward_dynamics", [&](const Dyn &d, const Eigen::VectorXd &q, const Eigen::VectorXd &dq, const Eigen::VectorXd &tau,
                                     const py::object &gravity, const py::object &external) {
                 const Eigen::Vector3d g = gravity.is_none() ? Eigen::Vector3d(Eigen::Vector3d::Zero()) : toVec3(gravity, "gravity");
                 if (external.is_none()) return d.forwardDynamics(q, dq, tau, g);
                 const Wrenches w = toWrenches(external.cast<Eigen::MatrixXd>());
                 return d.forwardDynamics(q, dq, tau, g, &w);
             }, py::arg("q"), py::arg("dq"), py::arg("tau"), py::arg("gravity") = py::none(), py::arg("external") = py::none())
        .def("kinetic_energy", &Dyn::kineticEnergy, py::arg("q"), py::arg("dq"))
        .def("potential_energy", [](const Dyn &d, const Eigen::VectorXd &q, const py::object &gravity) {
                 return d.potentialEnergy(q, toVec3(gravity, "gravity"));
             }, py::arg("q"), py::arg("gravity"));
}
