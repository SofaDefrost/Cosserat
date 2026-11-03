// Standalone pybind11 module `articulated_rigids` (no SOFA dependency): joints between two rigid bodies
// built on the Cosserat plugin Lie group code (SE3 product of exponentials).
//
// Pose convention: (x, y, z, qx, qy, qz, qw) -- same layout as a SOFA Rigid3d (quaternion xyzw).
// Twist convention (se3_exp / se3_log): (phi, rho), angular part first, as in SE3.h.
//
// Build with build_python_binding.sh.

#include <pybind11/pybind11.h>

#include <array>
#include <cmath>
#include <stdexcept>

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
}
