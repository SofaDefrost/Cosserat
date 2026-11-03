// Tests of the revolute joint / two articulated rigid bodies built on SE3.
// Expected values are computed independently of SE3::exp (Eigen::AngleAxis,
// closed-form planar kinematics, finite differences).

#include <gtest/gtest.h>
#include <liegroups/RevoluteJoint.h>

#include <Eigen/Geometry>
#include <cmath>

using namespace sofa::component::cosserat::liegroups;

namespace {

using Vec3 = Eigen::Vector3d;
using Vec6 = Eigen::Matrix<double, 6, 1>;

constexpr double kTol = 1e-12;

// Pose of a frame rotated by q about the line (point p, unit axis u), applied to T0:
//   x -> p + R(x - p), i.e. rotation R, translation p - R p. Built with Eigen only.
SE3d rotationAboutLine(const Vec3 &u, const Vec3 &p, double q) {
    const Eigen::Matrix3d R = Eigen::AngleAxisd(q, u.normalized()).toRotationMatrix();
    return SE3d(SO3d(R), p - R * p);
}

SE3d samplePose() {
    return SE3d(SO3d(0.7, Vec3(1, 2, -1)), Vec3(0.3, -0.4, 1.1));
}

double poseError(const SE3d &a, const SE3d &b) { return (a.matrix() - b.matrix()).norm(); }

} // namespace

TEST(RevoluteJointTest, ZeroAngleGivesRestPose) {
    const SE3d rest(SO3d(0.3, Vec3(0, 1, 0)), Vec3(0.5, 0, 0));
    const RevoluteJoint<double> joint(Vec3(0, 0, 1), Vec3(0.2, 0.1, 0), rest);
    EXPECT_LT(poseError(joint.relativePose(0.0), rest), kTol);
}

TEST(RevoluteJointTest, MatchesRotationAboutLineWithEigenReference) {
    // Axis not aligned with the frame axes, pivot off the origin.
    const Vec3 axis = Vec3(1, -2, 0.5).normalized();
    const Vec3 pivot(0.4, 0.2, -0.3);
    const SE3d rest(SO3d(0.5, Vec3(0, 0, 1)), Vec3(1.0, 0.0, 0.0));
    const RevoluteJoint<double> joint(axis, pivot, rest);

    for (double q : {-2.5, -1.0, -1e-3, 1e-7, 0.4, 1.3, 3.0}) {
        const SE3d expected = rotationAboutLine(axis, pivot, q) * rest;
        EXPECT_LT(poseError(joint.relativePose(q), expected), 1e-11) << "q = " << q;
    }
}

TEST(RevoluteJointTest, PivotPointIsFixedForBothBodies) {
    const Vec3 pivot(0.4, 0.2, -0.3);
    const SE3d rest(SO3d(0.5, Vec3(0, 0, 1)), Vec3(1.0, 0.0, 0.0));
    const RevoluteJoint<double> joint(Vec3(0, 1, 1), pivot, rest);
    const SE3d parent = samplePose();

    // The point of the child body that coincides with the pivot at q = 0 (child-frame coords).
    const Vec3 pivotInChild = rest.inverse().act(pivot);
    const Vec3 pivotWorld = parent.act(pivot);

    for (double q : {-1.2, 0.0, 0.9, 2.2}) {
        const Vec3 viaChild = joint.childPose(parent, q).act(pivotInChild);
        EXPECT_LT((viaChild - pivotWorld).norm(), kTol) << "q = " << q;
    }
}

TEST(RevoluteJointTest, RigidityDistancesAreConserved) {
    const RevoluteJoint<double> joint(Vec3(0, 0, 1), Vec3(1, 0, 0), SE3d(SO3d(), Vec3(1, 0, 0)));
    const SE3d parent = samplePose();
    const Vec3 pointA(0.3, 0.1, 0.2), pointB(-0.2, 0.5, 0.4); // two points fixed in the child frame

    const double d0 = (joint.childPose(parent, 0.0).act(pointA) - joint.childPose(parent, 0.0).act(pointB)).norm();
    for (double q : {0.5, 1.5, -2.0}) {
        const SE3d child = joint.childPose(parent, q);
        EXPECT_NEAR((child.act(pointA) - child.act(pointB)).norm(), d0, kTol);
    }
}

TEST(RevoluteJointTest, OneParameterSubgroupComposition) {
    // Exp(q1 ξ) Exp(q2 ξ) = Exp((q1+q2) ξ): angles add along a hinge.
    const RevoluteJoint<double> joint(Vec3(0.3, -1, 0.2), Vec3(0.1, 0.2, 0.3));
    const double q1 = 0.8, q2 = -1.9;
    const SE3d lhs = joint.relativePose(q1) * joint.relativePose(q2);
    EXPECT_LT(poseError(lhs, joint.relativePose(q1 + q2)), 1e-11);
}

TEST(RevoluteJointTest, JointAngleRoundTrip) {
    const SE3d rest(SO3d(0.5, Vec3(0, 0, 1)), Vec3(1.0, 0.0, 0.0));
    const RevoluteJoint<double> joint(Vec3(1, 1, -1), Vec3(0.4, 0.2, -0.3), rest);
    const SE3d parent = samplePose();

    for (double q : {-3.0, -1.1, -0.01, 0.0, 0.02, 0.7, 2.9}) {
        const SE3d child = joint.childPose(parent, q);
        EXPECT_NEAR(joint.jointAngle(parent, child), q, 1e-10) << "q = " << q;
        EXPECT_LT(joint.constraintResidual(parent.inverse() * child).norm(), 1e-10) << "q = " << q;
    }
}

TEST(RevoluteJointTest, ConstraintResidualDetectsNonHingeMotion) {
    const RevoluteJoint<double> joint(Vec3(0, 0, 1), Vec3(0, 0, 0));
    // Child displaced along the axis: not reachable by the hinge.
    const SE3d offAxis(SO3d(), Vec3(0, 0, 0.1));
    EXPECT_GT(joint.constraintResidual(offAxis).norm(), 0.05);
    // Child rotated about another axis.
    const SE3d wrongAxis(SO3d(0.3, Vec3(1, 0, 0)), Vec3::Zero());
    EXPECT_GT(joint.constraintResidual(wrongAxis).norm(), 0.1);
}

TEST(RevoluteJointTest, PlanarDoublePendulumClosedForm) {
    // Three bodies: world -> link1 -> link2, hinges about z at the end of each link.
    // Tip of link 2 (at distance l2 along its x axis) has a known analytic position.
    const double l1 = 0.8, l2 = 0.5;
    const RevoluteJoint<double> j1(Vec3(0, 0, 1), Vec3::Zero());                          // at the world origin
    const RevoluteJoint<double> j2(Vec3(0, 0, 1), Vec3(l1, 0, 0), SE3d(SO3d(), Vec3(l1, 0, 0))); // at the end of link 1

    for (double q1 : {-0.9, 0.0, 0.6, 2.0}) {
        for (double q2 : {-1.5, 0.3, 1.1}) {
            const SE3d link1 = j1.childPose(SE3d(), q1);
            const SE3d link2 = j2.childPose(link1, q2);
            const Vec3 tip = link2.act(Vec3(l2, 0, 0));
            const Vec3 expected(l1 * std::cos(q1) + l2 * std::cos(q1 + q2), l1 * std::sin(q1) + l2 * std::sin(q1 + q2), 0.0);
            EXPECT_LT((tip - expected).norm(), kTol) << "q1=" << q1 << " q2=" << q2;
            // Orientation of link 2: pure rotation of q1+q2 about z.
            EXPECT_NEAR(link2.rotation().log().z(), q1 + q2, 1e-11);
        }
    }
}

TEST(RevoluteJointTest, JointTwistMatchesFiniteDifference) {
    const SE3d rest(SO3d(0.5, Vec3(0, 0, 1)), Vec3(1.0, 0.0, 0.0));
    const RevoluteJoint<double> joint(Vec3(1, -2, 0.5), Vec3(0.4, 0.2, -0.3), rest);
    const double h = 1e-6;

    for (double q : {-1.0, 0.0, 0.8}) {
        // Central difference of the body twist: log(T(q-h)^-1 T(q+h)) / 2h
        const Vec6 fd = (joint.relativePose(q - h).inverse() * joint.relativePose(q + h)).log() / (2 * h);
        EXPECT_LT((fd - joint.jointTwistInChild(q)).norm(), 1e-6) << "q = " << q;
    }
}

TEST(RevoluteJointTest, ChildBodyTwistMatchesFiniteDifference) {
    // Parent moving with a constant body twist Vp, hinge angle q(t) = q0 + qdot t.
    // Child body twist = log(T_child(t)^-1 T_child(t+dt)) / dt must match childBodyTwist().
    const SE3d rest(SO3d(0.5, Vec3(0, 0, 1)), Vec3(1.0, 0.0, 0.0));
    const RevoluteJoint<double> joint(Vec3(0, 1, 1), Vec3(0.2, 0, 0.1), rest);
    Vec6 Vp;
    Vp << 0.2, -0.1, 0.3, 0.5, 0.1, -0.2; // [angular; linear], parent frame
    const SE3d parent0 = samplePose();
    const double q0 = 0.4, qdot = 1.3, h = 1e-6;

    auto childAt = [&](double t) {
        return joint.childPose(parent0 * SE3d::computeExp(t * Vp), q0 + qdot * t);
    };
    const Vec6 fd = (childAt(-h).inverse() * childAt(h)).log() / (2 * h);
    EXPECT_LT((fd - joint.childBodyTwist(Vp, q0, qdot)).norm(), 1e-6);
}

TEST(ArticulatedRigidPairTest, ChildFollowsParentAndAngle) {
    const RevoluteJoint<double> joint(Vec3(0, 0, 1), Vec3(1, 0, 0), SE3d(SO3d(), Vec3(1, 0, 0)));
    ArticulatedRigidPair<double> pair(joint, samplePose(), 0.0);

    // q = 0: child is the parent translated by the rest offset.
    EXPECT_LT(poseError(pair.childPose(), samplePose() * joint.restPose()), kTol);

    pair.setAngle(M_PI / 2);
    // Hinge pivot is invariant, world axis is the parent's z axis rotated.
    const Vec3 pivotInChild = joint.restPose().inverse().act(joint.pivot());
    EXPECT_LT((pair.childPose().act(pivotInChild) - pair.pivotWorld()).norm(), kTol);
    EXPECT_LT((pair.axisWorld() - samplePose().rotation().act(Vec3(0, 0, 1))).norm(), kTol);

    // A child-frame x axis ends up at +90 deg about the world axis relative to q = 0.
    pair.setAngle(0.0);
    const Vec3 ex0 = pair.childPose().rotation().act(Vec3(1, 0, 0));
    pair.setAngle(M_PI / 2);
    const Vec3 ex1 = pair.childPose().rotation().act(Vec3(1, 0, 0));
    const Vec3 expected = Eigen::AngleAxisd(M_PI / 2, pair.axisWorld()) * ex0;
    EXPECT_LT((ex1 - expected).norm(), kTol);

    // Moving the parent drags the child rigidly (relative pose independent of the parent pose).
    const SE3d rel = pair.parentPose().inverse() * pair.childPose();
    pair.setParentPose(SE3d(SO3d(-1.1, Vec3(0, 1, 0)), Vec3(5, 6, 7)));
    EXPECT_LT(poseError(pair.parentPose().inverse() * pair.childPose(), rel), kTol);
}
