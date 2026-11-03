// Tests of SphericalJoint (3 DoF) and TwoAxisJoint (2 DoF: condylar / saddle) on SE3.
// References are built with Eigen::AngleAxis only (independent of SE3::exp).

#include <gtest/gtest.h>
#include <liegroups/MultiDofJoints.h>

#include <Eigen/Geometry>
#include <Eigen/SVD>
#include <cmath>

using namespace sofa::component::cosserat::liegroups;

namespace {

using Vec3 = Eigen::Vector3d;
using Vec2 = Eigen::Vector2d;

// Rotation by q about the line (point p, axis u): x -> p + R (x - p).
SE3d rotationAboutLine(const Vec3 &u, const Vec3 &p, double q) {
    const Eigen::Matrix3d R = Eigen::AngleAxisd(q, u.normalized()).toRotationMatrix();
    return SE3d(SO3d(R), p - R * p);
}

SE3d rotationVectorAbout(const Vec3 &w, const Vec3 &p) {
    const double n = w.norm();
    const Eigen::Matrix3d R = n < 1e-300 ? Eigen::Matrix3d::Identity() : Eigen::AngleAxisd(n, w / n).toRotationMatrix();
    return SE3d(SO3d(R), p - R * p);
}

SE3d samplePose() { return SE3d(SO3d(0.7, Vec3(1, 2, -1)), Vec3(0.3, -0.4, 1.1)); }
SE3d sampleRest() { return SE3d(SO3d(0.5, Vec3(0, 0, 1)), Vec3(0.0, 0.0, 1.0)); }
double poseError(const SE3d &a, const SE3d &b) { return (a.matrix() - b.matrix()).norm(); }

// Singular values of the 6xN matrix of body-twist finite differences of the relative pose.
template <typename Joint, typename Q>
Eigen::VectorXd twistSingularValues(const Joint &joint, const Q &q, int n) {
    const double h = 1e-6;
    Eigen::MatrixXd J(6, n);
    for (int i = 0; i < n; ++i) {
        Q qp = q, qm = q;
        qp[i] += h;
        qm[i] -= h;
        J.col(i) = (joint.relativePose(qm).inverse() * joint.relativePose(qp)).log() / (2 * h);
    }
    return Eigen::JacobiSVD<Eigen::MatrixXd>(J).singularValues();
}

} // namespace

// ---------------------------------------------------------------- SphericalJoint

TEST(SphericalJointTest, DofAndRestPose) {
    const SphericalJoint<double> j(Vec3(0.1, 0.2, 0.3), sampleRest());
    EXPECT_EQ(j.dof(), 3);
    EXPECT_LT(poseError(j.relativePose(Vec3::Zero()), sampleRest()), 1e-14);
}

TEST(SphericalJointTest, MatchesEigenRotationAboutPoint) {
    const Vec3 pivot(0.4, -0.2, 0.9);
    const SphericalJoint<double> j(pivot, sampleRest());
    for (const Vec3 &w : {Vec3(0.3, -0.5, 0.8), Vec3(2.0, 1.0, -1.0), Vec3(1e-6, 0, 2e-6), Vec3(0, 0, -2.9), Vec3(-1.2, 0.4, 0.1)}) {
        EXPECT_LT(poseError(j.relativePose(w), rotationVectorAbout(w, pivot) * sampleRest()), 1e-11) << w.transpose();
    }
}

TEST(SphericalJointTest, PivotFixedForBothBodies) {
    const Vec3 pivot(0.4, -0.2, 0.9);
    const SE3d rest = sampleRest();
    const SphericalJoint<double> j(pivot, rest);
    const SE3d parent = samplePose();
    const Vec3 pivotInChild = rest.inverse().act(pivot);
    for (const Vec3 &w : {Vec3(0.3, -0.5, 0.8), Vec3(-2.0, 0.5, 1.0), Vec3(0, 1.0, 0)}) {
        EXPECT_LT((j.childPose(parent, w).act(pivotInChild) - parent.act(pivot)).norm(), 1e-12);
    }
}

TEST(SphericalJointTest, ThreeDegreesOfFreedom) {
    const SphericalJoint<double> j(Vec3(0.1, 0.2, 0.3), sampleRest());
    const Eigen::VectorXd sv = twistSingularValues(j, Vec3(0.4, -0.7, 0.9), 3);
    EXPECT_EQ(sv.size(), 3);
    EXPECT_GT(sv.minCoeff(), 0.3); // three independent motions
}

TEST(SphericalJointTest, JointAnglesRoundTripAndResidual) {
    const SphericalJoint<double> j(Vec3(0.4, -0.2, 0.9), sampleRest());
    for (const Vec3 &w : {Vec3(0.3, -0.5, 0.8), Vec3(0, 0, 0), Vec3(-1.0, 2.0, 0.5)}) {
        const SE3d rel = j.relativePose(w);
        EXPECT_LT((j.jointAngles(rel) - w).norm(), 1e-10);
        EXPECT_LT(j.constraintResidual(rel).norm(), 1e-10);
    }
    // Translating the child away from the pivot is not reachable.
    SE3d bad = j.relativePose(Vec3(0.3, 0, 0));
    bad.translation() += Vec3(0, 0, 0.2);
    EXPECT_GT(j.constraintResidual(bad).norm(), 0.1);
}

TEST(SphericalJointTest, RotationsAboutTheSameAxisAdd) {
    const SphericalJoint<double> j(Vec3(0.4, -0.2, 0.9));
    const Vec3 u = Vec3(1, 2, -1).normalized();
    EXPECT_LT(poseError(j.relativePose(0.5 * u) * j.relativePose(0.9 * u), j.relativePose(1.4 * u)), 1e-11);
}

// ---------------------------------------------------------------- TwoAxisJoint

TEST(TwoAxisJointTest, DofAndRestPose) {
    const TwoAxisJoint<double> j(Vec3(1, 0, 0), Vec3::Zero(), Vec3(0, 1, 0), Vec3::Zero(), sampleRest());
    EXPECT_EQ(j.dof(), 2);
    EXPECT_LT(poseError(j.relativePose(Vec2::Zero()), sampleRest()), 1e-14);
}

TEST(TwoAxisJointTest, CondylarMatchesEigenAndPointStaysFixed) {
    const Vec3 c(0.2, -0.3, 0.5); // common point of the two axes
    const Vec3 a1(1, 0, 0), a2(0, 1, 0);
    const SE3d rest = sampleRest();
    const TwoAxisJoint<double> j(a1, c, a2, c, rest);
    const SE3d parent = samplePose();
    const Vec3 cInChild = rest.inverse().act(c);

    for (const Vec2 &q : {Vec2(0.4, -0.7), Vec2(1.3, 0.9), Vec2(-2.0, 0.2), Vec2(0, 1.1)}) {
        const SE3d expected = rotationAboutLine(a1, c, q[0]) * rotationAboutLine(a2, c, q[1]) * rest;
        EXPECT_LT(poseError(j.relativePose(q), expected), 1e-11) << q.transpose();
        // Concurrent axes: the shared point is a fixed point of the whole articulation.
        EXPECT_LT((j.childPose(parent, q).act(cInChild) - parent.act(c)).norm(), 1e-12);
        // Rotation is exactly Rx(q1) Ry(q2) R_rest.
        const Eigen::Matrix3d R = (Eigen::AngleAxisd(q[0], a1) * Eigen::AngleAxisd(q[1], a2)).toRotationMatrix();
        EXPECT_LT((j.relativePose(q).rotation().matrix() - R * rest.rotation().matrix()).norm(), 1e-11);
    }
}

TEST(TwoAxisJointTest, SaddleOffsetAxesMatchEigenAndDifferFromCondylar) {
    // Axis 1 (x) through the origin, axis 2 (y) offset by L along the bone (z): perpendicular, not concurrent.
    const double L = 0.15;
    const Vec3 a1(1, 0, 0), a2(0, 1, 0), p1(0, 0, 0), p2(0, 0, L);
    const SE3d rest = sampleRest();
    const TwoAxisJoint<double> saddle(a1, p1, a2, p2, rest);
    const TwoAxisJoint<double> condylar(a1, p1, a2, p1, rest);

    const Vec2 q(0.6, -0.8);
    const SE3d expected = rotationAboutLine(a1, p1, q[0]) * rotationAboutLine(a2, p2, q[1]) * rest;
    EXPECT_LT(poseError(saddle.relativePose(q), expected), 1e-11);
    EXPECT_GT(poseError(saddle.relativePose(q), condylar.relativePose(q)), 1e-3);

    // With q2 = 0 only axis 1 acts: p1 is fixed. With q1 = 0 only axis 2 acts: p2 is fixed.
    const Vec3 p1InChild = rest.inverse().act(p1), p2InChild = rest.inverse().act(p2);
    EXPECT_LT((saddle.relativePose(Vec2(0.9, 0)).act(p1InChild) - p1).norm(), 1e-12);
    EXPECT_LT((saddle.relativePose(Vec2(0, -1.1)).act(p2InChild) - p2).norm(), 1e-12);
    // For q1 != 0 and q2 != 0 the offset axis does not keep p2 fixed any more.
    EXPECT_GT((saddle.relativePose(q).act(p2InChild) - p2).norm(), 1e-3);
}

TEST(TwoAxisJointTest, TwoDegreesOfFreedom) {
    const TwoAxisJoint<double> saddle(Vec3(1, 0, 0), Vec3(0, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 0.15), sampleRest());
    const TwoAxisJoint<double> condylar(Vec3(1, 0, 0), Vec3(0.1, 0, 0.2), Vec3(0, 1, 0), Vec3(0.1, 0, 0.2), sampleRest());
    for (const Vec2 &q : {Vec2(0.4, -0.7), Vec2(-1.0, 0.5)}) {
        EXPECT_GT(twistSingularValues(saddle, q, 2).minCoeff(), 0.3);
        EXPECT_GT(twistSingularValues(condylar, q, 2).minCoeff(), 0.3);
    }
    // No axial rotation: a pure rotation about the third axis (z) is not reachable (concurrent case).
    const Vec3 c(0, 0, 0);
    const TwoAxisJoint<double> k(Vec3(1, 0, 0), c, Vec3(0, 1, 0), c, SE3d());
    const SE3d axial = rotationAboutLine(Vec3(0, 0, 1), c, 0.5);
    double best = 1e9;
    for (double a = -1.5; a <= 1.5; a += 0.01)
        for (double b = -1.5; b <= 1.5; b += 0.01)
            best = std::min(best, poseError(k.relativePose(Vec2(a, b)), axial));
    EXPECT_GT(best, 0.2);
}

TEST(TwoAxisJointTest, ChildPoseFollowsParentRigidly) {
    const TwoAxisJoint<double> j(Vec3(1, 0, 0), Vec3(0, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 0.15), sampleRest());
    const Vec2 q(0.4, 0.9);
    const SE3d p1 = samplePose(), p2(SO3d(-1.1, Vec3(0, 1, 0)), Vec3(5, 6, 7));
    EXPECT_LT(poseError(p1.inverse() * j.childPose(p1, q), p2.inverse() * j.childPose(p2, q)), 1e-12);
}
