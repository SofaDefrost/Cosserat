// Tests of the SO3 members defined in SO3.inl: exp-map Jacobians, dexp/dlog, BCH,
// distance, interpolate and computeAd. Also checks that the SE3 Jacobians, which
// rely on SO3::dexp, link and agree with finite differences.

#include <gtest/gtest.h>
#include <liegroups/SE3.h>
#include <liegroups/SO3.h>

#include <cmath>
#include <random>

using namespace sofa::component::cosserat::liegroups;

namespace {

using Vec3 = Eigen::Vector3d;
using Mat3 = Eigen::Matrix3d;
using Rot = SO3d;

Vec3 randomRotationVector(std::mt19937_64 &rng, double angle) {
    std::uniform_real_distribution<double> u(-1.0, 1.0);
    return Vec3(u(rng), u(rng), u(rng)).normalized() * angle;
}

const double kAngles[] = {1e-9, 1e-6, 1e-3, 0.1, 1.0, 2.0, 3.0};

} // namespace

TEST(SO3Jacobians, RightAndLeftAreInverseOfTheirInverse) {
    std::mt19937_64 rng(1);
    for (double angle : kAngles) {
        const Vec3 w = randomRotationVector(rng, angle);
        EXPECT_LT((Rot::computeRightJacobian(w) * Rot::computeRightJacobianInverse(w) - Mat3::Identity()).norm(), 1e-12) << angle;
        EXPECT_LT((Rot::computeLeftJacobian(w) * Rot::computeLeftJacobianInverse(w) - Mat3::Identity()).norm(), 1e-12) << angle;
        EXPECT_LT((Rot::computeLeftJacobian(w) - Rot::computeRightJacobian(-w)).norm(), 1e-15) << angle;
    }
}

TEST(SO3Jacobians, MatchFiniteDifferencesOfExp) {
    // exp(w + d) = exp(w) exp(Jr(w) d) = exp(Jl(w) d) exp(w) to first order in d.
    std::mt19937_64 rng(2);
    for (double angle : kAngles) {
        const Vec3 w = randomRotationVector(rng, angle);
        const Vec3 d = randomRotationVector(rng, 1e-6);
        const Mat3 target = Rot::exp(w + d).matrix();
        const Mat3 viaRight = (Rot::exp(w) * Rot::exp(Rot::computeRightJacobian(w) * d)).matrix();
        const Mat3 viaLeft = (Rot::exp(Rot::computeLeftJacobian(w) * d) * Rot::exp(w)).matrix();
        EXPECT_LT((target - viaRight).norm(), 1e-11) << angle;
        EXPECT_LT((target - viaLeft).norm(), 1e-11) << angle;
    }
}

TEST(SO3Jacobians, SmallAngleHasNoCancellation) {
    // Jr(w) = I - [w]x/2 + O(|w|^2): the closed form (1-cos t)/t^2 is useless below ~1e-5.
    const Vec3 w(1e-9, -2e-9, 0.5e-9);
    const Mat3 expected = Mat3::Identity() - 0.5 * Rot::computeHat(w);
    EXPECT_LT((Rot::computeRightJacobian(w) - expected).norm(), 1e-17);
}

TEST(SO3Jacobians, DexpIsLeftJacobianAndDlogItsInverse) {
    std::mt19937_64 rng(3);
    for (double angle : kAngles) {
        const Vec3 w = randomRotationVector(rng, angle);
        EXPECT_LT((Rot::dexp(w) - Rot::computeLeftJacobian(w)).norm(), 1e-15) << angle;
        EXPECT_LT((Rot::exp(w).dlog() - Rot::computeLeftJacobianInverse(w)).norm(), 1e-9) << angle;
    }
}

TEST(SO3Jacobians, ComputeAdIsHatAndDoesNotRecurse) {
    const Vec3 v(0.3, -0.2, 0.7);
    EXPECT_LT((Rot::computeAd(v) - Rot::computeHat(v)).norm(), 1e-15);
    EXPECT_LT((Rot::ad(v) - Rot::computeHat(v)).norm(), 1e-15);
}

TEST(SO3Jacobians, BchThirdOrderErrorIsFourthOrder) {
    std::mt19937_64 rng(4);
    double previous = 0.0;
    for (double scale : {0.1, 0.05, 0.025}) {
        double worst = 0.0;
        for (int i = 0; i < 200; ++i) {
            const Vec3 v = randomRotationVector(rng, scale), w = randomRotationVector(rng, scale);
            const Vec3 exact = (Rot::exp(v) * Rot::exp(w)).log();
            worst = std::max(worst, (Rot::BCH(v, w, 3) - exact).norm());
            // The third-order term must improve on the second order one.
            EXPECT_LT((Rot::BCH(v, w, 3) - exact).norm(), (Rot::BCH(v, w, 2) - exact).norm() + 1e-15);
        }
        if (previous > 0.0) EXPECT_GT(previous / worst, 10.0) << "error should drop ~16x when halving the angles";
        previous = worst;
    }
}

TEST(SO3Jacobians, DistanceAndInterpolateAreGeodesic) {
    std::mt19937_64 rng(5);
    for (int i = 0; i < 100; ++i) {
        const Rot a = Rot::exp(randomRotationVector(rng, 1.5)), b = Rot::exp(randomRotationVector(rng, 1.2));
        const double t = 0.37, d = a.distance(b);
        EXPECT_NEAR(d, (a.inverse() * b).log().norm(), 1e-12);
        const Rot m = a.interpolate(b, t);
        EXPECT_NEAR(a.distance(m), t * d, 1e-12);
        EXPECT_NEAR(m.distance(b), (1.0 - t) * d, 1e-12);
    }
}

TEST(SE3Jacobians, ExpMapJacobiansLinkAndMatchFiniteDifferences) {
    std::mt19937_64 rng(6);
    for (double angle : {1e-6, 0.1, 1.0, 2.5}) {
        Eigen::Matrix<double, 6, 1> xi, d;
        xi << randomRotationVector(rng, angle), Vec3(0.4, -0.7, 0.2);
        d << randomRotationVector(rng, 1e-6), randomRotationVector(rng, 1e-6);
        const SE3d target = SE3d::exp(xi + d);
        const SE3d viaLeft = SE3d::exp(SE3d::computeLeftJacobian(xi) * d) * SE3d::exp(xi);
        const SE3d viaRight = SE3d::exp(xi) * SE3d::exp(SE3d::computeRightJacobian(xi) * d);
        EXPECT_LT((target.inverse() * viaLeft).log().norm(), 1e-11) << angle;
        EXPECT_LT((target.inverse() * viaRight).log().norm(), 1e-11) << angle;
    }
}
