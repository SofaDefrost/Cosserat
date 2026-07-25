#include <gtest/gtest.h>
#include <liegroups/CosseratBodyJacobian.h>
#include <liegroups/CosseratUncertaintyPropagator.h>
#include <liegroups/BezierSE3.h>
#include <liegroups/SE3.h>
#include <liegroups/Twist.h>
#include <liegroups/Wrench.h>
#include <cmath>

using namespace sofa::component::cosserat::liegroups;

using SE3d    = SE3<double>;
using Twistd  = Twist<double>;
using Wrenchd = Wrench<double>;
using Jacobian = CosseratBodyJacobian<double>;
using Propagator = CosseratUncertaintyPropagator<double>;
using GaussianSE3 = GaussianOnManifold<SE3d>;

// ─────────────────────────────────────────────────────────────────────────────
//  Helpers
// ─────────────────────────────────────────────────────────────────────────────

static SE3d makeSectionSE3(double kappa, double L) {
    Eigen::Matrix<double,6,1> strain = Eigen::Matrix<double,6,1>::Zero();
    strain[0] = kappa;  // torsion
    return SE3d::expCosserat(strain, L);
}

// Compute tangent-exp matrix via the same formula used in CosseratGeometryMapping
static Eigen::Matrix<double,6,6> computeTangExp(double L, const Eigen::Matrix<double,6,1>& strain) {
    // Reuse CosseratUncertaintyPropagator's private formula via a Section
    Propagator::Section sec;
    sec.strain = strain;
    sec.length = L;
    sec.strain_cov = Eigen::Matrix<double,6,6>::Zero();

    GaussianSE3 zero_state(SE3d::Identity(), Eigen::Matrix<double,6,6>::Zero());
    GaussianSE3 result = Propagator::propagateStep(zero_state, sec);
    // The tangent-exp is embedded in the covariance propagation. We instead
    // reconstruct it via finite differences of the propagation.
    // For unit tests, just check the formula via forward kinematics.
    (void)result;
    return Eigen::Matrix<double,6,6>::Identity() * L;  // placeholder
}

// ─────────────────────────────────────────────────────────────────────────────
//  CosseratBodyJacobian tests
// ─────────────────────────────────────────────────────────────────────────────

TEST(BodyJacobianTest, SingleStraightSection) {
    // Straight rod: strain = 0, local transform = Identity
    Eigen::Matrix<double,6,1> strain = Eigen::Matrix<double,6,1>::Zero();
    double L = 0.1;
    SE3d g = SE3d::expCosserat(strain, L);
    Eigen::Matrix<double,6,6> J_loc = Eigen::Matrix<double,6,6>::Identity() * L;

    Jacobian J(1);
    J.pushSection(g, J_loc);

    EXPECT_EQ(J.size(), 1);
}

TEST(BodyJacobianTest, ForwardPassIdentityrod) {
    // Straight rod with N sections, zero strain
    int N = 3;
    double L = 0.1;
    Eigen::Matrix<double,6,1> strain = Eigen::Matrix<double,6,1>::Zero();

    Jacobian J(N);
    for (int k = 0; k < N; ++k) {
        SE3d g = SE3d::expCosserat(strain, L);
        Eigen::Matrix<double,6,6> J_loc = Eigen::Matrix<double,6,6>::Identity() * L;
        J.pushSection(g, J_loc);
    }

    std::vector<Twistd> strain_rates(N, Twistd::Zero());
    auto twists = J.applyForward(strain_rates, Twistd::Zero());

    ASSERT_EQ(static_cast<int>(twists.size()), N + 1);
    for (const auto& xi : twists)
        EXPECT_TRUE(xi.isZero(1e-12));
}

TEST(BodyJacobianTest, ForwardPassWrongSizeThrows) {
    Jacobian J(2);
    SE3d g = SE3d::Identity();
    Eigen::Matrix<double,6,6> J_loc = Eigen::Matrix<double,6,6>::Identity();
    J.pushSection(g, J_loc);
    J.pushSection(g, J_loc);

    std::vector<Twistd> wrong_size(1, Twistd::Zero());  // wrong: should be 2
    EXPECT_THROW(J.applyForward(wrong_size), std::invalid_argument);
}

TEST(BodyJacobianTest, TransposeWrongSizeThrows) {
    Jacobian J(2);
    SE3d g = SE3d::Identity();
    Eigen::Matrix<double,6,6> J_loc = Eigen::Matrix<double,6,6>::Identity();
    J.pushSection(g, J_loc);
    J.pushSection(g, J_loc);

    std::vector<Wrenchd> wrong_size(2, Wrenchd::Zero());  // wrong: should be 3
    Wrenchd base_out;
    EXPECT_THROW(J.applyTranspose(wrong_size, base_out), std::invalid_argument);
}

TEST(BodyJacobianTest, VirtualPowerDualityForwardTranspose) {
    // Key property: <J * xi_dot, w> = <xi_dot, J^T * w>
    int N = 3;
    double L = 0.05;
    Eigen::Matrix<double,6,1> strain;
    strain << 0.1, 0.05, 0.0, 0.02, 0.0, 0.0;

    Jacobian J(N);
    for (int k = 0; k < N; ++k) {
        SE3d g = SE3d::expCosserat(strain, L);
        Eigen::Matrix<double,6,6> J_loc = Eigen::Matrix<double,6,6>::Identity() * L;
        J.pushSection(g, J_loc);
    }

    // Random strain rates and wrenches
    std::vector<Twistd> xi_dot(N);
    for (int k = 0; k < N; ++k)
        xi_dot[k] = Twistd(Eigen::Matrix<double,6,1>::Random());

    std::vector<Wrenchd> wrenches(N + 1);
    for (int k = 0; k <= N; ++k)
        wrenches[k] = Wrenchd(Eigen::Matrix<double,6,1>::Random());

    // Forward: twists
    auto twists = J.applyForward(xi_dot, Twistd::Zero());

    // Transpose: strain forces
    Wrenchd base_w;
    auto q = J.applyTranspose(wrenches, base_w);

    // Virtual power from the forward side: sum_k <w_k, eta_k>
    double power_fwd = 0.0;
    for (int k = 0; k <= N; ++k)
        power_fwd += wrenches[k].dot(twists[k]);

    // Virtual power from the transpose side: sum_k <q_k, xi_dot_k>
    double power_T = base_w.dot(Twistd::Zero());  // base strain rate = 0
    for (int k = 0; k < N; ++k)
        power_T += q[k].dot(xi_dot[k]);

    EXPECT_NEAR(power_fwd, power_T, 1e-10);
}

TEST(BodyJacobianTest, ExplicitJacobianBlockTriangular) {
    // For a rod with N sections, J_body(k) must have zeros in columns > k
    int N = 3;
    double L = 0.1;
    Eigen::Matrix<double,6,1> strain;
    strain << 0.1, 0.0, 0.05, 0.01, 0.0, 0.0;

    Jacobian J(N);
    for (int k = 0; k < N; ++k) {
        SE3d g = SE3d::expCosserat(strain, L);
        Eigen::Matrix<double,6,6> J_loc = Eigen::Matrix<double,6,6>::Identity() * L;
        J.pushSection(g, J_loc);
    }

    for (int k = 0; k < N; ++k) {
        Eigen::MatrixXd Jk = J.getJacobianAtSection(k);
        ASSERT_EQ(Jk.rows(), 6);
        ASSERT_EQ(Jk.cols(), 6 * N);

        // Columns after block k must be zero
        for (int j = k + 1; j < N; ++j)
            EXPECT_NEAR(Jk.block(0, j*6, 6, 6).norm(), 0.0, 1e-14)
                << "Non-zero at section k=" << k << " col block j=" << j;
    }
}

TEST(BodyJacobianTest, OutOfRangeJacobianThrows) {
    Jacobian J(1);
    J.pushSection(SE3d::Identity(), Eigen::Matrix<double,6,6>::Identity());
    EXPECT_THROW(J.getJacobianAtSection(-1), std::out_of_range);
    EXPECT_THROW(J.getJacobianAtSection(1),  std::out_of_range);
}

// ─────────────────────────────────────────────────────────────────────────────
//  CosseratUncertaintyPropagator tests
// ─────────────────────────────────────────────────────────────────────────────

TEST(UncertaintyPropagatorTest, ZeroNoiseStraightRod) {
    // Zero strain + zero noise → identity propagation, covariance unchanged
    Propagator::Section sec;
    sec.strain     = Eigen::Matrix<double,6,1>::Zero();
    sec.length     = 0.1;
    sec.strain_cov = Eigen::Matrix<double,6,6>::Zero();

    Eigen::Matrix<double,6,6> Sigma0 = Eigen::Matrix<double,6,6>::Identity();
    GaussianSE3 prior(SE3d::Identity(), Sigma0);
    GaussianSE3 post = Propagator::propagateStep(prior, sec);

    // Mean should be identity (zero strain step)
    EXPECT_TRUE(post.getMean().isApprox(SE3d::expCosserat(sec.strain, sec.length)));
    // Covariance should remain Identity (no noise, Ad_{I} = I)
    EXPECT_TRUE(post.getCovariance().isApprox(Sigma0, 1e-12));
}

TEST(UncertaintyPropagatorTest, NoisePropagatesForward) {
    // Non-zero strain noise must increase uncertainty
    Propagator::Section sec;
    sec.strain     = Eigen::Matrix<double,6,1>::Zero();
    sec.length     = 0.1;
    sec.strain_cov = 0.01 * Eigen::Matrix<double,6,6>::Identity();

    GaussianSE3 prior(SE3d::Identity(), Eigen::Matrix<double,6,6>::Zero());
    GaussianSE3 post = Propagator::propagateStep(prior, sec);

    // Covariance must be positive definite (all eigenvalues > 0)
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,6,6>> eig(post.getCovariance());
    EXPECT_GT(eig.eigenvalues().minCoeff(), 0.0);
}

TEST(UncertaintyPropagatorTest, PropagateAlongRodSize) {
    int N = 5;
    std::vector<Propagator::Section> sections(N);
    for (auto& sec : sections) {
        sec.strain     = Eigen::Matrix<double,6,1>::Zero();
        sec.length     = 0.1;
        sec.strain_cov = 0.001 * Eigen::Matrix<double,6,6>::Identity();
    }

    GaussianSE3 base(SE3d::Identity(), Eigen::Matrix<double,6,6>::Zero());
    auto results = Propagator::propagateAlongRod(base, sections);

    ASSERT_EQ(static_cast<int>(results.size()), N + 1);
}

TEST(UncertaintyPropagatorTest, CovarianceMonotonicallyIncreases) {
    // Constant noise at each step → trace of covariance must increase
    int N = 5;
    std::vector<Propagator::Section> sections(N);
    for (auto& sec : sections) {
        sec.strain     = Eigen::Matrix<double,6,1>::Zero();
        sec.length     = 0.05;
        sec.strain_cov = 0.001 * Eigen::Matrix<double,6,6>::Identity();
    }

    GaussianSE3 base(SE3d::Identity(), Eigen::Matrix<double,6,6>::Zero());
    auto results = Propagator::propagateAlongRod(base, sections);

    double prev_trace = results[0].getCovariance().trace();
    for (int k = 1; k <= N; ++k) {
        double cur_trace = results[k].getCovariance().trace();
        EXPECT_GT(cur_trace, prev_trace - 1e-12);
        prev_trace = cur_trace;
    }
}

TEST(UncertaintyPropagatorTest, ConfidenceRadiiPositive) {
    Propagator::Section sec;
    sec.strain     = Eigen::Matrix<double,6,1>::Zero();
    sec.length     = 0.1;
    sec.strain_cov = 0.01 * Eigen::Matrix<double,6,6>::Identity();

    GaussianSE3 prior(SE3d::Identity(), Eigen::Matrix<double,6,6>::Zero());
    GaussianSE3 post = Propagator::propagateStep(prior, sec);

    auto radii = Propagator::tipConfidenceRadii(post);
    EXPECT_GT(radii.minCoeff(), 0.0);
}

// ─────────────────────────────────────────────────────────────────────────────
//  BezierSE3 tests
// ─────────────────────────────────────────────────────────────────────────────

TEST(BezierSE3Test, Degree1IsGeodesicInterpolation) {
    SE3d g0 = SE3d::Identity();
    SE3d g1(SO3<double>(Eigen::AngleAxisd(M_PI/4, Eigen::Vector3d::UnitZ()).toRotationMatrix()),
             Eigen::Vector3d(1.0, 0.0, 0.0));

    BezierSE3d curve({g0, g1});
    EXPECT_EQ(curve.degree(), 1);

    // t=0 → g0, t=1 → g1
    EXPECT_TRUE(curve.evaluate(0.0).isApprox(g0, 1e-10));
    EXPECT_TRUE(curve.evaluate(1.0).isApprox(g1, 1e-10));

    // t=0.5 should be the midpoint along the geodesic
    SE3d mid_expected = g0.compose(SE3d::computeExp(0.5 * g0.computeInverse().compose(g1).log()));
    EXPECT_TRUE(curve.evaluate(0.5).isApprox(mid_expected, 1e-10));
}

TEST(BezierSE3Test, EndpointsReachedForHigherDegree) {
    // Cubic Bézier: endpoints are always the first and last control poses
    std::vector<SE3d> ctrl(4);
    ctrl[0] = SE3d::Identity();
    ctrl[1] = SE3d(SO3<double>(Eigen::Matrix3d::Identity()), Eigen::Vector3d(0.5, 0.0, 0.0));
    ctrl[2] = SE3d(SO3<double>(Eigen::Matrix3d::Identity()), Eigen::Vector3d(0.5, 0.5, 0.0));
    ctrl[3] = SE3d(SO3<double>(Eigen::AngleAxisd(M_PI/2, Eigen::Vector3d::UnitZ()).toRotationMatrix()),
                   Eigen::Vector3d(1.0, 1.0, 0.0));

    BezierSE3d curve(ctrl);
    EXPECT_EQ(curve.degree(), 3);
    EXPECT_TRUE(curve.evaluate(0.0).isApprox(ctrl[0], 1e-10));
    EXPECT_TRUE(curve.evaluate(1.0).isApprox(ctrl[3], 1e-10));
}

TEST(BezierSE3Test, SampleCountCorrect) {
    BezierSE3d curve({SE3d::Identity(), SE3d::Identity()});
    auto samples = curve.sample(10);
    EXPECT_EQ(static_cast<int>(samples.size()), 11);
}

TEST(BezierSE3Test, TooFewControlPosesThrows) {
    EXPECT_THROW(BezierSE3d({SE3d::Identity()}), std::invalid_argument);
}

TEST(BezierSE3Test, SplitReconstructsCurve) {
    // After splitting at t=0.5, evaluating each sub-curve at its own endpoints
    // should reproduce curve(0), curve(0.5), curve(1).
    std::vector<SE3d> ctrl = {
        SE3d::Identity(),
        SE3d(SO3<double>(Eigen::Matrix3d::Identity()), Eigen::Vector3d(1.0, 0.0, 0.0)),
        SE3d(SO3<double>(Eigen::AngleAxisd(M_PI/3, Eigen::Vector3d::UnitY()).toRotationMatrix()),
             Eigen::Vector3d(1.0, 1.0, 0.0))
    };
    BezierSE3d curve(ctrl);

    auto [left, right] = curve.split(0.5);

    EXPECT_TRUE(left.evaluate(0.0).isApprox(curve.evaluate(0.0), 1e-10));
    EXPECT_TRUE(left.evaluate(1.0).isApprox(curve.evaluate(0.5), 1e-10));
    EXPECT_TRUE(right.evaluate(0.0).isApprox(curve.evaluate(0.5), 1e-10));
    EXPECT_TRUE(right.evaluate(1.0).isApprox(curve.evaluate(1.0), 1e-10));
}

TEST(BezierSE3Test, PiecewiseCubicPathInterpolatesEndpoints) {
    std::vector<SE3d> waypoints = {
        SE3d::Identity(),
        SE3d(SO3<double>(Eigen::Matrix3d::Identity()), Eigen::Vector3d(0.5, 0.0, 0.0)),
        SE3d(SO3<double>(Eigen::AngleAxisd(M_PI/4, Eigen::Vector3d::UnitZ()).toRotationMatrix()),
             Eigen::Vector3d(1.0, 0.5, 0.0))
    };

    auto segments = makePiecewiseCubicPath<double>(waypoints);
    ASSERT_EQ(static_cast<int>(segments.size()), 2);

    // Each segment's endpoints match consecutive waypoints
    for (int i = 0; i < 2; ++i) {
        EXPECT_TRUE(segments[i].evaluate(0.0).isApprox(waypoints[i],     1e-10));
        EXPECT_TRUE(segments[i].evaluate(1.0).isApprox(waypoints[i + 1], 1e-10));
    }
}

TEST(BezierSE3Test, ArcLengthPositive) {
    BezierSE3d curve({
        SE3d::Identity(),
        SE3d(SO3<double>(Eigen::Matrix3d::Identity()), Eigen::Vector3d(1.0, 0.0, 0.0))
    });
    double L = curve.arcLength(100);
    EXPECT_GT(L, 0.0);
}
