#include <gtest/gtest.h>
#include <liegroups/Twist.h>
#include <liegroups/Wrench.h>
#include <liegroups/SE3.h>

using namespace sofa::component::cosserat::liegroups;

// ── Twist tests ───────────────────────────────────────────────────────────────

TEST(TwistTest, DefaultConstructorIsZero) {
    Twistd xi;
    EXPECT_TRUE(xi.isZero());
    EXPECT_TRUE(xi.angular().isZero());
    EXPECT_TRUE(xi.linear().isZero());
}

TEST(TwistTest, ConstructFromParts) {
    Eigen::Vector3d phi(1.0, 2.0, 3.0);
    Eigen::Vector3d rho(4.0, 5.0, 6.0);
    Twistd xi(phi, rho);

    EXPECT_TRUE(xi.angular().isApprox(phi));
    EXPECT_TRUE(xi.linear().isApprox(rho));
}

TEST(TwistTest, HatOperatorStructure) {
    // hat(xi) must be a 4x4 se(3) matrix:
    //   top-left 3x3 = skew(phi), top-right col = rho, bottom row = 0
    Eigen::Vector3d phi(0.0, 0.0, 1.0);  // rotation around Z
    Eigen::Vector3d rho(1.0, 0.0, 0.0);  // translation along X
    Twistd xi(phi, rho);

    Eigen::Matrix4d Xi = xi.hat();

    // skew(phi) for phi=[0,0,1]: [[0,-1,0],[1,0,0],[0,0,0]]
    EXPECT_DOUBLE_EQ(Xi(0, 1), -1.0);
    EXPECT_DOUBLE_EQ(Xi(1, 0),  1.0);
    EXPECT_DOUBLE_EQ(Xi(0, 2),  0.0);

    // linear part in last column
    EXPECT_DOUBLE_EQ(Xi(0, 3), rho[0]);
    EXPECT_DOUBLE_EQ(Xi(1, 3), rho[1]);
    EXPECT_DOUBLE_EQ(Xi(2, 3), rho[2]);

    // bottom row must be zero
    EXPECT_DOUBLE_EQ(Xi(3, 0), 0.0);
    EXPECT_DOUBLE_EQ(Xi(3, 3), 0.0);
}

TEST(TwistTest, HatIsAntisymmetricInAngularBlock) {
    Twistd xi(Eigen::Vector3d(1.0, 2.0, 3.0), Eigen::Vector3d(0.0, 0.0, 0.0));
    Eigen::Matrix4d Xi = xi.hat();
    // Angular block must be skew-symmetric
    Eigen::Matrix3d ang = Xi.topLeftCorner<3, 3>();
    EXPECT_TRUE(ang.isApprox(-ang.transpose()));
}

TEST(TwistTest, LieBracketAnticommutativity) {
    Twistd xi(Eigen::Vector3d(1.0, 0.0, 0.0), Eigen::Vector3d(0.0, 1.0, 0.0));
    Twistd eta(Eigen::Vector3d(0.0, 1.0, 0.0), Eigen::Vector3d(1.0, 0.0, 0.0));

    Twistd bracket_xi_eta  = xi.bracket(eta);
    Twistd bracket_eta_xi  = eta.bracket(xi);

    // [xi, eta] = -[eta, xi]
    EXPECT_TRUE(bracket_xi_eta.isApprox(-bracket_eta_xi));
}

TEST(TwistTest, LieBracketPureRotation) {
    // For pure rotations: [phi1, phi2] angular = phi1 x phi2, linear = 0
    Eigen::Vector3d e1(1.0, 0.0, 0.0), e2(0.0, 1.0, 0.0), e3(0.0, 0.0, 1.0);
    Twistd xi(e1, Eigen::Vector3d::Zero());
    Twistd eta(e2, Eigen::Vector3d::Zero());

    Twistd bracket = xi.bracket(eta);

    EXPECT_TRUE(bracket.angular().isApprox(e1.cross(e2)));  // = e3
    EXPECT_TRUE(bracket.linear().isZero());
}

TEST(TwistTest, Arithmetic) {
    Twistd a(Eigen::Vector3d(1, 2, 3), Eigen::Vector3d(4, 5, 6));
    Twistd b(Eigen::Vector3d(1, 0, 0), Eigen::Vector3d(0, 1, 0));

    Twistd sum = a + b;
    EXPECT_TRUE(sum.angular().isApprox(Eigen::Vector3d(2, 2, 3)));

    Twistd scaled = a * 2.0;
    EXPECT_TRUE(scaled.angular().isApprox(Eigen::Vector3d(2, 4, 6)));

    Twistd neg = -a;
    EXPECT_TRUE(neg.angular().isApprox(Eigen::Vector3d(-1, -2, -3)));
}

TEST(TwistTest, TransformByAdjoint) {
    // Identity transform: Ad_I = I, so twist is unchanged
    SE3d g = SE3d::Identity();
    Twistd xi(Eigen::Vector3d(1.0, 0.0, 0.0), Eigen::Vector3d(0.0, 1.0, 0.0));

    Twistd xi_transformed = xi.transformBy(g.adjoint());
    EXPECT_TRUE(xi_transformed.isApprox(xi));
}

TEST(TwistTest, TransformByPureRotation) {
    // Pure rotation by 90° around Z: e_x -> e_y for the angular part
    double angle = M_PI / 2.0;
    Eigen::AngleAxisd rot(angle, Eigen::Vector3d::UnitZ());
    SE3d g(SO3d(rot.toRotationMatrix()), Eigen::Vector3d::Zero());

    // Twist purely angular along X
    Twistd xi(Eigen::Vector3d(1.0, 0.0, 0.0), Eigen::Vector3d::Zero());
    Twistd xi_b = xi.transformBy(g.adjoint());

    // After 90° Z rotation, X axis maps to Y
    EXPECT_NEAR(xi_b.angular()[0], 0.0, 1e-10);
    EXPECT_NEAR(xi_b.angular()[1], 1.0, 1e-10);
    EXPECT_NEAR(xi_b.angular()[2], 0.0, 1e-10);
}

// ── Wrench tests ──────────────────────────────────────────────────────────────

TEST(WrenchTest, DefaultConstructorIsZero) {
    Wrenchd w;
    EXPECT_TRUE(w.isZero());
    EXPECT_TRUE(w.torque().isZero());
    EXPECT_TRUE(w.force().isZero());
}

TEST(WrenchTest, ConstructFromParts) {
    Eigen::Vector3d tau(1.0, 2.0, 3.0);
    Eigen::Vector3d f(4.0, 5.0, 6.0);
    Wrenchd w(tau, f);

    EXPECT_TRUE(w.torque().isApprox(tau));
    EXPECT_TRUE(w.force().isApprox(f));
}

TEST(WrenchTest, Arithmetic) {
    Wrenchd a(Eigen::Vector3d(1, 2, 3), Eigen::Vector3d(4, 5, 6));
    Wrenchd b(Eigen::Vector3d(0, 1, 0), Eigen::Vector3d(1, 0, 0));

    Wrenchd sum = a + b;
    EXPECT_TRUE(sum.torque().isApprox(Eigen::Vector3d(1, 3, 3)));

    Wrenchd scaled = a * 3.0;
    EXPECT_TRUE(scaled.force().isApprox(Eigen::Vector3d(12, 15, 18)));
}

// ── Duality pairing tests ─────────────────────────────────────────────────────

TEST(DualityTest, VirtualPowerOrthogonality) {
    // Power of a pure torque on a pure translation twist is zero
    Wrenchd w(Eigen::Vector3d(1.0, 0.0, 0.0), Eigen::Vector3d::Zero());  // pure torque
    Twistd xi(Eigen::Vector3d::Zero(), Eigen::Vector3d(1.0, 0.0, 0.0));  // pure translation

    EXPECT_DOUBLE_EQ(w.dot(xi), 0.0);
}

TEST(DualityTest, VirtualPowerValue) {
    // P = tau.phi + f.rho
    Wrenchd w(Eigen::Vector3d(1.0, 2.0, 3.0), Eigen::Vector3d(4.0, 5.0, 6.0));
    Twistd xi(Eigen::Vector3d(1.0, 1.0, 1.0), Eigen::Vector3d(1.0, 1.0, 1.0));

    // P = (1+2+3) + (4+5+6) = 6 + 15 = 21
    EXPECT_DOUBLE_EQ(w.dot(xi), 21.0);
    EXPECT_DOUBLE_EQ(virtualPower(w, xi), 21.0);
}

TEST(DualityTest, PowerPreservedUnderCoAdjointTransport) {
    // Key property: ⟨w_b, ξ_b⟩ = ⟨w_a, ξ_a⟩
    // where ξ_b = Ad * ξ_a  and  w_a = Ad^T * w_b
    SE3d g(SO3d(Eigen::AngleAxisd(0.5, Eigen::Vector3d::UnitY()).toRotationMatrix()),
           Eigen::Vector3d(1.0, 2.0, 3.0));

    Eigen::Matrix<double, 6, 6> Ad = g.adjoint();

    Wrenchd w_b(Eigen::Vector3d(1.0, 2.0, 0.5), Eigen::Vector3d(3.0, 1.0, 2.0));
    Twistd  xi_a(Eigen::Vector3d(0.1, 0.2, 0.3), Eigen::Vector3d(0.4, 0.5, 0.6));

    // Transport: ξ_b = Ad * ξ_a
    Twistd xi_b = xi_a.transformBy(Ad);

    // Co-adjoint transport: w_a = Ad^T * w_b
    Wrenchd w_a = w_b.transformBy(Ad);

    double power_a = w_a.dot(xi_a);
    double power_b = w_b.dot(xi_b);

    EXPECT_NEAR(power_a, power_b, 1e-10);
}

TEST(WrenchTest, TransformByIdentity) {
    SE3d g = SE3d::Identity();
    Wrenchd w(Eigen::Vector3d(1, 2, 3), Eigen::Vector3d(4, 5, 6));

    Wrenchd w_transported = w.transformBy(g.adjoint());
    EXPECT_TRUE(w_transported.isApprox(w));
}

// ── SmallAdjoint tests ────────────────────────────────────────────────────────

TEST(TwistTest, SmallAdjointConsistentWithBracket) {
    // ad_xi * eta = [xi, eta]
    Twistd xi(Eigen::Vector3d(1.0, 0.5, 0.3), Eigen::Vector3d(0.1, 0.2, 0.4));
    Twistd eta(Eigen::Vector3d(0.2, 0.7, 0.1), Eigen::Vector3d(0.3, 0.1, 0.5));

    Twistd bracket_direct = xi.bracket(eta);
    Twistd bracket_via_ad = Twistd(xi.smallAdjoint() * eta.data());

    EXPECT_TRUE(bracket_direct.isApprox(bracket_via_ad, 1e-10));
}
