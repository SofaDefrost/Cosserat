// Checks that the pieces of src/liegroups that are meant to follow the Cosserat strain layout
//     Xi = [phi (torsion, bending_y, bending_z) ; rho (elongation, shear_y, shear_z)]
// really do so: exp/log/adjoint/expCosserat/Twist/Wrench/body Jacobian, and SE3::ad, computeBCH and
// parallelTransport and the action Jacobians actionJacobians/actionJacobianGroup (all fixed from the old
// [rho ; phi] layout). computeTangent and dlog are deliberately not covered: they still return a different
// layout / object (see the audit notes).
//
// References: 4x4 matrix exponential/commutator (Taylor series), finite differences.

#include <gtest/gtest.h>
#include <liegroups/CosseratBodyJacobian.h>
#include <liegroups/SE3.h>
#include <liegroups/Twist.h>
#include <liegroups/Wrench.h>

#include <random>

using namespace sofa::component::cosserat::liegroups;

namespace {

using V3 = Eigen::Vector3d;
using V6 = Eigen::Matrix<double, 6, 1>;
using M4 = Eigen::Matrix4d;
using M6 = Eigen::Matrix<double, 6, 6>;

std::mt19937 rng(7);
V6 rnd6(double s) {
    std::normal_distribution<double> n(0, 1);
    V6 v;
    for (int i = 0; i < 6; i++) v[i] = s * n(rng);
    return v;
}

M4 hat(const V6 &x) {
    M4 H = M4::Zero();
    H.block<3, 3>(0, 0) = SO3d::computeHat(x.head<3>());
    H.block<3, 1>(0, 3) = x.tail<3>();
    return H;
}
V6 vee(const M4 &H) {
    V6 x;
    x << SO3d::computeVee(H.block<3, 3>(0, 0)), H.block<3, 1>(0, 3);
    return x;
}
M4 expm(const M4 &A) {
    int s = 0;
    double n = A.norm();
    while (n > 0.25) { n /= 2; s++; }
    const M4 A2 = A / std::pow(2.0, s);
    M4 T = M4::Identity(), R = M4::Identity();
    for (int k = 1; k < 25; k++) { T = T * A2 / k; R += T; }
    for (int i = 0; i < s; i++) R = R * R;
    return R;
}
V6 commutator(const V6 &a, const V6 &b) { return vee(hat(a) * hat(b) - hat(b) * hat(a)); }

M6 smallAd(const V6 &x) { // [[phi^ 0],[rho^ phi^]]
    M6 a = M6::Zero();
    a.block<3, 3>(0, 0) = SO3d::computeHat(x.head<3>());
    a.block<3, 3>(3, 3) = SO3d::computeHat(x.head<3>());
    a.block<3, 3>(3, 0) = SO3d::computeHat(x.tail<3>());
    return a;
}

// Tangent-exp matrix  L * Jl(L xi)  = int_0^L exp(s ad_xi) ds, closed form used by the plugin mappings
// (CosseratBeamGeometry::computeTangExpImplementation), valid for the FULL generator xi (with the e1 reference).
M6 tangExp(double L, const V6 &xi) {
    const M6 ad = smallAd(xi);
    const double th = xi.head<3>().norm();
    const M6 I = M6::Identity();
    if (th < 1e-9) return L * I + 0.5 * L * L * ad;
    const double Lt = L * th, c = std::cos(Lt), s = std::sin(Lt), t2 = th * th, t3 = t2 * th, t4 = t2 * t2, t5 = t4 * th;
    const double a1 = (4 - 4 * c - Lt * s) / (2 * t2), a2 = (4 * Lt + Lt * c - 5 * s) / (2 * t3);
    const double a3 = (2 - 2 * c - Lt * s) / (2 * t4), a4 = (2 * Lt + Lt * c - 3 * s) / (2 * t5);
    const M6 ad2 = ad * ad;
    return L * I + a1 * ad + a2 * ad2 + a3 * ad2 * ad + a4 * ad2 * ad2;
}

} // namespace

TEST(ComponentOrderTest, ExpAndLogFollowPhiRhoLayout) {
    for (int i = 0; i < 100; i++) {
        const V6 x = rnd6(0.8);
        const SE3d g = SE3d::computeExp(x);
        EXPECT_LT((g.matrix() - expm(hat(x))).norm(), 1e-12);
        EXPECT_LT((g.computeLog() - x).norm(), 1e-10);
    }
}

TEST(ComponentOrderTest, TwistHatBracketAndSmallAdjoint) {
    const V6 x = rnd6(1), y = rnd6(1);
    const Twistd tx(x), ty(y);
    EXPECT_LT((tx.hat() - hat(x)).norm(), 1e-14);
    EXPECT_LT((tx.bracket(ty).data() - commutator(x, y)).norm(), 1e-13);
    EXPECT_LT((tx.smallAdjoint() * y - commutator(x, y)).norm(), 1e-13);
}

TEST(ComponentOrderTest, AdjointConjugationAndDerivative) {
    for (int i = 0; i < 20; i++) {
        const SE3d g = SE3d::computeExp(rnd6(0.7));
        const V6 w = rnd6(1e-6);
        const V6 conj = (g * SE3d::computeExp(w) * g.inverse()).computeLog();
        EXPECT_LT((conj - g.computeAdjoint() * w).norm(), 1e-10);

        const V6 xi = rnd6(1);
        const double h = 1e-6;
        const M6 fd = (SE3d::computeExp(h * xi).computeAdjoint() - SE3d::computeExp(-h * xi).computeAdjoint()) / (2 * h);
        EXPECT_LT((fd - smallAd(xi)).norm(), 1e-8);
    }
}

TEST(ComponentOrderTest, CoAdjointAndWrenchAreDualToTwist) {
    const SE3d g = SE3d::computeExp(rnd6(0.8));
    EXPECT_LT((g.computeCoAdjoint() - g.computeAdjoint().inverse().transpose()).norm(), 1e-12);
    const V6 xa = rnd6(1), wb = rnd6(1); // Wrench = [torque ; force] pairs with Twist = [phi ; rho]
    const Twistd xi_a(xa), xi_b = xi_a.transformBy(g.computeAdjoint());
    const Wrenchd w_b(wb), w_a = w_b.transformBy(g.computeAdjoint());
    EXPECT_NEAR(w_b.dot(xi_b), w_a.dot(xi_a), 1e-12);
}

TEST(ComponentOrderTest, LeftJacobianInverseMatchesFiniteDifference) {
    const SE3d g = SE3d::computeExp(rnd6(0.8));
    const V6 l = g.computeLog();
    const double h = 1e-6;
    M6 fd;
    for (int i = 0; i < 6; i++) {
        V6 d = V6::Zero();
        d[i] = h;
        fd.col(i) = ((SE3d::computeExp(d) * g).computeLog() - (SE3d::computeExp(-d) * g).computeLog()) / (2 * h);
    }
    EXPECT_LT((SE3d::computeLeftJacobianInverse(l) - fd).norm(), 1e-7);
}

TEST(ComponentOrderTest, ExpCosseratIsExpOfStrainPlusReferenceE1) {
    for (int i = 0; i < 50; i++) {
        const V6 Xi = rnd6(1.5);
        const double L = 0.05 + 0.3 * std::abs(rnd6(1)[0]);
        V6 tot = Xi;
        tot[3] += 1.0; // reference elongation along x lives in rho_x (index 3)
        EXPECT_LT((SE3d::expCosserat(Xi, L).matrix() - SE3d::computeExp(L * tot).matrix()).norm(), 1e-12);
    }
}

TEST(ComponentOrderTest, ExpCosseratSemanticOfEachStrainComponent) {
    // Xi = {torsion, bending_y, bending_z, elongation, shear_y, shear_z}
    const double L = 0.2;
    auto pose = [&](int idx, double v) {
        V6 Xi = V6::Zero();
        Xi[idx] = v;
        return SE3d::expCosserat(Xi, L);
    };
    EXPECT_LT((pose(0, 2.0).rotation().log() - V3(0.4, 0, 0)).norm(), 1e-12);             // torsion: about x
    EXPECT_LT((pose(1, 2.0).rotation().log() - V3(0, 0.4, 0)).norm(), 1e-12);             // bending_y: about y
    EXPECT_LT((pose(2, 2.0).rotation().log() - V3(0, 0, 0.4)).norm(), 1e-12);             // bending_z: about z
    EXPECT_LT((pose(3, 0.5).translation() - V3(0.3, 0, 0)).norm(), 1e-12);                // elongation: (1+e) L along x
    EXPECT_LT((pose(4, 0.5).translation() - V3(0.2, 0.1, 0)).norm(), 1e-12);              // shear_y
    EXPECT_LT((pose(5, 0.5).translation() - V3(0.2, 0, 0.1)).norm(), 1e-12);              // shear_z
    // bending_z curves towards +y, bending_y towards -z (right-handed rotation of the x tangent)
    EXPECT_GT(pose(2, 2.0).translation().y(), 0);
    EXPECT_LT(pose(1, 2.0).translation().z(), 0);
}

TEST(ComponentOrderTest, BodyJacobianForwardMatchesFiniteDifference) {
    const int N = 4;
    const double L = 0.25;
    std::vector<V6> xi(N);
    for (auto &x : xi) x = rnd6(0.8);
    auto tip = [&](const std::vector<V6> &X) {
        SE3d g;
        for (int k = 0; k < N; k++) g = g * SE3d::expCosserat(X[k], L);
        return g;
    };
    CosseratBodyJacobian<double> J(N);
    for (int k = 0; k < N; k++) {
        V6 full = xi[k];
        full[3] += 1.0;
        J.pushSection(SE3d::expCosserat(xi[k], L), tangExp(L, full));
    }
    for (int dir = 0; dir < 6; dir++) {
        for (int sec = 0; sec < N; sec++) {
            std::vector<Twistd> rates(N, Twistd::Zero());
            rates[sec].data()[dir] = 1.0;
            const V6 pred = J.applyForward(rates).back().data();
            const double h = 1e-6;
            auto Xp = xi, Xm = xi;
            Xp[sec][dir] += h;
            Xm[sec][dir] -= h;
            const V6 fd = (tip(Xm).inverse() * tip(Xp)).computeLog() / (2 * h);
            EXPECT_LT((pred - fd).norm(), 1e-7) << "dir=" << dir << " sec=" << sec;
        }
    }
}

// ------------------------------------------------------------------ ad / BCH / parallelTransport

namespace {
V6 vee_commutator(const V6 &a, const V6 &b) { return commutator(a, b); }
} // namespace

TEST(ComponentOrderTest, AdIsTheAlgebraOperatorForPhiRhoLayout) {
    for (int i = 0; i < 20; i++) {
        const V6 v = rnd6(1), w = rnd6(1);
        const SE3d any = SE3d::computeExp(rnd6(1.0)), id;
        // ad_v = [[phi^ 0],[rho^ phi^]] whatever the group element it is called on (no R block)
        EXPECT_LT((any.ad(v) - smallAd(v)).norm(), 1e-14);
        EXPECT_LT((id.ad(v) - any.ad(v)).norm(), 1e-14);
        EXPECT_LT((SE3d::ad(v) - smallAd(v)).norm(), 1e-14);
        // ad_v w = [v, w] = vee([hat v, hat w])
        EXPECT_LT((SE3d::ad(v) * w - vee_commutator(v, w)).norm(), 1e-13);
        EXPECT_LT((SE3d::ad(v, w) - vee_commutator(v, w)).norm(), 1e-13);
        EXPECT_LT((SE3d::ad(v, w) + SE3d::ad(w, v)).norm(), 1e-13); // antisymmetry
        // overloads with a point: the point never entered the result; a pure point is a rho generator
        const V3 p(0.3, -0.2, 0.5);
        EXPECT_LT((any.ad(v, p) - smallAd(v)).norm(), 1e-14);
        EXPECT_LT((any.ad(p, v) - smallAd(v)).norm(), 1e-14);
        V6 pt = V6::Zero();
        pt.tail<3>() = p;
        EXPECT_LT((any.ad(p) - smallAd(pt)).norm(), 1e-14);
        // Jacobi identity [u,[v,w]] + [v,[w,u]] + [w,[u,v]] = 0
        const V6 u = rnd6(1);
        EXPECT_LT((SE3d::ad(u, SE3d::ad(v, w)) + SE3d::ad(v, SE3d::ad(w, u)) + SE3d::ad(w, SE3d::ad(u, v))).norm(), 1e-12);
        // d/dt Ad(exp(t v)) at t = 0 equals ad_v (finite differences)
        const double h = 1e-6;
        const M6 fd = (SE3d::computeExp(h * v).computeAdjoint() - SE3d::computeExp(-h * v).computeAdjoint()) / (2 * h);
        EXPECT_LT((fd - SE3d::ad(v)).norm(), 1e-8);
    }
}

TEST(ComponentOrderTest, BchIsSecondOrderForPhiRhoLayout) {
    // Reference: log(exp(a) exp(b)) with exp from the 4x4 matrix series (independent of SE3::computeExp).
    auto reference = [](const V6 &a, const V6 &b) { return SE3d(expm(hat(a)) * expm(hat(b))).computeLog(); };
    const V6 a0 = rnd6(1), b0 = rnd6(1);
    double prevErr = 0;
    for (double eps : {0.2, 0.1, 0.05, 0.025}) {
        const V6 a = eps * a0 / a0.norm(), b = eps * b0 / b0.norm();
        const V6 ref = reference(a, b);
        const double err = (SE3d::computeBCH(a, b) - ref).norm();
        // Second-order BCH: error O(eps^3), much smaller than the O(eps^2) error of a + b alone.
        EXPECT_LT(err, 0.5 * eps * eps * eps) << "eps=" << eps;
        EXPECT_LT(err, 0.1 * (a + b - ref).norm()) << "eps=" << eps;
        if (prevErr > 0) EXPECT_GT(prevErr / err, 6.0) << "eps=" << eps; // ~8 when eps is halved (3rd order)
        prevErr = err;
    }
    // Exact second-order term, component by component: [phi_a x phi_b ; phi_a x rho_b - phi_b x rho_a] / 2
    const V6 a = rnd6(1e-3), b = rnd6(1e-3);
    const V6 second = SE3d::computeBCH(a, b) - a - b;
    EXPECT_LT((second.head<3>() - 0.5 * a.head<3>().cross(b.head<3>())).norm(), 1e-15);
    EXPECT_LT((second.tail<3>() - 0.5 * (a.head<3>().cross(b.tail<3>()) - b.head<3>().cross(a.tail<3>()))).norm(), 1e-15);
}

TEST(ComponentOrderTest, ParallelTransportFollowsItsDefinition) {
    for (int i = 0; i < 20; i++) {
        const V6 u = rnd6(0.8), v = rnd6(1);
        const SE3d g(expm(hat(u))); // group element built with the matrix series
        const V6 expected = v - 0.5 * vee_commutator(g.computeLog(), v);
        EXPECT_LT((g.parallelTransport(v) - expected).norm(), 1e-13);
        // first-order accuracy of the definition v - 1/2 ad_u v  ==  (I - 1/2 ad_u) v
        EXPECT_LT((g.parallelTransport(v) - (M6::Identity() - 0.5 * smallAd(g.computeLog())) * v).norm(), 1e-13);
    }
    // pure rotation about x moves the y-angular component into z-angular (angular part first)
    V6 u = V6::Zero();
    u[0] = M_PI / 2;
    V6 v = V6::Zero();
    v[1] = 1.0;
    const V6 t = SE3d::computeExp(u).parallelTransport(v);
    EXPECT_NEAR(t[2], -M_PI / 4, 1e-14);
    EXPECT_NEAR(t[5], 0.0, 1e-14);
}

// ------------------------------------------------------------------ action Jacobians

namespace {
// y = (expm(hat(delta)) * G * [p;1]).head<3>(): left perturbation, exp from the 4x4 matrix series.
V3 leftPerturbedAction(const M4 &G, const V6 &delta, const V3 &p) {
    const Eigen::Vector4d ph(p.x(), p.y(), p.z(), 1.0);
    return (expm(hat(delta)) * G * ph).head<3>();
}
} // namespace

TEST(ComponentOrderTest, ActionJacobianGroupMatchesLeftPerturbationFiniteDifferences) {
    for (int i = 0; i < 20; i++) {
        const V6 xi0 = rnd6(0.9); // generic rotation AND non-zero translation
        const M4 G = expm(hat(xi0));
        const SE3d g(G);
        const V3 p(rnd6(1).head<3>());
        Eigen::Matrix<double, 3, 6> fd;
        const double h = 1e-6;
        for (int k = 0; k < 6; k++) {
            V6 d = V6::Zero();
            d[k] = h;
            fd.col(k) = (leftPerturbedAction(G, d, p) - leftPerturbedAction(G, -d, p)) / (2 * h);
        }
        const auto J = g.actionJacobianGroup(p);
        EXPECT_LT((J - fd).norm(), 1e-8);
        // layout [phi ; rho]: the rho columns (3..5) are the identity, the phi columns are -[g.p]x
        EXPECT_LT((J.block<3, 3>(0, 3) - Eigen::Matrix3d::Identity()).norm(), 1e-15);
        EXPECT_LT((J.block<3, 3>(0, 0) + SO3d::computeHat(g.act(p))).norm(), 1e-14);
        // the pair version returns the same group block, plus the point block R
        const auto [Jg, Jp] = g.actionJacobians(p);
        EXPECT_LT((Jg - J).norm(), 1e-15);
        EXPECT_LT((Jp - g.rotation().matrix()).norm(), 1e-15);
    }
}

TEST(ComponentOrderTest, ActionJacobianPointMatchesFiniteDifferences) {
    const SE3d g(expm(hat(rnd6(0.9))));
    const V3 p(0.4, -0.7, 1.1);
    Eigen::Matrix3d fd;
    const double h = 1e-6;
    for (int k = 0; k < 3; k++) {
        V3 d = V3::Zero();
        d[k] = h;
        fd.col(k) = (g.act(p + d) - g.act(p - d)) / (2 * h);
    }
    EXPECT_LT((g.actionJacobianPoint(p) - fd).norm(), 1e-8);
}

TEST(ComponentOrderTest, ActionJacobianIsNotTheRightPerturbationOne) {
    // Sanity: the left-perturbation Jacobian differs from the right one [-R [p]x | R] when t != 0 and R != I.
    const SE3d g(expm(hat(rnd6(0.9))));
    const V3 p(0.4, -0.7, 1.1);
    Eigen::Matrix<double, 3, 6> right;
    right << -g.rotation().matrix() * SO3d::computeHat(p), g.rotation().matrix();
    EXPECT_GT((g.actionJacobianGroup(p) - right).norm(), 1e-2);
}
