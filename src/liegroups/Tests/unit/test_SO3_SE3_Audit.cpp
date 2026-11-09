// Second audit pass over SO3 and SE3 (and the base-class defaults they inherit).
// Each test pins one defect found by reading SO3.h / SE3.h and measuring it; the commit that
// fixes it carries the numbers.

#include <gtest/gtest.h>
#include <liegroups/SE2.h>
#include <liegroups/SE23.h>
#include <liegroups/SE3.h>
#include <liegroups/SO3.h>

#include <random>

#include "LieGroupReference.h"

using namespace sofa::component::cosserat::liegroups;
using namespace liegroup_ref;

namespace {
constexpr double kPi = 3.14159265358979323846;

template<typename G>
int roundTripFailures(std::mt19937_64 &rng, int n) {
	int fails = 0;
	for (int i = 0; i < n; ++i) {
		const G a = G::computeRandom(rng), b = G::computeRandom(rng);
		// (a b) b^-1 == a up to rounding: the default tolerance must accept it.
		if (!((a * b) * b.inverse()).isApprox(a))
			++fails;
		if (!(a * a.inverse()).isApprox(G::Identity()))
			++fails;
	}
	return fails;
}
} // namespace

// The default tolerance of isApprox was the machine epsilon (2.2e-16): a rounding error of one ulp
// on a quaternion component or on a translation of norm ~1 made g * h * h^-1 != g
// (SO3 3.3 % of 20000 random pairs, SE3 58 %).
TEST(SO3SE3Audit, IsApproxDefaultToleranceAcceptsRoundingErrors) {
	std::mt19937_64 rng(1);
	EXPECT_EQ(roundTripFailures<SO3d>(rng, 5000), 0);
	EXPECT_EQ(roundTripFailures<SE3d>(rng, 5000), 0);
	EXPECT_EQ(roundTripFailures<SE2d>(rng, 5000), 0);
	EXPECT_EQ(roundTripFailures<SE23<double>>(rng, 5000), 0);
}

// ... without becoming vacuous: a 1e-9 rotation or translation offset is still detected.
TEST(SO3SE3Audit, IsApproxDefaultToleranceStillDiscriminates) {
	const SE3d g = SE3d::exp((SE3d::TangentVector() << 0.3, -0.2, 0.5, 1.0, 2.0, -1.0).finished());
	const SE3d dr = g * SE3d::exp((SE3d::TangentVector() << 1e-9, 0, 0, 0, 0, 0).finished());
	const SE3d dt = g * SE3d::exp((SE3d::TangentVector() << 0, 0, 0, 1e-9, 0, 0).finished());
	EXPECT_FALSE(g.isApprox(dr));
	EXPECT_FALSE(g.isApprox(dt));
	EXPECT_TRUE(g.isApprox(g));
}

// computeRandom drew a uniform angle in [-pi, pi] about a random axis: not the Haar measure
// (rotations cluster around the identity). For Haar, the angle has density (1 - cos t)/pi on
// [0, pi]: mean pi/2 + 2/pi = 2.2074 (the old generator gave 1.5714) and E[trace R] = 0
// (the old one gave 1).
TEST(SO3SE3Audit, RandomRotationsAreHaarDistributed) {
	std::mt19937_64 rng(2);
	const int n = 100000;
	double meanAngle = 0, meanTrace = 0, meanAngleSE3 = 0, meanTraceSE3 = 0;
	for (int i = 0; i < n; ++i) {
		const SO3d r = SO3d::computeRandom(rng);
		meanAngle += r.distance(SO3d::Identity()); // geodesic angle in [0, pi]
		meanTrace += r.matrix().trace();
		const SO3d rs = SE3d::computeRandom(rng).rotation();
		meanAngleSE3 += rs.distance(SO3d::Identity());
		meanTraceSE3 += rs.matrix().trace();
	}
	const double expected = kPi / 2 + 2 / kPi;
	EXPECT_NEAR(meanAngle / n, expected, 0.01);
	EXPECT_NEAR(meanTrace / n, 0.0, 0.02);
	EXPECT_NEAR(meanAngleSE3 / n, expected, 0.01);
	EXPECT_NEAR(meanTraceSE3 / n, 0.0, 0.02);
}

// randomUnitVector normalised a vector drawn in the cube [-1, 1]^3: directions along the cube
// diagonals are favoured. On the uniform sphere E|x| = 1/2 (cube-normalised: 0.516).
TEST(SO3SE3Audit, RandomUnitVectorIsUniformOnTheSphere) {
	std::mt19937_64 rng(3);
	const int n = 400000;
	double mx = 0, mxy = 0;
	for (int i = 0; i < n; ++i) {
		const Eigen::Vector3d u = Typesd::randomUnitVector<3>(rng);
		EXPECT_NEAR(u.norm(), 1.0, 1e-14);
		mx += std::abs(u.x());
		mxy += u.x() * u.y();
	}
	EXPECT_NEAR(mx / n, 0.5, 0.002);
	EXPECT_NEAR(mxy / n, 0.0, 0.002);
}

// SE3(Quaternion, t) went through quat.toRotationMatrix() and back: for a quaternion that is not
// unit (SO3(Quaternion) documents "will be normalized") the matrix is not a rotation and the
// stored quaternion had norm sqrt(2) for q = (1, 1, 0, 0), a 90 degree rotation about x was lost.
TEST(SO3SE3Audit, SE3FromNonUnitQuaternionIsNormalised) {
	const Eigen::Quaterniond q(1, 1, 0, 0);
	const Eigen::Vector3d t(0.1, -0.2, 0.3);
	const SE3d g(q, t);
	const SE3d expected(SO3d(q), t);
	EXPECT_NEAR(g.rotation().quaternion().norm(), 1.0, 1e-15);
	EXPECT_LT((g.matrix() - expected.matrix()).cwiseAbs().maxCoeff(), 1e-15);
	// and a unit quaternion is reproduced exactly (no matrix round trip any more)
	const Eigen::Quaterniond u = Eigen::Quaterniond(0.3, -0.5, 0.2, 0.7).normalized();
	EXPECT_LT((SE3d(u, t).rotation().quaternion().coeffs() - u.coeffs()).cwiseAbs().maxCoeff(), 3e-16);
}

// SO3(Matrix) kept the quaternion Eigen builds from a matrix as is: for a rotation matrix that is
// only orthonormal to 1e-6 (accumulated product, float input) the quaternion has norm 1 + O(1e-6),
// so isValid(), log() and matrix() of the result drifted; 1.01 R gave a quaternion of norm 1.0039.
TEST(SO3SE3Audit, SO3FromImperfectMatrixIsNormalised) {
	const Eigen::Matrix3d R = SO3d(0.7, Eigen::Vector3d(1, 2, 3)).matrix();
	for (double scale: {1.0 + 1e-6, 1.01, 0.99}) {
		const SO3d g(Eigen::Matrix3d(scale * R));
		EXPECT_NEAR(g.quaternion().norm(), 1.0, 1e-15) << "scale=" << scale;
		EXPECT_TRUE(g.isValid());
		EXPECT_LT((g.matrix() * g.matrix().transpose() - Eigen::Matrix3d::Identity()).cwiseAbs().maxCoeff(), 1e-14);
	}
	// a uniformly scaled rotation is still the same rotation
	EXPECT_LT((SO3d(Eigen::Matrix3d((1.0 + 1e-6) * R)).matrix() - R).cwiseAbs().maxCoeff(), 1e-5);
}
