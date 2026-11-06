// Audit of SE2 (and the SO2 / RealSpace Jacobian hooks) against long double references.
// SE2 tangent layout: xi = [rho_x ; rho_y ; theta].

#include <gtest/gtest.h>
#include <liegroups/RealSpace.h>
#include <liegroups/SE2.h>
#include <liegroups/SO2.h>

#include <random>

#include "LieGroupReference.h"

using namespace sofa::component::cosserat::liegroups;
using namespace liegroup_ref;

namespace {
using V3 = Eigen::Vector3d;
using M3 = Eigen::Matrix3d;
const double kPi = 3.14159265358979323846;
// theta may be negative; |theta| < pi because SO2 stores an angle in [-pi, pi).
const double kThetas[] = {1e-12, 1e-9, 1e-7, 1e-5, 1e-4, 2e-4, 1e-3, 1e-2, 0.1, 0.3, 1.0, 2.0, 3.0,
						  3.14, kPi - 1e-6, kPi - 1e-9, -1e-9, -1e-3, -1.0, -3.0, -kPi + 1e-6};

V3 randomXi(std::mt19937_64 &rng, double theta) {
	std::uniform_real_distribution<double> u(-1.0, 1.0);
	return V3(2 * u(rng), 2 * u(rng), theta);
}
/// SE2 element built from the long double exponential of xi.
SE2d referenceElement(const V3 &xi) {
	const MatL<3> T = expm<3>(hatSE2(xi.cast<LD>()));
	return SE2d(xi[2], Eigen::Vector2d(static_cast<double>(T(0, 2)), static_cast<double>(T(1, 2))));
}
MatL<3> adSE2(const V3 &x) {
	// [x, y] = hat^-1(hat x hat y - hat y hat x): rho part theta_x J rho_y - theta_y J rho_x
	MatL<3> A = MatL<3>::Zero();
	const LD th = x[2];
	A(0, 1) = -th;
	A(1, 0) = th;
	A(0, 2) = x[1];
	A(1, 2) = -static_cast<LD>(x[0]);
	return A;
}
} // namespace

TEST(SE2Audit, MatrixAndTangentConstructorsCompileAndAgree) {
	std::mt19937_64 rng(44);
	const V3 xi = randomXi(rng, 0.8);
	const SE2d g = SE2d::exp(xi);
	const SE2d fromMatrix(g.matrix());
	const SE2d fromTangent(xi);
	EXPECT_TRUE(fromMatrix.isApprox(g, 1e-12));
	EXPECT_TRUE(fromTangent.isApprox(g, 1e-12));
}

