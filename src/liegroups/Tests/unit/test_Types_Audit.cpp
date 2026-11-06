// Types<Scalar> carries its own SO3Log / SE3Exp / SE3Log (matrix based, layout [rho ; phi], i.e.
// linear first, unlike SE3 which is angular first). Checked against the long double references.

#include <gtest/gtest.h>
#include <liegroups/Types.h>

#include <random>

#include "LieGroupReference.h"

using namespace sofa::component::cosserat::liegroups;
using namespace liegroup_ref;

namespace {
using T = Typesd;
const double kPi = 3.14159265358979323846;
const double kAngles[] = {1e-9, 1e-6, 1e-5, 5e-5, 1e-4, 2e-4, 1e-3, 1e-2, 0.014, 0.1, 1.0, 3.0, 3.13, 3.14, kPi - 1e-6, kPi - 1e-9};
// The rotation matrix is rounded to double: near pi the axis is only known to ~1e-16 / sin(angle).
double tol(double angle) { return 5e-15 + (angle > 2.0 ? 3e-16 / std::max(1e-3, std::sin(angle)) : 0.0); }
} // namespace

TEST(TypesAudit, SO3LogInvertsTheRotationMatrix) {
	std::mt19937_64 rng(81);
	for (double angle : kAngles) {
		for (int i = 0; i < 5; ++i) {
			const Eigen::Matrix<double, 6, 1> x = randomSE3Tangent(rng, angle);
			const Eigen::Vector3d phi = x.head<3>();
			const Eigen::Matrix3d R = expm<4>(hatSE3(x.cast<LD>())).block<3, 3>(0, 0).cast<double>();
			EXPECT_LT((T::SO3Log(R) - phi).cwiseAbs().maxCoeff(), tol(angle))
					<< "angle " << angle;
			EXPECT_LT((T::logSO3(R) - T::skew3(phi)).cwiseAbs().maxCoeff(), tol(angle))
					<< "logSO3 angle " << angle;
		}
	}
}

TEST(TypesAudit, SE3ExpMatchesMatrixExponential) {
	std::mt19937_64 rng(82);
	for (double angle : kAngles) {
		for (int i = 0; i < 5; ++i) {
			const Eigen::Matrix<double, 6, 1> x = randomSE3Tangent(rng, angle); // [phi ; rho]
			Eigen::Matrix<double, 6, 1> y;
			y << x.tail<3>(), x.head<3>(); // Types layout [rho ; phi]
			EXPECT_LT(maxAbsDiff(T::SE3Exp(y), expm<4>(hatSE3(x.cast<LD>()))), 4e-15) << "angle " << angle;
		}
	}
}

TEST(TypesAudit, SE3LogInvertsExp) {
	std::mt19937_64 rng(83);
	for (double angle : kAngles) {
		for (int i = 0; i < 5; ++i) {
			const Eigen::Matrix<double, 6, 1> x = randomSE3Tangent(rng, angle);
			Eigen::Matrix<double, 6, 1> y;
			y << x.tail<3>(), x.head<3>();
			const Eigen::Matrix4d Tm = expm<4>(hatSE3(x.cast<LD>())).cast<double>();
			EXPECT_LT((T::SE3Log(Tm) - y).cwiseAbs().maxCoeff(), 10 * tol(angle))
					<< "angle " << angle;
		}
	}
}
