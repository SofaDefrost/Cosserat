// Audit tests of the generic helpers of LieGroupBase (LieGroupBase.inl) that are
// reached through SO3 / SE3: dexp_alternative, dexpInv_alternative, dlog, BCH and slerp.
// The references are the series of LieGroupReference.h (long double), not closed forms.

#include <gtest/gtest.h>
#include <liegroups/SE3.h>
#include <liegroups/SO3.h>

#include <random>

#include "LieGroupReference.h"

using namespace sofa::component::cosserat::liegroups;
using namespace liegroup_ref;

namespace {
using V6 = Eigen::Matrix<double, 6, 1>;
using M6 = Eigen::Matrix<double, 6, 6>;
const double kAngles[] = {1e-9, 1e-6, 1e-3, 0.1, 1.0, 2.0, 3.0};
} // namespace

TEST(LieGroupBaseGeneric, DexpAlternativeIsTheRightJacobian) {
	// The header says dexp_alternative is equivalent to dexp() == rightJacobian().
	std::mt19937_64 rng(11);
	for (double angle : kAngles) {
		const V6 x = randomSE3Tangent(rng, angle);
		const MatL<6> Jr = jacobianSeries<6>(-adSE3(x.cast<LD>()));
		EXPECT_LT(maxAbsDiff(SE3d::dexp_alternative(x), Jr), 1e-12) << "SE3 angle " << angle;
		const Eigen::Vector3d w = x.head<3>();
		const MatL<3> Jr3 = jacobianSeries<3>(-hat3(w.cast<LD>()));
		EXPECT_LT(maxAbsDiff(SO3d::dexp_alternative(w), Jr3), 1e-12) << "SO3 angle " << angle;
	}
}

TEST(LieGroupBaseGeneric, DexpInvAlternativeIsTheRightJacobianInverse) {
	std::mt19937_64 rng(12);
	for (double angle : kAngles) {
		const V6 x = randomSE3Tangent(rng, angle);
		const MatL<6> Jri = jacobianSeries<6>(-adSE3(x.cast<LD>())).inverse();
		EXPECT_LT(maxAbsDiff(SE3d::dexpInv_alternative(x), Jri), 1e-11) << "SE3 angle " << angle;
		const Eigen::Vector3d w = x.head<3>();
		const MatL<3> Jri3 = jacobianSeries<3>(-hat3(w.cast<LD>())).inverse();
		EXPECT_LT(maxAbsDiff(SO3d::dexpInv_alternative(w), Jri3), 1e-11) << "SO3 angle " << angle;
	}
}

TEST(LieGroupBaseGeneric, SE3DlogIsTheInverseOfDexpAtLog) {
	// dlog() is the differential of log at the current element, in the same
	// (right) convention as dexp()/dexpInv(): dlog() * dexp(log g) == I.
	std::mt19937_64 rng(13);
	for (double angle : kAngles) {
		const V6 x = randomSE3Tangent(rng, angle);
		const SE3d g = SE3d::exp(x);
		const V6 l = g.log();
		EXPECT_LT((g.dlog() * SE3d::dexp(l) - M6::Identity()).cwiseAbs().maxCoeff(), 1e-10) << angle;
		EXPECT_LT((g.dlog() - SE3d::dexpInv(l)).cwiseAbs().maxCoeff(), 1e-12) << angle;
	}
}

TEST(LieGroupBaseGeneric, SE3DlogMatchesFiniteDifferencesOfLog) {
	// log(g exp(d)) = log(g) + dlog(g) d + O(d^2)
	std::mt19937_64 rng(14);
	const V6 x = randomSE3Tangent(rng, 0.8);
	const SE3d g = SE3d::exp(x);
	const double h = 1e-6;
	M6 J;
	for (int i = 0; i < 6; ++i) {
		V6 d = V6::Zero();
		d[i] = h;
		V6 dm = -d;
		J.col(i) = ((g * SE3d::exp(d)).log() - (g * SE3d::exp(dm)).log()) / (2 * h);
	}
	EXPECT_LT((J - g.dlog()).cwiseAbs().maxCoeff(), 1e-7);
}

TEST(LieGroupBaseGeneric, BchOrdersConvergeOnSE3) {
	// log(exp(v) exp(w)) - BCH_n = O(s^(n+1)): halving the size must divide the error by ~2^(n+1).
	std::mt19937_64 rng(15);
	std::uniform_real_distribution<double> u(-1.0, 1.0);
	V6 v0, w0;
	for (int i = 0; i < 6; ++i) {
		v0[i] = u(rng);
		w0[i] = u(rng);
	}
	auto err = [&](double s, int order) {
		const V6 v = v0 * s, w = w0 * s;
		const V6 ref = (SE3d::exp(v) * SE3d::exp(w)).log();
		return (SE3d::BCH(v, w, order) - ref).norm();
	};
	const double s = 0.04;
	for (int order = 2; order <= 5; ++order) {
		const double e1 = err(s, order), e2 = err(s / 2, order);
		const double expected = std::pow(2.0, order + 1);
		EXPECT_GT(e1 / e2, 0.7 * expected) << "order " << order << " e1=" << e1 << " e2=" << e2;
		EXPECT_LT(e1 / e2, 1.5 * expected) << "order " << order << " e1=" << e1 << " e2=" << e2;
	}
	EXPECT_LT(err(0.05, 5), 1e-9);
}

TEST(LieGroupBaseGeneric, SlerpIsTheGeodesicAtParameterT) {
	const SO3d a = SO3d::exp(Eigen::Vector3d(0.3, 0.5, -0.2));
	const SO3d b = SO3d::exp(Eigen::Vector3d(-1.0, 0.4, 1.2));
	for (double t : {0.1, 0.25, 0.5, 0.9})
		EXPECT_LT((a.slerp(b, t).log() - a.interpolate(b, t).log()).norm(), 1e-12) << t;
	std::mt19937_64 rng(16);
	const SE3d g = SE3d::exp(randomSE3Tangent(rng, 0.7)), h = SE3d::exp(randomSE3Tangent(rng, 1.1));
	for (double t : {0.1, 0.25, 0.5, 0.9}) {
		// the geodesic g exp(t log(g^-1 h)) is at relative distance t from g
		const V6 d = (g.inverse() * g.slerp(h, t)).log();
		const V6 full = (g.inverse() * h).log();
		EXPECT_LT((d - t * full).norm(), 1e-10) << t;
	}
}
