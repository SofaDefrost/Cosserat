// Independent long double references for the Lie-group audit tests.
// Everything here is computed from the matrix definition of the groups (exp by
// scaling-and-squaring of the Taylor series, Jacobians by the series
// sum ad^k/(k+1)!), never from the closed forms used by the library, so a test
// that compares against it cannot pass on a wrong closed form.
#pragma once

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <random>

namespace liegroup_ref {

using LD = long double;
template<int N> using MatL = Eigen::Matrix<LD, N, N>;
template<int N> using VecL = Eigen::Matrix<LD, N, 1>;

/// exp(A) by Taylor series with scaling and squaring.
template<int N>
MatL<N> expm(const MatL<N> &A) {
	LD nrm = A.cwiseAbs().rowwise().sum().maxCoeff();
	int s = 0;
	while (nrm > 0.25L) {
		nrm /= 2;
		++s;
	}
	const MatL<N> B = A / std::ldexp(static_cast<LD>(1), s);
	MatL<N> R = MatL<N>::Identity(), T = MatL<N>::Identity();
	for (int k = 1; k < 30; ++k) {
		T = T * B / static_cast<LD>(k);
		R += T;
	}
	for (int i = 0; i < s; ++i)
		R = R * R;
	return R;
}

/// Left Jacobian of exp: sum_k ad^k / (k+1)!  (the right one is jacobianSeries(-ad)).
template<int N>
MatL<N> jacobianSeries(const MatL<N> &ad) {
	MatL<N> R = MatL<N>::Identity(), T = MatL<N>::Identity();
	for (int k = 1; k < 60; ++k) {
		T = T * ad / static_cast<LD>(k + 1);
		R += T;
	}
	return R;
}

inline MatL<3> hat3(const VecL<3> &v) {
	MatL<3> M;
	M << 0, -v[2], v[1], v[2], 0, -v[0], -v[1], v[0], 0;
	return M;
}

/// se(3) with the library layout xi = [phi; rho] (angular first).
inline MatL<4> hatSE3(const VecL<6> &x) {
	MatL<4> M = MatL<4>::Zero();
	M.block<3, 3>(0, 0) = hat3(x.head<3>());
	M.block<3, 1>(0, 3) = x.tail<3>();
	return M;
}
inline MatL<6> adSE3(const VecL<6> &x) {
	MatL<6> M = MatL<6>::Zero();
	M.block<3, 3>(0, 0) = hat3(x.head<3>());
	M.block<3, 3>(3, 3) = hat3(x.head<3>());
	M.block<3, 3>(3, 0) = hat3(x.tail<3>());
	return M;
}

/// se(2) with the library layout xi = [rho_x; rho_y; theta].
inline MatL<3> hatSE2(const VecL<3> &x) {
	MatL<3> M = MatL<3>::Zero();
	M(0, 1) = -x[2];
	M(1, 0) = x[2];
	M(0, 2) = x[0];
	M(1, 2) = x[1];
	return M;
}

template<typename A, typename B>
double maxAbsDiff(const A &a, const B &b) {
	return static_cast<double>((a.template cast<LD>() - b.template cast<LD>()).cwiseAbs().maxCoeff());
}

/// Random unit axis times angle plus a random vector, as a 6-vector [phi; rho].
template<typename Rng>
Eigen::Matrix<double, 6, 1> randomSE3Tangent(Rng &rng, double angle, double rhoScale = 2.0) {
	std::uniform_real_distribution<double> u(-1.0, 1.0);
	Eigen::Vector3d a(u(rng), u(rng), u(rng));
	a.normalize();
	Eigen::Matrix<double, 6, 1> x;
	x.head<3>() = a * angle;
	x.tail<3>() = Eigen::Vector3d(u(rng), u(rng), u(rng)) * rhoScale;
	return x;
}

} // namespace liegroup_ref
