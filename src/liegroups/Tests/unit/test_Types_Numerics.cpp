// Numerical audit of the scalar toolbox of Types.h used by SO3 / SE3 (sinc family, 1 - cos).
// References are long double Taylor series, never the closed forms of the library.

#include <gtest/gtest.h>
#include <liegroups/Types.h>

#include <cmath>

using namespace sofa::component::cosserat::liegroups;

namespace {
using LD = long double;

/// sum_k (-1)^k x^(2k) / (2k + m)!  : m = 2 gives (1-cos x)/x^2, m = 3 gives (x - sin x)/x^3.
LD seriesRef(LD x, int m) {
	LD term = 1, fact = 1;
	for (int i = 2; i <= m; ++i)
		fact *= i;
	LD sum = 0, sign = 1, p = 1;
	for (int k = 0; k < 30; ++k) {
		LD f = 1;
		for (int i = 2; i <= 2 * k + m; ++i)
			f *= i;
		sum += sign * p / f;
		p *= x * x;
		sign = -sign;
	}
	(void) term;
	(void) fact;
	return sum;
}

const double kXs[] = {1e-9, 1e-7, 1e-5, 5e-5, 1e-4, 1.01e-4, 2e-4, 1e-3, 1e-2, 0.1, 0.24, 0.26, 0.5, 1.0, 2.0, 3.0};

double relErr(double value, LD ref) { return static_cast<double>(std::abs((LD) value - ref) / std::abs(ref)); }
} // namespace

// (1 - cos x)/x^2 : the closed form cancels below 0.1 rad, the old switch at 1e-4 left 2.6e-9 just above it.
TEST(TypesNumerics, Sinc2RelativeError) {
	for (double x: kXs)
		EXPECT_LT(relErr(Typesd::sinc2(x), seriesRef(x, 2)), 5e-15) << "x=" << x;
}

// (x - sin x)/x^3 : same cancellation, 7.6e-9 relative error at 1.01e-4.
TEST(TypesNumerics, Sinc3RelativeError) {
	for (double x: kXs)
		EXPECT_LT(relErr(Typesd::sinc3(x), seriesRef(x, 3)), 5e-14) << "x=" << x;
}

// 1 - cos x was documented as numerically stable but computed as 1 - std::cos(x): 8e-8 relative error at 1e-5.
TEST(TypesNumerics, OneMinusCosRelativeError) {
	for (double x: kXs) {
		const LD ref = 2 * std::sin((LD) x / 2) * std::sin((LD) x / 2);
		EXPECT_LT(relErr(Typesd::oneMinusCos(x), ref), 5e-15) << "x=" << x;
	}
}

TEST(TypesNumerics, FloatPrecision) {
	for (float x: {1e-3f, 1e-2f, 0.1f, 0.3f, 0.9f}) {
		EXPECT_LT(relErr(Typesf::sinc2(x), seriesRef(x, 2)), 5e-7) << "x=" << x;
		EXPECT_LT(relErr(Typesf::sinc3(x), seriesRef(x, 3)), 5e-6) << "x=" << x;
		const LD ref = 2 * std::sin((LD) x / 2) * std::sin((LD) x / 2);
		EXPECT_LT(relErr(Typesf::oneMinusCos(x), ref), 5e-7) << "x=" << x;
	}
}
