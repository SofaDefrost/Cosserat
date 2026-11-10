/******************************************************************************
 *                 SOFA, Simulation Open-Framework Architecture * (c) 2006
 *INRIA, USTL, UJF, CNRS, MGH                     *
 *                                                                             *
 * This program is free software; you can redistribute it and/or modify it     *
 * under the terms of the GNU Lesser General Public License as published by    *
 * the Free Software Foundation; either version 2.1 of the License, or (at     *
 * your option) any later version.                                             *
 *                                                                             *
 * This program is distributed in the hope that it will be useful, but WITHOUT *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
 * for more details.                                                           *
 *                                                                             *
 * You should have received a copy of the GNU Lesser General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>.        *
 ******************************************************************************/

#include <Eigen/Core>
#include <liegroups/Rn.h>
#include <sofa/testing/BaseTest.h>

namespace sofa::component::cosserat::liegroups::testing {

	using namespace sofa::testing;

	/**
	 * @brief Test suite for Rn Lie group.
	 * This test fixture provides common types and a small epsilon for floating-point comparisons.
	 * It is templated to allow testing Rn with different dimensions.
	 * @tparam T The type of Rn to test (e.g., Rn<double, 3>).
	 */
	template<typename T>
	class RnTest : public BaseTest {
	protected:
		using Rn = T;
		using Scalar = typename Rn::Scalar;
		using Vector = typename Rn::Vector;
		using TangentVector = typename Rn::TangentVector;

		const Scalar eps = 1e-9;
	};

	using RnTypes = ::testing::Types<Rn<double, 1>, Rn<double, 2>, Rn<double, 3>>;
	TYPED_TEST_SUITE(RnTest, RnTypes);

	/**
	 * @brief Tests the constructors of the Rn class.
	 * Verifies that Rn objects can be constructed from default (zero vector) and from an Eigen vector.
	 */
	TYPED_TEST(RnTest, Constructors) {
		using Rn = typename TestFixture::Rn;
		using Scalar = typename TestFixture::Scalar;
		using Vector = typename TestFixture::Vector;

		Rn g1;
		EXPECT_TRUE(g1.computeIdentity().computeIsApprox(g1, this->eps));

		Vector v = Vector::Random();
		Rn g2(v);
		EXPECT_TRUE(g2.computeLog().isApprox(v, this->eps));
	}

	/**
	 * @brief Tests the identity element of the Rn group.
	 * Verifies that the `computeIdentity()` method returns a zero vector.
	 */
	TYPED_TEST(RnTest, Identity) {
		using Rn = typename TestFixture::Rn;
		using Scalar = typename TestFixture::Scalar;
		using Vector = typename TestFixture::Vector;

		Rn identity = Rn::computeIdentity();
		EXPECT_TRUE(identity.computeLog().isApprox(Vector::Zero(), this->eps));
	}

	/**
	 * @brief Tests the inverse operation of the Rn group.
	 * Verifies that the inverse of a vector is its negation, and that composing a vector with its inverse results in
	 * the identity element (zero vector).
	 */
	TYPED_TEST(RnTest, Inverse) {
		using Rn = typename TestFixture::Rn;
		using Scalar = typename TestFixture::Scalar;
		using Vector = typename TestFixture::Vector;

		Vector v = Vector::Random();
		Rn g(v);
		Rn inv_g = g.computeInverse();
		EXPECT_TRUE(inv_g.computeLog().isApprox(-v, this->eps));

		Rn composed = g.compose(inv_g);
		EXPECT_TRUE(composed.computeIdentity().computeIsApprox(composed, this->eps));
	}

	/**
	 * @brief Tests the exponential and logarithmic maps of the Rn group.
	 * Verifies that applying `computeExp` to a Lie algebra element and then `computeLog` to the resulting group element
	 * returns the original Lie algebra element. For Rn, both `exp` and `log` are identity functions.
	 */
	TYPED_TEST(RnTest, ExpLog) {
		using Rn = typename TestFixture::Rn;
		using Scalar = typename TestFixture::Scalar;
		using TangentVector = typename TestFixture::TangentVector;

		TangentVector v = TangentVector::Random();
		Rn g = Rn::computeExp(v);
		EXPECT_TRUE(g.computeLog().isApprox(v, this->eps));
	}

	/**
	 * @brief Tests the group action of Rn on a point.
	 * Verifies that the action is equivalent to vector addition.
	 */
	TYPED_TEST(RnTest, Action) {
		using Rn = typename TestFixture::Rn;
		using Scalar = typename TestFixture::Scalar;
		using Vector = typename TestFixture::Vector;

		Vector g_data = Vector::Random();
		Rn g(g_data);
		Vector point = Vector::Random();
		Vector transformed_point = g.computeAction(point);
		EXPECT_TRUE(transformed_point.isApprox(point + g_data, this->eps));
	}

	/**
	 * @brief Tests the approximate equality comparison for Rn elements.
	 * Verifies that two Rn elements are considered approximately equal within a given tolerance.
	 */
	TYPED_TEST(RnTest, IsApprox) {
		using Rn = typename TestFixture::Rn;
		using Scalar = typename TestFixture::Scalar;
		using Vector = typename TestFixture::Vector;

		Vector v = Vector::Random();
		Rn g1(v);
		Rn g2(v + Vector::Constant(this->eps / 2.0));
		Rn g3(v + Vector::Constant(this->eps * 2.0));

		EXPECT_TRUE(g1.computeIsApprox(g2, this->eps));
		EXPECT_FALSE(g1.computeIsApprox(g3, this->eps));
	}

	/**
	 * @brief Tests the hat and vee operators for Rn.
	 * Verifies that applying the hat operator and then the vee operator returns the original Lie algebra element.
	 */
	TYPED_TEST(RnTest, HatVee) {
		using Rn = typename TestFixture::Rn;
		using Scalar = typename TestFixture::Scalar;
		using Vector = typename TestFixture::Vector;
		using Matrix = typename Rn::Matrix;
		using TangentVector = typename Rn::TangentVector;

		TangentVector v = TangentVector::Random();
		Matrix hat_v = Rn::computeHat(v);
		TangentVector vee_hat_v = Rn::computeVee(hat_v);
		EXPECT_TRUE(vee_hat_v.isApprox(v, this->eps));
	}

	/**
	 * @brief Tests the adjoint representation of the Rn group.
	 * Verifies that the adjoint matrix for Rn is the zero matrix.
	 */
	TYPED_TEST(RnTest, Adjoint) {
		using Rn = typename TestFixture::Rn;
		using Scalar = typename TestFixture::Scalar;
		using TangentVector = typename TestFixture::TangentVector;
		using AdjointMatrix = typename Rn::AdjointMatrix;

		TangentVector v = TangentVector::Random();
		AdjointMatrix Ad = Rn::computeAd(v);
		EXPECT_TRUE(Ad.isApprox(AdjointMatrix::Zero(), this->eps));
	}

	/**
	 * @brief Tests the random element generation for Rn.
	 * Verifies that a randomly generated Rn element is valid.
	 */
	TYPED_TEST(RnTest, Random) {
		using Rn = typename TestFixture::Rn;
		std::mt19937 gen(0);
		Rn r = Rn::computeRandom(gen);
		EXPECT_TRUE(r.computeIsValid());
	}

	/**
	 * @brief Tests the stream output operator for Rn.
	 * Verifies that the `print` method produces non-empty output.
	 */
	TYPED_TEST(RnTest, Print) {
		using Rn = typename TestFixture::Rn;
		using Vector = typename TestFixture::Vector;

		Vector v = Vector::Random();
		Rn g(v);
		std::stringstream ss;
		g.print(ss);
		EXPECT_FALSE(ss.str().empty());
	}

	/**
	 * @brief Tests the validity check for Rn elements.
	 * Verifies that a valid Rn element is correctly identified as valid, and an invalid one (e.g., containing
	 * NaN) is identified as invalid.
	 */
	TYPED_TEST(RnTest, IsValid) {
		using Rn = typename TestFixture::Rn;
		using Vector = typename TestFixture::Vector;

		Rn g(Vector::Random());
		EXPECT_TRUE(g.computeIsValid());

		// Test with invalid data (e.g., NaN)
		Vector invalid_v = Vector::Constant(std::numeric_limits<typename Rn::Scalar>::quiet_NaN());
		Rn invalid_g(invalid_v);
		EXPECT_FALSE(invalid_g.computeIsValid());
	}

	/**
	 * @brief Tests the normalization of Rn elements.
	 * Verifies that `computeNormalize()` does not alter the Rn element, as no normalization is needed for
	 * Rn.
	 */
	TYPED_TEST(RnTest, Normalize) {
		using Rn = typename TestFixture::Rn;
		using Vector = typename TestFixture::Vector;

		Vector v = Vector::Random();
		Rn g(v);
		g.computeNormalize();
		EXPECT_TRUE(g.computeLog().isApprox(v, this->eps)); // Rn normalize does nothing
	}

	/**
	 * @brief Tests the squared distance calculation for Rn elements.
	 * Verifies that the squared distance is correctly computed as the squared Euclidean norm of the difference between
	 * the underlying vectors.
	 */
	TYPED_TEST(RnTest, SquaredDistance) {
		using Rn = typename TestFixture::Rn;
		using Scalar = typename TestFixture::Scalar;
		using Vector = typename TestFixture::Vector;

		Vector v1 = Vector::Random();
		Vector v2 = Vector::Random();
		Rn g1(v1);
		Rn g2(v2);

		Scalar expected_sq_dist = (v1 - v2).squaredNorm();
		EXPECT_NEAR(g1.squaredDistance(g2), expected_sq_dist, this->eps);
	}

	/**
	 * @brief Tests the interpolation function for Rn elements.
	 * Verifies that linear interpolation is correctly performed between two Rn elements.
	 */
	TYPED_TEST(RnTest, Interpolate) {
		using Rn = typename TestFixture::Rn;
		using Scalar = typename TestFixture::Scalar;
		using Vector = typename TestFixture::Vector;

		Vector v1 = Vector::Random();
		Vector v2 = Vector::Random();
		Rn g1(v1);
		Rn g2(v2);

		Rn interpolated = g1.interpolate(g2, 0.5);
		EXPECT_TRUE(interpolated.computeLog().isApprox((v1 + v2) / 2.0, this->eps));

		EXPECT_TRUE(g1.interpolate(g2, 0.0).computeIsApprox(g1, this->eps));
		EXPECT_TRUE(g1.interpolate(g2, 1.0).computeIsApprox(g2, this->eps));
	}

} // namespace sofa::component::cosserat::liegroups::testing
