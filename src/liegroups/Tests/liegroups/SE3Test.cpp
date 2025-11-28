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

#include <Cosserat/liegroups/SE3.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>
#include <sofa/testing/BaseTest.h>

namespace sofa::component::cosserat::liegroups::testing {

	using namespace sofa::testing;

	/**
	 * @brief Test suite for SE3 Lie group.
	 * This test fixture provides common types and a small epsilon for floating-point comparisons.
	 * @tparam T The type of SE3 to test (e.g., SE3<double>).
	 */
	template<typename T>
	class SE3Test : public BaseTest {
	protected:
		using SE3 = T;
		using Scalar = typename SE3::Scalar;
		using Vector3 = Eigen::Matrix<Scalar, 3, 1>;
		using TangentVector = typename SE3::TangentVector;
		using AdjointMatrix = typename SE3::AdjointMatrix;
		using Matrix = typename SE3::Matrix;
		using SO3 = typename SE3::SO3;

		const Scalar eps = 1e-9;

		/**
		 * @brief Set up method for the test fixture.
		 * Called before each test.
		 */
		void SetUp() override {}
		/**
		 * @brief Tear down method for the test fixture.
		 * Called after each test.
		 */
		void TearDown() override {}
	};

	using SE3Types = ::testing::Types<SE3<double>>;
	TYPED_TEST_SUITE(SE3Test, SE3Types);

	/**
	 * @brief Tests the constructors of the SE3 class.
	 * Verifies that SE3 objects can be constructed from default, and rotation and translation representations.
	 */
	TYPED_TEST(SE3Test, Constructors) {
		using SE3 = typename TestFixture::SE3;
		using Scalar = typename TestFixture::Scalar;
		using Vector3 = typename TestFixture::Vector3;
		using SO3 = typename TestFixture::SO3;

		SE3 g1;
		EXPECT_TRUE(g1.computeIdentity().computeIsApprox(g1, this->eps));

		SO3 rot(Scalar(M_PI / 2.0), Vector3::UnitZ());
		Vector3 trans = Vector3::Random();
		SE3 g2(rot, trans);
		EXPECT_TRUE(g2.rotation().computeIsApprox(rot, this->eps));
		EXPECT_TRUE(g2.translation().isApprox(trans, this->eps));
	}

	/**
	 * @brief Tests the identity element of the SE3 group.
	 * Verifies that the `computeIdentity()` method returns an SE3 element with identity rotation and zero translation.
	 */
	TYPED_TEST(SE3Test, Identity) {
		using SE3 = typename TestFixture::SE3;
		using Scalar = typename TestFixture::Scalar;
		using Vector3 = typename TestFixture::Vector3;
		using SO3 = typename TestFixture::SO3;

		SE3 identity = SE3::computeIdentity();
		EXPECT_TRUE(identity.rotation().computeIsApprox(SO3::computeIdentity(), this->eps));
		EXPECT_TRUE(identity.translation().isApprox(Vector3::Zero(), this->eps));
	}

	/**
	 * @brief Tests the inverse operation of the SE3 group.
	 * Verifies that composing an SE3 element with its inverse results in the identity element.
	 */
	TYPED_TEST(SE3Test, Inverse) {
		using SE3 = typename TestFixture::SE3;
		using Scalar = typename TestFixture::Scalar;
		using Vector3 = typename TestFixture::Vector3;
		using SO3 = typename TestFixture::SO3;

		SO3 rot(Scalar(0.1), Vector3(1.0, 0.0, 0.0));
		Vector3 trans(1.0, 2.0, 3.0);
		SE3 g(rot, trans);
		SE3 inv_g = g.computeInverse();

		SE3 composed = g.compose(inv_g);
		EXPECT_TRUE(composed.computeIdentity().computeIsApprox(composed, this->eps));
	}

	/**
	 * @brief Tests the exponential and logarithmic maps of the SE3 group.
	 * Verifies that applying `computeExp` to a Lie algebra element and then `computeLog` to the resulting group element
	 * returns the original Lie algebra element.
	 */
	TYPED_TEST(SE3Test, ExpLog) {
		using SE3 = typename TestFixture::SE3;
		using Scalar = typename TestFixture::Scalar;
		using TangentVector = typename TestFixture::TangentVector;

		TangentVector twist = TangentVector::Random();
		// Ensure the rotation part is not too large to avoid multiple coverings issues with log
		twist.template tail<3>() *= 0.1;

		SE3 g = SE3::computeExp(twist);
		TangentVector log_g = g.computeLog();
		EXPECT_TRUE(log_g.isApprox(twist, this->eps));
	}

	/**
	 * @brief Tests the group action of SE3 on a point.
	 * Verifies that a point is correctly transformed by an SE3 element.
	 */
	TYPED_TEST(SE3Test, Action) {
		using SE3 = typename TestFixture::SE3;
		using Scalar = typename TestFixture::Scalar;
		using Vector3 = typename TestFixture::Vector3;
		using SO3 = typename TestFixture::SO3;

		SO3 rot(Scalar(M_PI / 2.0), Vector3::UnitZ()); // 90 deg around Z
		Vector3 trans(1.0, 0.0, 0.0);
		SE3 g(rot, trans);

		Vector3 point(1.0, 0.0, 0.0);

		// Rotate (1,0,0) -> (0,1,0), then translate (+1,0,0) -> (1,1,0)
		Vector3 transformed_point = g.computeAction(point);

		EXPECT_NEAR(transformed_point[0], 1.0, this->eps);
		EXPECT_NEAR(transformed_point[1], 1.0, this->eps);
		EXPECT_NEAR(transformed_point[2], 0.0, this->eps);
	}

	/**
	 * @brief Tests the approximate equality comparison for SE3 elements.
	 * Verifies that two SE3 elements are considered approximately equal within a given tolerance.
	 */
	TYPED_TEST(SE3Test, IsApprox) {
		using SE3 = typename TestFixture::SE3;
		using Scalar = typename TestFixture::Scalar;
		using Vector3 = typename TestFixture::Vector3;
		using SO3 = typename TestFixture::SO3;

		SO3 rot(Scalar(0.1), Vector3::UnitX());
		Vector3 trans(1.0, 2.0, 3.0);
		SE3 g1(rot, trans);
		SE3 g2(rot, trans);

		EXPECT_TRUE(g1.computeIsApprox(g2, this->eps));
	}

	/**
	 * @brief Tests the hat and vee operators for SE3.
	 * Verifies that applying the hat operator and then the vee operator returns the original Lie algebra element.
	 */
	TYPED_TEST(SE3Test, HatVee) {
		using SE3 = typename TestFixture::SE3;
		using Scalar = typename TestFixture::Scalar;
		using TangentVector = typename TestFixture::TangentVector;
		using Matrix = typename SE3::Matrix; // Note: SE3::Matrix is usually 4x4, but Hat returns 4x4 for se(3)

		TangentVector twist = TangentVector::Random();

		auto hat_twist = SE3::computeHat(twist);
		TangentVector vee_hat_twist = SE3::computeVee(hat_twist);
		EXPECT_TRUE(vee_hat_twist.isApprox(twist, this->eps));
	}

	/**
	 * @brief Tests the adjoint representation of the SE3 group.
	 * Verifies that the adjoint matrix is not zero for a non-identity transformation.
	 */
	TYPED_TEST(SE3Test, Adjoint) {
		using SE3 = typename TestFixture::SE3;
		using Scalar = typename TestFixture::Scalar;
		using Vector3 = typename TestFixture::Vector3;
		using SO3 = typename TestFixture::SO3;
		using AdjointMatrix = typename TestFixture::AdjointMatrix;

		SO3 rot(Scalar(M_PI / 4.0), Vector3::UnitX());
		Vector3 trans(1.0, 2.0, 3.0);
		SE3 g(rot, trans);

		AdjointMatrix Ad_g = g.computeAdjoint();
		EXPECT_FALSE(Ad_g.isZero(this->eps));
	}

	/**
	 * @brief Tests the random element generation for SE3.
	 * Verifies that a randomly generated SE3 element is valid.
	 */
	TYPED_TEST(SE3Test, Random) {
		using SE3 = typename TestFixture::SE3;
		std::mt19937 gen(0);
		SE3 r = SE3::computeRandom(gen);
		EXPECT_TRUE(r.computeIsValid());
	}

	/**
	 * @brief Tests the stream output operator for SE3.
	 * Verifies that the `print` method produces non-empty output.
	 */
	TYPED_TEST(SE3Test, Print) {
		using SE3 = typename TestFixture::SE3;
		using Scalar = typename TestFixture::Scalar;
		using Vector3 = typename TestFixture::Vector3;
		using SO3 = typename TestFixture::SO3;

		SO3 rot(Scalar(0.1), Vector3::UnitX());
		Vector3 trans(1.0, 2.0, 3.0);
		SE3 g(rot, trans);
		std::stringstream ss;
		g.print(ss);
		EXPECT_FALSE(ss.str().empty());
	}

	/**
	 * @brief Tests the validity check for SE3 elements.
	 * Verifies that a valid SE3 element is correctly identified as valid.
	 */
	TYPED_TEST(SE3Test, IsValid) {
		using SE3 = typename TestFixture::SE3;
		using Scalar = typename TestFixture::Scalar;
		using Vector3 = typename TestFixture::Vector3;
		using SO3 = typename TestFixture::SO3;

		SO3 rot(Scalar(0.1), Vector3::UnitX());
		Vector3 trans(1.0, 2.0, 3.0);
		SE3 g(rot, trans);
		EXPECT_TRUE(g.computeIsValid());
	}

	/**
	 * @brief Tests the normalization of SE3 elements.
	 * Verifies that the rotation component is normalized after calling `computeNormalize()`.
	 */
	TYPED_TEST(SE3Test, Normalize) {
		using SE3 = typename TestFixture::SE3;
		using Scalar = typename TestFixture::Scalar;
		using Vector3 = typename TestFixture::Vector3;
		using SO3 = typename TestFixture::SO3;
		using Quaternion = typename SO3::Quaternion;

		// Create non-normalized rotation
		Quaternion q(0.1, 0.2, 0.3, 0.4);
		SO3 rot(q);
		Vector3 trans(1.0, 2.0, 3.0);

		SE3 g(rot, trans);
		g.computeNormalize();
		EXPECT_TRUE(g.rotation().computeIsValid());
	}

} // namespace sofa::component::cosserat::liegroups::testing
