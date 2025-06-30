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

#include <cmath>
#include <sofa/testing/BaseTest.h>
#include <Cosserat/liegroups/SO3.h>
#include <Eigen/Core>
#include <Eigen/Geometry>


namespace sofa::component::cosserat::liegroups::testing {

using namespace sofa::testing;

/**
 * @brief Test suite for SO3 Lie group.
 * This test fixture provides common types and a small epsilon for floating-point comparisons.
 * @tparam T The type of SO3 to test (e.g., SO3<double>).
 */
template <typename T>
class SO3Test : public BaseTest
{
protected:
    using SO3 = T;
    using Scalar = typename SO3::Scalar;
    using Vector = typename SO3::Vector;
    using Matrix = typename SO3::Matrix;
    using TangentVector = typename SO3::TangentVector;
    using AdjointMatrix = typename SO3::AdjointMatrix;
    using Quaternion = typename SO3::Quaternion;

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

using SO3Types = ::testing::Types<SO3<double>>;
TYPED_TEST_SUITE(SO3Test, SO3Types);

/**
 * @brief Tests the constructors of the SO3 class.
 * Verifies that SO3 objects can be constructed from default, angle-axis, quaternion, and matrix representations.
 */
TYPED_TEST(SO3Test, Constructors)
{
    using SO3 = typename TestFixture::SO3;
    using Scalar = typename TestFixture::Scalar;
    using Vector = typename TestFixture::Vector;
    using Quaternion = typename TestFixture::Quaternion;

    SO3 g1;
    EXPECT_TRUE(g1.computeIdentity().computeIsApprox(g1, this->eps));

    SO3 g2(Scalar(M_PI/2.0), Vector::UnitZ());
    EXPECT_TRUE(g2.computeIsValid());

    Quaternion q(0.707, 0.0, 0.0, 0.707); // 90 deg around Z
    SO3 g3(q);
    EXPECT_TRUE(g3.computeIsValid());

    Matrix m = g2.matrix();
    SO3 g4(m);
    EXPECT_TRUE(g4.computeIsValid());
}

/**
 * @brief Tests the identity element of the SO3 group.
 * Verifies that the `identity()` method returns a quaternion that is approximately the identity quaternion.
 */
TYPED_TEST(SO3Test, Identity)
{
    using SO3 = typename TestFixture::SO3;
    using Scalar = typename TestFixture::Scalar;
    using Quaternion = typename TestFixture::Quaternion;

    SO3 identity = SO3::computeIdentity();
    EXPECT_TRUE(identity.quaternion().isApprox(Quaternion::Identity(), this->eps));
}

/**
 * @brief Tests the inverse operation of the SO3 group.
 * Verifies that composing an SO3 element with its inverse results in the identity element.
 */
TYPED_TEST(SO3Test, Inverse)
{
    using SO3 = typename TestFixture::SO3;
    using Scalar = typename TestFixture::Scalar;
    using Vector = typename TestFixture::Vector;

    SO3 g(Scalar(M_PI/3.0), Vector::UnitX());
    SO3 inv_g = g.computeInverse();

    SO3 composed = g.compose(inv_g);
    EXPECT_TRUE(composed.computeIdentity().computeIsApprox(composed, this->eps));
}

/**
 * @brief Tests the exponential and logarithmic maps of the SO3 group.
 * Verifies that applying `exp` to a Lie algebra element and then `log` to the resulting group element returns the original Lie algebra element.
 */
TYPED_TEST(SO3Test, ExpLog)
{
    using SO3 = typename TestFixture::SO3;
    using Scalar = typename TestFixture::Scalar;
    using TangentVector = typename TestFixture::TangentVector;

    TangentVector omega;
    omega << 0.1, 0.2, 0.3; // Angular velocity

    SO3 g = SO3::computeExp(omega);
    TangentVector log_g = g.computeLog();
    EXPECT_TRUE(log_g.isApprox(omega, this->eps));
}

/**
 * @brief Tests the group action of SO3 on a 3D point.
 * Verifies that a point is correctly rotated by an SO3 element.
 */
TYPED_TEST(SO3Test, Action)
{
    using SO3 = typename TestFixture::SO3;
    using Scalar = typename TestFixture::Scalar;
    using Vector = typename TestFixture::Vector;

        SO3 g(Scalar(M_PI * 2.5)); // Angle > 2*PI
    Vector point(1.0, 0.0, 0.0);

    Vector transformed_point = g.computeAction(point);
    EXPECT_NEAR(transformed_point[0], 0.0, this->eps);
    EXPECT_NEAR(transformed_point[1], 1.0, this->eps);
    EXPECT_NEAR(transformed_point[2], 0.0, this->eps);
}

/**
 * @brief Tests the approximate equality comparison for SO3 elements.
 * Verifies that two SO3 elements are considered approximately equal within a given tolerance.
 */
TYPED_TEST(SO3Test, IsApprox)
{
    using SO3 = typename TestFixture::SO3;
    using Scalar = typename TestFixture::Scalar;
    using Vector = typename TestFixture::Vector;

    SO3 g1(Scalar(M_PI/4.0), Vector::UnitX());
    SO3 g2(Scalar(M_PI/4.0) + this->eps/2.0, Vector::UnitX());
    SO3 g3(Scalar(M_PI/4.0) + this->eps*2.0, Vector::UnitX());

    EXPECT_TRUE(g1.computeIsApprox(g2, this->eps));
    EXPECT_FALSE(g1.computeIsApprox(g3, this->eps));
}

/**
 * @brief Tests the hat and vee operators for SO3.
 * Verifies that applying the hat operator and then the vee operator returns the original Lie algebra element.
 */
TYPED_TEST(SO3Test, HatVee)
{
    using SO3 = typename TestFixture::SO3;
    using Scalar = typename TestFixture::Scalar;
    using TangentVector = typename TestFixture::TangentVector;
    using Matrix = typename SO3::Matrix;

    TangentVector omega;
    omega << 0.1, 0.2, 0.3;

    Matrix hat_omega = SO3::computeHat(omega);
    TangentVector vee_hat_omega = SO3::computeVee(hat_omega);
    EXPECT_TRUE(vee_hat_omega.isApprox(omega, this->eps));
}

/**
 * @brief Tests the adjoint representation of the SO3 group.
 * Verifies that the adjoint matrix is not zero for a non-identity rotation.
 */
TYPED_TEST(SO3Test, Adjoint)
{
    using SO3 = typename TestFixture::SO3;
    using Scalar = typename TestFixture::Scalar;
    using TangentVector = typename TestFixture::TangentVector;
    using AdjointMatrix = typename SO3::AdjointMatrix;

    SO3 g(Scalar(M_PI/4.0), Vector::UnitX());
    TangentVector omega;
    omega << 0.1, 0.2, 0.3;

    AdjointMatrix Ad_g = g.computeAdjoint();
    AdjointMatrix ad_omega = SO3::computeAd(omega);

    EXPECT_FALSE(Ad_g.isZero(this->eps));
    EXPECT_FALSE(ad_omega.isZero(this->eps));
}

/**
 * @brief Tests the random element generation for SO3.
 * Verifies that a randomly generated SO3 element is valid.
 */
TYPED_TEST(SO3Test, Random)
{
    using SO3 = typename TestFixture::SO3;
    std::mt19937 gen(0);
    SO3 r = SO3::computeRandom(gen);
    EXPECT_TRUE(r.computeIsValid());
}

/**
 * @brief Tests the stream output operator for SO3.
 * Verifies that the `print` method produces non-empty output.
 */
TYPED_TEST(SO3Test, Print)
{
    using SO3 = typename TestFixture::SO3;
    using Scalar = typename TestFixture::Scalar;
    using Vector = typename TestFixture::Vector;

    SO3 g(Scalar(M_PI/4.0), Vector::UnitX());
    std::stringstream ss;
    g.print(ss);
    EXPECT_FALSE(ss.str().empty());
}

/**
 * @brief Tests the validity check for SO3 elements.
 * Verifies that a valid SO3 element is correctly identified as valid.
 */
TYPED_TEST(SO3Test, IsValid)
{
    using SO3 = typename TestFixture::SO3;
    using Scalar = typename TestFixture::Scalar;
    using Vector = typename TestFixture::Vector;

    SO3 g(Scalar(M_PI/4.0), Vector::UnitX());
    EXPECT_TRUE(g.computeIsValid());
}

/**
 * @brief Tests the normalization of SO3 elements.
 * Verifies that a non-normalized quaternion is correctly normalized after calling `computeNormalize()`.
 */
TYPED_TEST(SO3Test, Normalize)
{
    using SO3 = typename TestFixture::SO3;
    using Scalar = typename TestFixture::Scalar;
    using Vector = typename TestFixture::Vector;
    using Quaternion = typename TestFixture::Quaternion;

    Quaternion q(0.1, 0.2, 0.3, 0.4); // Non-normalized quaternion
    SO3 g(q);
    g.computeNormalize();
    EXPECT_NEAR(g.quaternion().norm(), 1.0, this->eps);
}

} // namespace sofa::component::cosserat::liegroups::testing
