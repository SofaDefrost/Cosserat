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
#include <Cosserat/liegroups/SO2.h>
#include <Eigen/Core>
#include <Eigen/Geometry>

namespace sofa::component::cosserat::liegroups::testing {

using namespace sofa::testing;

/**
 * @brief Test suite for SO2 Lie group.
 * This test fixture provides common types and a small epsilon for floating-point comparisons.
 * @tparam T The type of SO2 to test (e.g., SO2<double>).
 */
template <typename T>
class SO2Test : public BaseTest
{
protected:
    using SO2 = T;
    using Scalar = typename SO2::Scalar;
    using Vector = typename SO2::Vector;
    using TangentVector = typename SO2::TangentVector;
    using AdjointMatrix = typename SO2::AdjointMatrix;
    using Complex = typename SO2::Complex;

    const Scalar eps = 1e-9;
    static constexpr Scalar M_PI_VAL = 3.14159265358979323846;

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

using SO2Types = ::testing::Types<SO2<double>>;
TYPED_TEST_SUITE(SO2Test, SO2Types);

/**
 * @brief Tests the constructors of the SO2 class.
 * Verifies that SO2 objects can be constructed from default and angle representations.
 */
TYPED_TEST(SO2Test, Constructors)
{
    using SO2 = typename TestFixture::SO2;
    using Scalar = typename TestFixture::Scalar;

    SO2 g1;
    EXPECT_TRUE(g1.computeIdentity().computeIsApprox(g1, this->eps));

    SO2 g2(Scalar(TestFixture::M_PI_VAL/2.0));
    EXPECT_NEAR(g2.angle(), Scalar(TestFixture::M_PI_VAL/2.0), this->eps);
}

/**
 * @brief Tests the identity element of the SO2 group.
 * Verifies that the `identity()` method returns an angle of 0.
 */
TYPED_TEST(SO2Test, Identity)
{
    using SO2 = typename TestFixture::SO2;
    using Scalar = typename TestFixture::Scalar;

    SO2 identity = SO2::computeIdentity();
    EXPECT_NEAR(identity.angle(), 0.0, this->eps);
}

/**
 * @brief Tests the inverse operation of the SO2 group.
 * Verifies that composing an SO2 element with its inverse results in the identity element.
 */
TYPED_TEST(SO2Test, Inverse)
{
    using SO2 = typename TestFixture::SO2;
    using Scalar = typename TestFixture::Scalar;

    SO2 g(Scalar(TestFixture::M_PI_VAL/3.0));
    SO2 inv_g = g.computeInverse();

    SO2 composed = g.compose(inv_g);
    EXPECT_TRUE(composed.computeIdentity().computeIsApprox(composed, this->eps));
}

/**
 * @brief Tests the exponential and logarithmic maps of the SO2 group.
 * Verifies that applying `exp` to a Lie algebra element and then `log` to the resulting group element returns the original Lie algebra element.
 */
TYPED_TEST(SO2Test, ExpLog)
{
    using SO2 = typename TestFixture::SO2;
    using Scalar = typename TestFixture::Scalar;
    using TangentVector = typename TestFixture::TangentVector;

    TangentVector omega;
    omega[0] = 0.5; // Angle

    SO2 g = SO2::computeExp(omega);
    TangentVector log_g = g.computeLog();
    EXPECT_NEAR(log_g[0], omega[0], this->eps);
}

/**
 * @brief Tests the group action of SO2 on a 2D point.
 * Verifies that a point is correctly rotated by an SO2 element.
 */
TYPED_TEST(SO2Test, Action)
{
    using SO2 = typename TestFixture::SO2;
    using Scalar = typename TestFixture::Scalar;
    using Vector = typename TestFixture::Vector;

    SO2 g(Scalar(TestFixture::M_PI_VAL/2.0)); // 90 degrees rotation
    Vector point;
    point << 1.0, 0.0;

    Vector transformed_point = g.computeAction(point);
    EXPECT_NEAR(transformed_point[0], 0.0, this->eps);
    EXPECT_NEAR(transformed_point[1], 1.0, this->eps);
}

/**
 * @brief Tests the approximate equality comparison for SO2 elements.
 * Verifies that two SO2 elements are considered approximately equal within a given tolerance.
 */
TYPED_TEST(SO2Test, IsApprox)
{
    using SO2 = typename TestFixture::SO2;
    using Scalar = typename TestFixture::Scalar;

    SO2 g1(Scalar(TestFixture::M_PI_VAL/4.0));
    SO2 g2(Scalar(TestFixture::M_PI_VAL/4.0) + this->eps/2.0);
    SO2 g3(Scalar(TestFixture::M_PI_VAL/4.0) + this->eps*2.0);

    EXPECT_TRUE(g1.computeIsApprox(g2, this->eps));
    EXPECT_FALSE(g1.computeIsApprox(g3, this->eps));
}

/**
 * @brief Tests the hat and vee operators for SO2.
 * Verifies that applying the hat operator and then the vee operator returns the original Lie algebra element.
 */
TYPED_TEST(SO2Test, HatVee)
{
    using SO2 = typename TestFixture::SO2;
    using Scalar = typename TestFixture::Scalar;
    using TangentVector = typename TestFixture::TangentVector;
    using Matrix = typename SO2::Matrix;

    TangentVector omega;
    omega[0] = 0.5;

    Matrix hat_omega = SO2::computeHat(omega);
    TangentVector vee_hat_omega = SO2::computeVee(hat_omega);
    EXPECT_NEAR(vee_hat_omega[0], omega[0], this->eps);
}

/**
 * @brief Tests the adjoint representation of the SO2 group.
 * Verifies that the adjoint matrix is identity for SO2 and zero for the Lie algebra element.
 */
TYPED_TEST(SO2Test, Adjoint)
{
    using SO2 = typename TestFixture::SO2;
    using Scalar = typename TestFixture::Scalar;
    using TangentVector = typename TestFixture::TangentVector;
    using AdjointMatrix = typename TestFixture::AdjointMatrix;

    SO2 g(Scalar(TestFixture::M_PI_VAL/4.0));
    TangentVector omega;
    omega[0] = 0.1;

    AdjointMatrix Ad_g = g.computeAdjoint();
    AdjointMatrix ad_omega = SO2::computeAd(omega);

    EXPECT_TRUE(Ad_g.isApprox(AdjointMatrix::Identity(), this->eps));
    EXPECT_TRUE(ad_omega.isApprox(AdjointMatrix::Zero(), this->eps));
}

/**
 * @brief Tests the random element generation for SO2.
 * Verifies that a randomly generated SO2 element is valid.
 */
TYPED_TEST(SO2Test, Random)
{
    using SO2 = typename TestFixture::SO2;
    std::mt19937 gen(0);
    SO2 r = SO2::computeRandom(gen);
    EXPECT_TRUE(r.computeIsValid());
}

/**
 * @brief Tests the stream output operator for SO2.
 * Verifies that the `print` method produces non-empty output.
 */
TYPED_TEST(SO2Test, Print)
{
    using SO2 = typename TestFixture::SO2;
    using Scalar = typename TestFixture::Scalar;

    SO2 g(Scalar(TestFixture::M_PI_VAL/4.0));
    std::stringstream ss;
    g.print(ss);
    EXPECT_FALSE(ss.str().empty());
}

/**
 * @brief Tests the validity check for SO2 elements.
 * Verifies that a valid SO2 element is correctly identified as valid.
 */
TYPED_TEST(SO2Test, IsValid)
{
    using SO2 = typename TestFixture::SO2;
    using Scalar = typename TestFixture::Scalar;

    SO2 g(Scalar(TestFixture::M_PI_VAL/4.0));
    EXPECT_TRUE(g.computeIsValid());
}

/**
 * @brief Tests the normalization of SO2 elements.
 * Verifies that an angle outside the [-π, π] range is correctly normalized after calling `computeNormalize()`.
 */
TYPED_TEST(SO2Test, Normalize)
{
    using SO2 = typename TestFixture::SO2;
    using Scalar = typename TestFixture::Scalar;

    SO2 g(Scalar(TestFixture::M_PI_VAL * 2.5)); // Angle > 2*PI
    g.computeNormalize();
    EXPECT_NEAR(g.angle(), Scalar(TestFixture::M_PI_VAL/2.0), this->eps);
}

} // namespace sofa::component::cosserat::liegroups::testing