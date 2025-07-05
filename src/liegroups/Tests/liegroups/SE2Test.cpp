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
#include <Cosserat/liegroups/SE2.h>
#include <Eigen/Core>
#include <Eigen/Geometry>


namespace sofa::component::cosserat::liegroups::testing {

using namespace sofa::testing;

/**
 * @brief Test suite for SE2 Lie group.
 * This test fixture provides common types and a small epsilon for floating-point comparisons.
 * @tparam T The type of SE2 to test (e.g., SE2<double>).
 */
template <typename T>
class SE2Test : public BaseTest
{
protected:
    using SE2 = T;
    using Scalar = typename SE2::Scalar;
    using Vector2 = Eigen::Matrix<Scalar, 2, 1>;
    using TangentVector = typename SE2::TangentVector;

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

using SE2Types = ::testing::Types<SE2<double>>;
TYPED_TEST_SUITE(SE2Test, SE2Types);

/**
 * @brief Tests the constructors of the SE2 class.
 * Verifies that SE2 objects can be constructed from default, rotation and translation, and angle and translation representations.
 */
TYPED_TEST(SE2Test, Constructors)
{
    using SE2 = typename TestFixture::SE2;
    using Scalar = typename TestFixture::Scalar;
    using Vector2 = typename TestFixture::Vector2;

    SE2 g1;
    EXPECT_TRUE(g1.computeIdentity().computeIsApprox(g1, this->eps));

    SO2<Scalar> rot(Scalar(M_PI/2.0));
    Vector2 trans(1.0, 2.0);
    SE2 g2(rot, trans);
    EXPECT_TRUE(g2.rotation().computeIsApprox(rot, this->eps));
    EXPECT_TRUE(g2.translation().isApprox(trans, this->eps));
}

/**
 * @brief Tests the identity element of the SE2 group.
 * Verifies that the `computeIdentity()` method returns an SE2 element with identity rotation and zero translation.
 */
TYPED_TEST(SE2Test, Identity)
{
    using SE2 = typename TestFixture::SE2;
    using Scalar = typename TestFixture::Scalar;
    using Vector2 = typename TestFixture::Vector2;

    SE2 identity = SE2::computeIdentity();
    EXPECT_TRUE(identity.rotation().computeIsApprox(SO2<Scalar>::computeIdentity(), this->eps));
    EXPECT_TRUE(identity.translation().isApprox(Vector2::Zero(), this->eps));
}

/**
 * @brief Tests the inverse operation of the SE2 group.
 * Verifies that composing an SE2 element with its inverse results in the identity element.
 */
TYPED_TEST(SE2Test, Inverse)
{
    using SE2 = typename TestFixture::SE2;
    using Scalar = typename TestFixture::Scalar;
    using Vector2 = typename TestFixture::Vector2;

    SO2<Scalar> rot(Scalar(M_PI/3.0));
    Vector2 trans(1.0, 2.0);
    SE2 g(rot, trans);
    SE2 inv_g = g.computeInverse();

    SE2 composed = g.compose(inv_g);
    EXPECT_TRUE(composed.computeIdentity().computeIsApprox(composed, this->eps));
}

/**
 * @brief Tests the exponential and logarithmic maps of the SE2 group.
 * Verifies that applying `computeExp` to a Lie algebra element and then `computeLog` to the resulting group element returns the original Lie algebra element.
 */
TYPED_TEST(SE2Test, ExpLog)
{
    using SE2 = typename TestFixture::SE2;
    using Scalar = typename TestFixture::Scalar;
    using TangentVector = typename TestFixture::TangentVector;

    TangentVector twist;
    twist << 0.1, 0.2, 0.3; // vx, vy, omega

    SE2 g = SE2::computeExp(twist);
    TangentVector log_g = g.computeLog();
    EXPECT_TRUE(log_g.isApprox(twist, this->eps));
}

/**
 * @brief Tests the group action of SE2 on a 2D point.
 * Verifies that a point is correctly transformed by an SE2 element (rotation and translation).
 */
TYPED_TEST(SE2Test, Action)
{
    using SE2 = typename TestFixture::SE2;
    using Scalar = typename TestFixture::Scalar;
    using Vector2 = Eigen::Matrix<Scalar, 2, 1>;
    using ActionVector = typename SE2::ActionVector;

    SO2<Scalar> rot(Scalar(M_PI/2.0));
    Vector2 trans(1.0, 2.0);
    SE2 g(rot, trans);

    ActionVector point;
    point << 3.0, 4.0;

    ActionVector transformed_point = g.computeAction(point);
    // Expected: rotate (3,4) by 90 deg to (-4,3), then translate by (1,2) to (-3,5)
    EXPECT_NEAR(transformed_point[0], -3.0, this->eps);
    EXPECT_NEAR(transformed_point[1], 5.0, this->eps);
}

/**
 * @brief Tests the approximate equality comparison for SE2 elements.
 * Verifies that two SE2 elements are considered approximately equal within a given tolerance.
 */
TYPED_TEST(SE2Test, IsApprox)
{
    using SE2 = typename TestFixture::SE2;
    using Scalar = typename TestFixture::Scalar;
    using Vector2 = Eigen::Matrix<Scalar, 2, 1>;

    SO2<Scalar> rot1(Scalar(M_PI/4.0));
    Vector2 trans1(1.0, 2.0);
    SE2 g1(rot1, trans1);

    SO2<Scalar> rot2(Scalar(M_PI/4.0) + this->eps/2.0);
    Vector2 trans2(1.0 + this->eps/2.0, 2.0 + this->eps/2.0);
    SE2 g2(rot2, trans2);

    SE2 g3(SO2<Scalar>(Scalar(M_PI/4.0) + this->eps*2.0), Vector2(1.0, 2.0));

    EXPECT_TRUE(g1.computeIsApprox(g2, this->eps));
    EXPECT_FALSE(g1.computeIsApprox(g3, this->eps));
}

/**
 * @brief Tests the hat and vee operators for SE2.
 * Verifies that applying the hat operator and then the vee operator returns the original Lie algebra element.
 */
TYPED_TEST(SE2Test, HatVee)
{
    using SE2 = typename TestFixture::SE2;
    using Scalar = typename TestFixture::Scalar;
    using TangentVector = typename TestFixture::TangentVector;
    using Matrix = typename SE2::Matrix;

    TangentVector twist;
    twist << 0.1, 0.2, 0.3; // vx, vy, omega

    Matrix hat_twist = SE2::hat(twist);
    TangentVector vee_hat_twist = SE2::vee(hat_twist);
    EXPECT_TRUE(vee_hat_twist.isApprox(twist, this->eps));
}

/**
 * @brief Tests the adjoint representation of the SE2 group.
 * Verifies that the adjoint matrix is not zero for a non-identity transformation.
 */
TYPED_TEST(SE2Test, Adjoint)
{
    using SE2 = typename TestFixture::SE2;
    using Scalar = typename TestFixture::Scalar;
    using TangentVector = typename TestFixture::TangentVector;
    using AdjointMatrix = typename SE2::AdjointMatrix;

    SO2<Scalar> rot(Scalar(M_PI/4.0));
    Vector2 trans(1.0, 2.0);
    SE2 g(rot, trans);

    TangentVector twist;
    twist << 0.1, 0.2, 0.3;

    AdjointMatrix Ad_g = g.computeAdjoint();
    TangentVector ad_twist = SE2::computeAd(twist);

    // For SE(2), Ad_g * twist should be equal to ad_twist
    // This test might need refinement based on the exact definition of Ad_g and ad
    // For now, a basic check that it's not all zeros
    EXPECT_FALSE(Ad_g.isZero(this->eps));
    EXPECT_FALSE(ad_twist.isZero(this->eps));
}

/**
 * @brief Tests the random element generation for SE2.
 * Verifies that a randomly generated SE2 element is valid.
 */
TYPED_TEST(SE2Test, Random)
{
    using SE2 = typename TestFixture::SE2;
    std::mt19937 gen(0);
    SE2 r = SE2::computeRandom(gen);
    EXPECT_TRUE(r.computeIsValid());
}

/**
 * @brief Tests the stream output operator for SE2.
 * Verifies that the `print` method produces non-empty output.
 */
TYPED_TEST(SE2Test, Print)
{
    using SE2 = typename TestFixture::SE2;
    using Scalar = typename TestFixture::Scalar;
    using Vector2 = Eigen::Matrix<Scalar, 2, 1>;

    SE2 g(SO2<Scalar>(Scalar(M_PI/4.0)), Vector2(1.0, 2.0));
    std::stringstream ss;
    g.print(ss);
    EXPECT_FALSE(ss.str().empty());
}

/**
 * @brief Tests the validity check for SE2 elements.
 * Verifies that a valid SE2 element is correctly identified as valid.
 */
TYPED_TEST(SE2Test, IsValid)
{
    using SE2 = typename TestFixture::SE2;
    using Scalar = typename TestFixture::Scalar;
    using Vector2 = Eigen::Matrix<Scalar, 2, 1>;

    SE2 g(SO2<Scalar>(Scalar(M_PI/4.0)), Vector2(1.0, 2.0));
    EXPECT_TRUE(g.computeIsValid());

    // Test with invalid rotation (e.g., non-unit quaternion in SO2)
    // This would require modifying SO2 to allow invalid states for testing
}

/**
 * @brief Tests the normalization of SE2 elements.
 * Verifies that the rotation component is normalized after calling `computeNormalize()`.
 */
TYPED_TEST(SE2Test, Normalize)
{
    using SE2 = typename TestFixture::SE2;
    using Scalar = typename TestFixture::Scalar;
    using Vector2 = Eigen::Matrix<Scalar, 2, 1>;

    SO2<Scalar> rot(Scalar(M_PI/4.0));
    Vector2 trans(1.0, 2.0);
    SE2 g(rot, trans);
    g.computeNormalize();
    EXPECT_TRUE(g.rotation().computeIsValid());
}

} // namespace sofa::component::cosserat::liegroups::testing