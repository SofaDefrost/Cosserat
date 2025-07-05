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

#include <sofa/testing/BaseTest.h>
#include <Cosserat/liegroups/RealSpace.h>
#include <Eigen/Core>

namespace sofa::component::cosserat::liegroups::testing {

using namespace sofa::testing;

/**
 * @brief Test suite for RealSpace Lie group.
 * This test fixture provides common types and a small epsilon for floating-point comparisons.
 * It is templated to allow testing RealSpace with different dimensions.
 * @tparam T The type of RealSpace to test (e.g., RealSpace<double, 3>).
 */
template <typename T>
class RealSpaceTest : public BaseTest
{
protected:
    using RealSpace = T;
    using Scalar = typename RealSpace::Scalar;
    using Vector = typename RealSpace::Vector;
    using TangentVector = typename RealSpace::TangentVector;

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

using RealSpaceTypes = ::testing::Types<RealSpace<double, 1>, RealSpace<double, 2>, RealSpace<double, 3>>;
TYPED_TEST_SUITE(RealSpaceTest, RealSpaceTypes);

/**
 * @brief Tests the constructors of the RealSpace class.
 * Verifies that RealSpace objects can be constructed from default (zero vector) and from an Eigen vector.
 */
TYPED_TEST(RealSpaceTest, Constructors)
{
    using RealSpace = typename TestFixture::RealSpace;
    using Scalar = typename TestFixture::Scalar;
    using Vector = typename TestFixture::Vector;

    RealSpace g1;
    EXPECT_TRUE(g1.computeIdentity().computeIsApprox(g1, this->eps));

    Vector v = Vector::Random();
    RealSpace g2(v);
    EXPECT_TRUE(g2.data().isApprox(v, this->eps));
}

/**
 * @brief Tests the identity element of the RealSpace group.
 * Verifies that the `computeIdentity()` method returns a zero vector.
 */
TYPED_TEST(RealSpaceTest, Identity)
{
    using RealSpace = typename TestFixture::RealSpace;
    using Scalar = typename TestFixture::Scalar;
    using Vector = typename TestFixture::Vector;

    RealSpace identity = RealSpace::computeIdentity();
    EXPECT_TRUE(identity.data().isApprox(Vector::Zero(), this->eps));
}

/**
 * @brief Tests the inverse operation of the RealSpace group.
 * Verifies that the inverse of a vector is its negation, and that composing a vector with its inverse results in the identity element (zero vector).
 */
TYPED_TEST(RealSpaceTest, Inverse)
{
    using RealSpace = typename TestFixture::RealSpace;
    using Scalar = typename TestFixture::Scalar;
    using Vector = typename TestFixture::Vector;

    Vector v = Vector::Random();
    RealSpace g(v);
    RealSpace inv_g = g.computeInverse();
    EXPECT_TRUE(inv_g.data().isApprox(-v, this->eps));

    RealSpace composed = g.compose(inv_g);
    EXPECT_TRUE(composed.computeIdentity().computeIsApprox(composed, this->eps));
}

/**
 * @brief Tests the exponential and logarithmic maps of the RealSpace group.
 * Verifies that applying `computeExp` to a Lie algebra element and then `computeLog` to the resulting group element returns the original Lie algebra element.
 * For RealSpace, both `exp` and `log` are identity functions.
 */
TYPED_TEST(RealSpaceTest, ExpLog)
{
    using RealSpace = typename TestFixture::RealSpace;
    using Scalar = typename TestFixture::Scalar;
    using TangentVector = typename TestFixture::TangentVector;

    TangentVector v = TangentVector::Random();
    RealSpace g = RealSpace::computeExp(v);
    EXPECT_TRUE(g.computeLog().isApprox(v, this->eps));
}

/**
 * @brief Tests the group action of RealSpace on a point.
 * Verifies that the action is equivalent to vector addition.
 */
TYPED_TEST(RealSpaceTest, Action)
{
    using RealSpace = typename TestFixture::RealSpace;
    using Scalar = typename TestFixture::Scalar;
    using Vector = typename TestFixture::Vector;

    Vector g_data = Vector::Random();
    RealSpace g(g_data);
    Vector point = Vector::Random();
    Vector transformed_point = g.computeAction(point);
    EXPECT_TRUE(transformed_point.isApprox(point + g_data, this->eps));
}

/**
 * @brief Tests the approximate equality comparison for RealSpace elements.
 * Verifies that two RealSpace elements are considered approximately equal within a given tolerance.
 */
TYPED_TEST(RealSpaceTest, IsApprox)
{
    using RealSpace = typename TestFixture::RealSpace;
    using Scalar = typename TestFixture::Scalar;
    using Vector = typename TestFixture::Vector;

    Vector v = Vector::Random();
    RealSpace g1(v);
    RealSpace g2(v + Vector::Constant(this->eps / 2.0));
    RealSpace g3(v + Vector::Constant(this->eps * 2.0));

    EXPECT_TRUE(g1.computeIsApprox(g2, this->eps));
    EXPECT_FALSE(g1.computeIsApprox(g3, this->eps));
}

/**
 * @brief Tests the hat and vee operators for RealSpace.
 * Verifies that applying the hat operator and then the vee operator returns the original Lie algebra element.
 */
TYPED_TEST(RealSpaceTest, HatVee)
{
    using RealSpace = typename TestFixture::RealSpace;
    using Scalar = typename TestFixture::Scalar;
    using Vector = typename TestFixture::Vector;
    using Matrix = typename RealSpace::Matrix;
    using TangentVector = typename RealSpace::TangentVector;

    TangentVector v = TangentVector::Random();
    Matrix hat_v = RealSpace::computeHat(v);
    TangentVector vee_hat_v = RealSpace::computeVee(hat_v);
    EXPECT_TRUE(vee_hat_v.isApprox(v, this->eps));
}

/**
 * @brief Tests the adjoint representation of the RealSpace group.
 * Verifies that the adjoint matrix for RealSpace is the zero matrix.
 */
TYPED_TEST(RealSpaceTest, Adjoint)
{
    using RealSpace = typename TestFixture::RealSpace;
    using Scalar = typename TestFixture::Scalar;
    using TangentVector = typename TestFixture::TangentVector;
    using AdjointMatrix = typename RealSpace::AdjointMatrix;

    TangentVector v = TangentVector::Random();
    AdjointMatrix Ad = RealSpace::computeAd(v);
    EXPECT_TRUE(Ad.isApprox(AdjointMatrix::Zero(), this->eps));
}

/**
 * @brief Tests the random element generation for RealSpace.
 * Verifies that a randomly generated RealSpace element is valid.
 */
TYPED_TEST(RealSpaceTest, Random)
{
    using RealSpace = typename TestFixture::RealSpace;
    std::mt19937 gen(0);
    RealSpace r = RealSpace::computeRandom(gen);
    EXPECT_TRUE(r.computeIsValid());
}

/**
 * @brief Tests the stream output operator for RealSpace.
 * Verifies that the `print` method produces non-empty output.
 */
TYPED_TEST(RealSpaceTest, Print)
{
    using RealSpace = typename TestFixture::RealSpace;
    using Vector = typename TestFixture::Vector;

    Vector v = Vector::Random();
    RealSpace g(v);
    std::stringstream ss;
    g.print(ss);
    EXPECT_FALSE(ss.str().empty());
}

/**
 * @brief Tests the validity check for RealSpace elements.
 * Verifies that a valid RealSpace element is correctly identified as valid, and an invalid one (e.g., containing NaN) is identified as invalid.
 */
TYPED_TEST(RealSpaceTest, IsValid)
{
    using RealSpace = typename TestFixture::RealSpace;
    using Vector = typename TestFixture::Vector;

    RealSpace g(Vector::Random());
    EXPECT_TRUE(g.computeIsValid());

    // Test with invalid data (e.g., NaN)
    Vector invalid_v = Vector::Constant(std::numeric_limits<typename RealSpace::Scalar>::quiet_NaN());
    RealSpace invalid_g(invalid_v);
    EXPECT_FALSE(invalid_g.computeIsValid());
}

/**
 * @brief Tests the normalization of RealSpace elements.
 * Verifies that `computeNormalize()` does not alter the RealSpace element, as no normalization is needed for RealSpace.
 */
TYPED_TEST(RealSpaceTest, Normalize)
{
    using RealSpace = typename TestFixture::RealSpace;
    using Vector = typename TestFixture::Vector;

    Vector v = Vector::Random();
    RealSpace g(v);
    g.computeNormalize();
    EXPECT_TRUE(g.data().isApprox(v, this->eps)); // RealSpace normalize does nothing
}

/**
 * @brief Tests the squared distance calculation for RealSpace elements.
 * Verifies that the squared distance is correctly computed as the squared Euclidean norm of the difference between the underlying vectors.
 */
TYPED_TEST(RealSpaceTest, SquaredDistance)
{
    using RealSpace = typename TestFixture::RealSpace;
    using Scalar = typename TestFixture::Scalar;
    using Vector = typename TestFixture::Vector;

    Vector v1 = Vector::Random();
    Vector v2 = Vector::Random();
    RealSpace g1(v1);
    RealSpace g2(v2);

    Scalar expected_sq_dist = (v1 - v2).squaredNorm();
    EXPECT_NEAR(g1.squaredDistance(g2), expected_sq_dist, this->eps);
}

/**
 * @brief Tests the interpolation function for RealSpace elements.
 * Verifies that linear interpolation is correctly performed between two RealSpace elements.
 */
TYPED_TEST(RealSpaceTest, Interpolate)
{
    using RealSpace = typename TestFixture::RealSpace;
    using Scalar = typename TestFixture::Scalar;
    using Vector = typename TestFixture::Vector;

    Vector v1 = Vector::Random();
    Vector v2 = Vector::Random();
    RealSpace g1(v1);
    RealSpace g2(v2);

    RealSpace interpolated = g1.interpolate(g2, 0.5);
    EXPECT_TRUE(interpolated.data().isApprox((v1 + v2) / 2.0, this->eps));

    EXPECT_TRUE(g1.interpolate(g2, 0.0).computeIsApprox(g1, this->eps));
    EXPECT_TRUE(g1.interpolate(g2, 1.0).computeIsApprox(g2, this->eps));
}

} // namespace sofa::component::cosserat::liegroups::testing
