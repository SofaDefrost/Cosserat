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
#include <Cosserat/liegroups/SGal3.h>
#include <Eigen/Core>
#include <Eigen/Geometry>

namespace sofa::component::cosserat::liegroups::testing {

using namespace sofa::testing;

/**
 * @brief Test suite for SGal3 Lie group.
 * This test fixture provides common types and a small epsilon for floating-point comparisons.
 * @tparam T The type of SGal3 to test (e.g., SGal3<double>).
 */
template <typename T>
class SGal3Test : public BaseTest
{
protected:
    using SGal3 = T;
    using Scalar = typename SGal3::Scalar;
    using Vector3 = Eigen::Matrix<Scalar, 3, 1>;
    using TangentVector = typename SGal3::TangentVector;

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

using SGal3Types = ::testing::Types<SGal3<double>>;
TYPED_TEST_SUITE(SGal3Test, SGal3Types);

/**
 * @brief Tests the constructors of the SGal3 class.
 * Verifies that SGal3 objects can be constructed from default, and pose, velocity, and time representations.
 */
TYPED_TEST(SGal3Test, Constructors)
{
    using SGal3 = typename TestFixture::SGal3;
    using Scalar = typename TestFixture::Scalar;
    using Vector3 = typename TestFixture::Vector3;

    SGal3 g1;
    EXPECT_TRUE(g1.computeIdentity().computeIsApprox(g1, this->eps));

    SE3<Scalar> pose;
    Vector3 vel = Vector3::Random();
    Scalar time = 1.0;
    SGal3 g2(pose, vel, time);
    EXPECT_TRUE(g2.pose().computeIsApprox(pose, this->eps));
    EXPECT_TRUE(g2.velocity().isApprox(vel, this->eps));
    EXPECT_NEAR(g2.time(), time, this->eps);
}

/**
 * @brief Tests the identity element of the SGal3 group.
 * Verifies that the `computeIdentity()` method returns an SGal3 element with identity pose, zero velocity, and zero time.
 */
TYPED_TEST(SGal3Test, Identity)
{
    using SGal3 = typename TestFixture::SGal3;
    using Scalar = typename TestFixture::Scalar;
    using Vector3 = typename TestFixture::Vector3;

    SGal3 identity = SGal3::computeIdentity();
    EXPECT_TRUE(identity.pose().computeIsApprox(SE3<Scalar>::computeIdentity(), this->eps));
    EXPECT_TRUE(identity.velocity().isApprox(Vector3::Zero(), this->eps));
    EXPECT_NEAR(identity.time(), 0.0, this->eps);
}

/**
 * @brief Tests the inverse operation of the SGal3 group.
 * Verifies that composing an SGal3 element with its inverse results in the identity element.
 */
TYPED_TEST(SGal3Test, Inverse)
{
    using SGal3 = typename TestFixture::SGal3;
    using Scalar = typename TestFixture::Scalar;
    using Vector3 = typename TestFixture::Vector3;

    SE3<Scalar> pose(SO3<Scalar>(0.1, Vector3(1,0,0)), Vector3(1,2,3));
    Vector3 vel(0.4, 0.5, 0.6);
    Scalar time = 1.0;
    SGal3 g(pose, vel, time);
    SGal3 inv_g = g.computeInverse();

    SGal3 composed = g.compose(inv_g);
    EXPECT_TRUE(composed.computeIdentity().computeIsApprox(composed, this->eps));
}

/**
 * @brief Tests the exponential and logarithmic maps of the SGal3 group.
 * Verifies that applying `computeExp` to a Lie algebra element and then `computeLog` to the resulting group element returns the original Lie algebra element.
 */
TYPED_TEST(SGal3Test, ExpLog)
{
    using SGal3 = typename TestFixture::SGal3;
    using Scalar = typename TestFixture::Scalar;
    using TangentVector = typename TestFixture::TangentVector;

    TangentVector twist = TangentVector::Random();
    SGal3 g = SGal3::computeExp(twist);
    TangentVector log_g = g.computeLog();
    EXPECT_TRUE(log_g.isApprox(twist, this->eps));
}

/**
 * @brief Tests the group action of SGal3 on a point-velocity-time tuple.
 * Verifies that a point-velocity-time tuple is correctly transformed by an SGal3 element.
 */
TYPED_TEST(SGal3Test, Action)
{
    using SGal3 = typename TestFixture::SGal3;
    using Scalar = typename TestFixture::Scalar;
    using Vector3 = Eigen::Matrix<Scalar, 3, 1>;
    using ActionVector = typename SGal3::ActionVector;

    SE3<Scalar> pose(SO3<Scalar>(0.1, Vector3(1,0,0)), Vector3(1,2,3));
    Vector3 vel(0.4, 0.5, 0.6);
    Scalar time = 1.0;
    SGal3 g(pose, vel, time);

    ActionVector point_vel_time;
    point_vel_time << 1.0, 2.0, 3.0, 0.1, 0.2, 0.3, 0.01, 0.02, 0.03, 0.5; // point, velocity, boost, time

    ActionVector transformed_point_vel_time = g.computeAction(point_vel_time);
    EXPECT_TRUE(transformed_point_vel_time.allFinite());
}

/**
 * @brief Tests the approximate equality comparison for SGal3 elements.
 * Verifies that two SGal3 elements are considered approximately equal within a given tolerance.
 */
TYPED_TEST(SGal3Test, IsApprox)
{
    using SGal3 = typename TestFixture::SGal3;
    using Scalar = typename TestFixture::Scalar;
    using Vector3 = Eigen::Matrix<Scalar, 3, 1>;

    SE3<Scalar> pose1(SO3<Scalar>(0.1, Vector3(1,0,0)), Vector3(1,2,3));
    Vector3 vel1(0.4, 0.5, 0.6);
    Scalar time1 = 1.0;
    SGal3 g1(pose1, vel1, time1);

                                SE3<Scalar> pose2(SO3<Scalar>(0.1, Vector3(1.0,0.0,0.0)), Vector3(1.0,2.0,3.0));
    Vector3 vel2(0.4, 0.5, 0.6);
    Scalar time2 = 1.0;
    SGal3 g2(pose2, vel2, time2);

    EXPECT_TRUE(g1.computeIsApprox(g2, this->eps));
}

/**
 * @brief Tests the random element generation for SGal3.
 * Verifies that a randomly generated SGal3 element is valid.
 */
TYPED_TEST(SGal3Test, Random)
{
    using SGal3 = typename TestFixture::SGal3;
    std::mt19937 gen(0);
    SGal3 r = SGal3::computeRandom(gen);
    EXPECT_TRUE(r.computeIsValid());
}

/**
 * @brief Tests the stream output operator for SGal3.
 * Verifies that the `print` method produces non-empty output.
 */
TYPED_TEST(SGal3Test, Print)
{
    using SGal3 = typename TestFixture::SGal3;
    using Scalar = typename TestFixture::Scalar;
    using Vector3 = Eigen::Matrix<Scalar, 3, 1>;

    SE3<Scalar> pose(SO3<Scalar>(0.1, Vector3(1,0,0)), Vector3(1,2,3));
    Vector3 vel(0.4, 0.5, 0.6);
    Scalar time = 1.0;
    SGal3 g(pose, vel, time);
    std::stringstream ss;
    g.print(ss);
    EXPECT_FALSE(ss.str().empty());
}

/**
 * @brief Tests the validity check for SGal3 elements.
 * Verifies that a valid SGal3 element is correctly identified as valid.
 */
TYPED_TEST(SGal3Test, IsValid)
{
    using SGal3 = typename TestFixture::SGal3;
    using Scalar = typename TestFixture::Scalar;
    using Vector3 = Eigen::Matrix<Scalar, 3, 1>;

    SE3<Scalar> pose(SO3<Scalar>(0.1, Vector3(1,0,0)), Vector3(1,2,3));
    Vector3 vel(0.4, 0.5, 0.6);
    Scalar time = 1.0;
    SGal3 g(pose, vel, time);
    EXPECT_TRUE(g.computeIsValid());
}

/**
 * @brief Tests the normalization of SGal3 elements.
 * Verifies that the pose component is normalized after calling `computeNormalize()`.
 */
TYPED_TEST(SGal3Test, Normalize)
{
    using SGal3 = typename TestFixture::SGal3;
    using Scalar = typename TestFixture::Scalar;
    using Vector3 = Eigen::Matrix<Scalar, 3, 1>;

    SE3<Scalar> pose(SO3<Scalar>(0.1, Vector3(1,0,0)), Vector3(1,2,3));
    Vector3 vel(0.4, 0.5, 0.6);
    Scalar time = 1.0;
    SGal3 g(pose, vel, time);
    g.computeNormalize();
    EXPECT_TRUE(g.pose().computeIsValid());
}

} // namespace sofa::component::cosserat::liegroups::testing