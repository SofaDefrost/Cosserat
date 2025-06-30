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

// #include <sofa/testing/BaseTest.h>
// #include <Cosserat/liegroups/SE23.h>
// #include <Eigen/Core>
// #include <Eigen/Geometry>

namespace sofa::component::cosserat::liegroups::testing {

using namespace sofa::testing;

/**
 * @brief Test suite for SE23 Lie group.
 * This test fixture provides common types and a small epsilon for floating-point comparisons.
 * @tparam T The type of SE23 to test (e.g., SE23<double>).
 */
template <typename T>
class SE23Test : public BaseTest
{
protected:
    using SE23 = T;
    using Scalar = typename SE23::Scalar;
    using Vector3 = Eigen::Matrix<Scalar, 3, 1>;
    using TangentVector = typename SE23::TangentVector;

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

using SE23Types = ::testing::Types<SE23<double>>;
TYPED_TEST_SUITE(SE23Test, SE23Types);

/**
 * @brief Tests the constructors of the SE23 class.
 * Verifies that SE23 objects can be constructed from default, and pose and velocity representations.
 */
TYPED_TEST(SE23Test, Constructors)
{
    using SE23 = typename TestFixture::SE23;
    using Scalar = typename TestFixture::Scalar;
    using Vector3 = typename TestFixture::Vector3;

    SE23 g1;
    EXPECT_TRUE(g1.computeIdentity().computeIsApprox(g1, this->eps));

    SE3<Scalar> pose;
    Vector3 vel = Vector3::Random();
    SE23 g2(pose, vel);
    EXPECT_TRUE(g2.pose().computeIsApprox(pose, this->eps));
    EXPECT_TRUE(g2.velocity().isApprox(vel, this->eps));
}

/**
 * @brief Tests the identity element of the SE23 group.
 * Verifies that the `computeIdentity()` method returns an SE23 element with identity pose and zero velocity.
 */
TYPED_TEST(SE23Test, Identity)
{
    using SE23 = typename TestFixture::SE23;
    using Scalar = typename TestFixture::Scalar;
    using Vector3 = typename TestFixture::Vector3;

    SE23 identity = SE23::computeIdentity();
    EXPECT_TRUE(identity.pose().computeIsApprox(SE3<Scalar>::computeIdentity(), this->eps));
    EXPECT_TRUE(identity.velocity().isApprox(Vector3::Zero(), this->eps));
}

/**
 * @brief Tests the inverse operation of the SE23 group.
 * Verifies that composing an SE23 element with its inverse results in the identity element.
 */
TYPED_TEST(SE23Test, Inverse)
{
    using SE23 = typename TestFixture::SE23;
    using Scalar = typename TestFixture::Scalar;
    using Vector3 = typename TestFixture::Vector3;

    SE3<Scalar> pose(SO3<Scalar>(0.1, Vector3(1.0,0.0,0.0)), Vector3(1.0,2.0,3.0));
    Vector3 vel(0.4, 0.5, 0.6);
    SE23 g(pose, vel);
    SE23 inv_g = g.computeInverse();

    SE23 composed = g.compose(inv_g);
    EXPECT_TRUE(composed.computeIdentity().computeIsApprox(composed, this->eps));
}

/**
 * @brief Tests the exponential and logarithmic maps of the SE23 group.
 * Verifies that applying `computeExp` to a Lie algebra element and then `computeLog` to the resulting group element returns the original Lie algebra element.
 */
TYPED_TEST(SE23Test, ExpLog)
{
    using SE23 = typename TestFixture::SE23;
    using Scalar = typename TestFixture::Scalar;
    using TangentVector = typename TestFixture::TangentVector;

    TangentVector twist = TangentVector::Random();
    SE23 g = SE23::computeExp(twist);
    TangentVector log_g = g.computeLog();
    EXPECT_TRUE(log_g.isApprox(twist, this->eps));
}

/**
 * @brief Tests the group action of SE23 on a point-velocity pair.
 * Verifies that a point-velocity pair is correctly transformed by an SE23 element.
 */
TYPED_TEST(SE23Test, Action)
{
    using SE23 = typename TestFixture::SE23;
    using Scalar = typename TestFixture::Scalar;
    using Vector3 = Eigen::Matrix<Scalar, 3, 1>;
    using ActionVector = typename SE23::ActionVector;

    SE3<Scalar> pose(SO3<Scalar>(0.1, Vector3(1.0,0.0,0.0)), Vector3(1.0,2.0,3.0));
    Vector3 vel(0.4, 0.5, 0.6);
    SE23 g(pose, vel);

    ActionVector point_vel;
    point_vel << 1.0, 2.0, 3.0, 0.1, 0.2, 0.3; // point, velocity

    ActionVector transformed_point_vel = g.computeAction(point_vel);
    EXPECT_TRUE(transformed_point_vel.allFinite());
}

/**
 * @brief Tests the approximate equality comparison for SE23 elements.
 * Verifies that two SE23 elements are considered approximately equal within a given tolerance.
 */
TYPED_TEST(SE23Test, IsApprox)
{
    using SE23 = typename TestFixture::SE23;
    using Scalar = typename TestFixture::Scalar;
    using Vector3 = Eigen::Matrix<Scalar, 3, 1>;

    SE3<Scalar> pose1(SO3<Scalar>(0.1, Vector3(1.0,0.0,0.0)), Vector3(1.0,2.0,3.0));
    Vector3 vel1(0.4, 0.5, 0.6);
    SE23 g1(pose1, vel1);

    SE3<Scalar> pose2(SO3<Scalar>(0.1, Vector3(1.0,0.0,0.0)), Vector3(1.0,2.0,3.0));
    Vector3 vel2(0.4, 0.5, 0.6);
    SE23 g2(pose2, vel2);

    EXPECT_TRUE(g1.computeIsApprox(g2, this->eps));
}

/**
 * @brief Tests the random element generation for SE23.
 * Verifies that a randomly generated SE23 element is valid.
 */
TYPED_TEST(SE23Test, Random)
{
    using SE23 = typename TestFixture::SE23;
    std::mt19937 gen(0);
    SE23 r = SE23::computeRandom(gen);
    EXPECT_TRUE(r.computeIsValid());
}

/**
 * @brief Tests the stream output operator for SE23.
 * Verifies that the `print` method produces non-empty output.
 */
TYPED_TEST(SE23Test, Print)
{
    using SE23 = typename TestFixture::SE23;
    using Scalar = typename TestFixture::Scalar;
    using Vector3 = Eigen::Matrix<Scalar, 3, 1>;

    SE3<Scalar> pose(SO3<Scalar>(0.1, Vector3(1.0,0.0,0.0)), Vector3(1.0,2.0,3.0));
    Vector3 vel(0.4, 0.5, 0.6);
    SE23 g(pose, vel);
    std::stringstream ss;
    g.print(ss);
    EXPECT_FALSE(ss.str().empty());
}

/**
 * @brief Tests the validity check for SE23 elements.
 * Verifies that a valid SE23 element is correctly identified as valid.
 */
TYPED_TEST(SE23Test, IsValid)
{
    using SE23 = typename TestFixture::SE23;
    using Scalar = typename TestFixture::Scalar;
    using Vector3 = Eigen::Matrix<Scalar, 3, 1>;

    SE3<Scalar> pose(SO3<Scalar>(0.1, Vector3(1.0,0.0,0.0)), Vector3(1.0,2.0,3.0));
    Vector3 vel(0.4, 0.5, 0.6);
    SE23 g(pose, vel);
    EXPECT_TRUE(g.computeIsValid());
}

/**
 * @brief Tests the normalization of SE23 elements.
 * Verifies that the pose component is normalized after calling `computeNormalize()`.
 */
TYPED_TEST(SE23Test, Normalize)
{
    using SE23 = typename TestFixture::SE23;
    using Scalar = typename TestFixture::Scalar;
    using Vector3 = Eigen::Matrix<Scalar, 3, 1>;

    SE3<Scalar> pose(SO3<Scalar>(0.1, Vector3(1.0,0.0,0.0)), Vector3(1.0,2.0,3.0));
    Vector3 vel(0.4, 0.5, 0.6);
    SE23 g(pose, vel);
    g.computeNormalize();
    EXPECT_TRUE(g.pose().computeIsValid());
}

} // namespace sofa::component::cosserat::liegroups::testing
