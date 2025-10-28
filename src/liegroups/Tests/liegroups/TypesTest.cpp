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
#include <Cosserat/liegroups/Types.h>
#include <Eigen/Core>
#include <cmath>


namespace sofa::component::cosserat::liegroups::testing {

using namespace sofa::testing;

/**
 * @brief Test suite for Types utility functions.
 * This test fixture provides common types and a small epsilon for floating-point comparisons.
 */
class TypesTest : public BaseTest
{
protected:
    using Typesd = Types<double>;
    using Vector2d = Eigen::Vector2d;
    using Vector3d = Eigen::Vector3d;
    using Matrix3d = Eigen::Matrix3d;

    const double pi = M_PI;
    const double eps = 1e-10;

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

/**
 * @brief Tests the epsilon and tolerance values defined in Types.
 * Verifies that epsilon and tolerance are positive and that tolerance is greater than epsilon.
 */
TEST_F(TypesTest, EpsilonAndTolerance)
{
    EXPECT_GT(Typesd::epsilon(), 0.0);
    EXPECT_GT(Typesd::tolerance(), 0.0);
    EXPECT_GT(Typesd::tolerance(), Typesd::epsilon());
}

/**
 * @brief Tests the `isZero` function in Types.
 * Verifies that values very close to zero are correctly identified as zero.
 */
TEST_F(TypesTest, IsZero)
{
    EXPECT_TRUE(Typesd::isZero(0.0));
    EXPECT_TRUE(Typesd::isZero(1e-12));
    EXPECT_TRUE(Typesd::isZero(-1e-12));
    EXPECT_FALSE(Typesd::isZero(0.1));
}

/**
 * @brief Tests the `isApprox` function in Types.
 * Verifies that two values are considered approximately equal within a given tolerance.
 */
TEST_F(TypesTest, IsApprox)
{
    EXPECT_TRUE(Typesd::isApprox(0.0, 0.0));
    EXPECT_TRUE(Typesd::isApprox(1.0, 1.0 + 1e-12));
    EXPECT_FALSE(Typesd::isApprox(0.0, 0.1));
}

/**
 * @brief Tests the `sinc` function in Types.
 * Verifies the correct calculation of sinc(x) for various inputs, including near zero.
 */
TEST_F(TypesTest, Sinc)
{
    EXPECT_NEAR(Typesd::sinc(0.0), 1.0, eps);
    EXPECT_NEAR(Typesd::sinc(pi/2), 2.0/pi, eps);
    EXPECT_NEAR(Typesd::sinc(pi), 0.0, eps);
}

/**
 * @brief Tests the `cosc` function in Types.
 * Verifies the correct calculation of cosc(x) for various inputs, including near zero.
 */
TEST_F(TypesTest, Cosc)
{
    EXPECT_NEAR(Typesd::cosc(0.0), 1.0, eps);
    EXPECT_NEAR(Typesd::cosc(pi/2), 0.0, eps);
    EXPECT_NEAR(Typesd::cosc(pi), -1.0/pi, eps);
}

/**
 * @brief Tests the `sinc2` function in Types.
 * Verifies the correct calculation of sinc2(x) for various inputs, including near zero.
 */
TEST_F(TypesTest, Sinc2)
{
    EXPECT_NEAR(Typesd::sinc2(0.0), 0.5, eps);
    EXPECT_NEAR(Typesd::sinc2(pi), 2.0/(pi*pi), eps);
}

/**
 * @brief Tests the `normalizeAngle` function in Types.
 * Verifies that angles are correctly normalized to the range [-π, π].
 */
TEST_F(TypesTest, NormalizeAngle)
{
    EXPECT_NEAR(Typesd::normalizeAngle(0.0), 0.0, eps);
    EXPECT_NEAR(Typesd::normalizeAngle(pi), pi, eps);
    EXPECT_NEAR(Typesd::normalizeAngle(2*pi), 0.0, eps);
    EXPECT_NEAR(Typesd::normalizeAngle(-pi/2), -pi/2, eps);
    EXPECT_NEAR(Typesd::normalizeAngle(3*pi/2), -pi/2, eps);
}

/**
 * @brief Tests the `skew3` and `unskew3` functions in Types.
 * Verifies that a 3D vector can be converted to a skew-symmetric matrix and back, and that the matrix is indeed skew-symmetric.
 */
TEST_F(TypesTest, SkewSymmetricMatrix)
{
    Vector3d v(1.0, 2.0, 3.0);
    Matrix3d skew_mat = Typesd::skew3(v);
    EXPECT_TRUE(Typesd::isSkewSymmetric(skew_mat));
    EXPECT_TRUE(Typesd::unskew3(skew_mat).isApprox(v, eps));
}

/**
 * @brief Tests the random number generation functions in Types.
 * Verifies that random scalars are within the expected range and random unit vectors have a norm of 1.
 */
TEST_F(TypesTest, RandomFunctions)
{
    std::mt19937 gen(0);
    Scalar random_scalar = Typesd::randomScalar(gen);
    EXPECT_GE(random_scalar, 0.0);
    EXPECT_LE(random_scalar, 1.0);

    Vector3d random_vec = Typesd::randomVector<3>(gen);
    EXPECT_TRUE(random_vec.allFinite());

    Vector3d random_unit_vec = Typesd::randomUnitVector<3>(gen);
    EXPECT_NEAR(random_unit_vec.norm(), 1.0, eps);
}

} // namespace sofa::component::cosserat::liegroups::testing