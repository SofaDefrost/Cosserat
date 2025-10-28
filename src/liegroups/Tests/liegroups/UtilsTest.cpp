/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                  *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
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
* along with this program. If not, see <http://www.gnu.org/licenses/\>.        *
******************************************************************************/

#include <cmath>
#include <sofa/testing/BaseTest.h>
#include <Cosserat/liegroups/Utils.h>
#include <Eigen/Core>
#include <cmath>

namespace sofa::component::cosserat::liegroups::testing {

using namespace sofa::testing;

/**
 * @brief Test suite for Lie group utilities.
 * Inherits from BaseTest to leverage SOFA's testing framework.
 */
class UtilsTest : public BaseTest
{
protected:
    using Utils = Types<double>;
    using Vector2d = Eigen::Vector2d;
    using Vector3d = Eigen::Vector3d;

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
 * @brief Tests the angle normalization function.
 * Verifies that angles are correctly normalized to the range [-π, π].
 */
TEST_F(UtilsTest, AngleNormalization)
{
    // Test normalization of angles within [-π, π]
    EXPECT_NEAR(Utils::normalizeAngle(0.0), 0.0, eps);
    EXPECT_NEAR(Utils::normalizeAngle(pi/2), pi/2, eps);
    EXPECT_NEAR(Utils::normalizeAngle(-pi/2), -pi/2, eps);
    EXPECT_NEAR(Utils::normalizeAngle(pi), pi, eps);
    EXPECT_NEAR(Utils::normalizeAngle(-pi), -pi, eps);
    
    // Test normalization of angles outside [-π, π]
    EXPECT_NEAR(Utils::normalizeAngle(3*pi/2), -pi/2, eps);
    EXPECT_NEAR(Utils::normalizeAngle(-3*pi/2), pi/2, eps);
    EXPECT_NEAR(Utils::normalizeAngle(2*pi), 0.0, eps);
    EXPECT_NEAR(Utils::normalizeAngle(4*pi), 0.0, eps);
    EXPECT_NEAR(Utils::normalizeAngle(-2*pi), 0.0, eps);
    
    // Test extreme cases
    EXPECT_NEAR(Utils::normalizeAngle(1000*pi), 0.0, eps);
    EXPECT_NEAR(Utils::normalizeAngle(1001*pi), pi, eps);
}

/**
 * @brief Tests the sinc function for numerical stability.
 * Verifies correct behavior for non-zero, small, and negative values.
 */
TEST_F(UtilsTest, Sinc)
{
    // Test non-zero values
    EXPECT_NEAR(Utils::sinc(pi/2), 2/pi, eps);
    EXPECT_NEAR(Utils::sinc(pi), 0.0, eps);
    EXPECT_NEAR(Utils::sinc(2*pi), 0.0, eps);
    
    // Test small values (near zero)
    EXPECT_NEAR(Utils::sinc(1e-10), 1.0, eps);
    EXPECT_NEAR(Utils::sinc(1e-8), 1.0, eps);
    EXPECT_NEAR(Utils::sinc(0.0), 1.0, eps);
    
    // Test negative values
    EXPECT_NEAR(Utils::sinc(-pi/2), -2/pi, eps);
}

/**
 * @brief Tests the oneMinusCos function for numerical stability.
 * Verifies correct behavior for non-zero, small, and negative values.
 */
TEST_F(UtilsTest, OneMinusCos)
{
    // Test non-zero values
    EXPECT_NEAR(Utils::oneMinusCos(pi/2), 1.0, eps);
    EXPECT_NEAR(Utils::oneMinusCos(pi), 2.0, eps);
    EXPECT_NEAR(Utils::oneMinusCos(2*pi), 0.0, eps);
    
    // Test small values (near zero)
    EXPECT_NEAR(Utils::oneMinusCos(1e-8), 5e-17, 1e-16);
    EXPECT_NEAR(Utils::oneMinusCos(1e-4), 5e-9, 1e-8);
    EXPECT_NEAR(Utils::oneMinusCos(0.0), 0.0, eps);
    
    // Test negative values
    EXPECT_NEAR(Utils::oneMinusCos(-pi/2), 1.0, eps);
    EXPECT_NEAR(Utils::oneMinusCos(-pi), 2.0, eps);
}

/**
 * @brief Tests the angle difference calculation.
 * Verifies correct differences, including cases with angle wrapping.
 */
TEST_F(UtilsTest, AngleDifference)
{
    // Test differences within [-π, π]
    EXPECT_NEAR(Utils::angleDifference(0.0, 0.0), 0.0, eps);
    EXPECT_NEAR(Utils::angleDifference(pi/4, 0.0), pi/4, eps);
    EXPECT_NEAR(Utils::angleDifference(0.0, pi/4), -pi/4, eps);
    
    // Test differences that wrap around
    EXPECT_NEAR(Utils::angleDifference(3*pi/4, -3*pi/4), 3*pi/2, eps);
    EXPECT_NEAR(Utils::angleDifference(-3*pi/4, 3*pi/4), -3*pi/2, eps);
    
    // Test extreme cases
    EXPECT_NEAR(Utils::angleDifference(pi, -pi), 0.0, eps);
}

/**
 * @brief Tests the angle distance calculation.
 * Verifies correct distances, including cases with angle wrapping.
 */
TEST_F(UtilsTest, AngleDistance)
{
    // Test distances within [-π, π]
    EXPECT_NEAR(Utils::angleDistance(0.0, 0.0), 0.0, eps);
    EXPECT_NEAR(Utils::angleDistance(pi/4, 0.0), pi/4, eps);
    EXPECT_NEAR(Utils::angleDistance(0.0, pi/4), pi/4, eps);
    
    // Test distances that wrap around
    EXPECT_NEAR(Utils::angleDistance(3*pi/4, -3*pi/4), 3*pi/2, eps);
    EXPECT_NEAR(Utils::angleDistance(-3*pi/4, 3*pi/4), 3*pi/2, eps);
    
    // Test extreme cases
    EXPECT_NEAR(Utils::angleDistance(pi, -pi), 0.0, eps);
}

/**
 * @brief Tests linear interpolation.
 * Verifies correct interpolation and extrapolation behavior.
 */
TEST_F(UtilsTest, LinearInterpolation)
{
    // Test standard cases
    EXPECT_NEAR(Utils::lerp(0.0, 1.0, 0.0), 0.0, eps);
    EXPECT_NEAR(Utils::lerp(0.0, 1.0, 1.0), 1.0, eps);
    EXPECT_NEAR(Utils::lerp(0.0, 1.0, 0.5), 0.5, eps);
    EXPECT_NEAR(Utils::lerp(-1.0, 1.0, 0.5), 0.0, eps);
    
    // Test extrapolation
    EXPECT_NEAR(Utils::lerp(0.0, 1.0, 2.0), 2.0, eps);
    EXPECT_NEAR(Utils::lerp(0.0, 1.0, -1.0), -1.0, eps);
}

/**
 * @brief Tests spherical linear interpolation (SLERP) for angles.
 * Verifies correct interpolation, including cases with angle wrapping.
 */
TEST_F(UtilsTest, SlerpAngle)
{
    // Test standard cases
    EXPECT_NEAR(Utils::slerpAngle(0.0, pi/2, 0.0), 0.0, eps);
    EXPECT_NEAR(Utils::slerpAngle(0.0, pi/2, 1.0), pi/2, eps);
    EXPECT_NEAR(Utils::slerpAngle(0.0, pi/2, 0.5), pi/4, eps);
    
    // Test wrapping
    EXPECT_NEAR(Utils::slerpAngle(-3*pi/4, 3*pi/4, 0.5), 0.0, eps);
    EXPECT_NEAR(Utils::slerpAngle(3*pi/4, -3*pi/4, 0.5), 0.0, eps);
    
    // Test extreme cases
    EXPECT_NEAR(Utils::slerpAngle(pi, -pi, 0.5), pi, eps);
}

/**
 * @brief Tests near-zero detection for angles.
 * Verifies that angles very close to zero are correctly identified.
 */
TEST_F(UtilsTest, NearZeroAngle)
{
    EXPECT_TRUE(Utils::isAngleNearZero(0.0));
    EXPECT_TRUE(Utils::isAngleNearZero(1e-12));
    EXPECT_TRUE(Utils::isAngleNearZero(-1e-12));
    
    EXPECT_FALSE(Utils::isAngleNearZero(0.1));
    EXPECT_FALSE(Utils::isAngleNearZero(-0.1));
}

/**
 * @brief Tests nearly equal detection for angles.
 * Verifies that angles that are approximately equal (considering wrapping) are correctly identified.
 */
TEST_F(UtilsTest, NearlyEqualAngles)
{
    EXPECT_TRUE(Utils::areAnglesNearlyEqual(0.0, 0.0));
    EXPECT_TRUE(Utils::areAnglesNearlyEqual(pi/4, pi/4 + 1e-12));
    EXPECT_TRUE(Utils::areAnglesNearlyEqual(pi, -pi));
    EXPECT_TRUE(Utils::areAnglesNearlyEqual(2*pi, 0.0));
    
    EXPECT_FALSE(Utils::areAnglesNearlyEqual(0.0, 0.1));
    EXPECT_FALSE(Utils::areAnglesNearlyEqual(pi/4, pi/2));
}

/**
 * @brief Tests safe vector normalization.
 * Verifies that non-zero, zero, and near-zero vectors are normalized correctly.
 */
TEST_F(UtilsTest, SafeNormalize)
{
    // Test non-zero vectors
    Vector2d v1(3.0, 4.0);
    Vector2d v1_norm = Utils::safeNormalize(v1);
    EXPECT_NEAR(v1_norm.norm(), 1.0, eps);
    EXPECT_NEAR(v1_norm[0], 0.6, eps);
    EXPECT_NEAR(v1_norm[1], 0.8, eps);
    
    // Test zero vector
    Vector2d v2(0.0, 0.0);
    Vector2d v2_norm = Utils::safeNormalize(v2);
    EXPECT_NEAR(v2_norm[0], 0.0, eps);
    EXPECT_NEAR(v2_norm[1], 0.0, eps);
    
    // Test near-zero vector
    Vector2d v3(1e-12, 1e-12);
    Vector2d v3_norm = Utils::safeNormalize(v3);
    EXPECT_NEAR(v3_norm[0], 0.0, eps);
    EXPECT_NEAR(v3_norm[1], 0.0, eps);
}

/**
 * @brief Tests vector projection.
 * Verifies correct projection onto various vectors, including zero vectors.
 */
TEST_F(UtilsTest, VectorProjection)
{
    // Standard projection
    Vector2d v(3.0, 3.0);
    Vector2d onto(1.0, 0.0);
    Vector2d proj = Utils::projectVector(v, onto);
    EXPECT_NEAR(proj[0], 3.0, eps);
    EXPECT_NEAR(proj[1], 0.0, eps);
    
    // Project onto zero vector
    Vector2d proj_zero = Utils::projectVector(v, Vector2d(0.0, 0.0));
    EXPECT_NEAR(proj_zero[0], 0.0, eps);
    EXPECT_NEAR(proj_zero[1], 0.0, eps);
    
    // Project onto self
    Vector2d proj_self = Utils::projectVector(v, v);
    EXPECT_NEAR(proj_self[0], v[0], eps);
    EXPECT_NEAR(proj_self[1], v[1], eps);
}



} // namespace sofa::component::cosserat::liegroups::testing