/******************************************************************************
 * Cosserat Plugin — test for the resolutions of CosseratSlidingConstraint    *
 * (separate translation unit: CosseratNeedleSlidingConstraint.inl brings a   *
 * BilateralConstraintResolution of SOFA into the same namespace)             *
 ******************************************************************************/

#include <Cosserat/constraint/CosseratSlidingConstraint.h>

#include <gtest/gtest.h>

#include <cmath>

// The resolutions of CosseratSlidingConstraint divided by the diagonal of W without checking it: a clamped point (no compliance)
// gave 0/0 = NaN and blew the simulation up.
TEST(CosseratSlidingConstraintResolutionTest, ZeroDiagonalGivesNoNaN)
{
    double w00 = 0.0;
    double* rows[1] = {&w00};
    double d[1] = {0.5};

    double force[1] = {0.0};
    sofa::component::constraintset::BilateralConstraintResolution bilateral;
    bilateral.resolution(0, rows, d, force, nullptr);
    EXPECT_TRUE(std::isfinite(force[0]));

    force[0] = 0.0;
    sofa::component::constraintset::UnilateralConstraintResolution unilateral;
    unilateral.resolution(0, rows, d, force, nullptr);
    EXPECT_TRUE(std::isfinite(force[0]));

    // and the usual case is unchanged
    w00 = 2.0;
    force[0] = 0.0;
    bilateral.resolution(0, rows, d, force, nullptr);
    EXPECT_DOUBLE_EQ(force[0], -0.25);
}
