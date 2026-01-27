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

#include <sofa/testing/BaseTest.h>
#include <Cosserat/liegroups/SO2.h>
#include <Eigen/Core>
#include <cmath>

namespace sofa::component::cosserat::liegroups::testing {

using namespace sofa::testing;

/**
 * Test suite for SO2 Lie group implementation
 */
class SO2Test : public BaseTest
{
protected:
    using SO2d = SO2<double>;
    using Vector2 = Eigen::Vector2d;
    using Matrix2 = Eigen::Matrix2d;

    const double pi = M_PI;
    const double eps = 1e-10;

    void SetUp() override {}
    void TearDown() override {}
};

/**
 * Test identity element properties
 */
TEST_F(SO2Test, Identity)
{
    SO2d id = SO2d::identity();
    EXPECT_NEAR(id.angle(), 0.0, eps);
    
    // Test right and left identity
    SO2d rot(pi/4);  // 45-degree rotation
    EXPECT_TRUE((rot * id).isApprox(rot));
    EXPECT_TRUE((id * rot).isApprox(rot));
    
    // Test identity matrix
    EXPECT_TRUE(id.matrix().isApprox(Matrix2::Identity()));
}

/**
 * Test group operation (rotation composition)
 */
TEST_F(SO2Test, GroupOperation)
{
    SO2d a(pi/4);    // 45 degrees
    SO2d b(pi/3);    // 60 degrees
    SO2d c = a * b;  // 105 degrees
    
    EXPECT_NEAR(c.angle(), pi/4 + pi/3, eps);
    
    // Test that composition matches matrix multiplication
    Matrix2 Ra = a.matrix();
    Matrix2 Rb = b.matrix();
    Matrix2 Rc = c.matrix();
    EXPECT_TRUE((Ra * Rb).isApprox(Rc));
}

/**
 * Test inverse operation
 */
TEST_F(SO2Test, Inverse)
{
    SO2d rot(pi/3);  // 60-degree rotation
    SO2d inv = rot.inverse();
    
    // Test that inverse rotation has opposite angle
    EXPECT_NEAR(inv.angle(), -pi/3, eps);
    
    // Test that rot * inv = inv * rot = identity
    EXPECT_TRUE((rot * inv).isApprox(SO2d::identity()));
    EXPECT_TRUE((inv * rot).isApprox(SO2d::identity()));
    
    // Test that inverse matches matrix inverse
    EXPECT_TRUE(inv.matrix().isApprox(rot.matrix().inverse()));
}

/**
 * Test exponential and logarithm maps
 */
TEST_F(SO2Test, ExpLog)
{
    // Test exp(log(g)) = g
    SO2d rot(pi/6);  // 30-degree rotation
    auto angle = rot.log();
    auto rot2 = SO2d().exp(angle);
    EXPECT_TRUE(rot.isApprox(rot2));
    
    // Test log(exp(w)) = w
    double w = pi/4;  // Angular velocity
    auto rot3 = SO2d().exp(Vector2::Constant(w));
    EXPECT_NEAR(rot3.log()[0], w, eps);
}

/**
 * Test adjoint representation
 */
TEST_F(SO2Test, Adjoint)
{
    SO2d rot(pi/4);  // 45-degree rotation
    
    // For SO(2), adjoint should always be identity
    EXPECT_TRUE(rot.adjoint().isApprox(Matrix2::Identity()));
}

/**
 * Test group action on points
 */
TEST_F(SO2Test, Action)
{
    SO2d rot(pi/2);  // 90-degree rotation
    Vector2 p(1.0, 0.0);  // Point on x-axis
    
    // 90-degree rotation should map (1,0) to (0,1)
    Vector2 q = rot.act(p);
    EXPECT_NEAR(q[0], 0.0, eps);
    EXPECT_NEAR(q[1], 1.0, eps);
    
    // Test that action matches matrix multiplication
    EXPECT_TRUE(q.isApprox(rot.matrix() * p));
}

/**
 * Test angle normalization
 */
TEST_F(SO2Test, AngleNormalization)
{
    // Test that angles are normalized to [-π, π]
    SO2d rot1(3*pi/2);  // 270 degrees
    SO2d rot2(-3*pi/2); // -270 degrees
    
    EXPECT_NEAR(rot1.angle(), -pi/2, eps);
    EXPECT_NEAR(rot2.angle(), pi/2, eps);
}

/**
 * Test interpolation
 */
TEST_F(SO2Test, Interpolation)
{
    SO2d start(0);          // 0 degrees
    SO2d end(pi/2);        // 90 degrees
    SO2d mid = slerp(start, end, 0.5);
    
    // Midpoint should be 45-degree rotation
    EXPECT_NEAR(mid.angle(), pi/4, eps);
}

} // namespace sofa::component::cosserat::liegroups::testing
