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
 * Test suite for RealSpace Lie group
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

    void SetUp() override {}
    void TearDown() override {}
};

using RealSpaceTypes = ::testing::Types<RealSpace<double, 1>, RealSpace<double, 2>, RealSpace<double, 3>>;
TYPED_TEST_SUITE(RealSpaceTest, RealSpaceTypes);

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

TYPED_TEST(RealSpaceTest, Identity)
{
    using RealSpace = typename TestFixture::RealSpace;
    using Scalar = typename TestFixture::Scalar;
    using Vector = typename TestFixture::Vector;

    RealSpace identity = RealSpace::computeIdentity();
    EXPECT_TRUE(identity.data().isApprox(Vector::Zero(), this->eps));
}

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

TYPED_TEST(RealSpaceTest, ExpLog)
{
    using RealSpace = typename TestFixture::RealSpace;
    using Scalar = typename TestFixture::Scalar;
    using TangentVector = typename TestFixture::TangentVector;

    TangentVector v = TangentVector::Random();
    RealSpace g = RealSpace::computeExp(v);
    EXPECT_TRUE(g.computeLog().isApprox(v, this->eps));
}

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

TYPED_TEST(RealSpaceTest, Adjoint)
{
    using RealSpace = typename TestFixture::RealSpace;
    using Scalar = typename TestFixture::Scalar;
    using TangentVector = typename RealSpace::TangentVector;
    using AdjointMatrix = typename RealSpace::AdjointMatrix;

    TangentVector v = TangentVector::Random();
    AdjointMatrix Ad = RealSpace::computeAd(v);
    EXPECT_TRUE(Ad.isApprox(AdjointMatrix::Zero(), this->eps));
}

TYPED_TEST(RealSpaceTest, Random)
{
    using RealSpace = typename TestFixture::RealSpace;
    std::mt19937 gen(0);
    RealSpace r = RealSpace::computeRandom(gen);
    EXPECT_TRUE(r.computeIsValid());
}

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

TYPED_TEST(RealSpaceTest, Normalize)
{
    using RealSpace = typename TestFixture::RealSpace;
    using Vector = typename TestFixture::Vector;

    Vector v = Vector::Random();
    RealSpace g(v);
    g.computeNormalize();
    EXPECT_TRUE(g.data().isApprox(v, this->eps)); // RealSpace normalize does nothing
}

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