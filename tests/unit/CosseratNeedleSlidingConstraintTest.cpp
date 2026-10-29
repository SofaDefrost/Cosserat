/******************************************************************************
 * Cosserat Plugin — test for CosseratNeedleSlidingConstraint                 *
 *                                                                            *
 * buildConstraintMatrix, getConstraintViolation and getConstraintResolution  *
 * must describe the same lines, in the same order, for every value of        *
 * `useDirections`: one line per constrained direction and per point, in the  *
 * order y, z (x is the sliding axis and is never constrained, whatever       *
 * `useDirections[0]`).  (Before the fix the violation was written at 2*i /   *
 * 2*i+1 whatever the directions: wrong when a single direction is used.)     *
 ******************************************************************************/

#include <Cosserat/constraint/CosseratNeedleSlidingConstraint.inl>

#include <gtest/gtest.h>

#include <iterator>

#include <sofa/component/statecontainer/MechanicalObject.h>
#include <sofa/core/ConstraintParams.h>
#include <sofa/core/ExecParams.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/linearalgebra/FullVector.h>
#include <sofa/simulation/Node.h>
#include <sofa/simulation/graph/DAGSimulation.h>

using sofa::component::constraintset::CosseratNeedleSlidingConstraint;
using sofa::component::statecontainer::MechanicalObject;
using sofa::defaulttype::Vec3Types;

class CosseratNeedleSlidingConstraintTest : public ::testing::Test
{
public:
    using Constraint = CosseratNeedleSlidingConstraint<Vec3Types>;
    sofa::simulation::Node::SPtr root;
    MechanicalObject<Vec3Types>::SPtr state;
    Constraint::SPtr constraint;

    void SetUp() override
    {
        root = sofa::simulation::getSimulation()->createNewNode("root");
        state = sofa::core::objectmodel::New<MechanicalObject<Vec3Types>>();
        root->addObject(state);
        state->resize(3);
        {
            auto x = sofa::helper::getWriteAccessor(*state->write(sofa::core::vec_id::write_access::position));
            x[0] = sofa::type::Vec3(0., 1., 2.);
            x[1] = sofa::type::Vec3(3., 4., 5.);
            x[2] = sofa::type::Vec3(6., 7., 8.);
        }
        constraint = sofa::core::objectmodel::New<Constraint>(state.get());
        root->addObject(constraint);
    }

    void TearDown() override
    {
        if (root) sofa::simulation::node::unload(root);
    }

    // Runs the three methods for the given directions and checks them against each other.
    void check(bool x, bool y, bool z)
    {
        constraint->findData("useDirections")->read(std::string(x ? "1 " : "0 ") + (y ? "1 " : "0 ") + (z ? "1" : "0"));
        constraint->init();

        sofa::core::ConstraintParams cparams(*sofa::core::execparams::defaultInstance());
        Vec3Types::MatrixDeriv matrixValue;
        sofa::Data<Vec3Types::MatrixDeriv> matrix(matrixValue, "matrix");
        unsigned int index = 7;  // the constraint does not start at line 0
        const auto& positions = state->read(sofa::core::vec_id::read_access::position);
        constraint->buildConstraintMatrix(&cparams, matrix, index, *positions);

        const unsigned int perPoint = (y ? 1u : 0u) + (z ? 1u : 0u);   // x is never constrained
        ASSERT_EQ(index, 7u + 3u * perPoint);

        // the violation is written on the lines of the matrix, in the same order
        sofa::linearalgebra::FullVector<double> violation(7 + 3 * perPoint);
        for (unsigned int i = 0; i < violation.size(); ++i) violation.set(i, -1000.);
        sofa::Data<Vec3Types::VecDeriv> velocities;
        constraint->getConstraintViolation(&cparams, &violation, *positions, velocities);

        const Vec3Types::MatrixDeriv& jacobian = matrix.getValue();
        std::vector<double> expected;
        for (unsigned int p = 0; p < 3; ++p)
            for (unsigned int c = 0; c < 3; ++c)
                if ((c == 1 && y) || (c == 2 && z)) expected.push_back(3. * p + c);
        ASSERT_EQ(expected.size(), 3u * perPoint);

        for (unsigned int line = 0; line < expected.size(); ++line)
        {
            EXPECT_DOUBLE_EQ(violation.element(7 + line), expected[line]) << "line " << line;
            // the line of the matrix has exactly one entry: 1 on the same coordinate of the same point
            auto row = jacobian.readLine(7 + line);
            ASSERT_EQ(std::distance(row.begin(), row.end()), 1) << "line " << line;
            const auto it = row.begin();
            EXPECT_EQ(static_cast<unsigned int>(it.index()), line / perPoint) << "line " << line;
            const double coordinate = expected[line] - 3. * (line / perPoint);
            EXPECT_DOUBLE_EQ(it.val()[static_cast<unsigned int>(coordinate)], 1.0) << "line " << line;
        }
        for (unsigned int line = 0; line < 7; ++line) EXPECT_DOUBLE_EQ(violation.element(line), -1000.);  // lines before are untouched

        std::vector<sofa::core::behavior::ConstraintResolution*> resolutions(7 + 3 * perPoint, nullptr);
        unsigned int offset = 7;
        constraint->getConstraintResolution(&cparams, resolutions, offset);
        EXPECT_EQ(offset, 7u + 3u * perPoint);
        for (auto* r : resolutions) delete r;
    }
};

TEST_F(CosseratNeedleSlidingConstraintTest, DefaultDirectionsAreYAndZ)
{
    check(false, true, true);
}

TEST_F(CosseratNeedleSlidingConstraintTest, OnlyZ)
{
    check(false, false, true);
}

TEST_F(CosseratNeedleSlidingConstraintTest, OnlyY)
{
    check(false, true, false);
}

TEST_F(CosseratNeedleSlidingConstraintTest, XIsTheSlidingAxisAndIsNeverConstrained)
{
    // x asked: no line for it, whatever the other directions
    check(true, false, false);
    check(true, true, true);
}
