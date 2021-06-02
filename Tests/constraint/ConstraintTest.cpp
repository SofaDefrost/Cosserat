//
// Created by younes on 02/06/2021.
//

#include <gtest/gtest.h>


#include "Constraint.h"
#include "../Example.h"

using testing::Test;

struct ConstraintTest
        : public ::testing::Test
{
    int* x;

    int GetX() const{
        return *x;
    }

    // Sets up the test fixture.
    virtual void SetUp() override {
        fprintf(stderr, "Starting up ! \n");
        x = new int(42);
    }

    // Tears down the test fixture.
    virtual void TearDown() override {
        fprintf(stderr, "Starting down ! \n");
        delete x;
    }
};

TEST_F(ConstraintTest,MAC)
{
    int y = 16;
    int sum = 100;
    int oldSum = sum + GetX() * y;
    int expectedNewSum = sum + GetX() * y;

    EXPECT_EQ(
            expectedNewSum,
            MAC(GetX(),y,sum)
    );

    EXPECT_EQ(
            expectedNewSum,
            oldSum
    );
}


TEST(ConstraintTest2, Square){
    int x = 5;
    int expectedSquare = x*x;
    EXPECT_EQ(
            expectedSquare,
            Square(x)
    );
}