//
// Created by younes on 31/05/2021.
//

#include <gtest/gtest.h>
#include "../Example.h"

TEST(ExempleTests,MAC)
{
    int x = 42;
    int y = 16;
    int sum = 100;
    int oldSum = sum;
    int expectedNewSum = oldSum + x * y;

    EXPECT_EQ(
            expectedNewSum,
            MAC(x,y,sum)
    );

    EXPECT_EQ(
            expectedNewSum,
            sum
    );
}

TEST(ExampleTests, Square){
    int x = 5;
    int expectedSquare = x*x;
    EXPECT_EQ(
          expectedSquare,
          Square(x)
          );
}