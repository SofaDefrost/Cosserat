//
// Created by younes on 31/05/2021.
//

#include <iostream>
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "ClassName.h"


using testing::Test;
using testing::Eq;

namespace{
    class classDeclaration : public testing::Test {
    public:
        ClassName obj;
        classDeclaration(){
            obj; //
        }
    };
}

int main(int argc, char *argv[])
{
    testing::InitGoogle(argc, argv);
    RUN_ALL_TESTS();
    return 0;
}