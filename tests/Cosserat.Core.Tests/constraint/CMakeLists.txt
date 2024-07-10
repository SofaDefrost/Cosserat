cmake_minimum_required(VERSION 3.1)

#[[project(Cosserat.test VERSION 1.0)]]
set(This ExampleTest)

set(SOURCE_FILES
        ExampleTest.cpp
    )

add_executable(${This} ${SOURCE_FILES})

target_link_libraries(${This} PUBLIC  
        SofaTest 
        gtest 
        gtest_main 
        Cosserat
        Example
)

add_test(
        NAME ${This} 
        COMMAND ${This}
)
