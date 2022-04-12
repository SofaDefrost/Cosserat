cmake_minimum_required(VERSION 3.1)

set(This Cosserat.Tests)

project(${This} C CXX)

#find_package(CosseratPlugin REQUIRED)
find_package(SofaTest REQUIRED)

enable_testing()

set(HEADER_FILES
        Example.h
        constraint/Constraint.h
        )
set(SOURCE_FILES
        Example.cpp
        constraint/ExampleTest.cpp
#        constraint/CosseratUnilateralInteractionConstraintTest.cpp
        forcefield/BeamHookeLawForceFieldTest.cpp
    )


add_executable(${This} ${SOURCE_FILES} ${HEADER_FILES})

target_link_libraries(${PROJECT_NAME} 
        SofaTest 
        gtest
        CosseratPlugin
)

target_include_directories(${This} 
        PUBLIC 
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../>"
)


add_test(
        NAME ${This} 
        COMMAND ${This}
)

#[[add_subdirectory(constraint)]]
