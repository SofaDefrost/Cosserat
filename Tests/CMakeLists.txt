cmake_minimum_required(VERSION 3.1)

set(This Cosserat_test)

project(${This} C CXX)

#find_package(Cosserat REQUIRED)
find_package(Sofa.Testing REQUIRED)

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
        Sofa.Testing
        Cosserat
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
