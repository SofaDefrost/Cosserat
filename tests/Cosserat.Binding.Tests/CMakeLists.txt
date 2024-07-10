cmake_minimum_required(VERSION 3.12)
project(Cosserat_core_test C CXX)

find_package(Sofa.Testing REQUIRED)

enable_testing()

set(HEADER_FILES
        constraint/Constraint.h
        )
set(SOURCE_FILES
        constraint/ExampleTest.cpp
        forcefield/BeamHookeLawForceFieldTest.cpp
    )


add_executable(${PROJECT_NAME} ${SOURCE_FILES} ${HEADER_FILES})

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
