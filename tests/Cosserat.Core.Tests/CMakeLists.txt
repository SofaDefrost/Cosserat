cmake_minimum_required(VERSION 3.12)
project(Cosserat.Core.Tests C CXX)

find_package(Sofa.Testing REQUIRED)

enable_testing()

set(HEADER_FILES)

set(SOURCE_FILES
        forcefield/BeamHookeLawForceFieldTest.cpp
    )

add_executable(${PROJECT_NAME} ${SOURCE_FILES} ${HEADER_FILES})

target_link_libraries(${PROJECT_NAME} 
        Sofa.Testing
        Cosserat
)

target_include_directories(${PROJECT_NAME}
        PUBLIC 
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../>"
)


add_test(
        NAME ${PROJECT_NAME}
        COMMAND ${PROJECT_NAME}
)

