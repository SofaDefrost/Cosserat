cmake_minimum_required(VERSION 3.12)
project(Cosserat_liegroups_tests VERSION 1.4)
find_package(Sofa.Testing REQUIRED)
sofa_find_package(Sofa.SimpleApi REQUIRED)

# Set project for liegroups tests


# Add test sources
set(LIEGROUPS_TEST_FILES
        integration/SO2Test.cpp
        integration/UtilsTest.cpp
        integration/LieGroupBaseTest.cpp
        integration/TypesTest.cpp
        integration/RealSpaceTest.cpp
        integration/SE2Test.cpp
        integration/SE23Test.cpp
    
        # Unit tests
        unit/test_jacobian_verification.cpp
        unit/test_HookeSerat_Phase1.cpp
        unit/test_Phase2.cpp
        unit/test_HookeSerat_Comprehensive.cpp
        unit/test_MultiSectionBeam.cpp
        unit/test_AdvancedFeatures.cpp
)

# If unit tests are enabled
message("=======> Buildin Cosserat lieGroup test", ${UNIT_TEST})

    # Create test executable
    message("=======> Buildin Cosserat lieGroup test" )
    add_executable(${PROJECT_NAME} ${LIEGROUPS_TEST_FILES})
    # Link against required libraries
    target_link_libraries(${PROJECT_NAME}
        Sofa.Testing
        Cosserat
        Sofa.SimpleApi
        #gtest
        #gtest_main
        #SofaGTestMain
    )

    # Include directories
    target_include_directories(${PROJECT_NAME} PRIVATE
        ${CMAKE_SOURCE_DIR}/src
    )
    add_definitions("-DCOSSETRAT_TEST_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/\"")

    # Add test
    add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})



# Add benchmarks subdirectory if benchmarks are enabled
if(COSSERAT_BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()

