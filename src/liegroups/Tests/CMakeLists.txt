cmake_minimum_required(VERSION 3.12)

# Set project for liegroups tests
project(Cosserat_liegroups_tests)

# Add test sources
set(LIEGROUPS_TEST_FILES
    liegroups/SO2Test.cpp
    liegroups/UtilsTest.cpp
)

# If unit tests are enabled
if(UNIT_TEST)
    # Create test executable
    add_executable(${PROJECT_NAME} ${LIEGROUPS_TEST_FILES})
    
    # Link against required libraries
    target_link_libraries(${PROJECT_NAME} 
        Cosserat
        gtest 
        gtest_main
        SofaGTestMain
    )
    
    # Include directories
    target_include_directories(${PROJECT_NAME} PRIVATE
        ${CMAKE_SOURCE_DIR}/src
    )
    
    # Add test
    add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
endif()

# Add benchmarks if enabled
if(COSSERAT_BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()
