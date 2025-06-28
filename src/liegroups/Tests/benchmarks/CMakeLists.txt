cmake_minimum_required(VERSION 3.12)

project(Cosserat_liegroups_benchmarks)

# Set benchmark source files
set(BENCHMARK_SOURCE_FILES
    liegroups/LieGroupsBenchmark.cpp
)

# Create benchmark executable
add_executable(${PROJECT_NAME} ${BENCHMARK_SOURCE_FILES})

# Link against required libraries
target_link_libraries(${PROJECT_NAME} 
    Cosserat
    benchmark::benchmark
    benchmark::benchmark_main
)

# Include directories
target_include_directories(${PROJECT_NAME} PRIVATE
    ${CMAKE_SOURCE_DIR}/src
)

# Add benchmark target
add_custom_target(run_liegroups_benchmarks
    COMMAND ${PROJECT_NAME}
    DEPENDS ${PROJECT_NAME}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Running Cosserat liegroups benchmarks"
)

