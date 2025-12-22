cmake_minimum_required(VERSION 3.12)
project(Cosserat_examples)

# Find Eigen
find_package(Eigen3 REQUIRED)

# Simple trajectory optimization example
add_executable(simple_trajectory_optimization
    simple_trajectory_optimization.cpp
)

target_include_directories(simple_trajectory_optimization PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${EIGEN3_INCLUDE_DIR}
)

target_link_libraries(simple_trajectory_optimization
    Eigen3::Eigen
)

# Set C++20 standard
target_compile_features(simple_trajectory_optimization PRIVATE cxx_std_20)

# Install target
install(TARGETS simple_trajectory_optimization
    RUNTIME DESTINATION bin/examples
)
