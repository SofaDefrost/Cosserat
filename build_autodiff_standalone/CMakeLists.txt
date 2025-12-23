cmake_minimum_required(VERSION 3.12)
project(CosseratAutodiffExamples)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Eigen
find_package(Eigen3 REQUIRED)

# Use autodiff as header-only (no need for find_package)
set(AUTODIFF_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../autodiff")
if(NOT EXISTS "${AUTODIFF_INCLUDE_DIR}/autodiff")
    message(FATAL_ERROR "autodiff headers not found at ${AUTODIFF_INCLUDE_DIR}")
endif()
message(STATUS "Using autodiff headers from: ${AUTODIFF_INCLUDE_DIR}")

# Forward mode example
add_executable(autodiff_forward_mode
    ../examples/autodiff_forward_mode.cpp
)

target_include_directories(autodiff_forward_mode PRIVATE
    ${CMAKE_SOURCE_DIR}/..
    ${EIGEN3_INCLUDE_DIR}
    ${AUTODIFF_INCLUDE_DIR}
)

target_link_libraries(autodiff_forward_mode
    Eigen3::Eigen
)

target_compile_definitions(autodiff_forward_mode PRIVATE
    COSSERAT_WITH_AUTODIFF
)

# Reverse mode example
add_executable(autodiff_reverse_mode
    ../examples/autodiff_reverse_mode.cpp
)

target_include_directories(autodiff_reverse_mode PRIVATE
    ${CMAKE_SOURCE_DIR}/..
    ${EIGEN3_INCLUDE_DIR}
    ${AUTODIFF_INCLUDE_DIR}
)

target_link_libraries(autodiff_reverse_mode
    Eigen3::Eigen
)

target_compile_definitions(autodiff_reverse_mode PRIVATE
    COSSERAT_WITH_AUTODIFF
)

message(STATUS "===========================================")
message(STATUS "Autodiff Examples Configuration Complete")
message(STATUS "===========================================")
