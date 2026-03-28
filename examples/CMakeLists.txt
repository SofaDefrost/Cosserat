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

# Autodiff examples (only if autodiff is enabled)
if(COSSERAT_WITH_AUTODIFF)
    message(STATUS "Building autodiff examples")
    
    # Forward mode example
    add_executable(autodiff_forward_mode
        autodiff_forward_mode.cpp
    )
    
    target_include_directories(autodiff_forward_mode PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${EIGEN3_INCLUDE_DIR}
    )
    
    target_link_libraries(autodiff_forward_mode
        Eigen3::Eigen
        autodiff::autodiff
    )
    
    target_compile_features(autodiff_forward_mode PRIVATE cxx_std_20)
    
    # Reverse mode example
    add_executable(autodiff_reverse_mode
        autodiff_reverse_mode.cpp
    )
    
    target_include_directories(autodiff_reverse_mode PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${EIGEN3_INCLUDE_DIR}
    )
    
    target_link_libraries(autodiff_reverse_mode
        Eigen3::Eigen
        autodiff::autodiff
    )
    
    target_compile_features(autodiff_reverse_mode PRIVATE cxx_std_20)
    
    # Install autodiff examples
    install(TARGETS autodiff_forward_mode autodiff_reverse_mode
        RUNTIME DESTINATION bin/examples
    )
else()
    message(STATUS "Autodiff examples disabled (COSSERAT_WITH_AUTODIFF=OFF)")
endif()
