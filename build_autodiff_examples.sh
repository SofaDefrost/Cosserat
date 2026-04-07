#!/bin/bash

# Build script for Cosserat autodiff examples
# This script builds the examples standalone without requiring full SOFA

set -e  # Exit on error

echo "════════════════════════════════════════════════════════════════"
echo "    Building Cosserat Autodiff Examples (Standalone)"
echo "════════════════════════════════════════════════════════════════"
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BUILD_DIR="${SCRIPT_DIR}/build_autodiff_standalone"

echo "Project directory: ${SCRIPT_DIR}"
echo "Build directory:   ${BUILD_DIR}"
echo ""

# Check if autodiff exists
AUTODIFF_DIR="${SCRIPT_DIR}/../autodiff"
if [ ! -d "${AUTODIFF_DIR}" ]; then
    echo -e "${RED}✗ Error: autodiff not found at ${AUTODIFF_DIR}${NC}"
    echo ""
    echo "Please ensure autodiff library is available at ../autodiff/"
    echo "Get it from: https://github.com/autodiff/autodiff"
    exit 1
fi

echo -e "${GREEN}✓ Found autodiff at ${AUTODIFF_DIR}${NC}"
echo ""

# Create build directory
mkdir -p "${BUILD_DIR}"
cd "${BUILD_DIR}"

echo "Configuring CMake..."
echo "────────────────────────────────────────────────────────────────"

# Create minimal CMakeLists.txt for standalone build
cat > CMakeLists.txt << 'EOF'
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
EOF

# Run CMake
cmake . || {
    echo -e "${RED}✗ CMake configuration failed${NC}"
    exit 1
}

echo ""
echo "Building examples..."
echo "────────────────────────────────────────────────────────────────"

# Build
make -j$(sysctl -n hw.ncpu) || {
    echo -e "${RED}✗ Build failed${NC}"
    exit 1
}

echo ""
echo "════════════════════════════════════════════════════════════════"
echo -e "${GREEN}✓ Build completed successfully!${NC}"
echo "════════════════════════════════════════════════════════════════"
echo ""
echo "Executables created:"
echo "  • autodiff_forward_mode"
echo "  • autodiff_reverse_mode"
echo ""
echo "Run examples:"
echo "  cd ${BUILD_DIR}"
echo "  ./autodiff_forward_mode"
echo "  ./autodiff_reverse_mode"
echo ""
