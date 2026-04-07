---
author: Yinoussa Adagolodjo
date: 2025-04-12
---

# Lie Groups Library for Cosserat Models

This library provides implementations of various Lie groups for use in Cosserat rod modeling and simulation. It is inspired by the [manif](https://github.com/artivis/manif) library but tailored specifically for the needs of the Cosserat plugin.

## Overview

Lie groups are mathematical structures that are both groups and differentiable manifolds. They are essential in robotics and computer vision for representing rigid body transformations and rotations. In the context of Cosserat rods, they allow us to represent the configuration of rod elements in a mathematically elegant and computationally efficient way.

This library implements the following Lie groups:

- **RealSpace**: Euclidean vector space ℝⁿ
- **SO(2)**: Special Orthogonal group in 2D (rotations in a plane)
- **SE(2)**: Special Euclidean group in 2D (rigid transformations in a plane)
- **SO(3)**: Special Orthogonal group in 3D (rotations in 3D space)
- **SE(3)**: Special Euclidean group in 3D (rigid transformations in 3D space)
- **Sim(3)**: Similarity transformations in 3D space (rotation + translation + scaling)
- **SE(2,3)**: Extended Special Euclidean group in 3D (rigid transformations with linear velocity)
- **SGal(3)**: Special Galilean group in 3D (Galilean transformations with time)

Additional utilities:

- **Bundle**: Product manifold for combining multiple Lie groups
- **GaussianOnManifold**: Gaussian distributions on Lie groups for uncertainty propagation

## Installation

The Lie groups library is part of the Cosserat plugin. Follow these steps to set up the development environment:

### Prerequisites

- **C++ Compiler**: C++14 or higher (GCC 7+, Clang 5+, MSVC 2017+)
- **CMake**: Version 3.10 or higher
- **Eigen**: Version 3.3 or higher
- **SOFA Framework**: Compatible version for the Cosserat plugin

### Building from Source

1. **Clone the repository**:
   ```bash
   git clone https://github.com/your-repo/cosserat-plugin.git
   cd cosserat-plugin
   ```

2. **Create a build directory**:
   ```bash
   mkdir build && cd build
   ```

3. **Configure with CMake**:
   ```bash
   cmake .. -DCMAKE_BUILD_TYPE=Release
   ```

   For development with debug symbols:
   ```bash
   cmake .. -DCMAKE_BUILD_TYPE=Debug
   ```

4. **Build the project**:
   ```bash
   make -j$(nproc)
   ```

5. **Install (optional)**:
   ```bash
   make install
   ```

### Integration with SOFA

To use the Lie groups in your SOFA-based Cosserat simulations:

1. **Include headers**:
   ```cpp
   #include <Cosserat/liegroups/SE3.h>
   #include <Cosserat/liegroups/SO3.h>
   // etc.
   ```

2. **Link libraries**:
   In your CMakeLists.txt:
   ```cmake
   find_package(Cosserat REQUIRED)
   target_link_libraries(your_target Cosserat::liegroups)
   ```

### Running Tests

After building, run the unit tests:
```bash
cd build
ctest --output-on-failure
```

Or run specific Lie groups tests:
```bash
./bin/test_liegroups
```

### Troubleshooting

#### Build Issues

- **Eigen not found**: Ensure Eigen is installed and CMake can find it. You may need to set `EIGEN3_INCLUDE_DIR`:
  ```bash
  cmake .. -DEIGEN3_INCLUDE_DIR=/path/to/eigen
  ```

- **SOFA compatibility**: Check that your SOFA version is compatible with the Cosserat plugin. Refer to the plugin documentation for version requirements.

- **Compilation errors**: Ensure your compiler supports C++14 features. Update to GCC 7+, Clang 5+, or MSVC 2017+ if needed.

- **CMake errors**: Clear the build directory and reconfigure:
  ```bash
  rm -rf build && mkdir build && cd build && cmake ..
  ```

#### Runtime Issues

- **Linking errors**: Ensure the Lie groups library is properly linked. Check your CMakeLists.txt includes:
  ```cmake
  target_link_libraries(your_target Cosserat::liegroups)
  ```

- **Memory issues**: Some operations allocate matrices. Ensure adequate stack/heap space for large problems.

#### Usage Issues

- **Unexpected transformation results**: Double-check composition order. `A.compose(B)` applies B after A.

- **Numerical instability**: For long transformation chains, consider periodic normalization:
  ```cpp
  pose = pose.normalize(); // If available
  ```

- **Performance problems**: Cache frequently used matrices and avoid repeated conversions between representations.

#### Common Mistakes

1. **Wrong composition order**: Remember that `A.compose(B)` means "B after A"
2. **Uninitialized matrices**: Always initialize Eigen matrices before use
3. **Dimension mismatches**: Check vector/matrix sizes match group dimensions
4. **Coordinate frame confusion**: Ensure all transformations use consistent coordinate frames

#### Getting Help

- Check the [detailed documentation](docs/) for your specific Lie group
- Look at [examples](../examples/) for usage patterns
- Run [unit tests](../tests/) to verify your setup
- Check [benchmarks](docs/benchmarks.md) for performance guidance

## Dependencies

- **Eigen**: For linear algebra operations and matrix computations
- **CMake**: For building the project and managing dependencies
- **SOFA**: Framework integration for Cosserat rod simulations

## Basic Usage

Here are some examples of how to use the library:

### RealSpace (Euclidean vector space)

```cpp
#include <Cosserat/liegroups/RealSpace.h>

// Create a 3D vector in RealSpace
Cosserat::RealSpace<double, 3> point(1.0, 2.0, 3.0);

// Create from Eigen vector
Eigen::Vector3d eigen_vec(4.0, 5.0, 6.0);
Cosserat::RealSpace<double, 3> another_point(eigen_vec);

// Compose points (vector addition)
auto result = point.compose(another_point);
```

### SO(2) (2D rotations)

```cpp
#include <Cosserat/liegroups/SO2.h>

// Create a rotation of 45 degrees
Cosserat::SO2<double> rotation(M_PI/4.0);

// Get the angle
double angle = rotation.angle();

// Convert to rotation matrix
Eigen::Matrix2d rot_mat = rotation.matrix();

// Compose rotations
Cosserat::SO2<double> another_rotation(M_PI/2.0);
auto composed = rotation.compose(another_rotation); // 135 degrees rotation

// Inverse rotation
auto inverse = rotation.inverse();
```

### SE(2) (2D rigid transformations)

```cpp
#include <Cosserat/liegroups/SE2.h>

// Create a transform: 45-degree rotation with translation (1,2)
Cosserat::SE2<double> transform(M_PI/4.0, 1.0, 2.0);

// Get components
double angle = transform.angle();
Eigen::Vector2d translation = transform.translation();

// Convert to homogeneous transformation matrix
Eigen::Matrix3d transform_mat = transform.matrix();

// Compose transformations
Cosserat::SE2<double> another_transform(M_PI/2.0, 3.0, 4.0);
auto composed = transform.compose(another_transform);

// Inverse transformation
auto inverse = transform.inverse();
```

## Detailed Documentation

For more detailed documentation, including mathematical foundations, implementation details, and advanced usage examples, see:

- [Mathematical Foundations](docs/math_foundations.md)
- [RealSpace Implementation](docs/realspace.md)
- [SO(2) Implementation](docs/so2.md)
- [SE(2) Implementation](docs/se2.md)
- [SO(3) Implementation](docs/so3.md)
- [SE(3) Implementation](docs/se3.md)
- [Sim(3) Implementation](docs/sim3.md)
- [SE(2,3) Implementation](docs/se23.md)
- [SGal(3) Implementation](docs/sgal3.md)
- [Bundle Implementation](docs/bundle.md)
- [Gaussian on Manifold](docs/gaussian_on_manifold.md)
- [Advanced Topics](docs/advanced_topics.md)
- [Usage Guide](docs/USAGE.md)
- [Performance Benchmarks](docs/benchmarks.md)
- [Comparison of Implementations](docs/comparison.md)
- [Dependency Tree](docs/dependency_tree.md)

## Benchmarking

The library includes benchmarks to measure the performance of various operations. You can run the benchmarks with:

```bash
cd build
make run_benchmarks
```

## Roadmap

For future development plans, see the [Roadmap](/src/Cosserat/liegroups/tasks.md/feature-lieAgebra.md).

This feature is inspire by the repo : https://github.com/artivis/manif and the devel branch.

## Todo see : [Roadmap](/src/Cosserat/liegroups/tasks.md/feature-lieAgebra.md)
