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

Future implementations will include:

- **SO(3)**: Special Orthogonal group in 3D (rotations in 3D space)
- **SE(3)**: Special Euclidean group in 3D (rigid transformations in 3D space)
- **Sim(3)**: Similarity transformations in 3D space

## Installation

The Lie groups library is part of the Cosserat plugin. Installation requirements:

- C++14 or higher
- Eigen 3.3 or higher
- Sofa Framework

## Dependencies

- Eigen: For linear algebra operations
- CMake: For building the project

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
- [Performance Benchmarks](docs/benchmarks.md)

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
