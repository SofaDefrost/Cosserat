# Quick Start Guide: Lie Groups for Cosserat Rods

This guide gets you up and running with the Lie groups library in under 15 minutes. We'll cover the essentials for common Cosserat rod applications.

## What You'll Learn

- Basic usage of SE(3) for 3D rigid body configurations
- Computing relative transformations (strains)
- Simple interpolation between configurations
- Integration with Eigen for linear algebra

## Prerequisites

- C++14 compiler
- Eigen 3.3+
- Basic understanding of 3D transformations

## Step 1: Setting Up Your First Transformation

```cpp
#include <Cosserat/liegroups/SE3.h>
#include <iostream>

int main() {
    // Create an SE(3) element: 90° rotation around Z + translation (1,2,3)
    Eigen::AngleAxisd rotation(M_PI/2, Eigen::Vector3d::UnitZ());
    Eigen::Vector3d translation(1.0, 2.0, 3.0);

    Cosserat::SE3<double> pose(rotation, translation);

    // Print the transformation matrix
    std::cout << "Transformation matrix:\n" << pose.matrix() << "\n";

    return 0;
}
```

**What this does**: Creates a rigid body pose combining rotation and translation.

## Step 2: Composing Transformations

```cpp
// Create two poses
Cosserat::SE3<double> pose1(
    Eigen::AngleAxisd(M_PI/4, Eigen::Vector3d::UnitZ()),
    Eigen::Vector3d(1.0, 0.0, 0.0)
);

Cosserat::SE3<double> pose2(
    Eigen::AngleAxisd(M_PI/4, Eigen::Vector3d::UnitX()),
    Eigen::Vector3d(0.0, 1.0, 0.0)
);

// Compose them: apply pose2 after pose1
Cosserat::SE3<double> composed = pose1.compose(pose2);

// Apply to a point
Eigen::Vector3d point(1.0, 1.0, 1.0);
Eigen::Vector3d transformed = composed.act(point);

std::cout << "Original point: " << point.transpose() << "\n";
std::cout << "Transformed point: " << transformed.transpose() << "\n";
```

**Key Concept**: `compose(A, B)` means "apply B after A" - order matters!

## Step 3: Computing Strains (Relative Transformations)

```cpp
// Two consecutive rod element poses
Cosserat::SE3<double> element_i(/* ... */);
Cosserat::SE3<double> element_ip1(/* ... */);

// Compute relative deformation (strain)
Cosserat::SE3<double> relative = element_i.inverse().compose(element_ip1);

// Convert to strain vector (Lie algebra)
Eigen::Matrix<double, 6, 1> strain = relative.log();

std::cout << "Strain vector (ω_x, ω_y, ω_z, v_x, v_y, v_z):\n" << strain.transpose() << "\n";
```

**Cosserat Application**: This strain vector represents the deformation between rod elements.

## Step 4: Interpolation Between Configurations

```cpp
Cosserat::SE3<double> start_pose(/* ... */);
Cosserat::SE3<double> end_pose(/* ... */);

// Convert to Lie algebra for linear interpolation
auto start_vec = start_pose.log();
auto end_vec = end_pose.log();

// Interpolate at t=0.5
double t = 0.5;
auto mid_vec = start_vec + t * (end_vec - start_vec);

// Convert back to group
Cosserat::SE3<double> mid_pose = Cosserat::SE3<double>::exp(mid_vec);

std::cout << "Interpolated pose at t=0.5:\n" << mid_pose.matrix() << "\n";
```

**Why this works**: The Lie algebra is a vector space where linear operations make sense.

## Step 5: Working with Velocities (SE_2(3))

```cpp
#include <Cosserat/liegroups/SE23.h>

// Create pose with velocity
Cosserat::SE3<double> pose(/* rotation and translation */);
Eigen::Vector3d velocity(0.1, 0.2, 0.3);

Cosserat::SE23<double> dynamic_pose(pose, velocity);

// Act on point-velocity pair (6D)
Eigen::Matrix<double, 6, 1> point_velocity;
point_velocity << 1.0, 1.0, 1.0, 0.0, 0.0, 0.0; // position + velocity

auto transformed = dynamic_pose.act(point_velocity);
```

**Use Case**: Dynamic Cosserat rods where velocity is part of the state.

## Step 6: Basic Optimization Example

```cpp
// Simple gradient descent to find pose that minimizes some objective
Cosserat::SE3<double> current_pose = Cosserat::SE3<double>::identity();
double learning_rate = 0.01;

for (int iter = 0; iter < 100; ++iter) {
    // Compute objective gradient in Lie algebra
    Eigen::Matrix<double, 6, 1> gradient = computeObjectiveGradient(current_pose);

    // Update in Lie algebra
    auto tangent = current_pose.log();
    tangent -= learning_rate * gradient;

    // Project back to group
    current_pose = Cosserat::SE3<double>::exp(tangent);
}
```

## Common Patterns and Best Practices

### 1. Always Check Your Coordinate Frames
```cpp
// Make sure transformations are in the same coordinate frame
Cosserat::SE3<double> world_to_body = /* ... */;
Cosserat::SE3<double> world_to_sensor = /* ... */;

// Transform sensor measurement to body frame
Cosserat::SE3<double> body_to_sensor = world_to_body.inverse().compose(world_to_sensor);
```

### 2. Use Log/Exp for Small Updates
```cpp
// For small perturbations, work in tangent space
Eigen::Matrix<double, 6, 1> small_update = 0.01 * Eigen::Matrix<double, 6, 1>::Random();
Cosserat::SE3<double> updated_pose = pose.compose(Cosserat::SE3<double>::exp(small_update));
```

### 3. Prefer Quaternions for SO(3)
```cpp
// For pure rotations, use SO(3) with quaternion representation
Eigen::Quaterniond quat(Eigen::AngleAxisd(M_PI/2, Eigen::Vector3d::UnitZ()));
Cosserat::SO3<double> rotation(quat);
```

## Next Steps

1. **Explore Advanced Groups**: Check out Sim(3) for scaling, SGal(3) for time evolution
2. **Read Detailed Docs**: Each group has comprehensive documentation
3. **Look at Examples**: Check the examples/ directory for complete simulations
4. **Performance**: See benchmarks.md for optimization tips

## Troubleshooting

- **Compilation errors**: Ensure Eigen headers are found
- **Runtime crashes**: Check matrix dimensions and initialization
- **Unexpected results**: Verify transformation composition order
- **Numerical issues**: Normalize quaternions periodically

This quick start covers the 80% of use cases for Cosserat rod modeling. For advanced topics, refer to the detailed documentation.