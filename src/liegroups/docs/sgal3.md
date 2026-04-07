# SGal(3) Implementation

## Overview

`SGal3<Scalar>` implements the Special Galilean group in 3D, which represents Galilean transformations in 3D space. SGal(3) is a 10-dimensional Lie group that combines a 6-dimensional SE(3) transformation (rotation and translation), a 3-dimensional linear velocity vector, and a 1-dimensional time parameter.

## Mathematical Properties

- **Dimension**: 10 (6 for SE(3) + 3 for linear velocity + 1 for time)
- **Group operation**: composition of Galilean transformations
- **Identity element**: zero rotation, zero translation, zero velocity, and zero time
- **Inverse**: inverse rotation, scaled and rotated negative translation, transformed velocity, and negated time
- **Lie algebra**: sgal(3), which can be represented as a 10D vector [ω, v, β, τ]
  where ω is the rotation part, v is the translation part, β is the boost (velocity change), and τ is the time rate
- **Exponential map**: converts from sgal(3) to SGal(3)
- **Logarithmic map**: converts from SGal(3) to sgal(3)

## Internal Representation

The `SGal3` class internally stores:
- A rotation component as an SO(3) element (typically as a quaternion)
- A translation component as a 3D vector
- A linear velocity component as a 3D vector
- A time parameter as a scalar

This representation ensures proper handling of the group structure and efficient computation of group operations.

## Implementation Details

The `SGal3` class is implemented as a template with one parameter:
- `Scalar`: The scalar type (typically `double` or `float`)

### Key Methods

```cpp
// Constructors
SGal3(); // Identity transformation
SGal3(const SE3<Scalar>& pose, const Vector3& velocity, const Scalar& time);
SGal3(const SO3<Scalar>& rotation, const Vector3& position, const Vector3& velocity, const Scalar& time);

// Group operations
SGal3 compose(const SGal3& other) const;
SGal3 inverse() const;

// Access to components
SE3<Scalar> pose() const;
Vector3 velocity() const;
Scalar time() const;

// Conversion to extended matrix
Eigen::Matrix<Scalar, 6, 6> extendedMatrix() const;

// Tangent space (Lie algebra) operations
Eigen::Matrix<Scalar, 10, 1> log() const;
static SGal3 exp(const Eigen::Matrix<Scalar, 10, 1>& tangent);

// Acting on extended points (position + velocity + time)
Eigen::Matrix<Scalar, 10, 1> act(const Eigen::Matrix<Scalar, 10, 1>& point_vel_time) const;

// Adjoint representation
Eigen::Matrix<Scalar, 10, 10> adjoint() const;
```

## Performance Characteristics

Based on benchmarks, SGal(3) operations have the following performance characteristics:

- **Composition**: O(1) time complexity - involves SE(3) composition, velocity transformation, and time addition
- **Inverse**: O(1) time complexity - computes inverse SE(3), transformed velocity, and negated time
- **Exponential map**: O(1) time complexity - converts from Lie algebra to the group
- **Logarithmic map**: O(1) time complexity - converts from the group to Lie algebra
- **Acting on points**: O(1) time complexity - applies the transformation to a 10D point-velocity-time vector

The actual performance depends on the hardware and compiler optimizations. The implementation leverages the SE(3) operations for the pose component.

## Example Usage

### Basic Operations

```cpp
#include <Cosserat/liegroups/SGal3.h>
#include <iostream>

int main() {
    // Create an SGal(3) element (90-degree rotation around z-axis, translation (1,2,3), velocity (0.1,0.2,0.3), time 1.0)
    Eigen::AngleAxisd rotation(M_PI/2, Eigen::Vector3d::UnitZ());
    Eigen::Vector3d translation(1.0, 2.0, 3.0);
    Eigen::Vector3d velocity(0.1, 0.2, 0.3);
    double time = 1.0;
    Cosserat::SGal3<double> transform(
        Cosserat::SO3<double>(rotation), 
        translation,
        velocity,
        time
    );

    // Get components
    Cosserat::SE3<double> pose = transform.pose();
    Eigen::Vector3d vel = transform.velocity();
    double t = transform.time();
    std::cout << "Pose: " << pose << "\n";
    std::cout << "Velocity: " << vel.transpose() << "\n";
    std::cout << "Time: " << t << "\n";

    // Convert to extended homogeneous transformation matrix
    Eigen::Matrix<double, 6, 6> mat = transform.extendedMatrix();
    std::cout << "Extended matrix:\n" << mat << "\n";

    // Create another transformation
    Cosserat::SGal3<double> another_transform(
        Cosserat::SO3<double>(Eigen::AngleAxisd(M_PI/4, Eigen::Vector3d::UnitX())),
        Eigen::Vector3d(4.0, 5.0, 6.0),
        Eigen::Vector3d(0.4, 0.5, 0.6),
        2.0
    );

    // Compose transformations
    auto composed = transform.compose(another_transform);

    // Inverse transformation
    auto inverse = transform.inverse();

    return 0;
}
```

### Tangent Space Operations

```cpp
#include <Cosserat/liegroups/SGal3.h>
#include <iostream>

int main() {
    // Create an SGal(3) element
    Eigen::AngleAxisd rotation(M_PI/3, Eigen::Vector3d::UnitZ());
    Eigen::Vector3d translation(1.0, 2.0, 3.0);
    Eigen::Vector3d velocity(0.1, 0.2, 0.3);
    double time = 1.0;
    Cosserat::SGal3<double> transform(
        Cosserat::SO3<double>(rotation), 
        translation,
        velocity,
        time
    );

    // Convert to Lie algebra (tangent space)
    Eigen::Matrix<double, 10, 1> tangent = transform.log();
    std::cout << "Tangent vector: " << tangent.transpose() << "\n";

    // Convert back from Lie algebra to SGal(3)
    auto recovered = Cosserat::SGal3<double>::exp(tangent);

    // Create directly from tangent vector
    Eigen::Matrix<double, 10, 1> new_tangent;
    new_tangent << 0.1, 0.2, 0.3, 1.0, 2.0, 3.0, 0.1, 0.2, 0.3, 0.5; // ω, v, β, τ
    auto new_transform = Cosserat::SGal3<double>::exp(new_tangent);

    return 0;
}
```

### Acting on Extended Points

```cpp
#include <Cosserat/liegroups/SGal3.h>
#include <iostream>

int main() {
    // Create a transformation
    Eigen::AngleAxisd rotation(M_PI/2, Eigen::Vector3d::UnitZ());
    Eigen::Vector3d translation(1.0, 2.0, 3.0);
    Eigen::Vector3d velocity(0.1, 0.2, 0.3);
    double time = 1.0;
    Cosserat::SGal3<double> transform(
        Cosserat::SO3<double>(rotation), 
        translation,
        velocity,
        time
    );

    // Transform a point-velocity-time tuple (10D vector: position + velocity + time)
    Eigen::Matrix<double, 10, 1> point_vel_time;
    point_vel_time << 1.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0; // position + velocity + boost + time
    Eigen::Matrix<double, 10, 1> transformed = transform.act(point_vel_time);
    std::cout << "Original point-velocity-time: " << point_vel_time.transpose() << "\n";
    std::cout << "Transformed point-velocity-time: " << transformed.transpose() << "\n";

    return 0;
}
```

## Applications

SGal(3) is particularly useful for:

- **Relativistic mechanics**: Modeling transformations that include time evolution
- **Multi-reference frame dynamics**: Systems where different parts move with different velocities
- **Galilean-invariant physics**: Simulations requiring Galilean transformation properties
- **Time-varying systems**: Dynamic systems where time is an explicit parameter

## Best Practices

1. **Use appropriate constructors** based on available data (full Galilean transform vs individual components)
2. **Cache extended matrices** when performing multiple operations with the same transformation
3. **For time integration**, use the exponential map with appropriate time steps
4. **When composing many transformations**, consider the computational cost of SE(3) operations
5. **For interpolation**, work in the Lie algebra space for smooth interpolation
6. **Remember that SGal(3) includes time evolution** - transformations affect both spatial and temporal coordinates

## Numerical Stability Considerations

- **Time parameter handling**: Ensure proper time evolution in transformations
- **Velocity transformations**: Ensure proper handling of velocity and boost components during composition
- **Small angle approximations**: Use specialized implementations for near-identity elements
- **Extended matrix operations**: Be aware of potential numerical drift in long transformation chains</content>
</xai:function_call/>
<xai:function_call name="write">
<parameter name="filePath">src/liegroups/docs/bundle.md