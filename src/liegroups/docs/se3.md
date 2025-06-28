# SE(3) Implementation

## Overview

`SE3<Scalar>` implements the Special Euclidean group in 3D, which represents rigid body transformations (rotation and translation) in 3D space. SE(3) is a 6-dimensional Lie group, combining a 3-dimensional rotation (SO(3)) with a 3-dimensional translation.

## Mathematical Properties

- **Dimension**: 6 (3 for rotation + 3 for translation)
- **Group operation**: composition of transformations
- **Identity element**: zero rotation and zero translation
- **Inverse**: inverse rotation and negated rotated translation
- **Lie algebra**: se(3), which can be represented as a 6D vector [ω, v]
  where ω is the rotation part and v is the translation part
- **Exponential map**: converts from se(3) to SE(3)
- **Logarithmic map**: converts from SE(3) to se(3)

## Internal Representation

The `SE3` class internally stores:
- A rotation component as an SO(3) element (typically as a quaternion)
- A translation component as a 3D vector

This representation ensures proper handling of the group structure and efficient computation of group operations.

## Implementation Details

The `SE3` class is implemented as a template with one parameter:
- `Scalar`: The scalar type (typically `double` or `float`)

### Key Methods

```cpp
// Constructors
SE3(); // Identity transformation
SE3(const SO3<Scalar>& rotation, const Vector3<Scalar>& translation);
SE3(const Matrix4<Scalar>& homogeneous_matrix);
SE3(const Quaternion<Scalar>& quat, const Vector3<Scalar>& translation);

// Group operations
SE3 compose(const SE3& other) const;
SE3 inverse() const;

// Access to components
SO3<Scalar> rotation() const;
Vector3<Scalar> translation() const;
Matrix4<Scalar> matrix() const; // Homogeneous transformation matrix

// Tangent space (Lie algebra) operations
Vector6<Scalar> log() const;
static SE3 exp(const Vector6<Scalar>& tangent);

// Acting on points
Vector3<Scalar> act(const Vector3<Scalar>& point) const;

// Adjoint representation
Matrix6<Scalar> adjoint() const;
```

## Performance Characteristics

Based on benchmarks, SE(3) operations have the following performance characteristics:

- **Composition**: O(1) time complexity - involves rotation composition and translation transformation
- **Inverse**: O(1) time complexity - computes inverse rotation and inverse transformed translation
- **Matrix conversion**: O(1) time complexity - creates a 4×4 homogeneous transformation matrix
- **Exponential map**: O(1) time complexity - converts from Lie algebra to the group
- **Logarithmic map**: O(1) time complexity - converts from the group to Lie algebra
- **Acting on points**: O(1) time complexity - applies the transformation to a 3D point

The actual performance depends on the hardware and compiler optimizations. Using quaternions for the rotational component provides better performance for most operations compared to rotation matrices.

## Example Usage

### Basic Operations

```cpp
#include <Cosserat/liegroups/SE3.h>
#include <iostream>

int main() {
    // Create an SE(3) element (90-degree rotation around z-axis with translation (1,2,3))
    Eigen::AngleAxisd rotation(M_PI/2, Eigen::Vector3d::UnitZ());
    Eigen::Vector3d translation(1.0, 2.0, 3.0);
    Cosserat::SE3<double> transform(Cosserat::SO3<double>(rotation), translation);
    
    // Get components
    Cosserat::SO3<double> rot = transform.rotation();
    Eigen::Vector3d trans = transform.translation();
    std::cout << "Rotation matrix:\n" << rot.matrix() << "\n";
    std::cout << "Translation: " << trans.transpose() << "\n";
    
    // Convert to homogeneous transformation matrix
    Eigen::Matrix4d mat = transform.matrix();
    std::cout << "Transformation matrix:\n" << mat << "\n";
    
    // Create another transformation
    Eigen::AngleAxisd another_rotation(M_PI/4, Eigen::Vector3d::UnitX());
    Eigen::Vector3d another_translation(4.0, 5.0, 6.0);
    Cosserat::SE3<double> another_transform(
        Cosserat::SO3<double>(another_rotation), 
        another_translation
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
#include <Cosserat/liegroups/SE3.h>
#include <iostream>

int main() {
    // Create an SE(3) element
    Eigen::AngleAxisd rotation(M_PI/3, Eigen::Vector3d::UnitZ());
    Eigen::Vector3d translation(1.0, 2.0, 3.0);
    Cosserat::SE3<double> transform(Cosserat::SO3<double>(rotation), translation);
    
    // Convert to Lie algebra (tangent space)
    Eigen::Matrix<double, 6, 1> tangent = transform.log();
    std::cout << "Tangent vector: " << tangent.transpose() << "\n";
    
    // Convert back from Lie algebra to SE(3)
    auto recovered = Cosserat::SE3<double>::exp(tangent);
    
    // Create directly from tangent vector
    Eigen::Matrix<double, 6, 1> new_tangent;
    new_tangent << 0.1, 0.2, 0.3, 1.0, 2.0, 3.0; // Small rotation with translation
    auto new_transform = Cosserat::SE3<double>::exp(new_tangent);
    
    return 0;
}
```

### Acting on Points

```cpp
#include <Cosserat/liegroups/SE3.h>
#include <iostream>
#include <vector>

int main() {
    // Create a transformation
    Eigen::AngleAxisd rotation(M_PI/2, Eigen::Vector3d::UnitZ());
    Eigen::Vector3d translation(1.0, 2.0, 3.0);
    Cosserat::SE3<double> transform(Cosserat::SO3<double>(rotation), translation);
    
    // Transform a single point
    Eigen::Vector3d point(1.0, 0.0, 0.0);
    Eigen::Vector3d transformed_point = transform.act(point);
    std::cout << "Original point: " << point.transpose() << "\n";
    std::cout << "Transformed point: " << transformed_point.transpose() << "\n";
    
    // Transform multiple points
    std::vector<Eigen::Vector3d> points = {
        Eigen::Vector3d(1.0, 0.0, 0.0),
        Eigen::Vector3d(0.0, 1.0, 0.0),
        Eigen::Vector3d(0.0, 0.0, 1.0)
    };
    
    std::vector<Eigen::Vector3d> transformed_points;
    for (const auto& p : points) {
        transformed_points.push_back(transform.act(p));
    }
    
    // Check that applying the transformation twice is equivalent to composing
    Eigen::Vector3d twice_transformed = transform.act(transform.act(point));
    Eigen::Vector3d composed_transform = transform.compose(transform).act(point);
    
    return 0;
}
```

### Interpolation

```cpp
#include <Cosserat/liegroups/SE3.h>
#include <iostream>

int main() {
    // Create two transformations
    Eigen::AngleAxisd start_rot(0, Eigen::Vector3d::UnitZ());
    Eigen::Vector3d start_trans(0.0, 0.0, 0.0);
    Cosserat::SE3<double> start(Cosserat::SO3<double>(start_rot), start_trans);
    
    Eigen::AngleAxisd end_rot(M_PI/2, Eigen::Vector3d::UnitZ());
    Eigen::Vector3d end_trans(1.0, 2.0, 3.0);
    Cosserat::SE3<double> end(Cosserat::SO3<double>(end_rot), end_trans);
    
    // Interpolate using the exponential map (screw linear interpolation)
    Eigen::Matrix<double, 6, 1> start_tangent = start.log();
    Eigen::Matrix<double, 6, 1> end_tangent = end.log();
    Eigen::Matrix<double, 6, 1> mid_tangent = start_tangent + 0.5 * (end_tangent - start_tangent);
    Cosserat::SE3<double> mid = Cosserat::SE3<double>::exp(mid_tangent);
    
    // Create a sequence of interpolated poses
    std::vector<Cosserat::SE3<double>> path;
    int steps = 10;
    for (int i = 0; i <= steps; ++i) {
        double t = static_cast<double>(i) / steps;
        Eigen::Matrix<double, 6, 1> interp_tangent = 
            start_tangent + t * (end_tangent - start_tangent);
        path.push_back(Cosserat::SE3<double>::exp(interp_tangent));
    }
    
    return 0;
}
```

## Best Practices

1. **Choose the right representation** for your specific use case:
   - Use quaternions for the rotation component for better numerical stability
   - Store translation separately rather than using 4x4 matrices for better performance
   
2. **Optimize composition operations** when transforming many points:
   - Compose transformations first, then apply the result once
   - This is much more efficient than applying each transformation sequentially
   
3. **For interpolation between poses**:
   - Use screw linear interpolation (SLERP in the Lie algebra) for smooth and natural interpolation
   - Avoid linear interpolation of matrices or Euler angles
   
4. **For derivatives and velocities**:
   - Work in the tangent space (Lie algebra) where velocities have a natural representation
   - Convert back to the group only when needed
   
5. **Numerical stability**:
   - Normalize quaternions periodically to prevent numerical drift
   - For exponential map with small rotation, use specialized implementations
   
6. **Adjoint representation**:
   - Use the adjoint to transform velocities between different coordinate frames
   - This is more efficient than explicitly transforming velocity vectors
   
7. **Avoid repeated conversions**:
   - When possible, stay within a single representation
   - Convert only when necessary for specific operations

8. **Cache expensive operations**:
   - Matrix computations can be expensive, so cache matrices when they'll be reused
   - Consider caching log() results for frequently accessed transforms

9. **Remember that SE(3) is not commutative**:
   - The order of composition matters (A*B ≠ B*A)
   - Be careful about the order of transformations in your code

10. **Use the right tool for the job**:
    - For pure rotations, use SO(3) instead of SE(3)
    - For pure translations, use RealSpace instead of SE(3)
    - For similarity transforms (rotation + translation + scaling), use Sim(3)

## Numerical Stability Considerations

- **Quaternion normalization**: Periodically renormalize quaternions to maintain unit length
- **Small angle approximations**: For small rotations, use specialized implementations of exp/log
- **Near-singularity handling**: Add checks for divide-by-zero conditions in log() calculations
- **Composition of many transformations**: Be aware of error accumulation in long chains of transformations
- **Interpolation path selection**: Ensure interpolation takes the shortest path, especially for rotations
- **Exponential map implementation**: Use a robust implementation that handles all cases correctly
- **Conversion between representations**: Be careful about numerical issues at singularities
- **Double cover handling**: For interpolation, ensure quaternions are on the same hemisphere

