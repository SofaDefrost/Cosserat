# SE(2) Implementation

## Overview

`SE2<Scalar>` implements the Special Euclidean group in 2D, which represents rigid body transformations (rotation and translation) in a plane. SE(2) is a 3-dimensional Lie group, combining a 1-dimensional rotation (SO(2)) with a 2-dimensional translation.

## Mathematical Properties

- **Dimension**: 3 (1 for rotation + 2 for translation)
- **Group operation**: composition of transformations
- **Identity element**: zero rotation and zero translation
- **Inverse**: inverse rotation and negated rotated translation
- **Lie algebra**: se(2), which can be represented as a 3D vector [θ, x, y]
- **Exponential map**: converts from se(2) to SE(2)
- **Logarithmic map**: converts from SE(2) to se(2)

## Internal Representation

The `SE2` class internally stores:
- A rotation component as an `SO(2)` element
- A translation component as a 2D vector

This representation ensures proper handling of the group structure and efficient computation of group operations.

## Implementation Details

The `SE2` class is implemented as a template with one parameter:
- `Scalar`: The scalar type (typically `double` or `float`)

### Key Methods

```cpp
// Constructors
SE2(); // Identity transformation
SE2(Scalar angle, Scalar x, Scalar y); // From angle and translation
SE2(const SO2<Scalar>& rotation, const Eigen::Matrix<Scalar, 2, 1>& translation);
SE2(const Eigen::Matrix<Scalar, 3, 3>& homogeneous_matrix);

// Group operations
SE2 compose(const SE2& other) const;
SE2 inverse() const;

// Access to components
SO2<Scalar> rotation() const;
Eigen::Matrix<Scalar, 2, 1> translation() const;
Scalar angle() const;

// Conversion to matrix representation
Eigen::Matrix<Scalar, 3, 3> matrix() const; // Homogeneous transformation matrix

// Tangent space (Lie algebra) operations
Eigen::Matrix<Scalar, 3, 1> log() const;
static SE2 exp(const Eigen::Matrix<Scalar, 3, 1>& tangent);

// Acting on points
Eigen::Matrix<Scalar, 2, 1> act(const Eigen::Matrix<Scalar, 2, 1>& point) const;
```

## Performance Characteristics

Based on our benchmarks, SE(2) operations have the following performance characteristics:

- **Composition**: O(1) time complexity - involves rotation composition and translation addition
- **Inverse**: O(1) time complexity - computes inverse rotation and negated rotated translation
- **Matrix conversion**: O(1) time complexity - creates a 3×3 homogeneous transformation matrix
- **Exponential map**: O(1) time complexity - converts from Lie algebra to the group
- **Logarithmic map**: O(1) time complexity - converts from the group to Lie algebra
- **Acting on points**: O(1) time complexity - applies the transformation to a 2D point

The actual performance depends on the hardware and compiler optimizations, but these operations are typically very fast due to their low dimensionality.

## Example Usage

### Basic Operations

```cpp
#include <Cosserat/liegroups/SE2.h>
#include <iostream>

int main() {
    // Create an SE(2) element (45-degree rotation with translation (1,2))
    Cosserat::SE2<double> transform(M_PI/4.0, 1.0, 2.0);
    
    // Get components
    double angle = transform.angle();
    Eigen::Vector2d translation = transform.translation();
    std::cout << "Angle: " << angle << " radians\n";
    std::cout << "Translation: [" << translation.transpose() << "]\n";
    
    // Convert to homogeneous transformation matrix
    Eigen::Matrix3d mat = transform.matrix();
    std::cout << "Matrix:\n" << mat << "\n";
    
    // Create another transformation
    Cosserat::SE2<double> another_transform(M_PI/2.0, 3.0, 4.0);
    
    // Compose transformations
    auto composed = transform.compose(another_transform);
    std::cout << "Composed angle: " << composed.angle() << " radians\n";
    std::cout << "Composed translation: [" << composed.translation().transpose() << "]\n";
    
    // Inverse transformation
    auto inverse = transform.inverse();
    std::cout << "Inverse angle: " << inverse.angle() << " radians\n";
    std::cout << "Inverse translation: [" << inverse.translation().transpose() << "]\n";
    
    return 0;
}
```

### Tangent Space Operations

```cpp
#include <Cosserat/liegroups/SE2.h>
#include <iostream>

int main() {
    // Create an SE(2) element
    Cosserat::SE2<double> transform(M_PI/4.0, 1.0, 2.0);
    
    // Convert to Lie algebra (tangent space)
    Eigen::Vector3d tangent = transform.log();
    std::cout << "Tangent vector: [" << tangent.transpose() << "]\n";
    
    // Convert back from Lie algebra to SE(2)
    auto recovered = Cosserat::SE2<double>::exp(tangent);
    std::cout << "Recovered angle: " << recovered.angle() << " radians\n";
    std::cout << "Recovered translation: [" << recovered.translation().transpose() << "]\n";
    
    // Create directly from tangent vector
    Eigen::Vector3d new_tangent;
    new_tangent << M_PI/6.0, 3.0, 4.0;
    auto new_transform = Cosserat::SE2<double>::exp(new_tangent);
    
    return 0;
}
```

### Acting on Points

```cpp
#include <Cosserat/liegroups/SE2.h>
#include <iostream>
#include <vector>

int main() {
    // Create an SE(2) element
    Cosserat::SE2<double> transform(M_PI/4.0, 1.0, 2.0);
    
    // Transform a single point
    Eigen::Vector2d point(3.0, 4.0);
    Eigen::Vector2d transformed_point = transform.act(point);
    std::cout << "Original point: [" << point.transpose() << "]\n";
    std::cout << "Transformed point: [" << transformed_point.transpose() << "]\n";
    
    // Transform multiple points
    std::vector<Eigen::Vector2d> points = {
        Eigen::Vector2d(1.0, 0.0),
        Eigen::Vector2d(0.0, 1.0),
        Eigen::Vector2d(1.0, 1.0)
    };
    
    std::vector<Eigen::Vector2d> transformed_points;
    for (const auto& p : points) {
        transformed_points.push_back(transform.act(p));
    }
    
    // Check that applying the transformation twice is equivalent to composing
    Eigen::Vector2d twice_transformed = transform.act(transform.act(point));
    Eigen::Vector2d composed_transform = transform.compose(transform).act(point);
    
    return 0;
}
```

## Best Practices

1. **Use the most appropriate constructor** for your use case to avoid unnecessary conversions.
2. **Avoid repeated matrix construction** when performing multiple operations with the same transformation.
3. **Use the compose() method rather than manual matrix multiplication** for better semantics and potentially better performance.
4. **When transforming many points**, extract the rotation matrix and translation vector once outside the loop for better performance.
5. **For interpolation between poses**, convert to the Lie algebra (log), perform the interpolation there, and then convert back (exp).
6. **When working with velocities**, use the Lie algebra representation which directly corresponds to angular and linear velocities.
7. **For small displacements**, the exponential map can be approximated, but use the full implementation for general cases.
8. **Remember that SE(2) is not commutative** - the order of composition matters.

