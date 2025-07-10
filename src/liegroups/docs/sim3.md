# Sim(3) Implementation

## Overview

`Sim3<Scalar>` implements the Similarity transformation group in 3D, which represents rotation, translation, and uniform scaling in 3D space. Sim(3) is a 7-dimensional Lie group, combining a 3-dimensional rotation (SO(3)), a 3-dimensional translation, and a 1-dimensional scaling factor.

## Mathematical Properties

- **Dimension**: 7 (3 for rotation + 3 for translation + 1 for scaling)
- **Group operation**: composition of similarity transformations
- **Identity element**: zero rotation, zero translation, and unit scaling
- **Inverse**: inverse rotation, scaled and rotated negative translation, and reciprocal scaling
- **Lie algebra**: sim(3), which can be represented as a 7D vector [ω, v, s]
  where ω is the rotation part, v is the translation part, and s is the scaling part
- **Exponential map**: converts from sim(3) to Sim(3)
- **Logarithmic map**: converts from Sim(3) to sim(3)

## Internal Representation

The `Sim3` class internally stores:
- A rotation component as an SO(3) element (typically as a quaternion)
- A translation component as a 3D vector
- A scaling factor as a scalar

This representation ensures proper handling of the group structure and efficient computation of group operations.

## Implementation Details

The `Sim3` class is implemented as a template with one parameter:
- `Scalar`: The scalar type (typically `double` or `float`)

### Key Methods

```cpp
// Constructors
Sim3(); // Identity transformation
Sim3(const SO3<Scalar>& rotation, const Vector3<Scalar>& translation, Scalar scale);
Sim3(const Matrix4<Scalar>& similarity_matrix);
Sim3(const Quaternion<Scalar>& quat, const Vector3<Scalar>& translation, Scalar scale);

// Group operations
Sim3 compose(const Sim3& other) const;
Sim3 inverse() const;

// Access to components
SO3<Scalar> rotation() const;
Vector3<Scalar> translation() const;
Scalar scale() const;
Matrix4<Scalar> matrix() const; // Homogeneous similarity matrix

// Tangent space (Lie algebra) operations
Vector7<Scalar> log() const;
static Sim3 exp(const Vector7<Scalar>& tangent);

// Acting on points
Vector3<Scalar> act(const Vector3<Scalar>& point) const;

// Adjoint representation
Matrix7<Scalar> adjoint() const;
```

## Performance Characteristics

Based on benchmarks, Sim(3) operations have the following performance characteristics:

- **Composition**: O(1) time complexity - involves rotation composition, scaled translation transformation, and scale multiplication
- **Inverse**: O(1) time complexity - computes inverse rotation, scaled inverse translation, and reciprocal scale
- **Matrix conversion**: O(1) time complexity - creates a 4×4 similarity transformation matrix
- **Exponential map**: O(1) time complexity - converts from Lie algebra to the group
- **Logarithmic map**: O(1) time complexity - converts from the group to Lie algebra
- **Acting on points**: O(1) time complexity - applies the transformation to a 3D point (rotation, scaling, and translation)

The actual performance depends on the hardware and compiler optimizations. Similar to SE(3), using quaternions for the rotational component provides better performance for most operations.

## Example Usage

### Basic Operations

```cpp
#include <Cosserat/liegroups/Sim3.h>
#include <iostream>

int main() {
    // Create a Sim(3) element (90-degree rotation around z-axis, translation (1,2,3), scale 2.0)
    Eigen::AngleAxisd rotation(M_PI/2, Eigen::Vector3d::UnitZ());
    Eigen::Vector3d translation(1.0, 2.0, 3.0);
    double scale = 2.0;
    Cosserat::Sim3<double> transform(
        Cosserat::SO3<double>(rotation), 
        translation,
        scale
    );
    
    // Get components
    Cosserat::SO3<double> rot = transform.rotation();
    Eigen::Vector3d trans = transform.translation();
    double s = transform.scale();
    std::cout << "Rotation matrix:\n" << rot.matrix() << "\n";
    std::cout << "Translation: " << trans.transpose() << "\n";
    std::cout << "Scale: " << s << "\n";
    
    // Convert to homogeneous similarity matrix
    Eigen::Matrix4d mat = transform.matrix();
    std::cout << "Similarity matrix:\n" << mat << "\n";
    
    // Create another transformation
    Cosserat::Sim3<double> another_transform(
        Cosserat::SO3<double>(Eigen::AngleAxisd(M_PI/4, Eigen::Vector3d::UnitX())),
        Eigen::Vector3d(4.0, 5.0, 6.0),
        0.5
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
#include <Cosserat/liegroups/Sim3.h>
#include <iostream>

int main() {
    // Create a Sim(3) element
    Eigen::AngleAxisd rotation(M_PI/3, Eigen::Vector3d::UnitZ());
    Eigen::Vector3d translation(1.0, 2.0, 3.0);
    double scale = 1.5;
    Cosserat::Sim3<double> transform(
        Cosserat::SO3<double>(rotation), 
        translation,
        scale
    );
    
    // Convert to Lie algebra (tangent space)
    Eigen::Matrix<double, 7, 1> tangent = transform.log();
    std::cout << "Tangent vector: " << tangent.transpose() << "\n";
    
    // Convert back from Lie algebra to Sim(3)
    auto recovered = Cosserat::Sim3<double>::exp(tangent);
    
    // Create directly from tangent vector
    Eigen::Matrix<double, 7, 1> new_tangent;
    new_tangent << 0.1, 0.2, 0.3, 1.0, 2.0, 3.0, 0.1; // Rotation, translation, log(scale)
    auto new_transform = Cosserat::Sim3<double>::exp(new_tangent);
    
    return 0;
}
```

