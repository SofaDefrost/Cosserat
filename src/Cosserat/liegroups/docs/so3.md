# SO(3) Implementation

## Overview

`SO3<Scalar>` implements the Special Orthogonal group in 3D, which represents rotations in 3D space. SO(3) is a 3-dimensional Lie group, where the dimension corresponds to the degrees of freedom for 3D rotations.

## Mathematical Foundations

### Group Structure

SO(3) is the group of all 3×3 orthogonal matrices with determinant 1. Mathematically:

SO(3) = {R ∈ ℝ³ˣ³ | R^T R = I, det(R) = 1}

Key properties:
- **Group operation**: Matrix multiplication (composition of rotations)
- **Identity element**: 3×3 identity matrix (no rotation)
- **Inverse element**: Transpose of the rotation matrix (R^T = R^(-1))
- **Not commutative**: In general, R₁R₂ ≠ R₂R₁ (order of rotations matters)
- **Compact**: All elements are bounded (closed and bounded subset of ℝ³ˣ³)
- **Connected**: Any rotation can be continuously deformed to any other

### Lie Algebra

The Lie algebra of SO(3), denoted so(3), consists of all 3×3 skew-symmetric matrices:

so(3) = {ω̂ ∈ ℝ³ˣ³ | ω̂^T = -ω̂}

A general element of so(3) has the form:

```
ω̂ = [  0  -ω₃  ω₂ ]
     [ ω₃   0  -ω₁ ]
     [-ω₂   ω₁   0  ]
```

where ω = [ω₁, ω₂, ω₃] ∈ ℝ³ is the corresponding angular velocity vector.

The "hat" operator (^) maps a 3D vector to a skew-symmetric matrix:
- ω̂ = hat(ω)

The inverse "vee" operator (∨) maps a skew-symmetric matrix to a 3D vector:
- ω = vee(ω̂)

### Exponential and Logarithmic Maps

#### Exponential Map: so(3) → SO(3)

The exponential map converts an element of the Lie algebra (skew-symmetric matrix or rotation vector) to a rotation matrix:

R = exp(ω̂) = I + sin(θ)/θ · ω̂ + (1-cos(θ))/θ² · ω̂²

where θ = ‖ω‖ is the magnitude of the rotation vector.

This is known as Rodrigues' formula. For small rotations, numerical approximations are used to avoid division by near-zero angles.

#### Logarithmic Map: SO(3) → so(3)

The logarithmic map converts a rotation matrix to an element of the Lie algebra:

ω̂ = log(R)

The rotation vector ω can be computed as:

```
θ = arccos((trace(R) - 1)/2)
ω = θ/(2sin(θ)) · [R₃₂-R₂₃, R₁₃-R₃₁, R₂₁-R₁₂]^T
```

Special care is needed for rotations near the identity (θ ≈ 0) and for 180° rotations.

### Quaternion Representation

Unit quaternions provide a compact and numerically stable representation of rotations.

A quaternion q = [w, x, y, z] represents the rotation:
- Around axis [x, y, z]/‖[x, y, z]‖
- By angle 2·arccos(w)

The relationship between a rotation vector ω = θ·n (where n is a unit vector) and a quaternion is:

```
q = [cos(θ/2), sin(θ/2)·n]
```

Key properties:
- The quaternion must have unit length (‖q‖ = 1)
- Quaternions q and -q represent the same rotation (double cover)
- Composition of rotations corresponds to quaternion multiplication

## Implementation Details

### Internal Representation

The `SO3` class supports several internal representations:

- **Quaternion representation** (primary): Most efficient for composition and inversion
- **Rotation matrix**: Used for acting on points and for interfacing with other systems
- **Axis-angle**: Used for certain operations and for intuitive construction

The implementation uses the quaternion as the primary storage format due to its numerical stability and efficiency.

### Key Methods

```cpp
// Constructors
SO3(); // Identity rotation
SO3(const Quaternion<Scalar>& quat); // From quaternion
SO3(const Matrix3<Scalar>& rot_matrix); // From rotation matrix
SO3(const AngleAxis<Scalar>& angle_axis); // From angle-axis

// Group operations
SO3 compose(const SO3& other) const;
SO3 inverse() const;

// Access to components
Quaternion<Scalar> quaternion() const;
Matrix3<Scalar> matrix() const;
AngleAxis<Scalar> angleAxis() const;

// Tangent space (Lie algebra) operations
Vector3<Scalar> log() const;
static SO3 exp(const Vector3<Scalar>& tangent);

// Acting on points
Vector3<Scalar> act(const Vector3<Scalar>& point) const;

// Adjoint representation
Matrix3<Scalar> adjoint() const;
```

### Performance Considerations

Operation performance with quaternion representation:

- **Composition**: O(1), 16 multiplications + 12 additions
- **Inverse**: O(1), constant time (quaternion conjugate)
- **Acting on points**: O(1), slightly slower than using matrices directly
- **Conversion to matrix**: O(1), but more expensive
- **Logarithmic map**: O(1), but involves transcendental functions
- **Exponential map**: O(1), transcendental functions + normalization

Compared to matrix representation, quaternions are:
- More memory efficient (4 vs. 9 elements)
- Faster for composition and inversion
- Slower for point transformation (unless batched)
- More numerically stable under repeated operations

### Numerical Stability

The implementation includes special handling for numerical stability:

- **Quaternion normalization**: Applied periodically to prevent drift
- **Small angle approximations**: For exp/log near the identity
- **Near-singular cases**: Special handling for 180° rotations
- **Double cover handling**: Ensuring consistent choice between q and -q

## Usage Examples

### Basic Operations

```cpp
#include <Cosserat/liegroups/SO3.h>
#include <iostream>

int main() {
    // Creating rotations in different ways
    
    // 1. From axis-angle representation
    Eigen::Vector3d axis(1.0, 1.0, 1.0);
    axis.normalize();
    double angle = M_PI/3;  // 60 degrees
    Cosserat::SO3<double> rotation1(Eigen::AngleAxisd(angle, axis));
    
    // 2. From rotation matrix
    Eigen::Matrix3d R;
    R = Eigen::AngleAxisd(M_PI/2, Eigen::Vector3d::UnitZ());
    Cosserat::SO3<double> rotation2(R);
    
    // 3. From quaternion
    Eigen::Quaterniond q;
    q = Eigen::AngleAxisd(M_PI/4, Eigen::Vector3d::UnitX());
    Cosserat::SO3<double> rotation3(q);
    
    // Compose rotations (apply rotation2 after rotation1)
    auto composed = rotation1.compose(rotation2);
    
    // Inverse rotation
    auto inverse = rotation1.inverse();
    
    // Convert to different representations
    Eigen::Matrix3d rot_mat = rotation1.matrix();
    Eigen::Quaterniond quat = rotation1.quaternion();
    Eigen::AngleAxisd aa = rotation1.angleAxis();
    
    // Accessing properties
    std::cout << "Rotation matrix:\n" << rot_mat << "\n";
    std::cout << "Quaternion: " << quat.coeffs().transpose() << "\n";
    std::cout << "Axis: " << aa.axis().transpose() << ", angle: " << aa.angle() << "\n";
    
    // Rotate a point
    Eigen::Vector3d point(1.0, 0.0, 0.0);
    Eigen::Vector3d rotated_point = rotation1.act(point);
    std::cout << "Original point: " << point.transpose() << "\n";
    std::cout << "Rotated point: " << rotated_point.transpose() << "\n";
    
    return 0;
}
```

### Interpolation

```cpp
#include <Cosserat/liegroups/SO3.h>
#include <vector>

// Spherical Linear Interpolation (SLERP)
Cosserat::SO3<double> slerp(
    const Cosserat::SO3<double>& start,
    const Cosserat::SO3<double>& end,
    double t
) {
    // Get quaternions
    Eigen::Quaterniond q_start = start.quaternion();
    Eigen::Quaterniond q_end = end.quaternion();
    
    // Ensure shortest path (handle double cover)
    if (q_start.dot(q_end) < 0) {
        q_end.coeffs() = -q_end.coeffs();
    }
    
    // Use Eigen's SLERP
    Eigen::Quaterniond q_interp = q_start.slerp(t, q_end);
    return Cosserat::SO3<double>(q_interp);
}

// Create a smooth rotation path
std::vector<Cosserat::SO3<double>> createRotationPath(
    const Cosserat::SO3<double>& start,
    const Cosserat::SO3<double>& end,
    int steps
) {
    std::vector<Cosserat::SO3<double>> path;
    for (int i = 0; i <= steps; ++i) {
        double t = static_cast<double>(i) / steps;
        path.push_back(slerp(start, end, t));
    }
    return path;
}
```

### Integration with Dynamics

```cpp
#include <Cosserat/liegroups/SO3.h>

// Angular velocity integration (discrete)
Cosserat::SO3<double> integrateRotation(
    const Cosserat::SO3<double>& R,
    const Eigen::Vector3d& omega,
    double dt
) {
    // Convert body angular velocity to a rotation increment
    Eigen::Vector3d delta_rot = omega * dt;
    
    // Apply the increment using the exponential map
    Cosserat::SO3<double> delta_R = Cosserat::SO3<double>::exp(delta_rot);
    
    // Apply to current rotation (left multiplication for body frame)
    return R.compose(delta_R);
}

// Rigid body dynamics example
void rigidBodyStep(
    Cosserat::SO3<double>& orientation,
    Eigen::Vector3d& angular_velocity,
    const Eigen::Matrix3d& inertia_tensor,
    const Eigen::Vector3d& torque,
    double dt
) {
    // Angular momentum
    Eigen::Vector3d angular_momentum = inertia_tensor * angular_velocity;
    
    // Update angular momentum with torque
    angular_momentum += torque * dt;
    
    // Update angular velocity
    angular_velocity = inertia_tensor.inverse() * angular_momentum;
    
    // Update orientation
    orientation = integrateRotation(orientation, angular_velocity, dt);
}
```

### Common Applications

```cpp
#include <Cosserat/liegroups/SO3.h>

// Camera orientation tracking
class Camera {
private:
    Cosserat::SO3<double> orientation;
    Eigen::Vector3d position;
    
public:
    // Update from IMU measurements
    void updateFromIMU(const Eigen::Vector3d& gyro, double dt) {
        orientation = integrateRotation(orientation, gyro, dt);
    }
    
    // Get view matrix for rendering
    Eigen::Matrix4d getViewMatrix() const {
        Eigen::Matrix4d view = Eigen::Matrix4d::Identity();
        
        // Set rotation part
        view.block<3,3>(0,0) = orientation.inverse().matrix();
        
        // Set translation part
        view.block<3,1>(0,3) = -orientation.inverse().act(position);
        
        return view;
    }
};

// Attitude control for spacecraft
class Spacecraft {
private:
    Cosserat::SO3<double> orientation;
    Eigen::Vector3d angular_velocity;
    Eigen::Matrix3d inertia;
    
public:
    // Compute control torque to reach target orientation
    Eigen::Vector3d computeControlTorque(const Cosserat::SO3<double>& target) {
        // Error in orientation (in the Lie algebra)
        Eigen::Vector3d error_vector = 
            (orientation.inverse().compose(target)).log();
        
        // PD controller
        double Kp = 1.0;  // Proportional gain
        double Kd = 0.5;  // Derivative gain
        
        return Kp * error_vector - Kd * angular_velocity;
    }
};
```

## Best Practices and Edge Cases

### Handling Singularities

1. **Logarithmic Map Singularities**:
   - For rotations near the identity, use Taylor series approximation
   - For 180° rotations, find a valid rotation axis by checking for non-zero elements in R - R^T

```cpp
// Robust logarithmic map implementation
Eigen::Vector3d robustLog(const Eigen::Matrix3d& R) {
    // Check if near identity
    double trace = R.trace();
    if (std::abs(trace - 3.0) < 1e-10) {
        // Near identity - use first-order approximation
        return Eigen::Vector3d(
            R(2,1) - R(1,2),
            R(0,2) - R(2,0),
            R(1,0) - R(0,1)
        ) * 0.5;
    }
    
    // Check if near 180° rotation
    if (std::abs(trace + 1.0) < 1e-10) {
        // Find axis by checking non-zero elements in R + I
        // [Implementation omitted for brevity]
    }
    
    // Standard case
    double theta = std::acos((trace - 1.0) / 2.0);
    return Eigen::Vector3d(
        R(2,1) - R(1,2),
        R(0,2) - R(2,0),
        R(1,0) - R(0,1)
    ) * (theta / (2.0 * std::sin(theta)));
}
```

2. **Exponential Map Stability**:
   - For small rotation angles, use Taylor series approximation
   - Handle zero-magnitude rotation vectors

```cpp
// Robust exponential map implementation
Eigen::Matrix3d robustExp(const Eigen::Vector3d& omega) {
    double theta = omega.norm();
    
    // Check if near-zero rotation
    if (theta < 1e-10) {
        // Use

# SO(3) Implementation

## Overview

`SO3<Scalar>` implements the Special Orthogonal group in 3D, which represents rotations in 3D space. SO(3) is a 3-dimensional Lie group, where the dimension corresponds to the degrees of freedom for 3D rotations.

## Mathematical Properties

- **Dimension**: 3 (rotations around the x, y, and z axes)
- **Group operation**: composition of rotations (matrix multiplication)
- **Identity element**: rotation by 0 radians (identity matrix)
- **Inverse**: transpose of the rotation matrix (inverse rotation)
- **Lie algebra**: so(3), which is the space of 3×3 skew-symmetric matrices
- **Exponential map**: converts from so(3) to SO(3)
- **Logarithmic map**: converts from SO(3) to so(3)

## Internal Representation

The `SO3` class can be internally represented in several ways:
- Rotation matrix (3×3 orthogonal matrix with determinant 1)
- Unit quaternion (more compact and numerically stable)
- Angle-axis representation (for certain operations)

Our implementation primarily uses unit quaternions for storage, with conversion methods to other representations as needed.

## Implementation Details

The `SO3` class is implemented as a template with one parameter:
- `Scalar`: The scalar type (typically `double` or `float`)

### Key Methods

```cpp
// Constructors
SO3(); // Identity rotation
SO3(const Quaternion<Scalar>& quat); // From quaternion
SO3(const Matrix3<Scalar>& rot_matrix); // From rotation matrix
SO3(const AngleAxis<Scalar>& angle_axis); // From angle-axis

// Group operations
SO3 compose(const SO3& other) const;
SO3 inverse() const;

// Access to components
Quaternion<Scalar> quaternion() const;
Matrix3<Scalar> matrix() const;
AngleAxis<Scalar> angleAxis() const;

// Tangent space (Lie algebra) operations
Vector3<Scalar> log() const;
static SO3 exp(const Vector3<Scalar>& tangent);

// Acting on points
Vector3<Scalar> act(const Vector3<Scalar>& point) const;

// Adjoint representation
Matrix3<Scalar> adjoint() const;
```

## Mathematical Background

### Rotation Matrix

A 3×3 rotation matrix R must satisfy:
- Orthogonality: R^T R = I (identity matrix)
- Proper rotation: det(R) = 1

### Lie Algebra

The Lie algebra so(3) consists of 3×3 skew-symmetric matrices of the form:
```
S = [  0  -w3   w2 ]
    [  w3   0  -w1 ]
    [ -w2   w1   0 ]
```

where [w1, w2, w3] is the axis-angle representation scaled by the angle. This vector can be thought of as the angular velocity vector.

### Exponential Map

The exponential map from so(3) to SO(3) can be computed using Rodrigues' formula:

For a rotation vector ω = θ·a (where a is a unit vector and θ is the angle):
```
exp(ω) = I + sin(θ)/θ · [ω]× + (1-cos(θ))/θ² · [ω]×²
```

where [ω]× is the skew-symmetric matrix formed from ω.

For small rotations, numerical approximations are used to avoid division by near-zero angles.

### Logarithmic Map

The logarithmic map from SO(3) to so(3) extracts the rotation vector from a rotation matrix or quaternion:

From a rotation matrix R with trace tr(R):
```
θ = acos((tr(R) - 1)/2)
ω = θ/(2*sin(θ)) · [R32-R23, R13-R31, R21-R12]^T
```

Again, special handling is required for small angles to ensure numerical stability.

## Performance Characteristics

Based on benchmarks, SO(3) operations have the following performance characteristics:

- **Composition**: O(1) time complexity - involves quaternion multiplication
- **Inverse**: O(1) time complexity - involves quaternion conjugation
- **Matrix conversion**: O(1) time complexity - creates a 3×3 rotation matrix
- **Exponential map**: O(1) time complexity - converts from axis-angle to quaternion
- **Logarithmic map**: O(1) time complexity - converts from quaternion to axis-angle
- **Acting on points**: O(1) time complexity - applies the rotation to a 3D point

For most operations, the quaternion representation offers better performance and numerical stability compared to rotation matrices, especially for composition and inversion.

## Example Usage

### Basic Operations

```cpp
#include <Cosserat/liegroups/SO3.h>
#include <iostream>

int main() {
    // Create an SO(3) element (90-degree rotation around z-axis)
    Eigen::AngleAxisd angle_axis(M_PI/2, Eigen::Vector3d::UnitZ());
    Cosserat::SO3<double> rotation(angle_axis);
    
    // Get representations
    Eigen::Matrix3d rot_mat = rotation.matrix();
    Eigen::Quaterniond quat = rotation.quaternion();
    
    std::cout << "Rotation matrix:\n" << rot_mat << "\n";
    std::cout << "Quaternion: " << quat.coeffs().transpose() << "\n";
    
    // Create another rotation (45-degree rotation around x-axis)
    Eigen::AngleAxisd another_angle_axis(M_PI/4, Eigen::Vector3d::UnitX());
    Cosserat::SO3<double> another_rotation(another_angle_axis);
    
    // Compose rotations
    auto composed = rotation.compose(another_rotation);
    
    // Inverse rotation
    auto inverse = rotation.inverse();
    
    return 0;
}
```

### Tangent Space Operations

```cpp
#include <Cosserat/liegroups/SO3.h>
#include <iostream>

int main() {
    // Create a rotation from axis-angle
    Eigen::Vector3d axis(1.0, 1.0, 1.0);
    axis.normalize();
    double angle = M_PI/3;
    Eigen::AngleAxisd angle_axis(angle, axis);
    Cosserat::SO3<double> rotation(angle_axis);
    
    // Convert to Lie algebra (tangent space)
    Eigen::Vector3d tangent = rotation.log();
    std::cout << "Tangent vector: " << tangent.transpose() << "\n";
    
    // Convert back from Lie algebra to SO(3)
    auto recovered = Cosserat::SO3<double>::exp(tangent);
    
    // Create directly from tangent vector
    Eigen::Vector3d new_tangent(0.1, 0.2, 0.3); // Small rotation
    auto new_rotation = Cosserat::SO3<double>::exp(new_tangent);
    
    return 0;
}
```

### Acting on Points

```cpp
#include <Cosserat/liegroups/SO3.h>
#include <iostream>
#include <vector>

int main() {
    // Create a rotation
    Eigen::AngleAxisd angle_axis(M_PI/2, Eigen::Vector3d::UnitZ());
    Cosserat::SO3<double> rotation(angle_axis);
    
    // Rotate a single point
    Eigen::Vector3d point(1.0, 0.0, 0.0);
    Eigen::Vector3d rotated_point = rotation.act(point);
    std::cout << "Original point: " << point.transpose() << "\n";
    std::cout << "Rotated point: " << rotated_point.transpose() << "\n";
    
    // Rotate multiple points
    std::vector<Eigen::Vector3d> points = {
        Eigen::Vector3d(1.0, 0.0, 0.0),
        Eigen::Vector3d(0.0, 1.0, 0.0),
        Eigen::Vector3d(0.0, 0.0, 1.0)
    };
    
    std::vector<Eigen::Vector3d> rotated_points;
    for (const auto& p : points) {
        rotated_points.push_back(rotation.act(p));
    }
    
    return 0;
}
```

### Interpolation

```cpp
#include <Cosserat/liegroups/SO3.h>
#include <iostream>

int main() {
    // Create two rotations
    Cosserat::SO3<double> start(Eigen::AngleAxisd(0, Eigen::Vector3d::UnitZ()));
    Cosserat::SO3<double> end(Eigen::AngleAxisd(M_PI/2, Eigen::Vector3d::UnitZ()));
    
    // Interpolate between them (spherical linear interpolation)
    Cosserat::SO3<double> mid = Cosserat::slerp(start, end, 0.5);
    
    // Interpolate using the exponential map
    Eigen::Vector3d start_tangent = start.log();
    Eigen::Vector3d end_tangent = end.log();
    Eigen::Vector3d mid_tangent = start_tangent + 0.5 * (end_tangent - start_tangent);
    Cosserat::SO3<double> mid2 = Cosserat::SO3<double>::exp(mid_tangent);
    
    return 0;
}
```

## Best Practices

1. **Use quaternion representation** for most operations due to better numerical stability and performance.
2. **Normalize quaternions regularly** to prevent numerical drift.
3. **When performing many rotations**, compose them first and then apply the result, rather than applying each rotation individually.
4. **For interpolation between rotations**, use spherical linear interpolation (SLERP) rather than linear interpolation of matrices or Euler angles.
5. **Be aware of the double cover** - two quaternions can represent the same rotation (q and -q).
6. **For small rotations**, the exponential map can be approximated, but use the full implementation for general cases.
7. **When converting between representations**, be aware of numerical issues at singularities (e.g., Euler angles gimbal lock).
8. **For time-varying rotations**, work in the tangent space for derivatives.
9. **Use the appropriate constructor** based on your input data to avoid unnecessary conversions.
10. **Remember that rotations are not commutative** - the order of composition matters.

## Numerical Stability Considerations

- **Quaternion normalization**: Ensure quaternions stay unit length by normalizing periodically
- **Logarithm near identity**: Use specialized implementations for small rotations
- **Avoiding gimbal lock**: Work with quaternions or axis-angle rather than Euler angles
- **Double cover handling**: Ensure that interpolation takes the shortest path
- **Composition of many rotations**: Reorthogonalize occasionally to prevent drift

# SO(3) Implementation

## Overview

`SO3<Scalar>` implements the Special Orthogonal group in 3D, which represents rotations in 3D space. SO(3) is a 3-dimensional Lie group, where the dimension corresponds to the degrees of freedom for 3D rotations.

## Mathematical Properties

- **Dimension**: 3 (rotations around the x, y, and z axes)
- **Group operation**: composition of rotations (matrix multiplication)
- **Identity element**: rotation by 0 radians (identity matrix)
- **Inverse**: transpose of the rotation matrix (inverse rotation)
- **Lie algebra**: so(3), which is the space of 3×3 skew-symmetric matrices
- **Exponential map**: converts from so(3) to SO(3)
- **Logarithmic map**: converts from SO(3) to so(3)

## Internal Representation

The `SO3` class can be internally represented in several ways:
- Rotation matrix (3×3 orthogonal matrix with determinant 1)
- Unit quaternion (more compact and numerically stable)
- Angle-axis representation (for certain operations)

Our implementation primarily uses unit quaternions for storage, with conversion methods to other representations as needed.

## Implementation Details

The `SO3` class is implemented as a template with one parameter:
- `Scalar`: The scalar type (typically `double` or `float`)

### Key Methods

```cpp
// Constructors
SO3(); // Identity rotation
SO3(const Quaternion<Scalar>& quat); // From quaternion
SO3(const Matrix3<Scalar>& rot_matrix); // From rotation matrix
SO3(const AngleAxis<Scalar>& angle_axis); // From angle-axis

// Group operations
SO3 compose(const SO3& other) const;
SO3 inverse() const;

// Access to components
Quaternion<Scalar> quaternion() const;
Matrix3<Scalar> matrix() const;
AngleAxis<Scalar> angleAxis() const;

// Tangent space (Lie algebra) operations
Vector3<Scalar> log() const;
static SO3 exp(const Vector3<Scalar>& tangent);

// Acting on points
Vector3<Scalar> act(const Vector3<Scalar>& point) const;

// Adjoint representation
Matrix3<Scalar> adjoint() const;
```

## Performance Characteristics

Based on benchmarks, SO(3) operations have the following performance characteristics:

- **Composition**: O(1) time complexity - involves quaternion multiplication
- **Inverse**: O(1) time complexity - involves quaternion conjugation
- **Matrix conversion**: O(1) time complexity - creates a 3×3 rotation matrix
- **Exponential map**: O(1) time complexity - converts from axis-angle to quaternion
- **Logarithmic map**: O(1) time complexity - converts from quaternion to axis-angle
- **Acting on points**: O(1) time complexity - applies the rotation to a 3D point

For most operations, the quaternion representation offers better performance and numerical stability compared to rotation matrices, especially for composition and inversion.

## Example Usage

### Basic Operations

```cpp
#include <Cosserat/liegroups/SO3.h>
#include <iostream>

int main() {
    // Create an SO(3) element (90-degree rotation around z-axis)
    Eigen::AngleAxisd angle_axis(M_PI/2, Eigen::Vector3d::UnitZ());
    Cosserat::SO3<double> rotation(angle_axis);
    
    // Get representations
    Eigen::Matrix3d rot_mat = rotation.matrix();
    Eigen::Quaterniond quat = rotation.quaternion();
    
    std::cout << "Rotation matrix:\n" << rot_mat << "\n";
    std::cout << "Quaternion: " << quat.coeffs().transpose() << "\n";
    
    // Create another rotation (45-degree rotation around x-axis)
    Eigen::AngleAxisd another_angle_axis(M_PI/4, Eigen::Vector3d::UnitX());
    Cosserat::SO3<double> another_rotation(another_angle_axis);
    
    // Compose rotations
    auto composed = rotation.compose(another_rotation);
    
    // Inverse rotation
    auto inverse = rotation.inverse();
    
    return 0;
}
```

### Tangent Space Operations

```cpp
#include <Cosserat/liegroups/SO3.h>
#include <iostream>

int main() {
    // Create a rotation from axis-angle
    Eigen::Vector3d axis(1.0, 1.0, 1.0);
    axis.normalize();
    double angle = M_PI/3;
    Eigen::AngleAxisd angle_axis(angle, axis);
    Cosserat::SO3<double> rotation(angle_axis);
    
    // Convert to Lie algebra (tangent space)
    Eigen::Vector3d tangent = rotation.log();
    std::cout << "Tangent vector: " << tangent.transpose() << "\n";
    
    // Convert back from Lie algebra to SO(3)
    auto recovered = Cosserat::SO3<double>::exp(tangent);
    
    // Create directly from tangent vector
    Eigen::Vector3d new_tangent(0.1, 0.2, 0.3); // Small rotation
    auto new_rotation = Cosserat::SO3<double>::exp(new_tangent);
    
    return 0;
}
```

### Acting on Points

```cpp
#include <Cosserat/liegroups/SO3.h>
#include <iostream>
#include <vector>

int main() {
    // Create a rotation
    Eigen::AngleAxisd angle_axis(M_PI/2, Eigen::Vector3d::UnitZ());
    Cosserat::SO3<double> rotation(angle_axis);
    
    // Rotate a single point
    Eigen::Vector3d point(1.0, 0.0, 0.0);
    Eigen::Vector3d rotated_point = rotation.act(point);
    std::cout << "Original point: " << point.transpose() << "\n";
    std::cout << "Rotated point: " << rotated_point.transpose() << "\n";
    
    // Rotate multiple points
    std::vector<Eigen::Vector3d> points = {
        Eigen::Vector3d(1.0, 0.0, 0.0),
        Eigen::Vector3d(0.0, 1.0, 0.0),
        Eigen::Vector3d(0.0, 0.0, 1.0)
    };
    
    std::vector<Eigen::Vector3d> rotated_points;
    for (const auto& p : points) {
        rotated_points.push_back(rotation.act(p));
    }
    
    return 0;
}
```

### Interpolation

```cpp
#include <Cosserat/liegroups/SO3.h>
#include <iostream>

int main() {
    // Create two rotations
    Cosserat::SO3<double> start(Eigen::AngleAxisd(0, Eigen::Vector3d::UnitZ()));
    Cosserat::SO3<double> end(Eigen::AngleAxisd(M_PI/2, Eigen::Vector3d::UnitZ()));
    
    // Interpolate between them (spherical linear interpolation)
    Cosserat::SO3<double> mid = Cosserat::slerp(start, end, 0.5);
    
    // Interpolate using the exponential map
    Eigen::Vector3d start_tangent = start.log();
    Eigen::Vector3d end_tangent = end.log();
    Eigen::Vector3d mid_tangent = start_tangent + 0.5 * (end_tangent - start_tangent);
    Cosserat::SO3<double> mid2 = Cosserat::SO3<double>::exp(mid_tangent);
    
    return 0;
}
```

## Best Practices

1. **Use quaternion representation** for most operations due to better numerical stability and performance.
2. **Normalize quaternions regularly** to prevent numerical drift.
3. **When performing many rotations**, compose them first and then apply the result, rather than applying each rotation individually.
4. **For interpolation between rotations**, use spherical linear interpolation (SLERP) rather than linear interpolation of matrices or Euler angles.
5. **Be aware of the double cover** - two quaternions can represent the same rotation (q and -q).
6. **For small rotations**, the exponential map can be approximated, but use the full implementation for general cases.
7. **When converting between representations**, be aware of numerical issues at singularities (e.g., Euler angles gimbal lock).
8. **For time-varying rotations**, work in the tangent space for derivatives.
9. **Use the appropriate constructor** based on your input data to avoid unnecessary conversions.
10. **Remember that rotations are not commutative** - the order of composition matters.

## Numerical Stability Considerations

- **Quaternion normalization**: Ensure quaternions stay unit length by normalizing periodically
- **Logarithm near identity**: Use specialized implementations for small rotations
- **Avoiding gimbal lock**: Work with quaternions or axis-angle rather than Euler angles
- **Double cover handling**: Ensure that interpolation takes the shortest path
- **Composition of many rotations**: Reorthogonalize occasionally to prevent drift

