# Comparison of Lie Group Implementations

This document compares the various Lie group implementations in the Cosserat plugin, highlighting their features, performance characteristics, and appropriate use cases.

## Feature Comparison

| Feature | RealSpace | SO(2) | SO(3) | SE(2) | SE(3) | Sim(3) |
|---------|-----------|-------|-------|-------|-------|--------|
| **Dimension** | n (templated) | 1 | 3 | 3 | 6 | 7 |
| **Represents** | Vectors | 2D rotation | 3D rotation | 2D rigid transform | 3D rigid transform | 3D similarity transform |
| **Group operation** | Addition | Rotation composition | Rotation composition | Rigid motion composition | Rigid motion composition | Similarity composition |
| **Internal representation** | Vector | Angle or complex | Quaternion | Angle + vector | Quaternion + vector | Quaternion + vector + scale |
| **Has rotation component** | No | Yes | Yes | Yes | Yes | Yes |
| **Has translation component** | Yes (represents position) | No | No | Yes | Yes | Yes |
| **Has scale component** | No | No | No | No | No | Yes |
| **Commutative** | Yes | Yes | No | No | No | No |
| **Primary application** | Points, vectors | 2D rotations | 3D rotations | 2D mechanics | 3D mechanics | Computer vision |

## Lie Algebra Properties

| Property | RealSpace | SO(2) | SO(3) | SE(2) | SE(3) | Sim(3) |
|----------|-----------|-------|-------|-------|-------|--------|
| **Algebra dimension** | n | 1 | 3 | 3 | 6 | 7 |
| **Represents** | Vectors | Angular velocity | Angular velocity | Twist (ang.+lin. vel) | Twist | Twist + scaling |
| **Exponential map complexity** | Trivial | Simple | Medium | Medium | Complex | Complex |
| **Logarithmic map complexity** | Trivial | Simple | Medium | Medium | Complex | Complex |
| **Primary application** | Velocity | Angular velocity | Angular velocity | 2D body velocity | 3D body velocity | Scale-velocity |

## Performance Comparison

The following table shows approximate relative performance for common operations (normalized to the fastest implementation, lower is better):

| Operation | RealSpace | SO(2) | SO(3) | SE(2) | SE(3) | Sim(3) |
|-----------|-----------|-------|-------|-------|-------|--------|
| **Composition** | 1.0 | 1.2 | 2.5 | 3.0 | 5.0 | 5.5 |
| **Inverse** | 1.0 | 1.1 | 2.0 | 2.2 | 3.5 | 4.0 |
| **Log** | 1.0 | 2.0 | 4.0 | 5.0 | 10.0 | 12.0 |
| **Exp** | 1.0 | 1.5 | 3.5 | 4.5 | 9.0 | 11.0 |
| **Acting on point** | 1.0 | 1.2 | 1.8 | 2.0 | 2.2 | 2.5 |
| **Memory footprint** | n | 1 | 4 | 3 | 7 | 8 |

Note: These numbers are approximate and can vary based on hardware, compiler optimizations, and the specific data being processed.

## Use Case Recommendations

### When to use RealSpace

- When working with simple vectors or points
- For displacements in Euclidean space
- When performance is critical
- When the dynamics do not involve rotations

Example: Position vectors, linear velocities, forces

```cpp
// Representing a 3D position
Cosserat::RealSpace<double, 3> position(1.0, 2.0, 3.0);

// Representing a velocity vector
Cosserat::RealSpace<double, 3> velocity(-0.5, 0.0, 2.0);

// Simple vector addition (displacement)
auto new_position = position.compose(velocity);

// Direct access to vector components
double x = position[0];
double y = position[1];
double z = position[2];
```

### When to use SO(2)

- For 2D rotations
- When working with planar mechanisms
- For representing orientations in 2D space
- For angular velocities in a plane

Example: 2D robot orientation, pendulum angle, compass heading

```cpp
// Representing a 45-degree rotation
Cosserat::SO2<double> rotation(M_PI/4);

// Rotating a 2D point
Eigen::Vector2d point(1.0, 0.0);
Eigen::Vector2d rotated_point = rotation.act(point);

// Sequential rotations
Cosserat::SO2<double> first_rotation(M_PI/6);   // 30 degrees
Cosserat::SO2<double> second_rotation(M_PI/3);  // 60 degrees
Cosserat::SO2<double> combined = first_rotation.compose(second_rotation);
```

### When to use SO(3)

- For 3D rotations
- When representing orientations in 3D space
- For spacecraft attitude control
- For camera orientation

Example: Robot joint orientation, IMU attitude, camera rotation

```cpp
// Creating a rotation around an arbitrary axis
Eigen::Vector3d axis(1.0, 1.0, 1.0);
axis.normalize();
double angle = M_PI/3;  // 60 degrees
Cosserat::SO3<double> rotation(Eigen::AngleAxisd(angle, axis));

// Rotating a 3D point
Eigen::Vector3d point(1.0, 0.0, 0.0);
Eigen::Vector3d rotated_point = rotation.act(point);

// Converting to different representations
Eigen::Matrix3d rot_matrix = rotation.matrix();
Eigen::Quaterniond quaternion = rotation.quaternion();
```

### When to use SE(2)

- For 2D rigid body transformations
- For planar robot kinematics
- For 2D SLAM (Simultaneous Localization and Mapping)
- For 2D path planning

Example: Planar robot pose, 2D object transformation

```cpp
// Creating a 2D pose (90-degree rotation + translation (1,2))
Cosserat::SE2<double> pose(M_PI/2, 1.0, 2.0);

// Transforming a 2D point
Eigen::Vector2d point(3.0, 4.0);
Eigen::Vector2d transformed_point = pose.act(point);

// Composing transformations (e.g., for following a path)
Cosserat::SE2<double> step1(M_PI/4, 1.0, 0.0);  // 45-degree turn, move 1 unit forward
Cosserat::SE2<double> step2(0.0, 2.0, 0.0);     // Move 2 units forward
Cosserat::SE2<double> combined_path = step1.compose(step2);
```

### When to use SE(3)

- For 3D rigid body transformations
- For robot kinematics and dynamics
- For 3D SLAM and computer vision
- For physics simulations

Example: Robot pose, camera transformation, articulated body

```cpp
// Creating a 3D pose (rotation around Z + translation)
Eigen::AngleAxisd rotation(M_PI/2, Eigen::Vector3d::UnitZ());
Eigen::Vector3d translation(1.0, 2.0, 3.0);
Cosserat::SE3<double> pose(Cosserat::SO3<double>(rotation), translation);

// Transforming a 3D point
Eigen::Vector3d point(1.0, 0.0, 0.0);
Eigen::Vector3d transformed_point = pose.act(point);

// Robot forward kinematics example (simplified)
Cosserat::SE3<double> base_to_joint1 = getJoint1Transform();
Cosserat::SE3<double> joint1_to_joint2 = getJoint2Transform();
Cosserat::SE3<double> joint2_to_endEffector = getEndEffectorTransform();

// Computing end effector position in base frame
Cosserat::SE3<double> base_to_endEffector = 
    base_to_joint1.compose(joint1_to_joint2.compose(joint2_to_endEffector));
```

### When to use Sim(3)

- For camera calibration with unknown scale
- For monocular SLAM
- For multi-scale registration
- For morphing and animation

Example: Camera calibration, model alignment with scaling

```cpp
// Creating a similarity transform (rotation + translation + scaling)
Eigen::AngleAxisd rotation(M_PI/4, Eigen::Vector3d::UnitZ());
Eigen::Vector3d translation(1.0, 2.0, 3.0);
double scale = 2.5;
Cosserat::Sim3<double> transform(
    Cosserat::SO3<double>(rotation),
    translation,
    scale
);

// Transforming a point (rotate, scale, translate)
Eigen::Vector3d point(1.0, 0.0, 0.0);
Eigen::Vector3d transformed_point = transform.act(point);

// Aligning two datasets with different scales
Cosserat::Sim3<double> alignment = findOptimalAlignment(source_points, target_points);
std::vector<Eigen::Vector3d> aligned_points;
for (const auto& p : source_points) {
    aligned_points.push_back(alignment.act(p));
}
```

## Implementation Trade-offs

When implementing Lie groups, several important trade-offs must be considered:

### Representation Choice

| Representation | Advantages | Disadvantages |
|----------------|------------|--------------|
| **Rotation Matrix** | - Direct geometric interpretation<br>- Easy to visualize<br>- Simple composition (just matrix multiply) | - 9 parameters for 3 DOF rotation<br>- Numerical drift (losing orthogonality)<br>- More memory usage |
| **Quaternion** | - Compact (4 parameters for 3 DOF)<br>- Numerically stable<br>- Efficient composition | - Less intuitive<br>- Double cover (q = -q)<br>- Requires normalization |
| **Angle-Axis** | - Minimal representation for SO(3)<br>- Direct physical interpretation | - Singularity at zero angle<br>- Less efficient for composition |
| **Euler Angles** | - Intuitive for humans<br>- Minimal representation | - Gimbal lock<br>- Order-dependent<br>- Poor computational properties |

### Storage vs. Computation

1. **Storage Efficiency**:
   - Storing minimal representations (e.g., quaternion for SO(3)) saves memory
   - Particularly important for large datasets or memory-constrained environments

2. **Computational Efficiency**:
   - Caching frequently accessed representations (matrices, quaternions)
   - Pre-computing components for frequent operations

3. **Numerical Precision**:
   - Higher precision requires more memory
   - Double precision is typically recommended for most applications

### Template Parameters

1. **Scalar Type**:
   - `float`: Faster, less memory, but lower precision
   - `double`: Better precision, standard for most applications
   - `long double`: Highest precision, but slower and more memory-intensive

2. **Dimension**:
   - Fixed dimension: Better performance, compile-time checking
   - Dynamic dimension: More flexibility, runtime cost

### Inheritance vs. Composition

1. **Inheritance Approach**:
   - Useful for algorithms generic across different Lie groups
   - Enables polymorphism for heterogeneous collections
   - May have virtual function call overhead

2. **Composition Approach**:
   - More direct control over implementation
   - Can be more efficient (no virtual calls)
   - Potentially less code reuse

### Optimization Considerations

1. **Expression Templates**:
   - Can improve performance by avoiding temporary objects
   - Increases compile time and code complexity

2. **SIMD Optimization**:
   - Significant performance improvements, especially for batch operations
   - May require platform-specific code or intrinsics

3. **Memory Layout**:
   - Cache-friendly organization for batch operations
   - Trade-off between math clarity and optimization

### API Design

1. **Method Naming**:
   - `compose()` vs. operator `*` for group operation
   - `inverse()` vs. operator `-` for inverse element
   - Consistency with mathematical notation vs. programming conventions

2. **Error Handling**:
   - Assertions vs. exceptions vs. error returns
   - Performance impact of error checking in critical paths

