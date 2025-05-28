# Lie Groups Usage Guide

This guide provides detailed examples and best practices for using the Lie groups implementation in Cosserat mechanics applications.

## Table of Contents
1. Common Usage Patterns
2. Advanced Applications
3. Best Practices
4. Edge Cases and Pitfalls
5. Integration with SOFA
6. Performance Optimization

## 1. Common Usage Patterns

### Material Frame Construction
```cpp
// Create a material frame from position and directors
SE3d makeMaterialFrame(const Vector3& position,
                      const Vector3& tangent,
                      const Vector3& normal) {
    // Align z-axis with tangent
    Vector3 z_axis(0, 0, 1);
    Vector3 rot_axis = z_axis.cross(tangent);
    double angle = std::acos(z_axis.dot(tangent));
    SO3d R(angle, rot_axis.normalized());

    // Additional rotation to align normal
    Vector3 current_normal = R.act(Vector3(1, 0, 0));
    Vector3 target_normal = normal - normal.dot(tangent) * tangent;
    double twist = std::acos(current_normal.dot(target_normal.normalized()));
    SO3d twist_rot(twist, tangent);
    
    return SE3d(twist_rot * R, position);
}
```

### Strain and Curvature
```cpp
// Compute strain between two frames
Vector6 computeStrain(const SE3d& frame1, const SE3d& frame2, double ds) {
    SE3d relative = frame1.inverse() * frame2;
    return relative.log() / ds;  // Normalize by segment length
}

// Extract curvature from strain
Vector3 getCurvature(const Vector6& strain) {
    return strain.tail<3>();  // Angular velocity components
}
```

### Multi-Body Systems
```cpp
// Define a system with multiple bodies
using MultiBody = Bundle<SE3d, SE3d, SE3d>;

// Create and update configuration
MultiBody system(body1, body2, body3);
system.get<0>() = newBody1State;
system.get<1>() = newBody2State;
system.get<2>() = newBody3State;

// Compute relative motions
SE3d relative01 = system.get<0>().inverse() * system.get<1>();
SE3d relative12 = system.get<1>().inverse() * system.get<2>();
```

## 2. Advanced Applications

### Cosserat Rod Dynamics
```cpp
// Define rod state including velocity
using RodState = Bundle<SE23d, RealSpace3d>;  // Configuration + strain

// Update rod configuration
void updateRodState(RodState& state, double dt) {
    // Get current values
    const auto& config = state.get<0>();
    const auto& strain = state.get<1>();
    
    // Create twist from velocity
    Vector6 twist = config.velocity();
    
    // Update configuration using exponential map
    SE3d delta = SE3d().exp(dt * twist);
    SE3d new_pose = config.pose() * delta;
    
    // Update state
    state = RodState(
        SE23d(new_pose, config.velocity()),
        strain
    );
}
```

### Time-Varying Trajectories
```cpp
// Create smooth interpolation
std::vector<SE3d> interpolateTrajectory(
    const SE3d& start,
    const SE3d& end,
    int steps
) {
    std::vector<SE3d> trajectory;
    trajectory.reserve(steps);
    
    for (int i = 0; i < steps; ++i) {
        double t = static_cast<double>(i) / (steps - 1);
        trajectory.push_back(interpolate(start, end, t));
    }
    
    return trajectory;
}
```

## 3. Best Practices

### Memory Management
- Prefer stack allocation for small fixed-size groups
- Use references for large composite groups
- Cache frequently used values (e.g., rotation matrices)

### Numerical Stability
```cpp
// Handle small angles in exponential map
if (angle < eps) {
    // Use small angle approximation
    return identity() + hat(omega);
} else {
    // Use full Rodriguez formula
    return computeRodriguez(angle, axis);
}
```

### Type Safety
```cpp
// Use type aliases for clarity
using StrainVector = RealSpace<double, 6>;
using StressVector = RealSpace<double, 6>;
using RodSection = Bundle<SE3d, StrainVector, StressVector>;
```

## 4. Edge Cases and Pitfalls

### Singularities
- Handle gimbal lock in Euler angle conversions
- Check for zero denominators in logarithm maps
- Handle parallel vectors in cross products

### Numerical Issues
```cpp
// Normalize quaternions periodically
void normalizeRotation(SO3d& rotation) {
    rotation = SO3d(rotation.quaternion().normalized());
}

// Handle angle wrapping
double normalizeAngle(double angle) {
    return std::fmod(angle + pi, 2*pi) - pi;
}
```

## 5. Integration with SOFA

### State Vectors
```cpp
// Convert to/from SOFA state vectors
void toSOFAVector(const SE3d& transform, Vector& state) {
    // Position
    state.segment<3>(0) = transform.translation();
    // Orientation (as quaternion)
    const auto& q = transform.rotation().quaternion();
    state.segment<4>(3) << q.w(), q.x(), q.y(), q.z();
}
```

### Mappings
```cpp
// Example mapping between frames
void computeJacobian(const SE3d& frame, Matrix& J) {
    // Fill Jacobian blocks
    J.block<3,3>(0,0) = Matrix3::Identity();
    J.block<3,3>(0,3) = -SO3d::hat(frame.translation());
    J.block<3,3>(3,3) = Matrix3::Identity();
}
```

## 6. Performance Optimization

### Caching Strategies
```cpp
class CachedTransform {
    SE3d transform;
    mutable Matrix4d matrix;
    mutable bool matrix_valid = false;

public:
    const Matrix4d& getMatrix() const {
        if (!matrix_valid) {
            matrix = transform.matrix();
            matrix_valid = true;
        }
        return matrix;
    }
    
    void setTransform(const SE3d& new_transform) {
        transform = new_transform;
        matrix_valid = false;
    }
};
```

### Parallel Operations
```cpp
// Parallel strain computation for rod segments
void computeStrains(
    const std::vector<SE3d>& frames,
    std::vector<Vector6>& strains
) {
    const size_t n = frames.size() - 1;
    strains.resize(n);
    
    #pragma omp parallel for
    for (size_t i = 0; i < n; ++i) {
        strains[i] = computeStrain(frames[i], frames[i+1]);
    }
}
```

### Memory Layout
```cpp
// Optimize for cache coherency
struct RodSegment {
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    SE3d frame;
    Vector6 strain;
    Vector6 stress;
};
```

## Additional Resources

- See benchmark results in `LieGroupBenchmark.cpp`
- Check integration tests in `LieGroupIntegrationTest.cpp`
- Refer to the mathematical background in README.md

