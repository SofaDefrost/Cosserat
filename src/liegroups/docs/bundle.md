# Bundle Implementation

## Overview

`Bundle<Groups...>` implements a product manifold (bundle) of multiple Lie groups, allowing them to be treated as a single composite Lie group. The Bundle class enables working with combinations of different Lie groups in a unified framework, maintaining the product structure while providing all necessary group operations.

## Mathematical Properties

- **Dimension**: Sum of dimensions of all component groups
- **Group operation**: Component-wise composition (g₁h₁, g₂h₂, ..., gₙhₙ)
- **Identity element**: Tuple of identity elements (e₁, e₂, ..., eₙ)
- **Inverse**: Tuple of component inverses (g₁⁻¹, g₂⁻¹, ..., gₙ⁻¹)
- **Lie algebra**: Direct sum of component Lie algebras (𝔤₁ ⊕ 𝔤₂ ⊕ ... ⊕ 𝔤ₙ)
- **Exponential map**: Component-wise exponential maps
- **Logarithmic map**: Component-wise logarithmic maps

## Implementation Details

The `Bundle` class is implemented as a variadic template with multiple Lie group parameters:
- `Groups...`: The Lie groups to bundle together (must all have the same scalar type)

### Key Methods

```cpp
// Constructors
Bundle(); // Identity bundle
Bundle(const Groups&... groups); // From individual group elements
Bundle(const GroupTuple& groups); // From tuple of groups
Bundle(const TangentVector& algebra_element); // From Lie algebra vector

// Group operations
Bundle operator*(const Bundle& other) const; // Composition
Bundle& operator*=(const Bundle& other); // In-place composition
Bundle inverse() const; // Inverse

// Lie algebra operations
Bundle exp(const TangentVector& algebra_element) const; // Exponential map
TangentVector log() const; // Logarithmic map
AdjointMatrix adjoint() const; // Adjoint representation

// Group actions
Vector act(const Vector& point) const; // Group action on point
template<int N> Eigen::Matrix<Scalar, Eigen::Dynamic, N> act(const Eigen::Matrix<Scalar, Eigen::Dynamic, N>& points) const; // Batch action

// Utility functions
bool isApprox(const Bundle& other, const Scalar& eps = Types<Scalar>::epsilon()) const;
Bundle interpolate(const Bundle& other, const Scalar& t) const; // Linear interpolation
template<typename Generator> static Bundle random(Generator& gen, const Scalar& scale = Scalar(1));

// Accessors
template<std::size_t I> const auto& get() const; // Access individual group (const)
template<std::size_t I> auto& get(); // Access individual group (mutable)
template<std::size_t I> void set(const GroupType<I>& group); // Set individual group
const GroupTuple& groups() const; // Get underlying tuple
int algebraDimension() const; // Get Lie algebra dimension
int actionDimension() const; // Get action space dimension
```

## Performance Characteristics

- **Composition**: O(sum of component group complexities)
- **Inverse**: O(sum of component group complexities)
- **Exponential/Logarithmic maps**: O(sum of component group complexities)
- **Acting on points**: O(sum of component action complexities)
- **Memory footprint**: Sum of component memory footprints

The Bundle class is designed to have minimal overhead beyond the component groups themselves.

## Example Usage

### Basic Bundle Creation

```cpp
#include <Cosserat/liegroups/Bundle.h>
#include <iostream>

// Create a bundle of SE(3) pose and 6D velocity
using RigidBodyState = Cosserat::Bundle<Cosserat::SE3<double>, Cosserat::RealSpace<double, 6>>;

int main() {
    // Create identity bundle
    RigidBodyState identity;
    std::cout << "Identity: " << identity << "\n";

    // Create from components
    Cosserat::SE3<double> pose(Eigen::AngleAxisd(M_PI/4, Eigen::Vector3d::UnitZ()),
                              Eigen::Vector3d(1.0, 2.0, 3.0));
    Cosserat::RealSpace<double, 6> velocity;
    velocity.coeffs() << 0.1, 0.2, 0.3, 0.0, 0.0, 0.0; // Linear and angular velocity

    RigidBodyState state(pose, velocity);
    std::cout << "State: " << state << "\n";

    return 0;
}
```

### Group Operations

```cpp
#include <Cosserat/liegroups/Bundle.h>

using RigidBodyState = Cosserat::Bundle<Cosserat::SE3<double>, Cosserat::RealSpace<double, 6>>;

int main() {
    // Create two states
    RigidBodyState state1(
        Cosserat::SE3<double>(Eigen::AngleAxisd(M_PI/4, Eigen::Vector3d::UnitZ()),
                             Eigen::Vector3d(1.0, 0.0, 0.0)),
        Cosserat::RealSpace<double, 6>(Eigen::VectorXd::Random(6))
    );

    RigidBodyState state2(
        Cosserat::SE3<double>(Eigen::AngleAxisd(M_PI/6, Eigen::Vector3d::UnitY()),
                             Eigen::Vector3d(0.0, 1.0, 0.0)),
        Cosserat::RealSpace<double, 6>(Eigen::VectorXd::Random(6))
    );

    // Compose states (relative transformation)
    RigidBodyState composed = state1 * state2;
    std::cout << "Composed state: " << composed << "\n";

    // Inverse
    RigidBodyState inv = state1.inverse();
    std::cout << "Inverse of state1: " << inv << "\n";

    // Check that state1 * inv ≈ identity
    RigidBodyState should_be_identity = state1 * inv;
    std::cout << "state1 * inv ≈ identity: " << should_be_identity.isApprox(RigidBodyState()) << "\n";

    return 0;
}
```

### Accessing Components

```cpp
#include <Cosserat/liegroups/Bundle.h>

using RigidBodyState = Cosserat::Bundle<Cosserat::SE3<double>, Cosserat::RealSpace<double, 6>>;

int main() {
    RigidBodyState state(
        Cosserat::SE3<double>(Eigen::AngleAxisd(M_PI/3, Eigen::Vector3d::UnitZ()),
                             Eigen::Vector3d(1.0, 2.0, 3.0)),
        Cosserat::RealSpace<double, 6>(Eigen::VectorXd::Random(6))
    );

    // Access components (const)
    const auto& pose = state.get<0>(); // SE(3) pose
    const auto& velocity = state.get<1>(); // 6D velocity

    std::cout << "Pose: " << pose << "\n";
    std::cout << "Velocity: " << velocity.coeffs().transpose() << "\n";

    // Modify components
    Cosserat::SE3<double> new_pose(Eigen::AngleAxisd(M_PI/2, Eigen::Vector3d::UnitX()),
                                  Eigen::Vector3d(4.0, 5.0, 6.0));
    state.set<0>(new_pose);

    Cosserat::RealSpace<double, 6> new_velocity;
    new_velocity.coeffs() << 1.0, 0.0, 0.0, 0.0, 0.0, 1.0;
    state.set<1>(new_velocity);

    std::cout << "Modified state: " << state << "\n";

    return 0;
}
```

### Lie Algebra Operations

```cpp
#include <Cosserat/liegroups/Bundle.h>

using RigidBodyState = Cosserat::Bundle<Cosserat::SE3<double>, Cosserat::RealSpace<double, 6>>;

int main() {
    RigidBodyState state(
        Cosserat::SE3<double>(Eigen::AngleAxisd(0.1, Eigen::Vector3d::UnitZ()),
                             Eigen::Vector3d(0.1, 0.2, 0.3)),
        Cosserat::RealSpace<double, 6>(Eigen::VectorXd::Random(6) * 0.1)
    );

    // Convert to Lie algebra
    auto algebra_element = state.log();
    std::cout << "Lie algebra element: " << algebra_element.transpose() << "\n";

    // Convert back to group
    auto recovered = state.exp(algebra_element);
    std::cout << "Recovered state ≈ original: " << recovered.isApprox(state) << "\n";

    // Create state from Lie algebra element
    Eigen::VectorXd new_algebra(12); // 6 (SE3) + 6 (velocity)
    new_algebra << 0.1, 0.0, 0.0, 0.0, 0.0, 0.1,  // Small SE(3) perturbation
                   0.01, 0.02, 0.03, 0.0, 0.0, 0.01; // Small velocity
    RigidBodyState new_state(new_algebra);

    return 0;
}
```

### Group Actions

```cpp
#include <Cosserat/liegroups/Bundle.h>

using RigidBodyState = Cosserat::Bundle<Cosserat::SE3<double>, Cosserat::RealSpace<double, 6>>;

int main() {
    RigidBodyState transform(
        Cosserat::SE3<double>(Eigen::AngleAxisd(M_PI/4, Eigen::Vector3d::UnitZ()),
                             Eigen::Vector3d(1.0, 0.0, 0.0)),
        Cosserat::RealSpace<double, 6>::zero() // No velocity for pure transformation
    );

    // The action dimension depends on the component groups
    // For this bundle: SE(3) acts on 3D points, RealSpace<6> acts on 6D vectors
    // Total action dimension is computed at runtime
    int action_dim = transform.actionDimension();
    std::cout << "Action dimension: " << action_dim << "\n";

    // Create a point in the action space
    Eigen::VectorXd point(action_dim);
    point << 1.0, 0.0, 0.0,  // 3D point for SE(3)
            0.1, 0.0, 0.0, 0.0, 0.0, 0.1; // 6D vector for RealSpace<6>

    // Apply the group action
    Eigen::VectorXd transformed_point = transform.act(point);
    std::cout << "Original point: " << point.transpose() << "\n";
    std::cout << "Transformed point: " << transformed_point.transpose() << "\n";

    return 0;
}
```

## Applications

Bundles are particularly useful for:

- **Multi-body systems**: Representing states of articulated systems
- **Extended configurations**: Combining pose with velocity or other state variables
- **Sensor fusion**: Combining measurements from different coordinate frames
- **Optimization**: Working with high-dimensional configuration spaces
- **Control systems**: Representing full system states including dynamics

## Type Aliases

The library provides convenient type aliases for common bundles:

```cpp
// Rigid body with velocity
template<typename Scalar>
using SE3_Velocity = Bundle<SE3<Scalar>, RealSpace<Scalar, 6>>;

// 2D rigid body with velocity
template<typename Scalar>
using SE2_Velocity = Bundle<SE2<Scalar>, RealSpace<Scalar, 3>>;

// Robot with multiple joints
template<typename Scalar, int N>
using SE3_Joints = Bundle<SE3<Scalar>, RealSpace<Scalar, N>>;

// Float and double versions
template<typename... Groups>
using Bundlef = Bundle<Groups...>;

template<typename... Groups>
using Bundled = Bundle<Groups...>;
```

## Best Practices

1. **Choose appropriate component groups** based on your physical system
2. **Use compile-time indices** for accessing components to enable optimization
3. **Cache bundle objects** when performing multiple operations
4. **Consider memory layout** for cache efficiency in performance-critical code
5. **Use the appropriate constructor** based on how you have the component data
6. **Validate dimensions** when working with Lie algebra elements or action vectors
7. **Remember that bundles maintain product structure** - operations are component-wise

## Performance Considerations

- **Template instantiation**: Each unique bundle type requires separate compilation
- **Memory layout**: Components are stored in a tuple, which may not be cache-optimal for all access patterns
- **Runtime action dimensions**: Some properties are computed at runtime due to varying action dimensions
- **Composition complexity**: Scales with the complexity of component group operations

## Error Handling

The Bundle class includes runtime checks for:
- **Algebra element dimensions**: Ensures input vectors match the expected Lie algebra dimension
- **Action input dimensions**: Validates that points have the correct dimension for the action space
- **Interpolation parameters**: Checks that interpolation parameter t is in [0,1]

Invalid inputs throw `std::invalid_argument` exceptions with descriptive error messages.</content>
</xai:function_call/>
<xai:function_call name="write">
<parameter name="filePath">src/liegroups/docs/gaussian_on_manifold.md