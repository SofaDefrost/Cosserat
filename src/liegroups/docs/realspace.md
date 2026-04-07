# RealSpace Implementation

## Overview

`RealSpace<Scalar, Dim>` implements the Euclidean vector space ℝⁿ as a Lie group. While it may seem trivial compared to other Lie groups, implementing it within the same framework provides consistency and allows for generic programming across different Lie groups.

In RealSpace, the group operation is vector addition, and the inverse operation is negation.

## Mathematical Properties

- **Dimension**: n (templated as `Dim`)
- **Group operation**: vector addition
- **Identity element**: zero vector
- **Inverse**: negation
- **Lie algebra**: also ℝⁿ
- **Exponential map**: identity function
- **Logarithmic map**: identity function

## Implementation Details

The `RealSpace` class is implemented as a template with two parameters:
- `Scalar`: The scalar type (typically `double` or `float`)
- `Dim`: The dimension (a positive integer or `Eigen::Dynamic` for runtime-sized vectors)

The internal representation uses an Eigen vector (`Eigen::Matrix<Scalar, Dim, 1>`).

### Key Methods

```cpp
// Constructor from coordinates
template <typename... Args>
RealSpace(Args&&... args);

// Constructor from Eigen vector
RealSpace(const VectorType& data);

// Group operations
RealSpace compose(const RealSpace& other) const;
RealSpace inverse() const;

// Access to internal representation
const VectorType& coeffs() const;
VectorType& coeffs();

// Tangent space (Lie algebra) operations
VectorType log() const;
static RealSpace exp(const VectorType& tangent);

// Accessors for individual components
Scalar operator[](int index) const;
Scalar& operator[](int index);
```

## Performance Characteristics

Based on benchmarks, RealSpace operations scale as follows:

- **Composition (addition)**: O(n) time complexity, where n is the dimension
- **Inverse (negation)**: O(n) time complexity
- **Exponential/Logarithmic maps**: O(1) time complexity (identity functions)

Performance is primarily limited by memory bandwidth for large dimensions.

## Example Usage

### Basic Operations

```cpp
#include <Cosserat/liegroups/RealSpace.h>
#include <iostream>

using Vec3 = Cosserat::RealSpace<double, 3>;

int main() {
    // Create points in 3D space
    Vec3 a(1.0, 2.0, 3.0);
    Vec3 b(4.0, 5.0, 6.0);
    
    // Composition (addition)
    Vec3 c = a.compose(b);
    std::cout << "a + b = [" << c.coeffs().transpose() << "]\n";
    
    // Inverse (negation)
    Vec3 a_inv = a.inverse();
    std::cout << "-a = [" << a_inv.coeffs().transpose() << "]\n";
    
    // Identity element
    Vec3 identity = Vec3::Identity();
    std::cout << "identity = [" << identity.coeffs().transpose() << "]\n";
    
    // Component access
    std::cout << "a[0] = " << a[0] << ", a[1] = " << a[1] << ", a[2] = " << a[2] << "\n";
    
    // Tangent space operations (trivial for RealSpace)
    Eigen::Vector3d log_a = a.log();
    Vec3 exp_log_a = Vec3::exp(log_a);
    std::cout << "exp(log(a)) = [" << exp_log_a.coeffs().transpose() << "]\n";
    
    return 0;
}
```

### Integration with Eigen

```cpp
#include <Cosserat/liegroups/RealSpace.h>
#include <Eigen/Dense>

using Vec3 = Cosserat::RealSpace<double, 3>;

int main() {
    // Create from Eigen vector
    Eigen::Vector3d eigen_vec(1.0, 2.0, 3.0);
    Vec3 point(eigen_vec);
    
    // Convert to Eigen vector
    Eigen::Vector3d converted = point.coeffs();
    
    // Use with Eigen operations
    Eigen::Vector3d doubled = 2.0 * point.coeffs();
    Vec3 doubled_point(doubled);
    
    return 0;
}
```

### Dynamic-sized Vectors

```cpp
#include <Cosserat/liegroups/RealSpace.h>

using VecX = Cosserat::RealSpace<double, Eigen::Dynamic>;

int main() {
    // Create a 5D vector
    VecX point(Eigen::VectorXd::Random(5));
    
    // Resize (only possible with dynamic-sized vectors)
    point.coeffs().resize(10);
    point.coeffs().setRandom();
    
    return 0;
}
```

## Best Practices

1. **Use fixed-size vectors when dimension is known at compile time** for better performance.
2. **Avoid unnecessary copies** by using references when passing RealSpace objects.
3. **Prefer direct access to coeffs()** when performing multiple operations on the internal vector.
4. **Use the compose() and inverse() methods** instead of direct arithmetic for consistency with other Lie groups.
5. **When using with Eigen functions**, access the internal representation using coeffs().

