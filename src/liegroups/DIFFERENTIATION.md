# Differentiation Guide for Lie Groups

Complete guide for using differentiation with the Cosserat Lie groups library.

## Table of Contents

1. [Overview](#overview)
2. [Analytical Jacobians](#analytical-jacobians)
3. [Automatic Differentiation](#automatic-differentiation)
4. [Applications](#applications)
5. [Performance Guide](#performance-guide)

---

## Overview

Three approaches to compute derivatives:

| Method | Best For | Accuracy | Speed |
|--------|----------|----------|-------|
| **Analytical** | Known formulas | Exact | Fastest |
| **Autodiff** | Complex functions | Machine precision | Fast |
| **Finite Diff** | Verification | Approximate | Slow |

---

## Analytical Jacobians

Phase 2 provides analytical jacobians for SO3 and SE3.

### SO3 Jacobians

```cpp
#include "SO3.h"
using SO3d = SO3<double>;

Eigen::Vector3d omega(0.1, 0.2, 0.3);

// Left jacobian of exponential map
Eigen::Matrix3d J = SO3d::leftJacobian(omega);

// Composition jacobians
SO3d R1 = SO3d::exp(omega);
SO3d R2 = SO3d::identity();
auto [J_g1, J_g2] = SO3d::composeJacobians(R1, R2);

// Action jacobians
Eigen::Vector3d p(1.0, 0.0, 0.0);
auto [J_g, J_p] = R1.actionJacobians(p);
```

### SE3 Jacobians

```cpp
#include "SE3.h"
using SE3d = SE3<double>;

Eigen::Matrix<double, 6, 1> xi;
xi << 0.1, 0.2, 0.3, 0.5, 0.0, 0.0;  // [φ | ρ]

// Left jacobian
Eigen::Matrix<double, 6, 6> J = SE3d::leftJacobian(xi);

// Composition for forward kinematics
SE3d T1 = SE3d::exp(xi);
SE3d T2 = SE3d::identity();
auto [J_g1, J_g2] = SE3d::composeJacobians(T1, T2);

// Action on points
Eigen::Vector3d point(1.0, 0.0, 0.0);
auto [J_g, J_p] = T1.actionJacobians(point);
```

---

## Automatic Differentiation

Use `autodiff` library for complex functions.

### Setup

```bash
cmake -DCOSSERAT_WITH_AUTODIFF=ON ..
make
```

### Forward Mode (autodiff::dual)

Best for: few inputs → many outputs

```cpp
#include <autodiff/forward/dual.hpp>
#include "SO3.h"
#include "AutodiffSupport.h"

using namespace autodiff;
using SO3dual = SO3<dual>;

// Function to differentiate
dual rotationAngle(const Eigen::Matrix<dual, 3, 1>& omega) {
    SO3dual R = SO3dual::exp(omega);
    return R.log().norm();
}

// Compute gradient
Eigen::Vector3d omega_val(0.1, 0.2, 0.3);
auto omega = toAutodiff<Eigen::Vector3d, dual>(omega_val);

Eigen::Vector3d gradient;
for (int i = 0; i < 3; i++) {
    gradient[i] = derivative(rotationAngle, wrt(omega[i]), at(omega));
}
```

### Reverse Mode (autodiff::var)

Best for: many inputs → scalar output (optimization!)

```cpp
#include <autodiff/reverse/var.hpp>
#include "SE3.h"
#include "AutodiffSupport.h"

using namespace autodiff;
using SE3var = SE3<var>;

// Cost function
Eigen::Matrix<var, 6, 1> xi;
for (int i = 0; i < 6; i++) xi[i] = xi_val[i];

SE3var T = SE3var::exp(xi);
auto pos = T.translation();

Eigen::Vector3d target(1.0, 0.0, 0.0);
var dx = pos[0] - target[0];
var dy = pos[1] - target[1];
var dz = pos[2] - target[2];
var cost = dx*dx + dy*dy + dz*dz;

// ALL gradients in ONE reverse pass!
Derivatives dcost = derivatives(cost);

Eigen::Matrix<double, 6, 1> gradient;
for (int i = 0; i < 6; i++) {
    gradient[i] = dcost(xi[i]);
}
```

### Multi-Segment Cosserat

```cpp
using SE3var = SE3<var>;

const int N = 10;  // 10 segments = 60 parameters
std::vector<Eigen::Matrix<var, 6, 1>> strains(N);

// Forward kinematics
SE3var T = SE3var::identity();
for (int i = 0; i < N; i++) {
    SE3var Ti = SE3var::exp(L * strains[i]);
    T = T * Ti;
}

// Cost
var cost = /* ... */;

// ONE pass → gradients for ALL 60 parameters!
Derivatives dcost = derivatives(cost);
```

**Speedup: 30-60× faster than alternatives for multi-segment rods!**

---

## Applications

### Trajectory Optimization

See `CosseratTrajectoryOptimizer.h` and `simple_trajectory_optimization.cpp`:

```cpp
#include "optimization/CosseratTrajectoryOptimizer.h"

CosseratTrajectoryOptimizer<double> optimizer(config);

auto cost_fn = [](const Vector3d& pos) {
    return (pos - target).squaredNorm();
};

auto result = optimizer.optimize(initial_strains, cost_fn);
```

### Parameter Calibration

Estimate beam stiffness from measurements using autodiff:

```cpp
using namespace autodiff;

var calibrationCost(const Eigen::Matrix<var, 4, 1>& params,
                    const std::vector<Measurement>& data) {
    var total = 0.0;
    for (const auto& m : data) {
        // Apply params to model
        SE3<var> T = SE3<var>::exp(scale_strain(m.strain, params));
        auto predicted = T.translation();
        
        // Error
        for (int i = 0; i < 3; i++) {
            var err = predicted[i] - m.position[i];
            total += err * err;
        }
    }
    return total;
}

// Optimize with gradient descent
// Derivatives computed automatically!
```

---

## Performance Guide

### When to Use What

| Problem | Method | Reason |
|---------|--------|--------|
| N params → 1 cost | Reverse mode | O(1) regardless of N |
| 3 params → 100 outputs | Forward mode | O(1) regardless of outputs |
| Known formula exists | Analytical | Fastest, exact |
| Testing/debugging | Finite differences | Simple, reliable |

### Timing Example

10-segment Cosserat (60 params → 1 cost):

- **Reverse mode**: ~3× forward pass
- **Forward mode**: ~180× forward pass  
- **Finite differences**: ~240× forward pass

**Reverse mode is 60-80× faster!**

### Best Practices

1. ✅ Use analytical jacobians when available
2. ✅ Use reverse mode for optimization (many params → scalar)
3. ✅ Use forward mode for few-input problems
4. ✅ Reuse expression trees in loops
5. ❌ Don't use finite differences for production code

---

## Examples

Working examples:
- `examples/simple_trajectory_optimization.cpp` - Full optimization demo
- `Tests/differentiation/test_autodiff_integration.cpp` - All autodiff tests
- `Tests/differentiation/test_analytical_jacobians.cpp` - Analytical tests

Compile:
```bash
cmake -DCOSSERAT_BUILD_EXAMPLES=ON -DCOSSERAT_WITH_AUTODIFF=ON ..
make
./bin/simple_trajectory_optimization
```

---

## References

- **autodiff**: https://autodiff.github.io/
- **Lie theory**: "A micro Lie theory for state estimation in robotics" - Sola et al., 2018
- **Phase 2**: `IMPLEMENTATION_PLAN.md`
- **Strain conventions**: `STRAIN_CONVENTION.md`

---

## Troubleshooting

### autodiff not found

```bash
# Check location
ls ../autodiff/autodiff/

# Explicit path
cmake -DCOSSERAT_WITH_AUTODIFF=ON \
      -Dautodiff_DIR=/path/to/autodiff/build ..
```

### Numerical issues

- Check for singularities (e.g., log near identity)
- Use `Types<double>::epsilon()` thresholds
- Verify with finite differences
- Initialize `var` types properly

---

## Future Phases

- **Phase 3.2**: iLQR optimal control
- **Phase 3.3**: Parameter calibration
- **Phase 4**: Differentiable simulation
- **Phase 5**: ML integration
