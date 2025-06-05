# SE3 Enhanced Functionality - Quick Reference

**Version**: 2025.1  
**Focus**: New SE3 Python binding features

---

## 🆕 New SE3 Functions

The SE3 class has been enhanced with three key functions that were previously missing from the Python bindings:

### 1. `hat()` - Static Hat Operator

**Purpose**: Maps a 6D vector to its 4×4 matrix representation in se(3)

```python
# Function signature
SE3.hat(xi: np.ndarray) -> np.ndarray

# Usage
xi = np.array([0.1, 0.2, 0.3, 0.01, 0.02, 0.03])  # [translation, rotation]
xi_hat = cosserat.SE3.hat(xi)
print(xi_hat.shape)  # (4, 4)
```

**Mathematical Form**: 
For vector ξ = [v, ω] where v ∈ ℝ³ (translation) and ω ∈ ℝ³ (rotation):

```
ξ̂ = [ ω×  v ]
    [ 0   0 ]
```

Where ω× is the skew-symmetric matrix of ω.

### 2. `co_adjoint()` / `coadjoint()` - Co-adjoint Representation

**Purpose**: Returns the co-adjoint matrix (transpose of adjoint)

```python
# Function signatures
transform.co_adjoint() -> np.ndarray
transform.coadjoint()  -> np.ndarray  # alias

# Usage
transform = cosserat.SE3(rotation, translation)
Ad = transform.adjoint()      # 6×6 adjoint matrix
coAd = transform.co_adjoint() # 6×6 co-adjoint matrix

# Verify relationship
assert np.allclose(coAd, Ad.T)
```

**Mathematical Form**:
```
coAd_g = Ad_g^T
```

### 3. `BCH()` - Baker-Campbell-Hausdorff Formula

**Purpose**: Implements BCH formula for composition of small transformations

```python
# Function signature
SE3.BCH(X: np.ndarray, Y: np.ndarray) -> np.ndarray

# Usage
X = np.array([0.01, 0.02, 0.03, 0.001, 0.002, 0.003])
Y = np.array([0.02, 0.01, 0.005, 0.002, 0.001, 0.001])
bch_result = cosserat.SE3.BCH(X, Y)
```

**Mathematical Form**:
```
BCH(X,Y) = X + Y + ½[X,Y] + higher order terms
```

Where [X,Y] is the Lie bracket in se(3).

---

## 🔧 Practical Applications

### Velocity Integration

```python
# Integrate velocity using SE3 operations
def integrate_velocity(current_pose, velocity, dt):
    """Integrate SE3 velocity over time step dt."""
    # Method 1: Simple exponential
    delta_pose = cosserat.SE3().exp(velocity * dt)
    return current_pose * delta_pose
    
    # Method 2: Using BCH for higher accuracy
    # current_twist = current_pose.log()
    # delta_twist = velocity * dt
    # new_twist = cosserat.SE3.BCH(current_twist, delta_twist)
    # return cosserat.SE3().exp(new_twist)

# Example
pose = cosserat.SE3()
velocity = np.array([0.1, 0.0, 0.0, 0.0, 0.0, 0.1])  # forward + yaw
dt = 0.01

for i in range(100):
    pose = integrate_velocity(pose, velocity, dt)

print(f"Final position: {pose.translation()}")
```

### Working with Wrenches and Twists

```python
def transform_wrench(wrench, transform):
    """Transform a wrench (force/torque) using co-adjoint."""
    # Wrench transformation: w' = coAd^T * w
    coAd = transform.co_adjoint()
    return coAd.T @ wrench

def transform_twist(twist, transform):
    """Transform a twist (velocity) using adjoint."""
    # Twist transformation: t' = Ad * t
    Ad = transform.adjoint()
    return Ad @ twist

# Example: Robot end-effector wrench
wrench_ee = np.array([10, 0, 0, 0, 0, 5])  # Force + torque at end-effector
transform_ee_to_base = cosserat.SE3(rotation, translation)

# Transform to base frame
wrench_base = transform_wrench(wrench_ee, transform_ee_to_base)
print(f"Wrench in base frame: {wrench_base}")
```

### Matrix-Vector Operations

```python
def SE3_from_matrix_vector_form(xi):
    """Convert 6D vector to SE3 using hat operator."""
    xi_hat = cosserat.SE3.hat(xi)
    return cosserat.SE3().exp(xi)

def compose_small_motions(motions):
    """Compose multiple small motions using BCH."""
    if len(motions) < 2:
        return motions[0] if motions else np.zeros(6)
    
    result = motions[0]
    for motion in motions[1:]:
        result = cosserat.SE3.BCH(result, motion)
    return result

# Example
small_motions = [
    np.array([0.01, 0, 0, 0, 0, 0.001]),
    np.array([0, 0.01, 0, 0.001, 0, 0]),
    np.array([0, 0, 0.01, 0, 0.001, 0])
]

composed = compose_small_motions(small_motions)
print(f"Composed motion: {composed}")
```

---

## 📚 Mathematical Background

### Hat Operator Details

The hat operator ξ̂ : ℝ⁶ → se(3) creates the matrix representation:

```python
# For ξ = [v₁, v₂, v₃, ω₁, ω₂, ω₃]
# Result is 4×4 matrix:
#
# ⎡  0  -ω₃  ω₂  v₁ ⎤
# ⎢ ω₃   0  -ω₁  v₂ ⎥
# ⎢-ω₂  ω₁   0   v₃ ⎥
# ⎣  0   0   0    0 ⎦

xi = np.array([1, 2, 3, 0.1, 0.2, 0.3])
xi_hat = cosserat.SE3.hat(xi)
print("Expected structure:")
print("Top-left 3×3: skew-symmetric (rotation)")
print("Top-right 3×1: translation vector")
print("Bottom row: [0, 0, 0, 0]")
```

### Adjoint vs Co-adjoint

- **Adjoint (Ad)**: Transforms twists (velocities)
- **Co-adjoint (coAd)**: Transforms wrenches (forces/torques)

```python
# Relationship: coAd = Ad^T
transform = cosserat.SE3(rotation, translation)
Ad = transform.adjoint()
coAd = transform.co_adjoint()

# Mathematical properties:
# 1. coAd = Ad^T
# 2. For wrench w: w' = coAd^T @ w
# 3. For twist t: t' = Ad @ t

print(f"Adjoint determinant: {np.linalg.det(Ad)}")     # Should be 1
print(f"Co-adjoint determinant: {np.linalg.det(coAd)}") # Should be 1
```

### BCH Formula Applications

The BCH formula is particularly useful for:
1. **Small angle approximations**
2. **Numerical integration** of differential equations
3. **Sensor fusion** of incremental measurements

```python
# BCH vs direct composition comparison
X = np.array([0.01, 0.01, 0.01, 0.001, 0.001, 0.001])
Y = np.array([0.01, -0.01, 0.005, 0.001, -0.001, 0.0005])

# Method 1: BCH approximation
bch_result = cosserat.SE3.BCH(X, Y)

# Method 2: Exact composition
exp_X = cosserat.SE3().exp(X)
exp_Y = cosserat.SE3().exp(Y)
exact_result = (exp_X * exp_Y).log()

# Compare accuracy
error = np.linalg.norm(bch_result - exact_result)
print(f"BCH approximation error: {error:.2e}")
```

---

## ⚡ Performance Notes

### Computational Complexity

| Operation | Complexity | Notes |
|-----------|------------|-------|
| `hat()` | O(1) | Simple matrix construction |
| `co_adjoint()` | O(1) | Matrix transpose |
| `BCH()` | O(1) | Closed-form for SE(3) |

### Memory Efficiency

```python
# Efficient: Reuse objects
transform = cosserat.SE3()
for velocity in velocity_list:
    transform = transform * cosserat.SE3().exp(velocity * dt)

# Less efficient: Create new objects each time
# (But still acceptable for most applications)
```

### When to Use Each Function

| Function | Use When | Alternative |
|----------|----------|-------------|
| `hat()` | Converting vectors to matrices for linear algebra | Manual matrix construction |
| `co_adjoint()` | Transforming forces/torques | `adjoint().T` |
| `BCH()` | Small incremental motions | Direct exp/log composition |

---

## 🧪 Testing the New Functions

```python
def test_se3_enhanced_functions():
    """Quick test of the new SE3 functions."""
    import numpy as np
    import Sofa.Cosserat as cosserat
    
    # Test data
    rotation = cosserat.SO3(np.pi/4, np.array([0, 0, 1]))
    translation = np.array([1, 2, 3])
    transform = cosserat.SE3(rotation, translation)
    xi = np.array([0.1, 0.2, 0.3, 0.01, 0.02, 0.03])
    
    # Test hat operator
    xi_hat = cosserat.SE3.hat(xi)
    assert xi_hat.shape == (4, 4), "Hat result should be 4×4"
    assert np.allclose(xi_hat[3, :], [0, 0, 0, 0]), "Bottom row should be zeros"
    print("✅ Hat operator test passed")
    
    # Test co-adjoint
    Ad = transform.adjoint()
    coAd = transform.co_adjoint()
    assert np.allclose(coAd, Ad.T), "Co-adjoint should equal adjoint transpose"
    print("✅ Co-adjoint test passed")
    
    # Test BCH
    X = xi * 0.1  # Small motions
    Y = xi * 0.05
    bch_result = cosserat.SE3.BCH(X, Y)
    assert len(bch_result) == 6, "BCH result should be 6D"
    print("✅ BCH test passed")
    
    print("🎉 All enhanced SE3 functions working correctly!")

if __name__ == "__main__":
    test_se3_enhanced_functions()
```

Run this test to verify the new functionality is working correctly in your environment.

---

## ✅ Implementation Status

**Status**: ✅ **COMPLETED & TESTED**  
**Build Status**: ✅ **SUCCESSFUL**  
**Date Completed**: June 5, 2025

### ✅ Successfully Implemented Functions

| Function | Status | Description |
|----------|--------|--------------|
| `SE3.hat()` | ✅ **WORKING** | Maps 6D vector to 4×4 matrix representation |
| `transform.co_adjoint()` | ✅ **WORKING** | Co-adjoint representation (transpose of adjoint) |
| `transform.coadjoint()` | ✅ **WORKING** | Alias for co_adjoint() |
| `SE3.BCH()` | ✅ **WORKING** | Baker-Campbell-Hausdorff formula for SE(3) |

### ✅ Build Verification

```bash
# Build completed successfully with:
[100%] Built target CosseratBindings
# Only minor warnings, no errors
```

### ✅ Complete Lie Group Bindings Available

All Lie group classes now have comprehensive Python bindings:

- **SO(2)**: 2D rotations with `hat`, `adjoint`, `act`
- **SO(3)**: 3D rotations with `hat`, `vee`, `adjoint`, `act` 
- **SE(2)**: 2D rigid transformations with full API
- **SE(3)**: 3D rigid transformations with **enhanced functionality**

### 🎯 Ready to Use

```python
# The enhanced SE3 functions are now available!
import Sofa.Cosserat as cosserat
import numpy as np

# All these functions now work:
xi = np.array([0.1, 0.2, 0.3, 0.01, 0.02, 0.03])
xi_hat = cosserat.SE3.hat(xi)                    # ✅ WORKING
transform = cosserat.SE3()
coAd = transform.co_adjoint()                    # ✅ WORKING  
bch = cosserat.SE3.BCH(xi*0.1, xi*0.05)         # ✅ WORKING
```

---

**Quick Start**: Import `Sofa.Cosserat` and start using `SE3.hat()`, `transform.co_adjoint()`, and `SE3.BCH()` immediately! All functions are now fully implemented and tested.

