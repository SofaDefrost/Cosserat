# Step 1: Jacobian Verification - Test Documentation

## Overview

This test suite verifies the correctness of the `computeTangExpImplementation()` method in `CosseratGeometryMapping`.

## Test File

**Location:** `tests/test_jacobian_verification.cpp`

## Test Cases

### 1. SmallStrain
- **Purpose:** Verify behavior near zero strain
- **Expected:** Should use first-order approximation (≈ curv_abs * I)
- **Tolerance:** 1e-4

### 2. ZeroStrain  
- **Purpose:** Verify exact behavior at zero
- **Expected:** Exactly curv_abs * I
- **Tolerance:** 1e-6

### 3. ModerateStrain
- **Purpose:** Test typical use case
- **Checks:** Non-zero, non-identity, all finite values

### 4. LargeStrain
- **Purpose:** Stress test with large angular strain
- **Checks:** Numerical stability, finite values

### 5. NumericalAccuracy
- **Purpose:** Compare with finite difference approximation
- **Method:** Numerical Jacobian via finite differences
- **Tolerance:** 1e-4 max error

### 6. SymmetryProperties
- **Purpose:** Verify consistent computation with symmetric strain
- **Checks:** Finite values, structural properties

### 7. ScalingProperties
- **Purpose:** Verify behavior under length scaling
- **Checks:** Results change appropriately with curv_abs

### 8. Consistency
- **Purpose:** Verify deterministic computation
- **Method:** Multiple calls with same inputs
- **Tolerance:** 1e-15 (machine precision)

### 9. PerformanceBenchmark
- **Purpose:** Measure computation speed
- **Iterations:** 10,000
- **Expected:** < 10 microseconds per call

## Running the Tests

### Build Tests
```bash
cd /path/to/build
cmake -DUNIT_TEST=ON ..
make Cosserat_test
```

### Run Tests
```bash
./tests/Cosserat_test --gtest_filter="JacobianVerificationTest.*"
```

### Run Specific Test
```bash
./tests/Cosserat_test --gtest_filter="JacobianVerificationTest.NumericalAccuracy"
```

## Success Criteria

All tests must pass for Step 1 to be considered complete:

- [ ] SmallStrain passes
- [ ] ZeroStrain passes  
- [ ] ModerateStrain passes
- [ ] LargeStrain passes
- [ ] NumericalAccuracy passes (< 1e-4 error)
- [ ] SymmetryProperties passes
- [ ] ScalingProperties passes
- [ ] Consistency passes
- [ ] PerformanceBenchmark passes (< 10μs)

## Next Steps

After all tests pass:
1. Document any findings
2. Update implementation if needed
3. Move to Step 2: Input Validation
