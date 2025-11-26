# HookeSeratBaseMapping - Incremental Improvement Plan

**Branch:** `feature/hookeserat-incremental-improvements`  
**Status:** Step 0 - Baseline Established  
**Last Updated:** 2025-11-26

---

## Overview

This document outlines the incremental improvement plan for `HookeSeratBaseMapping`. After removing untested advanced features, we're rebuilding with **proper validation at each step**.

---

## Current Baseline (Step 0) ✅

### What Works Now
- Core Lie groups integration (SE3/SO3)
- Basic `SectionInfo` with transformations
- Basic `FrameInfo` for beam frames  
- `computeTangExpImplementation()` for Jacobians
- Simple strain management with `Vector3` + `TangentVector`

### What Was Removed
All advanced features removed due to lack of validation:
- BeamStateEstimator, StrainState Bundle, BeamTopology
- JacobianStats, ML integration, validation methods
- Interpolation methods

**Reason:** Not tested step-by-step → bugs → removed → will re-add incrementally

---

## Improvement Steps

### Step 1: Verify Jacobian Implementation 🔴 **NEXT**

**Goal:** Ensure `computeTangExpImplementation()` is correct

**Tasks:**
1. Compare with `SE3::rightJacobian()` if available
2. Add unit tests (small/moderate/large strains)
3. Verify numerical accuracy (< 1e-6)
4. Document implementation choice

**Success Criteria:**
- [ ] Tests pass for all strain ranges
- [ ] Numerical accuracy < 1e-6
- [ ] Performance acceptable
- [ ] Documented

**Time:** 1-2 days

---

### Step 2: Add Input Validation 🟡

**Goal:** Prevent invalid inputs

**Tasks:**
1. Validate `setStrain()` - size & finite values
2. Validate `setLength()` - positive
3. Validate `setTransformation()` - valid SE3
4. Add clear error messages
5. Unit tests for edge cases

**Success Criteria:**
- [ ] Invalid inputs caught
- [ ] Clear error messages
- [ ] No performance regression
- [ ] Tests pass

**Time:** 1 day

---

### Step 3: Type-Safe Strain (Bundle) 🟡

**Goal:** Use `Bundle<SO3, RealSpace>` for type safety

**Tasks:**
1. Add `StrainState` type alias
2. Add member to `SectionInfo`
3. Keep legacy for compatibility
4. Add conversion methods
5. Test both interfaces

**Success Criteria:**
- [ ] Bundle compiles
- [ ] Conversions work
- [ ] Legacy works
- [ ] No regression
- [ ] Tests pass

**Time:** 2-3 days

---

### Step 4: Basic Interpolation 🟢

**Goal:** Add `lerp()` for trajectories

**Tasks:**
1. Implement `SectionInfo::lerp()`
2. Test boundaries (t=0, 0.5, 1)
3. Test invalid inputs
4. Document with examples

**Success Criteria:**
- [ ] Smooth interpolation
- [ ] Boundaries work
- [ ] Tests pass
- [ ] Example provided

**Time:** 1 day

---

### Step 5: Distance Metrics 🟢

**Goal:** Measure configuration similarity

**Tasks:**
1. Implement `distanceTo()` using SE3
2. Implement `strainDistance()`
3. Test properties (non-negative, identity, triangle inequality)
4. Document

**Success Criteria:**
- [ ] Metrics work
- [ ] Properties verified
- [ ] Tests pass

**Time:** 1 day

---

### Step 6: Geometry Validation 🟢

**Goal:** Detect invalid configurations

**Tasks:**
1. Implement `validateSectionProperties()`
2. Implement `checkInterSectionContinuity()`
3. Test valid/invalid configs
4. Document rules

**Success Criteria:**
- [ ] Detects invalid configs
- [ ] No false positives
- [ ] Tests pass

**Time:** 2 days

---

## Testing Strategy

### Per Step
1. **Unit tests** - Feature in isolation
2. **Integration tests** - With existing code
3. **Regression tests** - Nothing broke
4. **Performance tests** - If relevant

### Test Files
```
tests/
  ├── test_jacobian.cpp
  ├── test_validation.cpp
  ├── test_strain_bundle.cpp
  ├── test_interpolation.cpp
  ├── test_distance.cpp
  └── test_geometry_validation.cpp
```

---

## Success Metrics

**Per Step:**
- ✅ Compiles without errors/warnings
- ✅ All tests pass
- ✅ Documentation updated
- ✅ Code reviewed

**Overall:**
- Stable (no crashes)
- Correct (tests pass)
- Fast (no regression)
- Clear (documented)
- Useful (features used)

---

## Next Actions

1. ✅ Create branch `feature/hookeserat-incremental-improvements`
2. ✅ Update documentation
3. 🔴 **START:** Step 1 - Verify Jacobian
4. 🔴 Create test file
5. 🔴 Implement & run tests

---

## Notes

- Update this doc as we progress
- Mark completed steps with ✅
- Document issues/deviations
- Capture lessons learned
