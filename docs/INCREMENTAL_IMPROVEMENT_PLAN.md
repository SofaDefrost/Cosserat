# CosseratBeamGeometry - Incremental Improvement Plan

**Status:** Steps 2 to 7 open (Step 1 done)
**Last updated:** 2026-11-11, checked against `src/Cosserat/mapping/CosseratBeamGeometry.h`

The per-section and per-frame data (`SectionInfo`, `FrameInfo`, `computeTangExpImplementation`) now live in the
`CosseratBeamGeometry<TStrain>` mixin, used by `CosseratGeometryMapping`. Each remaining step is validated before
the next one starts.

---

## Remaining steps

### Step 2: Complete input validation 🟡

`SectionInfo::setLength()` already rejects non-positive values. Still missing:

1. `setStrain()`: check the size and that the values are finite
2. `setTransformation()`: check that the SE3 is valid
3. Clear error messages
4. Unit tests for edge cases

**Done when:** invalid inputs are caught, no performance regression, tests pass.

---

### Step 3: Type-safe strain (Bundle) 🟡

No `StrainState` type in `SectionInfo` yet.

1. Add a `StrainState` alias (`Bundle<SO3, RealSpace>`)
2. Add it as a member of `SectionInfo`, keep the legacy interface
3. Add conversion methods
4. Test both interfaces

**Done when:** the Bundle compiles, conversions work, the legacy interface is unchanged, tests pass.

---

### Step 4: Basic interpolation 🟢

No `lerp()` on `SectionInfo`.

1. Implement `SectionInfo::lerp()`
2. Test boundaries (t = 0, 0.5, 1) and invalid inputs
3. Document with an example

**Done when:** interpolation is smooth, boundaries work, tests pass.

---

### Step 5: Strain distance and tests 🟢

`SectionInfo::distanceTo()` exists (SE3 distance). Still missing:

1. `strainDistance()`
2. Tests of the metric properties (non-negative, identity, triangle inequality) for both distances

---

### Step 6: Geometry validation 🟢

`validateSectionProperties()` and `validateJacobianAccuracy()` exist. Still missing:

1. `checkInterSectionContinuity()`
2. Tests on valid and invalid configurations, no false positives
3. Document the rules

---

### Step 7: Lie-group features (optional) 🟢

1. Replace the hand-written trigonometric series of `computeTangExpImplementation()` by the Lie-group Jacobian
   (`SE3::rightJacobian`), if `test_jacobian_verification.cpp` shows the same values to 1e-6
2. Uncertainty propagation on sections (`GaussianOnManifold<SE3>`, `CosseratUncertaintyPropagator`)
3. BCH formula for the composition of sections

---

## Tests

Only `test_jacobian_verification.cpp` exists. To create, one per step:

```
test_validation.cpp
test_strain_bundle.cpp
test_interpolation.cpp
test_distance.cpp
test_geometry_validation.cpp
```

Per step: unit tests (feature in isolation), integration tests (with existing code), regression tests, and
performance tests where relevant.

## Success criteria (per step)

- Compiles without errors or warnings
- All tests pass
- Documentation updated
- Code reviewed

## Next action

Start with Step 2 (`setStrain()` and `setTransformation()` validation), then create `test_validation.cpp`.
