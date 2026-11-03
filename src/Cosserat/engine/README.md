# Cosserat engines

Components of `src/Cosserat/engine/`. This directory used to document a `GeometricStiffnessEngine`
(with `AdaptiveGeometricStiffnessEngine` and `ThreadedGeometricStiffnessEngine` variants); those
classes **do not exist in this repository** (no source file, not registered), so that description was
removed. Geometric stiffness is still a TODO for the staggered (Painless) beam.

## Staggered ("Painless") discretisation

Beam with N+1 node positions (R³) and N segment orientations (SO(3), at edge midpoints), after
Romanyà-Serrasolsas, Casafranca, Otaduy, *Painless Differentiable Rotation Dynamics*, ACM ToG /
SIGGRAPH 2025 (DOI 10.1145/3730944).

| Component | File | Role |
|---|---|---|
| `CosseratIntrinsicState` | `CosseratIntrinsicState.h/.cpp` | `BaseObject` (deliberately **not** a `MechanicalState`) holding `positions` (N+1 Vec3d) and `orientations` stored as **rotation vectors ω = log(R)** (identity = `0 0 0`). Rest lengths are taken from the initial positions at `init()`. Also provides the cached strains `Γ_i = R_iᵀ(x_{i+1}−x_i)/h_i − e₁` and `Ω_i = log(R_{i−1}ᵀR_i)/h̃_i`, and `getInverseLieJacobian(ω)` = J_r⁻¹(ω). |
| `CosseratTopologyBuilder` | `CosseratTopologyBuilder.h/.cpp` | Builds a **straight** beam (through `CosseratBeamGeometryEngine`) and writes it into the linked state at `init()`, overwriting any positions/orientations given in the scene. Computes EA, GA, EIy, EIz, GJ as outputs; they are **not** linked to the force field (pass them explicitly, see `tutorials/advanced/staggered/_common.py`). |
| `CosseratBeamGeometryEngine` | `CosseratBeamGeometryEngine.h/.cpp` | Static helpers: `buildStraightConfiguration`, `buildConfiguration` (N frames; rotation applied *before* the step; "curvatures" are per-segment rotation increments, not rates), `buildCurvedConfiguration` (N+1 frames; rotation applied *after* the step — a diverging convention, do not mix), `calculateStiffnessMatrix` (circular or rectangular section). |
| `PainlessBeamForceField` | `PainlessBeamForceField.h/.cpp` | Elastic forces on nodes and torques on segments (body frame), exact negative gradient of the discrete energy (checked against finite differences). Results are written to the `nodalForces` / `segmentTorques` Data on `AnimateBeginEvent`; an external (Python) integrator must run on `AnimateEndEvent`. `addDForce` / `addKToMatrix` are stubs and `computeDForces` is **not** the complete tangent stiffness (cross position↔rotation blocks, shear coupling term and geometric stiffness are missing) — see `IMPLICIT_ROADMAP.md`. `printEvery` (default 0) enables a periodic summary; always pass EA/GA/EIy/EIz/GJ explicitly (the defaults are deliberately wrong). |

The matching visualisation component is `mapping/StaggeredCosseratMapping` (Rigid3d frames at nodes and
segment midpoints + drawing; not a SOFA `Mapping`). Example scenes: `tutorials/advanced/staggered/`.

## Other engines (not Painless)

| Component | File | Notes |
|---|---|---|
| `PointsManager` | `PointsManager.h/.inl/.cpp` | Keyboard-driven management of point sets on a point-set topology (registered in `initCosserat.cpp`). |
| `ProjectionEngine` | `ProjectionEngine.h/.inl/.cpp` | `DataEngine` working with constraint data (registered in `initCosserat.cpp`; see the header). |
| `CosseratState` | `CosseratState.h/.inl/.cpp` | `TemperatureState` template derived from `MechanicalObject`; **not** part of the build (not in `CMakeLists.txt`, not registered). |
| `engineQuat` | `engineQuat.h/.cpp` | Empty placeholder class; **not** part of the build. |

## Related documentation

- `../forcefield/README.md` — PCS force fields
- `../mapping/README.md` — mappings
- `../../liegroups/` — Lie group library (`[φ ; ρ]` ordering: angular part first)
