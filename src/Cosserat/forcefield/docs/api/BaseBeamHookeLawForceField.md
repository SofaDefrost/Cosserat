# BaseBeamHookeLawForceField API Documentation

## Overview
Base class implementation for beam elements using Hooke's law in the Cosserat model. Supports both uniform and variant cross-sections.

## Class Features

### Cross Section Types
- Circular (tube with external and internal radius)
- Rectangular (with lengthY and lengthZ dimensions)

### Configuration Options
- Uniform or variant sections
- Material properties (Young's modulus, Poisson ratio)
- Direct inertia parameters
- Multi-threading support

## Public Interface

### Constructor Parameters
- `crossSectionShape`: Shape of the cross-section ("circular" or "rectangular")
- `youngModulus`: Material stiffness
- `poissonRatio`: Material compressibility
- `length`: List of beam section lengths

### Main Methods
- `init()`: Initialize force field and validate parameters
- `addForce()`: Compute and add forces
- `getRadius()`: Get beam external radius
- `getPotentialEnergy()`: Calculate strain energy

### Configuration Data
- `d_variantSections`: Enable variable section properties
- `d_useInertiaParams`: Use direct inertia parameters
- `d_useMultiThreading`: Enable parallel computation

## Implementation Details

See implementation documentation for internal details about:
- Stiffness matrix computation
- Force calculation methods
- Section property validation
- Multi-threading implementation

## Usage Example

```cpp
// Create a beam with uniform circular cross-section
auto* beam = new BaseBeamHookeLawForceField();
beam->d_crossSectionShape = "circular";
beam->d_radius = 0.01;  // 10mm radius
beam->d_youngModulus = 1e9;  // 1 GPa
beam->d_poissonRatio = 0.3;
beam->init();
```

## See Also
- Integration_Plan.md for Lie Group implementation details
- BeamHookeLawForceField for standard implementation
- BeamHookeLawForceFieldRigid for rigid body variant
