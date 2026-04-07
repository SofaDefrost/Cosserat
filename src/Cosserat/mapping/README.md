# Cosserat Mappings

## Overview

The mapping module provides coordinate transformations between different representations of Cosserat rod elements. These mappings are essential for transferring mechanical states (positions, velocities, forces) between different parameterizations of the rod configuration.

## Key Features

- **DiscreteCosseratMapping**: Maps between local Cosserat coordinates (twist and bending) and world-frame rigid transformations
- **AdaptiveBeamMapping**: Maps between beam centerline and frame representations with adaptive discretization
- **CosseratNonLinearMapping2D**: Non-linear mapping for 2D Cosserat models
- **DifferenceMultiMapping**: Computes differences between mechanical states, useful for tracking relative displacements

## Usage Examples

### DiscreteCosseratMapping Example

```python
# Example of DiscreteCosseratMapping usage from a scene
cosseratInSofaFrameNode.addObject(
    "DiscreteCosseratMapping",
    curv_abs_input=[0, 10, 20, 30],  # section curve abscissa
    curv_abs_output=[0, 10, 20, 30],  # frames curve abscissa
    name="cosseratMapping",
    input1="@../cosseratCoordinate/cosserat_state",  # Vec3d local coordinates
    input2="@../rigid_base/cosserat_base_mo",        # Rigid3d base frame
    output="@./FramesMO",                           # Rigid3d output frames
    debug=0,
    radius=1.0,
)
```

### Using Mappings in the CosseratBase Prefab

```python
# From the CosseratBase prefab class
def _addCosseratFrame(self, framesF, curv_abs_inputS, curv_abs_outputF):
    cosseratInSofaFrameNode = self.rigidBaseNode.addChild("cosseratInSofaFrameNode")
    self.cosseratCoordinateNode.addChild(cosseratInSofaFrameNode)
    
    framesMO = cosseratInSofaFrameNode.addObject(
        "MechanicalObject",
        template="Rigid3d",
        name="FramesMO",
        position=framesF
    )

    cosseratInSofaFrameNode.addObject(
        "UniformMass", totalMass=self.beam_mass, showAxisSizeFactor="0"
    )

    cosseratInSofaFrameNode.addObject(
        "DiscreteCosseratMapping",
        curv_abs_input=curv_abs_inputS,
        curv_abs_output=curv_abs_outputF,
        name="cosseratMapping",
        input1=self.cosseratCoordinateNode.cosseratCoordinateMO.getLinkPath(),
        input2=self.rigidBaseNode.RigidBaseMO.getLinkPath(),
        output=framesMO.getLinkPath(),
        debug=0,
        radius=self.radius,
    )
    return cosseratInSofaFrameNode
```

## API Documentation

### DiscreteCosseratMapping

Maps between local Cosserat coordinates (torsion and bending) and rigid frames in world coordinates.

**Template Parameters**:
- `TIn1`: First input type, typically Vec3d (local Cosserat coordinates)
- `TIn2`: Second input type, typically Rigid3d (base frame)
- `TOut`: Output type, typically Rigid3d (frames along the rod)

**Data Fields**:
- `input1`: Input MechanicalObject for local Cosserat coordinates (Vec3d)
- `input2`: Input MechanicalObject for the base frame (Rigid3d)
- `output`: Output MechanicalObject for the resulting frames (Rigid3d)
- `curv_abs_input`: Curvilinear abscissa for input sections
- `curv_abs_output`: Curvilinear abscissa for output frames
- `radius`: Radius of the rod (for visualization)
- `debug`: Enable debug visualization and logging

### AdaptiveBeamMapping

Maps beam centerline points to frames with adaptive discretization based on curvature.

**Template Parameters**:
- `TIn`: Input type (typically Vec3d)
- `TOut`: Output type (typically Rigid3d)

**Data Fields**:
- `input`: Input MechanicalObject (centerline points)
- `output`: Output MechanicalObject (frames)
- `adaptiveDiscretization`: Enable dynamic refinement
- `maxError`: Maximum allowed error for adaptive discretization

### CosseratNonLinearMapping2D

Non-linear mapping for 2D Cosserat elements with large deformations.

**Template Parameters**:
- `TIn`: Input type
- `TOut`: Output type

**Data Fields**:
- `input`: Input MechanicalObject
- `output`: Output MechanicalObject
- `useQuat`: Use quaternions for rotation representation

## Integration with Other Components

Mappings are a crucial part of the Cosserat rod simulation pipeline:

1. **Force Fields**: Mappings transform forces computed by `BeamHookeLawForceField` in material coordinates to spatial coordinates
2. **Constraints**: Allow constraints to be applied in the most convenient coordinate system
3. **Visualization**: Enable rendering of rods with proper geometry and orientation
4. **Multi-model Integration**: Connect Cosserat models to other SOFA components or physics models

### Integration Example with Force Field

```python
def createScene(rootNode):
    # Setup base node and coordinate node
    rigid_base = _add_rigid_base(rootNode)
    bending_node = _add_cosserat_state(rootNode, bending_states, list_sections_length)
    
    # Create the mapping node that connects them
    cosserat_frame_node = _add_cosserat_frame(
        rigid_base,
        bending_node,
        cosserat_G_frames,
        section_curv_abs,
        frames_curv_abs,
        beam_radius,
    )
    
    # Now you can add constraints, rendering, or other components to any of these nodes
    return rootNode
```

## Related Documentation

- [Lie Groups Library](../liegroups/Readme.md) - Mathematical foundation for transformations
- [Force Fields](../forcefield/README.md) - Physical models that work with these mappings
- [Python Base Class](../../examples/python3/cosserat/cosserat.py) - Prefab for easy construction of Cosserat rods

