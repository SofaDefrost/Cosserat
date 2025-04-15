# Cosserat Geometric Stiffness Engines

## Overview

The engine module provides specialized computational components for efficient calculation of geometric stiffness matrices in Cosserat rod simulations. Geometric stiffness (also known as stress stiffness or initial stress stiffness) arises from the change in orientation of forces due to deformation and is essential for accurate modeling of slender structures under large deformations.

In Cosserat rod theory, accurate computation of geometric stiffness is crucial for:
- Capturing correct buckling behavior
- Ensuring stability in large deformation scenarios
- Modeling pre-stressed structures
- Achieving convergence in highly nonlinear simulations

## Available Engines

### GeometricStiffnessEngine

The primary engine for computing geometric stiffness terms efficiently. This component accelerates the calculation of the geometric stiffness matrix, which would otherwise be computationally expensive, especially for systems with many degrees of freedom.

### Specialized Variants

- **AdaptiveGeometricStiffnessEngine**: Dynamically adjusts computation precision based on deformation state
- **ThreadedGeometricStiffnessEngine**: Exploits multi-threading for parallel computation of stiffness terms

## Usage Examples

### Basic Usage with BeamHookeLawForceField

```python
# Example of using GeometricStiffnessEngine with a beam force field
node.addObject('BeamHookeLawForceField',
              youngModulus=1e6,
              poissonRatio=0.3,
              radius=0.01)

# Add the geometric stiffness engine
node.addObject('GeometricStiffnessEngine',
              forceField='@BeamHookeLawForceField',
              computeGlobalMatrix=True,
              debugLevel=0)
```

### Advanced Configuration

```python
# Example with advanced configuration
node.addObject('GeometricStiffnessEngine',
              forceField='@BeamHookeLawForceField',
              computeGlobalMatrix=True,
              useMultiThreading=True,
              threadCount=8,
              updateStiffnessMatrix=True,
              debugLevel=1)
```

### Integration in a Complete Simulation

```python
def createScene(rootNode):
    # Setup solver
    rootNode.addObject('EulerImplicitSolver', rayleighStiffness=0.1, rayleighMass=0.1)
    rootNode.addObject('SparseLDLSolver', template='CompressedRowSparseMatrixd')
    
    # Create Cosserat rod elements
    beam = rootNode.addChild('beam')
    beam.addObject('MechanicalObject', template='Vec3d')
    
    # Add force field
    ff = beam.addObject('BeamHookeLawForceField',
                       youngModulus=1e6,
                       poissonRatio=0.3,
                       radius=0.01)
    
    # Add geometric stiffness engine
    beam.addObject('GeometricStiffnessEngine',
                  forceField=ff.getLinkPath(),
                  computeGlobalMatrix=True)
    
    return rootNode
```

## API Documentation

### GeometricStiffnessEngine

Computes and manages geometric stiffness matrices for Cosserat rod elements.

**Template Parameters**:
- `DataTypes`: The mechanical state data type (typically Vec3d or Rigid3d)

**Data Fields**:
- `forceField`: Link to the associated force field (e.g., BeamHookeLawForceField)
- `computeGlobalMatrix`: If true, assembles the global stiffness matrix instead of just local matrices
- `useMultiThreading`: Enables parallel computation of the stiffness matrix
- `threadCount`: Number of threads to use when multi-threading is enabled
- `updateStiffnessMatrix`: Determine when to update the stiffness matrix (always, or only when needed)
- `debugLevel`: Level of debugging information to output (0-3)

**Methods**:
- `init()`: Initializes the engine and establishes connections with associated components
- `reinit()`: Reinitializes the engine, typically called after a change in configuration
- `computeGeometricStiffness()`: Core method that computes the geometric stiffness
- `handleEvent()`: Processes system events such as time step changes

## Integration Guidelines

The geometric stiffness engines are designed to work seamlessly with other components of the Cosserat plugin:

1. **Force Fields**: Connect the engine to a compatible force field like BeamHookeLawForceField
2. **Solvers**: Ensure your solver can utilize the geometric stiffness contribution
3. **Scene Graph**: Place the engine in the same node as the force field

### Recommended Integration Steps

1. Create your Cosserat rod elements with appropriate mappings
2. Add a compatible force field (BeamHookeLawForceField)
3. Add the GeometricStiffnessEngine, linking it to the force field
4. Configure the solver to handle the additional stiffness terms

## Performance Considerations

Geometric stiffness computation can be computationally intensive. Here are guidelines to optimize performance:

### Multi-threading

For large systems (>50 elements), enable multi-threading:
```python
engine.useMultiThreading = True
engine.threadCount = 8  # Adjust based on your CPU cores
```

### Selective Updates

If your simulation includes phases with minimal deformation:
```python
engine.updateStiffnessMatrix = False  # Manual control
# Later when needed:
engine.updateStiffnessMatrix = True   # Update once
```

### Memory Optimization

For very large systems, consider:
```python
engine.storeGlobalMatrix = False  # Compute on demand instead of storing
```

### Benchmarks

Typical performance improvements:
- Single-threaded vs. multi-threaded (8 cores): 4-6x speedup
- Adaptive computation: 2-3x speedup in scenes with localized deformation
- Memory reduction with on-demand computation: Up to 70% for large systems

## Related Documentation

- [Force Fields](../forcefield/README.md) - Force fields that use geometric stiffness
- [Mappings](../mapping/README.md) - Transformations that may affect geometric stiffness computation
- [Lie Groups Library](../liegroups/Readme.md) - Mathematical foundation for rod kinematics
- [Performance Benchmarks](../../docs/text/benchmarks.md) - Detailed performance analysis

