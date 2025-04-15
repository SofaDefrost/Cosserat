# Cosserat Constraints

## Overview

The constraint module provides specialized constraint implementations for Cosserat rod elements. These constraints enable the definition of boundary conditions, contact handling, and internal constraints that are particularly important for accurate simulation of slender structures like beams, cables, and tubes.

## Key Features

- **UnilateralPlaneConstraint**: Constrains rod elements to remain on one side of a plane, with contact handling
- **Fixed-end Constraints**: Immobilize rod ends or specific points
- **Contact Handling**: Specialized contact models for rod-surface interaction
- **Self-collision**: Detection and handling of self-intersections in complex rod configurations

## Usage Examples

### Using Fixed Constraints

```python
# Example of constraining the base of a Cosserat rod
rigid_base_node.addObject(
    "RestShapeSpringsForceField",
    name="spring",
    stiffness=1e8,
    angularStiffness=1.0e8,
    external_points=0,
    points=0,
    template="Rigid3d"
)
```

### Unilateral Plane Constraint Example

```python
# Example of adding a unilateral plane constraint
node.addObject(
    "UnilateralPlaneConstraint",
    plane=[0, 1, 0, -10],  # Plane equation: ax + by + cz + d = 0
    indices=[0, 1, 2, 3],  # Indices of constrained points
    activateDebugOutput=False
)
```

### Collision Model for Rod-Environment Interaction

```python
# From the CosseratBase prefab class
def addCollisionModel(self):
    tab_edges = generate_edge_list(self.frames3D)
    return addEdgeCollision(self.cosseratFrame, self.frames3D, tab_edges)

def _addPointCollisionModel(self, nodeName="CollisionPoints"):
    tab_edges = generate_edge_list(self.frames3D)
    return addPointsCollision(
        self.cosseratFrame, self.frames3D, tab_edges, nodeName
    )
```

## API Documentation

### UnilateralPlaneConstraint

Constrains points to remain on one side of a plane, handling contact and collision responses.

**Template Parameters**:
- `DataTypes`: The type of the constrained points (typically Vec3d or Rigid3d)

**Data Fields**:
- `plane`: Plane equation coefficients (a, b, c, d) where ax + by + cz + d = 0
- `indices`: Indices of points to be constrained
- `activateDebugOutput`: Enable visualization and logging of constraint forces
- `bilateral`: If true, constrains points to stay exactly on the plane rather than on one side

### RestShapeSpringsForceField

While technically a force field, it's commonly used to constrain rod endpoints in Cosserat simulations.

**Template Parameters**:
- `DataTypes`: Type of points to constrain (Rigid3d for rod ends)

**Data Fields**:
- `stiffness`: Translational stiffness of the constraint
- `angularStiffness`: Rotational stiffness (important for Cosserat rods)
- `external_points`: Indices of external reference points
- `points`: Indices of points to be constrained
- `template`: Template type for the constraint

## Integration with Other Components

Constraints work closely with other components in the Cosserat plugin:

1. **Force Fields**: Constraints provide boundary conditions for rod mechanics
2. **Mappings**: Constraints can be applied in different coordinate systems
3. **Collision Models**: Constraints work with collision detection for contact handling
4. **Solvers**: Constraints are resolved by the constraint solver, which needs to be configured in the scene

### Integration Example with Solver

```python
def createScene(rootNode):
    # Configure animation and solver components for constraints
    rootNode.addObject("FreeMotionAnimationLoop")
    rootNode.addObject("GenericConstraintSolver", tolerance=1e-5, maxIterations=5e2)
    
    solverNode = rootNode.addChild("solverNode")
    solverNode.addObject(
        "EulerImplicitSolver", rayleighStiffness="0.2", rayleighMass="0.1"
    )
    solverNode.addObject(
        "SparseLDLSolver", name="solver", template="CompressedRowSparseMatrixd"
    )
    solverNode.addObject("GenericConstraintCorrection")
    
    # Create the Cosserat model
    cosserat = solverNode.addChild(CosseratBase(parent=solverNode, params=Params))
    
    # Fix the base
    cosserat.rigidBaseNode.addObject(
        "RestShapeSpringsForceField",
        name="spring",
        stiffness=1e8,
        angularStiffness=1.0e8,
        external_points=0,
        points=0,
        template="Rigid3d"
    )
    
    # Add collision model for constraints
    cosserat.addCollisionModel()
    
    return rootNode
```

## Related Documentation

- [Lie Groups Library](../liegroups/Readme.md) - Mathematical foundation for rigid body constraints
- [Force Fields](../forcefield/README.md) - Physical models that are constrained
- [Mappings](../mapping/README.md) - Transformations between different constraint representations
- [Python Base Class](../../examples/python3/cosserat/cosserat.py) - Prefab for easy construction of constrained Cosserat rods

