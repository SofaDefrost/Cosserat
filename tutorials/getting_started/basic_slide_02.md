---
title: Basic Slide 02
---

## Let's the beam fall under the influence of gravity.

```python
    # Add gravity
    root_node.gravity = [0, -9.81, 0]  # Add gravity!

    # Configure time integration and solver
    solver_node = root_node.addChild("solver")

    # Damping parameter for dynamics
    v_damping_param: float = 8.e-1 
    
    solver_node.addObject(
        "EulerImplicitSolver",
        rayleighStiffness="0.0",
        rayleighMass="0.0",
        vdamping=v_damping_param,  # Damping parameter for dynamics
    )
    solver_node.addObject("SparseLDLSolver", name="solver")
```

---

## Explanation of the Code

- The `gravity` parameter is set to `[0, -9.81, 0]`, which simulates the effect of gravity acting downwards in the y-direction.
- The `EulerImplicitSolver` is used to integrate the dynamics of the beam.
- The `rayleighStiffness` and `rayleighMass` parameters are set to `0.0`, meaning that no additional stiffness or mass damping is applied to the system.
- The `vdamping` parameter is set to a value of 0.8, which introduces damping to the system, helping to stabilize the simulation and reduce oscillations.
- The `SparseLDLSolver` is used to solve the linear system of equations that arise during the simulation.

---

