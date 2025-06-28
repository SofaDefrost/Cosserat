---
title: "Tutorial 02: Dynamic Simulation with Gravity"
date: 2025-06-20
---

# Tutorial 02: Dynamic Simulation with Gravity

## Adding Dynamics to Our Beam

In this tutorial, we'll add dynamic behavior to our Cosserat beam by:
1. Adding gravity
2. Configuring time integration
3. Adding mass to our beam
4. Comparing different beam discretizations under gravity

---

## Setting Up the Solver

To enable dynamic simulation, we need to add proper solvers:

```python
# Add gravity
root_node.gravity = [0, -9.81, 0]  # Standard gravity in the Y direction

# Configure time integration and solver
solver_node = root_node.addChild("solver")

# Damping parameter for stable dynamics
v_damping_param: float = 8.e-1 

solver_node.addObject(
    "EulerImplicitSolver",
    rayleighStiffness="0.0",
    rayleighMass="0.0",
    vdamping=v_damping_param,  # Damping helps stabilize the simulation
)
solver_node.addObject("SparseLDLSolver", name="solver")
```

---

## Understanding the Solver Components

### EulerImplicitSolver
- Performs implicit integration of the dynamics
- More stable than explicit methods, especially for stiff systems
- `vdamping` parameter adds artificial damping to stabilize the simulation
- `rayleighStiffness` and `rayleighMass` are set to 0, meaning no Rayleigh damping

### SparseLDLSolver
- Efficiently solves the linear system of equations that arise during simulation
- Uses a sparse matrix representation for better performance
- LDL decomposition is a variant of Cholesky factorization suited for these problems

---

## Adding Mass to the Beam

For dynamics to work, we need to add mass to our beam:

```python
# Create cosserat frame with mass (important for dynamics!)
frame_node = create_cosserat_frame(
    base_node, bending_node, beam_geometry, beam_mass=5.0
)
```

The mass is uniformly distributed across all frames of the beam.

---

## Experiment: Different Numbers of Sections

In this tutorial, we'll compare two beams:

### Beam 1: 30 sections, 30 frames
- High physical accuracy (30 sections)
- Smooth visualization (30 frames)
- More computationally intensive

### Beam 2: 3 sections, 30 frames
- Lower physical accuracy (3 sections)
- Same visual smoothness (30 frames)
- More efficient computation

---

## Expected Results

When we run the simulation:

1. Both beams will fall due to gravity
2. The first beam (30 sections) will bend more realistically, with a smooth curve
3. The second beam (3 sections) will bend in a more segmented way

This demonstrates that while frames affect visualization, the number of sections affects the physical behavior.

---

## Running the Example

Run the example with:

```bash
runSofa tutorial_02_basic_beam.py
```

Press the "Animate" button to start the simulation and observe how the beams fall under gravity.

---

## Key Takeaways

1. **Time Integration**: Proper solver setup is crucial for stable dynamics
2. **Mass Distribution**: Mass must be added for dynamic simulation
3. **Discretization Effects**: 
   - Number of sections affects physical accuracy
   - For realistic bending under loads, more sections are needed

In the next tutorial, we'll explore how to apply external forces to our beam.
