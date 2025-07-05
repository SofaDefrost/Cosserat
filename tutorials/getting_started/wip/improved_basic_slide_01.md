---
title: "Tutorial 01: Understanding Beam Geometry and Discretization"
date: 2025-06-20
---

# Tutorial 01: Understanding Beam Geometry and Discretization

## Key Concepts

In this tutorial, we'll explore:
1. The effect of different beam discretization parameters
2. How sections (physics) and frames (visualization) are related
3. How to create multiple beams with different parameters

---

## Sections vs. Frames

In Cosserat rod modeling, we make a distinction between:

- **Sections**: Discrete segments where physics (strain) is calculated
- **Frames**: Points along the beam where we visualize the rod

The number of sections affects physical accuracy, while the number of frames affects visual smoothness.

---

## Example Code Structure

```python
# Define beam geometry parameters
beam_geometry_params = BeamGeometryParameters(
    beam_length=30.0,  # Total beam length
    nb_section=3,      # Number of sections for physics
    nb_frames=12       # Number of frames for visualization
)

# Create geometry object
beam_geometry = CosseratGeometry(beam_geometry_params)

# Create beam components
base_node = create_rigid_base(root_node, node_name="rigid_base1")
bending_node = create_cosserat_state(root_node, beam_geometry, 
                                   node_name="cos_coord1")
create_cosserat_frame(base_node, bending_node, beam_geometry, 
                     node_name="cosserat_frame_node1")
```

---

## Comparing Different Discretizations

In this tutorial, we'll create two beams with the same physics (number of sections) but different visualizations (number of frames):

### Beam 1: 3 sections, 3 frames
- Physical accuracy: Lower (3 sections)
- Visual smoothness: Lower (3 frames)
- Computation: Faster

### Beam 2: 3 sections, 12 frames
- Physical accuracy: Same as Beam 1 (3 sections)
- Visual smoothness: Higher (12 frames)
- Computation: Slightly slower (only visualization is affected)

---

## Visual Comparison

![Beam Comparison](../../docs/images/example_beam_comparison.png)

Even though both beams have the same physical behavior (same number of sections), the second beam appears smoother due to the higher number of frames.

---

## Running the Example

Open and run `tutorial_01_basic_beam.py` to see both beams side by side:

```bash
runSofa tutorial_01_basic_beam.py
```

Try experimenting with different combinations of sections and frames to see how they affect performance and appearance.

---

## Understanding the Code

Let's look at the key parts of the code:

```python
# First beam: 3 sections, 3 frames
beam_geometry_params = BeamGeometryParameters(
    beam_length=30.0, nb_section=3, nb_frames=3
)
beam_geometry1 = CosseratGeometry(beam_geometry_params)

# Second beam: 3 sections, 12 frames
beam_geometry_params2 = BeamGeometryParameters(
    beam_length=30.0, nb_section=3, nb_frames=12
)
beam_geometry2 = CosseratGeometry(beam_geometry_params2)
```

The physical behavior is identical, but the visual representation is different.

---

## Key Takeaways

1. **Physics vs. Visualization**: Sections control physics, frames control visualization
2. **Efficiency**: Choose the right balance between accuracy and performance
3. **CosseratGeometry**: The helper class makes it easy to experiment with different configurations

Next, we'll add dynamics to our beams by introducing gravity and time integration.
