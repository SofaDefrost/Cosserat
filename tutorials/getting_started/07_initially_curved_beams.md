---
title: Tutorial 07 - Initially Curved Beams
date: 2025-06-28
---

# Tutorial 07: Modeling Initially Curved Beams

Not all objects are straight in their resting state. This tutorial explains how to model a Cosserat beam that has an intrinsic, stress-free curvature.

[[./07_initially_curved_beams.py]]

## The Goal: A Curved Arch

We will create a beam that forms a semi-circular arch at rest. When gravity is applied, it will deform from this curved shape, not from a straight line.

## The Key: `bendingState`

The `CosseratGeometry` calculates the initial positions of the frames (`frames`) based on the initial curvature (`bendingState`). So far, we have mostly used a `bendingState` of all zeros, resulting in a straight beam.

To create a curved beam, we must provide a non-zero `bendingState` that describes the curvature along the beam's length.

```python
# --- Define Initial Curvature ---
# To create a curved beam, we define a non-zero resting curvature.
# Here, we'll create a simple circular arch by applying a constant
# curvature around the Z-axis for each section.
num_sections = beam_geometry.get_number_of_sections()
# The total curvature will define a semi-circle (pi radians)
total_angle = np.pi
curvature_per_section = total_angle / beam_geometry.get_beam_length()

custom_bending_states = [[0, 0, curvature_per_section] for _ in range(num_sections)]
```

### Code Explanation

- **`bendingState`**: This is a list of 3D vectors, where each vector `[κ₁, κ₂, κ₃]` represents the material curvature for one section of the beam.
- **`total_angle`**: We want our beam to form a semi-circle, which spans an angle of π radians.
- **`curvature_per_section`**: We calculate the constant curvature required at each section to achieve the total desired angle over the entire beam length. We apply this curvature around the Z-axis (`κ₃`).
- **`custom_bending_states`**: We create a list where every section has the same calculated curvature. This results in a uniform, circular arch.

When you pass this `custom_bending_states` to the `_add_cosserat_state` function, the `DiscreteCosseratMapping` will automatically compute the correct curved initial positions for the visual frames.
