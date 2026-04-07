---
title: Tutorial 05 - Constraints
date: 2025-06-28
---

# Tutorial 05: Constraints and Boundary Conditions

In this tutorial, we'll learn how to apply constraints to a Cosserat beam. Constraints are essential for building complex models, as they allow you to fix parts of your model in place or connect them to other objects.

[[./05_constraints_and_boundary_conditions.py]]

## The Goal: A Simple Bridge

We will create a simple bridge by creating a single beam and fixing *both* of its ends. This will show how the beam deforms under its own weight when supported at both ends.

## Key Component: `FixedConstraint`

The `FixedConstraint` component is used to lock specific degrees of freedom of a `MechanicalObject`. In our case, we will use it to fix the position and orientation of the beam's tip.

```python
# --- CONSTRAINT ---
# Fix the tip of the beam to create a bridge
tip_frame_index = beam_geometry.get_number_of_frames() - 1

# Add a FixedConstraint to the last frame of the beam.
# This will lock its position and orientation.
frame_node.addObject(
    "FixedConstraint",
    name="bridgeConstraint",
    indices=[tip_frame_index], # Index of the frame to constrain
)
```

### Code Explanation

- `tip_frame_index`: We get the index of the last frame of the beam. Remember that indices are 0-based.
- `frame_node.addObject("FixedConstraint", ...)`: We add the `FixedConstraint` to the same node that contains the beam's frames.
- `indices=[tip_frame_index]`: This is the most important parameter. It tells the constraint which point(s) to affect. Here, we are only constraining the very last frame.

By default, `FixedConstraint` locks all 6 degrees of freedom (3 translations, 3 rotations). You can selectively constrain specific axes if needed, but this is the simplest way to create a fixed connection.

Now, when you run the simulation, you will see the beam start straight and then sag in the middle due to gravity, forming a bridge shape.
