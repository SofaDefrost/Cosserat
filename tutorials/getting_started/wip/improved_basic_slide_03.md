---
title: "Tutorial 03: Applied Forces and Interactions"
date: 2025-06-20
---

# Tutorial 03: Applied Forces and Interactions

In this tutorial, we'll learn how to apply different types of forces to our Cosserat beam.

---

## Types of Forces in SOFA

SOFA provides several ways to apply forces to objects:

1. **ConstantForceField**: Applies a constant force to specified points
2. **Controller-based forces**: Custom forces computed and applied each frame
3. **Interaction forces**: Forces resulting from interactions with other objects

---

## Force Type 1: Constant Force

The simplest way to apply a force is using a `ConstantForceField`:

```python
# Apply a constant force to the tip of the beam
tip_frame_index = beam_geometry.get_number_of_frames() - 1  # Last frame
applied_force = [-10.0, 0.0, 0.0, 0.0, 0.0, 0.0]  # Force in -X direction

frame_node.addObject(
    "ConstantForceField",
    name="tipForce",
    indices=[tip_frame_index],
    forces=[applied_force],
    showArrowSize=1,  # Visualize the force
)
```

This applies a constant force in the negative X direction to the last frame of the beam.

---

## Force Type 2: Controller-Based Forces

For more complex forces, we can use a controller that updates forces every frame:

```python
# Define force components
def compute_force(self):
    with self.forceNode.forces.writeable() as force:
        vec = [0., 0., 0., 0., self.forceCoeff / sqrt(2), self.forceCoeff / sqrt(2)]
        for i, v in enumerate(vec):
            force[0][i] = v
```

This example from our `ForceController` class applies a force that increases over time in the Y and Z angular directions.

---

## Force Type 3: Orientation-Based Forces

We can also apply forces that are based on the current orientation of the beam:

```python
def compute_orthogonal_force(self):
    # Get the last frame position and orientation
    position = self.frames.position[self.nb_frames]  
    orientation = Quat(
        position[3], position[4], position[5], position[6]
    )
    
    # Calculate force orthogonal to the beam's axis
    with self.forceNode.forces.writeable() as force:
        vec = orientation.rotate([0.0, self.forceCoeff * 5.0e-2, 0.0])
        for count in range(3):
            force[0][count] = vec[count]
```

This applies a force that remains perpendicular to the beam's axis as it bends.

---

## Force Type 4: End-Effector Control

For more complex control, we can create an end-effector that interacts with the beam:

```python
# Create a rigid body to control the end effector
tip_controller = root_node.addChild('tip_controller')
controller_state = tip_controller.addObject(
    'MechanicalObject', 
    template='Rigid3d', 
    name="controlEndEffector",
    showObjectScale=0.3, 
    position=[beam_length, 0, 0, 0, 0, 0, 1],
    showObject=True
)
```

This creates a rigid body at the tip position that can be used to manipulate the beam.

---

## Interacting with Forces

Our `ForceController` also includes keyboard interaction:

```python
def onKeypressedEvent(self, event):
    key = event["key"]
    if key == "+":
        self.applyForce = True
    elif key == "-":
        self.applyForce = False
```

This allows us to increase or decrease the applied force during simulation.

---

## Running the Example

When you run `tutorial_03_basic_beam.py`, you'll see:

1. A beam under gravity
2. A force applied to the tip
3. The ability to modify the force using keyboard input

Try pressing '+' to increase the force and '-' to decrease it!

---

## Key Takeaways

1. **Multiple Force Types**: SOFA offers several ways to apply forces
2. **Dynamic Forces**: Controllers allow for time-varying forces
3. **Interactive Simulation**: Keyboard and mouse events can modify forces in real-time
4. **Visualization**: Forces can be visualized in the scene for better understanding

In the next tutorial, we'll explore more advanced interactions and constraints.
