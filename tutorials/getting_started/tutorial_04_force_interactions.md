---
title: Basic Slide 04
date: 2025-06-17
---

# Add force Interactions

In this tutorial, we will explore three different ways to apply forces to the beam.

[[./04_force_interactions.py]]

## First kind of Force: Constant Force

This is the simplest way to apply a force. We use a `ConstantForceField` to apply a constant force to a specific frame of the beam.

```python
def compute_force(self):
        with self.forceNode.forces.writeable() as force:
            vec = [0., 0., 0., 0., self.forceCoeff / sqrt(2), self.forceCoeff / sqrt(2)]
            for i, v in enumerate(vec):
                force[0][i] = v
```

### Code Explanation

- `forceNode.forces.writeable()`: This gives us access to the forces that will be applied to the beam.
- `vec`: This is a 6D vector representing the force and torque to be applied. The first 3 components are the force in x, y, and z, and the last 3 are the torque around x, y, and z.
- `force[0][i] = v`: This sets the force for the first frame in the `forceNode`.

## Second kind of Force: Proportional Force

In this case, the applied force is proportional to the displacement of the end of the beam. This is a simple form of feedback control.

```python
def compute_force(self):
        with self.forceNode.forces.writeable() as force:
            # get the displacement of the end of the beam
            displacement = self.frame_mo.position.value[self.beam_geometry.nb_frames][0]
            
            vec = [0., 0., 0., 0., -self.forceCoeff * displacement, 0.]
            for i, v in enumerate(vec):
                force[0][i] = v
```

### Code Explanation

- `self.frame_mo.position.value[self.beam_geometry.nb_frames][0]`: This gets the x-displacement of the last frame of the beam.
- `vec = [0., 0., 0., 0., -self.forceCoeff * displacement, 0.]`: The applied torque around the y-axis is proportional to the displacement in the x-direction.

## Third kind of Force: Force from a Controller

In this case, we use a separate rigid body to control the end of the beam. The force is applied to the end of the beam to make it follow the controller.

```python
def compute_force(self):
        with self.forceNode.forces.writeable() as force:
            # get the position of the controller
            controller_pos = self.controller_state.position.value[0]
            
            # get the position of the end of the beam
            end_effector_pos = self.frame_mo.position.value[self.beam_geometry.nb_frames]
            
            # compute the error
            error = controller_pos - end_effector_pos
            
            vec = [self.forceCoeff * error[0], self.forceCoeff * error[1], self.forceCoeff * error[2], 0., 0., 0.]
            for i, v in enumerate(vec):
                force[0][i] = v
```

### Code Explanation

- `self.controller_state.position.value[0]`: This gets the position of the controller.
- `self.frame_mo.position.value[self.beam_geometry.nb_frames]`: This gets the position of the end of the beam.
- `error = controller_pos - end_effector_pos`: This computes the difference between the desired position (the controller) and the actual position (the end of the beam).
- `vec = [self.forceCoeff * error[0], self.forceCoeff * error[1], self.forceCoeff * error[2], 0., 0., 0.]`: The applied force is proportional to the error, which will drive the end of the beam towards the controller.

