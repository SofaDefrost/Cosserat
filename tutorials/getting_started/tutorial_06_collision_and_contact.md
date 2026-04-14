---
title: Tutorial 06 - Collision and Contact
date: 2025-06-28
---

# Tutorial 06: Collision and Contact

This tutorial covers how to make your Cosserat beam interact with other objects in the scene through collision and contact.

[[./06_collision_and_contact.py]]

## The Goal: A Falling Beam

We will simulate a beam falling under gravity and colliding with a solid floor. This will introduce the core components required for any simulation involving contact.

## The Collision Pipeline

To enable collision in SOFA, you must first set up a collision pipeline at the root of your scene. This pipeline defines the different phases of collision detection.

```python
# --- Collision Pipeline ---
root_node.addObject("DefaultPipeline", name="CollisionPipeline")
root_node.addObject("BruteForceBroadPhase")
root_node.addObject("BVHNarrowPhase")
root_node.addObject("DefaultContactManager", name="ContactManager", response="FrictionContact", responseParams="mu=0.6")
root_node.addObject("PenaltyContactForceField", name="PenaltyContactForceField", response="default", stiffness=10000, friction=0.6)
```

- **DefaultPipeline**: The standard pipeline for collision.
- **BruteForceBroadPhase** / **BVHNarrowPhase**: These components handle the detection of potential collisions (broad phase) and the precise calculation of intersections (narrow phase).
- **DefaultContactManager**: Manages the contacts that are detected.
- **PenaltyContactForceField**: This is the component that generates the forces to prevent objects from penetrating each other. `stiffness` controls how hard the contact is, and `friction` defines the friction coefficient.

## Collision Models

Next, you need to give your objects a physical shape for the collision detection system to see. This is done by adding **Collision Models**.

### For the Beam:

```python
# Add a collision model to the beam
beam_collision_node = frame_node.addChild("collision")
beam_collision_node.addObject("MechanicalObject", name="collisionMO", template="Vec3d", position=beam_geometry.frames_as_vec3())
beam_collision_node.addObject("LineCollisionModel", name="beamCollisionModel", proximity=0.5)
beam_collision_node.addObject("BarycentricMapping")
```

- We create a sub-node to hold the collision components.
- We use a `LineCollisionModel`, which is suitable for representing a slender beam.
- A `BarycentricMapping` is used to link the collision model back to the beam's main mechanical object.

### For the Floor:

```python
# --- Floor ---
floor_node = root_node.addChild("floor")
# ...
floor_collision_node.addObject("MeshObjLoader", name="loader", filename="mesh/floor.obj", scale=20)
# ...
floor_collision_node.addObject("TriangleCollisionModel", name="floorCollisionModel")
# ...
floor_collision_node.addObject("RigidMapping")
```

- For the floor, we load a 3D model from an `.obj` file.
- We use a `TriangleCollisionModel` because the floor is represented by a mesh of triangles.
- A `RigidMapping` is used to link the collision model to the floor's rigid body mechanical object.

Now, when you run the simulation, the beam will fall and realistically bounce and rest on the floor.
