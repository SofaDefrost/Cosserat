---
title: Tutorial 08 - Hybrid Modeling
date: 2025-06-28
---

# Tutorial 08: Hybrid Modeling - Combining Cosserat and FEM

This tutorial showcases one of the most powerful features of SOFA: the ability to combine different modeling techniques in a single simulation. We will attach a 3D deformable object, modeled using the Finite Element Method (FEM), to the end of our 1D Cosserat beam.

[[./08_hybrid_modeling.py]]

## The Goal: A Beam with a Soft Gripper

We will model a flexible arm (the Cosserat beam) with a soft, squishy gripper (the FEM object) attached to its tip. This is a common pattern in soft robotics.

## Creating the FEM Object

First, we create a new node in the scene to house our FEM model.

```python
# --- FEM Gripper ---
fem_node = root_node.addChild("fem_gripper")
# ...
fem_node.addObject("MeshVTKLoader", name="loader", filename="mesh/liver.vtk")
fem_node.addObject("TetrahedronSetTopologyContainer", name="topology", src="@loader")
fem_node.addObject("MechanicalObject", name="femMO", template="Vec3d", dx=20, dy=0, dz=0)
fem_node.addObject("TetrahedronFEMForceField", name="femForceField", youngModulus=1000, poissonRatio=0.4)
# ...
```

- We load a volumetric mesh from a `.vtk` file.
- We use a `TetrahedronSetTopologyContainer` because our mesh is composed of tetrahedra.
- A `MechanicalObject` stores the positions of the mesh nodes.
- A `TetrahedronFEMForceField` computes the internal forces that make the object behave like a deformable solid.

## The Magic: `BarycentricMapping`

Now, we need to connect the two models. The `BarycentricMapping` is the key component for this.

```python
# --- Connection ---
connection_node = fem_node.addChild("connection")
connection_node.addObject("BarycentricMapping", name="mapping",
                          map_from=frame_node.FramesMO.getLinkPath(),
                          map_to=fem_node.femMO.getLinkPath(),
                          use_rigid=True,
                          rigid_indices=[beam_geometry.get_number_of_frames() - 1])
```

- **`map_from`**: This points to the "parent" or "master" object, which is our Cosserat beam's `FramesMO`.
- **`map_to`**: This points to the "child" or "slave" object, which is our FEM gripper's `MechanicalObject`.
- **`use_rigid=True`**: This tells the mapping to perform a rigid mapping, which is what we want when attaching an object to a specific frame.
- **`rigid_indices`**: This specifies which frame of the master object the child object should be attached to. We choose the last frame of the beam.

This mapping essentially makes the entire FEM object follow the motion of the beam's tip. Now you have a complex, hybrid robot that combines the efficiency of a 1D model with the detail of a 3D model.
