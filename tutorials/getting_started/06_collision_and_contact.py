# -*- coding: utf-8 -*-
"""
Tutorial 06: Collision and Contact
===================================

This tutorial demonstrates how to enable collision detection and contact
response for a Cosserat beam.

We will create a beam that falls under gravity and collides with a static plane.

Key concepts:
- Setting up a collision pipeline in SOFA.
- Adding collision models to objects.
- Using `PenaltyContactForceField` for contact response.
"""

import os
import sys

# Add the python package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "python"))

from python.cosserat import BeamGeometryParameters, CosseratGeometry

from introduction_and_setup import (_add_cosserat_frame, _add_cosserat_state,
                                    _add_rigid_base, add_mini_header)

v_damping_param: float = 0.4

def createScene(root_node):
    """Create a Cosserat beam scene with collision."""
    add_mini_header(root_node)
    root_node.gravity = [0, -9.81, 0]

    # --- Collision Pipeline ---
    root_node.addObject("DefaultPipeline", name="CollisionPipeline")
    root_node.addObject("BruteForceBroadPhase")
    root_node.addObject("BVHNarrowPhase")
    root_node.addObject("DefaultContactManager", name="ContactManager", response="FrictionContact", responseParams="mu=0.6")
    root_node.addObject("PenaltyContactForceField", name="PenaltyContactForceField", response="default", stiffness=10000, friction=0.6)


    # --- Solver ---
    solver_node = root_node.addChild("solver")
    solver_node.addObject(
        "EulerImplicitSolver",
        rayleighStiffness="0.0",
        rayleighMass="0.0",
        vdamping=v_damping_param,
    )
    solver_node.addObject("SparseLDLSolver", name="solver")

    # --- Beam ---
    beam_geometry_params = BeamGeometryParameters(
        beam_length=20.0,
        nb_section=20,
        nb_frames=20,
    )
    beam_geometry = CosseratGeometry(beam_geometry_params)

    base_node = _add_rigid_base(solver_node, positions=[0, 20, 0, 0, 0, 0, 1]) # Start the beam higher
    bending_node = _add_cosserat_state(solver_node, beam_geometry)
    frame_node = _add_cosserat_frame(
        base_node, bending_node, beam_geometry, beam_mass=5.0
    )
    # Add a collision model to the beam
    beam_collision_node = frame_node.addChild("collision")
    beam_collision_node.addObject("MechanicalObject", name="collisionMO", template="Vec3d", position=beam_geometry.frames_as_vec3())
    beam_collision_node.addObject("LineCollisionModel", name="beamCollisionModel", proximity=0.5)
    beam_collision_node.addObject("BarycentricMapping")


    # --- Floor ---
    floor_node = root_node.addChild("floor")
    floor_node.addObject("MechanicalObject", name="floorMO", template="Rigid3d", position=[0, 0, 0, 0, 0, 0, 1], velocity=[0,0,0,0,0,0])
    floor_collision_node = floor_node.addChild("collision")
    floor_collision_node.addObject("MeshObjLoader", name="loader", filename="mesh/floor.obj", scale=20)
    floor_collision_node.addObject("MeshTopology", src="@loader")
    floor_collision_node.addObject("MechanicalObject", name="floorCollisionMO", src="@loader")
    floor_collision_node.addObject("TriangleCollisionModel", name="floorCollisionModel")
    floor_collision_node.addObject("LineCollisionModel", name="floorLineCollisionModel")
    floor_collision_node.addObject("PointCollisionModel", name="floorPointCollisionModel")
    floor_collision_node.addObject("RigidMapping")

    print("✨ Created a beam that will fall and collide with the floor.")

    return root_node

