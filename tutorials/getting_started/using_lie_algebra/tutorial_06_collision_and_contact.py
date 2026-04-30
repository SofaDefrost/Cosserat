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
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "python"))

from cosserat import BeamGeometryParameters, CosseratGeometry

from introduction_and_setup import (_add_cosserat_frame_v2, _add_cosserat_state_v2,
                                    _add_rigid_base, add_mini_header)

v_damping_param: float = 0.4

def createScene(root_node):
    """Create a Cosserat beam scene with collision."""
    add_mini_header(root_node)

    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.Collision.Detection.Algorithm') # Needed to use components [BVHNarrowPhase,BruteForceBroadPhase,CollisionPipeline]  
    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.Collision.Detection.Intersection') # Needed to use components [MinProximityIntersection]  
    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.Collision.Response.Contact') # Needed to use components [RuleBasedContactManager]  
    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.LinearSolver.Direct') # Needed to use components [SparseLDLSolver]  
    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.ODESolver.Backward') # Needed to use components [EulerImplicitSolver] 
    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.Collision.Geometry') # Needed to use components [LineCollisionModel]  
    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.Mapping.Linear') # Needed to use components [IdentityMapping]  
    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.Topology.Container.Dynamic') # Needed to use components [EdgeSetTopologyContainer]    
    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.IO.Mesh') # Needed to use components [MeshOBJLoader]  
    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.Topology.Container.Constant') # Needed to use components [MeshTopology]  
    root_node.addObject('RequiredPlugin', pluginName='Sofa.GL.Component.Rendering3D') # Needed to use components [OglModel]
    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.Mapping.NonLinear') # Needed to use components [RigidMapping]  
    
    root_node.gravity = [0, -9.81, 0]

    # --- Collision Pipeline ---
    root_node.addObject("CollisionPipeline")
    root_node.addObject("BruteForceBroadPhase")
    root_node.addObject("BVHNarrowPhase")
    root_node.addObject("MinProximityIntersection", name="Proximity", alarmDistance=2.0, contactDistance=1.0)
    root_node.addObject("RuleBasedContactManager", name="Response", response="PenalityContactForceField")


    # --- Solver ---
    solver_node = root_node.addChild("solver")
    solver_node.addObject(
        "EulerImplicitSolver",
        rayleighStiffness="0.0",
        rayleighMass="0.0",
        vdamping=v_damping_param,
    )
    solver_node.addObject("SparseLDLSolver", name="solver", template="CompressedRowSparseMatrixd")

    # --- Beam ---
    beam_geometry_params = BeamGeometryParameters(
        beam_length=20.0,
        nb_section=20,
        nb_frames=20,
    )
    beam_geometry = CosseratGeometry(beam_geometry_params)

    # base_node = _add_rigid_base(solver_node, positions=[0, 20, 0, 0, 0, 0, 1]) # Start the beam higher
    # on ne peut pas utiliser la fct _add_rigid_base dans laquelle on fixe la base de la tige
    base_node = solver_node.addChild("rigid_base")
    base_node.addObject("MechanicalObject", template="Rigid3d", name="cosserat_base_mo", position=[0, 0, 0, 0, 0, 0, 1.], showObject=True, showObjectScale=0.5)


    bending_node = _add_cosserat_state_v2(solver_node, beam_geometry)
    frame_node = _add_cosserat_frame_v2(
        base_node, bending_node, beam_geometry, beam_mass=10.0
    )
    # Add a collision model to the beam
    beam_collision_node = frame_node.addChild("collision")
    beam_collision_node.addObject("EdgeSetTopologyContainer", position=beam_geometry.frames, edges=beam_geometry.get_edge_list())
    beam_collision_node.addObject("MechanicalObject", name="collisionMO", template="Vec3d", position=beam_geometry.frames)
    beam_collision_node.addObject("LineCollisionModel", name="beamCollisionModel", contactStiffness=10)
    beam_collision_node.addObject("PointCollisionModel", contactStiffness=10)

    beam_collision_node.addObject("IdentityMapping")


    # --- Floor ---
    floor_node = root_node.addChild("floor")
    floor_node.addObject("MeshOBJLoader", name="loader", filename="mesh/floor.obj", scale=0.55, translation=[0, -10, 0.], rotation=[0, 0., 0.])
    floor_node.addObject("MechanicalObject", name="floorMO", template="Rigid3d", position=[0, 0, 0, 0, 0, 0, 1], velocity=[0,0,0,0,0,0])
    floor_visual_node = floor_node.addChild("VisualModel")
    floor_visual_node.addObject("OglModel", name="Visual", src="@../loader", color=[0.0, 1.0, 0.0])
    floor_visual_node.addObject("IdentityMapping", input="@../floorMO", output="@Visual")


    
    floor_collision_node = floor_node.addChild("collision")
    floor_collision_node.addObject("MeshTopology", src="@../loader")
    floor_collision_node.addObject("MechanicalObject", name="floorCollisionMO", src="@../loader")
    floor_collision_node.addObject("TriangleCollisionModel", name="floorCollisionModel", contactStiffness=200)
    floor_collision_node.addObject("LineCollisionModel", name="floorLineCollisionModel", contactStiffness=200)
    floor_collision_node.addObject("PointCollisionModel", name="floorPointCollisionModel", contactStiffness=200)
    floor_collision_node.addObject("RigidMapping")

    print("✨ Created a beam that will fall and collide with the floor.")

    return root_node

