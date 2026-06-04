# -*- coding: utf-8 -*-
"""
Test de la fonction applyJT des contraintes de Frames2StrainCosseratMapping
"""

import os
import sys

# Add the python package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "python"))

from cosserat import BeamGeometryParameters, CosseratGeometry

from basic_functions import (_add_cosserat_frame, _add_cosserat_state,
                                    _add_rigid_base, add_mini_header)

v_damping_param: float =  8e-1  # Damping parameter for dynamics
beam_length :float = 10
nb_section : int = 15
beam_mass: float = 1.0

def createScene(root_node):
    """Create a Cosserat beam scene with forces and dynamics."""
    # Configure scene with time integration
    add_mini_header(root_node)
    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.Collision.Geometry') # Needed to use components [LineCollisionModel]  
    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.Mapping.Linear') # Needed to use components [IdentityMapping]  
    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.Topology.Container.Dynamic') # Needed to use components [EdgeSetTopologyContainer]    
    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.IO.Mesh') # Needed to use components [MeshOBJLoader]  
    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.Topology.Container.Constant') # Needed to use components [MeshTopology]  
    root_node.addObject('RequiredPlugin', pluginName='Sofa.GL.Component.Rendering3D') # Needed to use components [OglModel]
    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.Mapping.NonLinear') # Needed to use components [RigidMapping]  
    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.Collision.Detection.Algorithm') # Needed to use components [BVHNarrowPhase,BruteForceBroadPhase,CollisionPipeline]  
    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.Collision.Detection.Intersection') # Needed to use components [MinProximityIntersection]  
    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.Collision.Response.Contact') # Needed to use components [RuleBasedContactManager]  
    
    # Add gravity
    root_node.gravity = [0, -9.81, 0]  # Add gravity!

    # --- Collision Pipeline ---
    root_node.addObject("CollisionPipeline")
    root_node.addObject("BruteForceBroadPhase")
    root_node.addObject("BVHNarrowPhase")
    root_node.addObject("MinProximityIntersection", name="Proximity", alarmDistance=2.0, contactDistance=1.0)
    root_node.addObject("RuleBasedContactManager", name="Response", response="PenalityContactForceField")


    # root_node.dt = 0.01
    # Configure time integration and solver
    solver = root_node.addChild("solver")

    solver.addObject(
        "EulerImplicitSolver",
        firstOrder="0",
        rayleighStiffness="0.1",
        rayleighMass="0.1",
        vdamping=v_damping_param,  # Damping parameter for dynamics
    )
    solver.addObject("CGLinearSolver", iterations=1000, tolerance=1e-10, threshold=1e-10)

    # --- Beam ---
    beam_geometry_params = BeamGeometryParameters(
        beam_length=beam_length,
        nb_section=nb_section,
        nb_frames=nb_section,
    )
    beam_geometry = CosseratGeometry(beam_geometry_params)

        ## base node
    rigid_base = _add_rigid_base(solver, node_name="rigid_base")
    
    ## frame node
    frame_node = _add_cosserat_frame(solver, beam_geometry, beam_mass=beam_mass, constraint=False)

    ## bending node
    strain_node = _add_cosserat_state(rigid_base, frame_node, beam_geometry)
 
 
    # Add a collision model to the beam
    beam_collision_node = frame_node.addChild("collision")
    beam_collision_node.addObject("EdgeSetTopologyContainer", position=beam_geometry.frames, edges=beam_geometry.get_edge_list())
    beam_collision_node.addObject("MechanicalObject", name="collisionMO", template="Vec3d", position=beam_geometry.frames)
    beam_collision_node.addObject("LineCollisionModel", name="beamCollisionModel", contactStiffness=10)
    beam_collision_node.addObject("PointCollisionModel", contactStiffness=10)

    beam_collision_node.addObject("IdentityMapping")

    # --- Floor ---
    floor_node = root_node.addChild("floor")
    floor_node.addObject("MeshOBJLoader", name="loader", filename="mesh/floor.obj", scale=0.4, translation=[0, -5, 0.], rotation=[0, 0., 0.])
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


    return root_node