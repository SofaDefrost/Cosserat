# -*- coding: utf-8 -*-
"""
Tutorial 08: Hybrid Modeling - Combining Cosserat and FEM
=========================================================

This tutorial demonstrates a powerful feature of SOFA: combining different
physical models in a single simulation. We will attach a 3D deformable
(FEM) object to the tip of a 1D Cosserat beam.

Key concepts:
- Creating a simple FEM model.
- Using `RigidMapping` to connect two different models.
"""

import os
import sys

# Add the python package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "python"))

from cosserat import BeamGeometryParameters, CosseratGeometry

from introduction_and_setup import (_add_cosserat_frame_v2, _add_cosserat_state_v2,
                                    _add_rigid_base, add_mini_header)

def createScene(root_node):
    """Create a scene with a Cosserat beam and an attached FEM object."""
    add_mini_header(root_node)
    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.IO.Mesh') # Needed to use components [MeshVTKLoader]  
    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.LinearSolver.Direct') # Needed to use components [SparseLDLSolver]  
    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.ODESolver.Backward') # Needed to use components [EulerImplicitSolver]  
    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.SolidMechanics.FEM.Elastic') # Needed to use components [TetrahedronFEMForceField]  
    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.Topology.Container.Dynamic') # Needed to use components [TetrahedronSetTopologyContainer]     
    root_node.addObject('RequiredPlugin', pluginName='Sofa.GL.Component.Rendering3D') # Needed to use components [OglModel]  
    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.Mapping.Linear') # Needed to use components [IdentityMapping]      
    root_node.addObject('RequiredPlugin', name='Sofa.Component.Mapping.NonLinear') # Needed to use components [RigidMapping]      
    
    root_node.gravity = [0, -9.81, 0]

    # --- Solver ---
    solver_node = root_node.addChild("solver")
    solver_node.addObject("EulerImplicitSolver", rayleighStiffness="0.0", rayleighMass="0.0", vdamping=0.5)
    solver_node.addObject("SparseLDLSolver", name="solver", template="CompressedRowSparseMatrixd")

    # --- Beam ---
    beam_geometry_params = BeamGeometryParameters(beam_length=20.0, nb_section=10, nb_frames=10)
    beam_geometry = CosseratGeometry(beam_geometry_params)

    base_node = _add_rigid_base(solver_node)
    bending_node = _add_cosserat_state_v2(solver_node, beam_geometry)
    frame_node = _add_cosserat_frame_v2(base_node, bending_node, beam_geometry, beam_mass=2.0)

    # --- FEM Gripper ---
    fem_node = frame_node.addChild("fem_gripper") # fem_node child node of frame_node (Essential for the futur Mapping)
    fem_node.addObject("EulerImplicitSolver", name="fem_solver")
    fem_node.addObject("SparseLDLSolver", name="fem_ldl_solver", template="CompressedRowSparseMatrixd")
    fem_node.addObject("MeshGmshLoader", name="loader", filename="mesh/liver.msh", scale=1.5) # from Gmsh in order to use TetrahedronFEMForceField
    fem_node.addObject("TetrahedronSetTopologyContainer", name="topology", src="@loader")
    fem_node.addObject("MechanicalObject", name="femMO", template="Vec3d", dx=20, dy=0, dz=0)
    fem_node.addObject("TetrahedronFEMForceField", name="femForceField", youngModulus=1000, poissonRatio=0.4)
    fem_node.addObject("UniformMass", totalMass=1.0)


    # # --- Connection ---
    # map the last frame of the beam to the liver (fem) points
    fem_node.addObject("RigidMapping", name="mapping",
                              input=frame_node.FramesMO.getLinkPath(),
                              output=fem_node.femMO.getLinkPath(),
                              index=beam_geometry.get_number_of_frames() - 1)

    fem_visual_node = fem_node.addChild("fem_visual")
    fem_visual_node.addObject("OglModel", name="VisualModel", src="@../loader", color="red", texturename="textures/liver-texture-square.png")
    fem_visual_node.addObject("IdentityMapping")


    print("✨ Created a hybrid model: FEM gripper attached to a Cosserat beam.")

    return root_node

