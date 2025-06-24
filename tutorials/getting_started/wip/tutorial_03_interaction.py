# -*- coding: utf-8 -*-
"""
Tutorial 03: Interactive Cosserat Beam
=====================================

This tutorial demonstrates how to use the CosseratBase prefab class for creating
an interactive Cosserat beam with collision detection and springs.

Key concepts:
- CosseratBase: High-level prefab class for complete beam setup
- BeamPhysicsParameters: Defines material properties
- Parameters: Combines geometry and physics parameters
- Interactive simulation with forces

This shows the highest-level API - CosseratBase handles everything automatically!
"""

__authors__ = "adagolodjo"
__contact__ = "adagolodjo@protonmail.com"
__version__ = "1.0.0"
__copyright__ = "(c) 2021,Inria"
__date__ = "October, 26 2021"

import sys
import os
# Add the python package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'python'))

from cosserat import addHeader, addSolverNode, Parameters
from cosserat import BeamPhysicsParameters, BeamGeometryParameters
from cosserat import CosseratBase

# Define beam geometry using the new clean approach
geoParams = BeamGeometryParameters(
    beam_length=30.0,
    nb_section=32,
    nb_frames=32,
    build_collision_model=0
)

# Define physics parameters 
physicsParams = BeamPhysicsParameters(
    beam_mass=0.3,
    young_modulus=1.0e3,
    poisson_ratio=0.38,
    beam_radius=1.0,
    beam_length=30.0
)

# Combine parameters
Params = Parameters(beam_geo_params=geoParams, beam_physics_params=physicsParams)

print(f"🎮 Setting up interactive beam with:")
print(f"   - Length: {geoParams.beam_length}")
print(f"   - Sections: {geoParams.nb_section}")
print(f"   - Young's modulus: {physicsParams.young_modulus}")
print(f"   - Mass: {physicsParams.beam_mass}")

def createScene(root_node):
    """Create an interactive Cosserat beam scene using CosseratBase prefab."""
    # Setup scene with solver
    addHeader(root_node)
    root_node.gravity = [0, -9.81, 0.]
    solver_node = addSolverNode(root_node, name="solver_node")

    print(f"🎯 Creating CosseratBase with {Params.beam_geo_params.nb_section} sections...")
    
    # === HIGHEST LEVEL API: CosseratBase handles everything! ===
    # The CosseratBase prefab automatically:
    # - Creates the rigid base
    # - Sets up cosserat coordinates
    # - Creates frames and mappings
    # - Handles all the geometry calculations
    beam = solver_node.addChild(CosseratBase(parent=solver_node, params=Params))
    
    # Access the rigid base node and add a spring force field for attachment
    # Note: rigid_base_node is the property name in the new structure
    beam.rigid_base_node.addObject(
        "RestShapeSpringsForceField",
        name="spring",
        stiffness=1e8,
        angularStiffness=1.0e8,
        external_points=0,
        points=0,
        template="Rigid3d"
    )
    
    # Optional: Add a force at the tip for interactivity
    tip_force = [-5.0, 0.0, 0.0, 0, 0, 0]  # Gentle force in -X direction
    last_frame = Params.beam_geo_params.nb_frames
    
    beam.cosserat_frame.addObject(
        'ConstantForceField',
        name='interactiveForce', 
        indices=[last_frame],
        forces=[tip_force],
        showArrowSize=0.1
    )
    
    print(f"✨ CosseratBase beam created successfully!")
    print(f"   - Base spring stiffness: 1e8")
    print(f"   - Applied tip force: {tip_force[:3]}")
    print(f"   - Interactive simulation ready")

    return root_node
