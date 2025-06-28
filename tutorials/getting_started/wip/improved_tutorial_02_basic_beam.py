# -*- coding: utf-8 -*-
"""
Tutorial 02: Dynamic Simulation with Gravity
===========================================

This tutorial builds on the previous tutorials by adding:
- Gravity forces
- Time integration for dynamics
- Mass distribution on the beam
- Comparison of different beam discretizations under gravity

Key concepts:
- Setting up solvers for dynamic simulation
- Adding mass to the beam for proper dynamics
- Understanding how discretization affects physical behavior
"""

import os
import sys

# Add the python package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "python"))

from cosserat import BeamGeometryParameters, CosseratGeometry

# Import helper functions from the first tutorial
from tutorials.getting_started.wip.improved_tutorial_00_basic_beam import (add_required_plugins,
                                                                           create_rigid_base,
                                                                           create_cosserat_state,
                                                                           create_cosserat_frame)

# Damping parameter for dynamics - helps stabilize the simulation
v_damping_param: float = 8.e-1

def create_solver_node(parent_node, name="solver"):
    """
    Create a solver node with time integration for dynamics.
    
    Parameters:
        parent_node: The SOFA node to attach this to
        name: Name for the node
        
    Returns:
        The created solver node
    """
    solver_node = parent_node.addChild(name)
    
    # Add an implicit Euler solver for time integration
    solver_node.addObject(
        "EulerImplicitSolver",
        firstOrder="0",              # Second-order integration
        rayleighStiffness="0.0",     # No stiffness-based damping
        rayleighMass="0.0",          # No mass-based damping
        vdamping=v_damping_param,    # Velocity damping for stability
    )
    
    # Add a sparse LDL solver for efficient linear system solving
    solver_node.addObject("SparseLDLSolver", name="solver")
    
    return solver_node

def createScene(root_node):
    """
    Create a scene with two beams under gravity, using different discretizations.
    """
    # Configure scene with all required plugins
    add_required_plugins(root_node)

    # Add gravity! This will cause the beams to fall and bend
    root_node.gravity = [0, -9.81, 0]

    # =========================================================================
    # First beam: 30 sections, 30 frames (high physical accuracy)
    # =========================================================================
    
    # Create solver node for the first beam
    solver_node1 = create_solver_node(root_node, name="solver_beam1")
    
    # Define beam geometry parameters with high discretization
    beam_geometry_params1 = BeamGeometryParameters(
        beam_length=30.0,  # Total beam length
        nb_section=30,     # High number of sections for accurate physics
        nb_frames=30,      # Same number of frames for visualization
    )

    # Create geometry object
    beam_geometry1 = CosseratGeometry(beam_geometry_params1)

    print(f"🚀 Created high-resolution beam with:")
    print(f"   - Length: {beam_geometry1.get_beam_length()}")
    print(f"   - Sections: {beam_geometry1.get_number_of_sections()}")
    print(f"   - Frames: {beam_geometry1.get_number_of_frames()}")
    print(f"   - Mass will be distributed across frames")

    # Create rigid base at the left side
    base_node1 = create_rigid_base(
        solver_node1, 
        positions=[-15, 0, 0, 0, 0, 0, 1]  # Position on the left
    )

    # Initialize with straight beam (no initial bending)
    custom_bending_states1 = []
    for i in range(beam_geometry1.get_number_of_sections()):
        custom_bending_states1.append([0, 0.0, 0.0])

    # Create cosserat state using geometry
    bending_node1 = create_cosserat_state(
        solver_node1, 
        beam_geometry1, 
        node_name="cosserat_states1",
        custom_bending_states=custom_bending_states1
    )

    # Create cosserat frame with mass (important for dynamics!)
    frame_node1 = create_cosserat_frame(
        base_node1, 
        bending_node1, 
        beam_geometry1, 
        node_name="cosserat_frames1",
        beam_mass=5.0  # Add mass to the beam
    )

    # =========================================================================
    # Second beam: 3 sections, 30 frames (lower physical accuracy, same visual quality)
    # =========================================================================
    
    # Create solver node for the second beam
    solver_node2 = create_solver_node(root_node, name="solver_beam2")
    
    # Define beam geometry parameters with lower physical discretization
    beam_geometry_params2 = BeamGeometryParameters(
        beam_length=30.0,  # Same beam length
        nb_section=3,      # Only 3 sections - less physical accuracy
        nb_frames=30,      # Same number of frames for visualization
    )

    # Create geometry object
    beam_geometry2 = CosseratGeometry(beam_geometry_params2)

    print(f"🚀 Created low-resolution beam with:")
    print(f"   - Length: {beam_geometry2.get_beam_length()}")
    print(f"   - Sections: {beam_geometry2.get_number_of_sections()}")
    print(f"   - Frames: {beam_geometry2.get_number_of_frames()}")
    print(f"   - Mass will be distributed across frames")

    # Create rigid base at the right side
    base_node2 = create_rigid_base(
        solver_node2, 
        node_name="rigid_base2",
        positions=[15, 0, 0, 0, 0, 0, 1]  # Position on the right
    )

    # Initialize with straight beam (no initial bending)
    custom_bending_states2 = []
    for i in range(beam_geometry2.get_number_of_sections()):
        custom_bending_states2.append([0, 0.0, 0.0])

    # Create cosserat state using geometry
    bending_node2 = create_cosserat_state(
        solver_node2, 
        beam_geometry2, 
        node_name="cosserat_states2",
        custom_bending_states=custom_bending_states2
    )

    # Create cosserat frame with same mass as the first beam
    frame_node2 = create_cosserat_frame(
        base_node2, 
        bending_node2, 
        beam_geometry2, 
        node_name="cosserat_frames2",
        beam_mass=5.0  # Same mass as first beam
    )

    # When you run this scene and press 'Animate', you'll see:
    # 1. Both beams fall under gravity
    # 2. The first beam (30 sections) shows smooth, realistic bending
    # 3. The second beam (3 sections) shows more segmented bending
    # This demonstrates that physical accuracy depends on the number of sections!

    return root_node
