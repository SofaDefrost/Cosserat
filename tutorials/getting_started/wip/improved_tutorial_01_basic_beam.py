# -*- coding: utf-8 -*-
"""
Tutorial 01: Understanding Beam Geometry and Discretization
=========================================================

This tutorial demonstrates how the number of sections and frames affects 
beam representation. We create two beams with identical physical parameters
but different discretization to show the difference between:

1. Sections: The physical discretization (affects accuracy of physics)
2. Frames: The visual discretization (affects smoothness of visualization)

Key concepts:
- BeamGeometryParameters: Configuring beam discretization
- Creating multiple beams with different parameters
- Visualizing the impact of discretization on beam appearance
"""

import os
import sys

# Add the python package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'python'))

from cosserat import (BeamGeometryParameters, BeamPhysicsBaseParameters,
                      CosseratGeometry)

# Import helper functions from the first tutorial
from tutorials.getting_started.wip.improved_tutorial_00_basic_beam import (add_required_plugins,
                                                                           create_rigid_base,
                                                                           create_cosserat_state,
                                                                           create_cosserat_frame)

def createScene(root_node):
    """
    Create a scene with two beams that have different discretization parameters.
    
    Both beams have the same physical parameters and number of sections,
    but different numbers of frames for visualization.
    """
    # Load required plugins
    add_required_plugins(root_node)
    
    # No gravity in this simple example
    root_node.gravity = [0, 0.0, 0]

    # =========================================================================
    # First beam: 3 sections, 3 frames (same number for both)
    # =========================================================================
    beam_geometry_params = BeamGeometryParameters(
        beam_length=30.0,  # Total beam length
        nb_section=3,      # Number of sections for physics
        nb_frames=3        # Number of frames for visualization (matches sections)
    )

    # Create geometry object - this automatically calculates all needed values
    beam_geometry1 = CosseratGeometry(beam_geometry_params)

    print("✨ Created first beam with:")
    print(f"   - Length: {beam_geometry1.get_beam_length()}")
    print(f"   - Sections: {beam_geometry1.get_number_of_sections()}")
    print(f"   - Frames: {beam_geometry1.get_number_of_frames()}")
    print(f"   - Section lengths: {beam_geometry1.section_lengths}")

    # Create rigid base for first beam (offset to the left)
    base_node1 = create_rigid_base(
        root_node, 
        node_name="rigid_base1", 
        positions=[-15, 0, 0, 0, 0, 0, 1]  # Offset to the left
    )

    # Define custom bending states for this tutorial (slight bend in z-direction)
    custom_bending_states = [
        [0.0, 0.0, 0.1],  # Section 1: bend around z-axis
        [0.0, 0.0, 0.1],  # Section 2: bend around z-axis
        [0.0, 0.0, 0.1]   # Section 3: bend around z-axis
    ]

    # Create cosserat state using the geometry object
    bending_node1 = create_cosserat_state(
        root_node, 
        beam_geometry1, 
        node_name="cos_coord1",
        custom_bending_states=custom_bending_states
    )

    # Create cosserat frame using the geometry object
    create_cosserat_frame(
        base_node1, 
        bending_node1, 
        beam_geometry1, 
        node_name="cosserat_frame_node1",
        beam_mass=0.0
    )

    # =========================================================================
    # Second beam: 3 sections, 12 frames (more frames for smoother visualization)
    # =========================================================================
    beam_geometry_params2 = BeamGeometryParameters(
        beam_length=30.0,  # Same beam length
        nb_section=3,      # Same number of sections (physics is identical)
        nb_frames=12       # More frames for smoother visualization
    )

    # Create second geometry object
    beam_geometry2 = CosseratGeometry(beam_geometry_params2)
    
    print("✨ Created second beam with:")
    print(f"   - Length: {beam_geometry2.get_beam_length()}")
    print(f"   - Sections: {beam_geometry2.get_number_of_sections()}")
    print(f"   - Frames: {beam_geometry2.get_number_of_frames()}")
    print(f"   - Section lengths: {beam_geometry2.section_lengths}")

    # Create rigid base for second beam (offset to the right)
    base_node2 = create_rigid_base(
        root_node, 
        node_name="rigid_base2", 
        positions=[15, 0, 0, 0, 0, 0, 1]  # Offset to the right
    )
    
    # Create cosserat state for the second beam (same bending states)
    bending_node2 = create_cosserat_state(
        root_node, 
        beam_geometry2, 
        node_name="cos_coord2",
        custom_bending_states=custom_bending_states
    )
    
    # Create cosserat frame for the second beam
    create_cosserat_frame(
        base_node2, 
        bending_node2, 
        beam_geometry2, 
        node_name="cosserat_frame_node2",
        beam_mass=0.0
    )

    # Note: Both beams have identical physical behavior (3 sections with the same
    # bending states), but the second beam looks smoother because it has more frames
    # for visualization.

    return root_node
