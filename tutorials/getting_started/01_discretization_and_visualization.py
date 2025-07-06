# -*- coding: utf-8 -*-
"""
Tutorial 01: Basic Cosserat Beam
===============================

This tutorial demonstrates how to create a basic Cosserat beam using the new
CosseratGeometry class. This approach is much cleaner than manually calculating
geometry parameters.

Key concepts:
- BeamGeometryParameters: Defines beam dimensions and discretization
- CosseratGeometry: Automatically calculates beam geometry
- Clean, reusable beam creation functions
"""

import importlib.util
import os
import sys

# Add the python package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'python'))

from introduction_and_setup import _add_cosserat_frame, _add_cosserat_state, _add_rigid_base, add_mini_header

from cosserat import (BeamGeometryParameters, BeamPhysicsBaseParameters,
                      CosseratGeometry)

def createScene(root_node):
    """Create a basic Cosserat beam scene using the new CosseratGeometry class."""
    # Load required plugins
    add_mini_header(
        root_node
    )
    root_node.gravity = [0, 0.0, 0]

    # === NEW APPROACH: Use CosseratGeometry ===
    # Define beam geometry parameters
    beam_geometry_params = BeamGeometryParameters(
        beam_length=30.0,  # Total beam length
        nb_section=3,      # Number of sections for physics
        nb_frames=3        # Number of frames for visualization
    )

    # Create geometry object - this automatically calculates all the geometry!
    beam_geometry1 = CosseratGeometry(beam_geometry_params)

    print("✨ Created beam with:")
    print(f"   - Length: {beam_geometry1.get_beam_length()}")
    print(f"   - Sections: {beam_geometry1.get_number_of_sections()}")
    print(f"   - Frames: {beam_geometry1.get_number_of_frames()}")
    print(f"   - Section lengths: {beam_geometry1.section_lengths}")

    # Create rigid base
    base_node1 = _add_rigid_base(root_node, node_name="rigid_base1")

    # Custom bending states for this tutorial (slight bend)
    custom_bending_states = [
        [0.0, 0.0, 0.1],  # Section 1: slight bend in y and z
        [0.0, 0.0, 0.1],  # Section 2: slight bend in y and z
        [0.0, 0.0, 0.1]   # Section 3: slight bend in y and z
    ]

    # Create cosserat state using the geometry object
    bending_node = _add_cosserat_state(root_node, beam_geometry1, node_name="cos_coord1",
                                       custom_bending_states=custom_bending_states)

    # Create cosserat frame using the geometry object
    _add_cosserat_frame(base_node1, bending_node, beam_geometry1, node_name="cosserat_in_Sofa_frame_node2",
                        beam_mass=0.0)

    #---------
    # ----------- Create a second beam with different parameters -----------
    # Define second beam geometry parameters
    ##############
    beam_geometry_params2 = BeamGeometryParameters(
        beam_length=30.0,  # Total beam length
        nb_section=3,  # Number of sections for physics
        nb_frames=12  # Number of frames for visualization
    )

    # Create second geometry object
    # Note: We use more frames (12) than sections (3). This is a common
    # practice to get a smooth visual representation of the beam while
    # keeping the physics simulation efficient with fewer sections.
    beam_geometry2 = CosseratGeometry(beam_geometry_params2)
    print("✨ Created second beam with:")
    print(f"   - Length: {beam_geometry2.get_beam_length()}")
    print(f"   - Sections: {beam_geometry2.get_number_of_sections()}")
    print(f"   - Frames: {beam_geometry2.get_number_of_frames()}")
    print(f"   - Section lengths: {beam_geometry2.section_lengths}")

    # Create rigid base
    base_node2 = _add_rigid_base(root_node, node_name="rigid_base2")
    # Create cosserat state for the second beam
    bending_node2 = _add_cosserat_state(root_node, beam_geometry2, node_name="cos_coord2",
                                        custom_bending_states=custom_bending_states)
    # Create cosserat frame for the second beam
    _add_cosserat_frame(base_node2, bending_node2, beam_geometry2, node_name="cosserat_in_Sofa_frame_node2" , beam_mass=0.0)

    return root_node
