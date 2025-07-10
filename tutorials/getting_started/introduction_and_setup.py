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

import os
import sys

# Add the python package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "python"))

from cosserat import (BeamGeometryParameters, BeamPhysicsBaseParameters,
                      CosseratGeometry)

# Global parameters
stiffness_param: float = 1.0e10
beam_radius: float = 1.0


def _add_rigid_base(p_node, node_name="rigid_base", positions=None):
    """Create a rigid base node for the beam."""
    if positions is None:
        positions = [0, 0, 0, 0, 0, 0, 1]
    rigid_base_node = p_node.addChild(node_name)
    rigid_base_node.addObject(
        "MechanicalObject",
        template="Rigid3d",
        name="cosserat_base_mo",
        position=positions,
        showObject=True,
        showObjectScale="0.1",
    )
    rigid_base_node.addObject(
        "RestShapeSpringsForceField",
        name="spring",
        stiffness=stiffness_param,
        angularStiffness=stiffness_param,
        external_points="0",
        mstate="@cosserat_base_mo",
        points="0",
        template="Rigid3d",
    )
    return rigid_base_node


def _add_cosserat_state(p_node, geometry: CosseratGeometry, node_name="cosserat_coordinates", custom_bending_states=None):
    """Create the cosserat coordinate node using CosseratGeometry."""
    cosserat_coordinate_node = p_node.addChild(node_name)

    # Use geometry data or custom bending states
    bending_states = (
        custom_bending_states if custom_bending_states else geometry.bendingState
    )

    cosserat_coordinate_node.addObject(
        "MechanicalObject",
        template="Vec3d",
        name="cosserat_state",
        position=bending_states,
    )
    cosserat_coordinate_node.addObject(
        "BeamHookeLawForceField",
        crossSectionShape="circular",
        length=geometry.section_lengths,  # Use geometry data
        radius=2.0,
        youngModulus=1.0e3,
        poissonRatio=0.4,
    )
    return cosserat_coordinate_node


def _add_cosserat_frame(
    p_node, bending_node, geometry: CosseratGeometry, node_name="cosserat_in_Sofa_frame_node", beam_mass=0.0
):
    """Create the cosserat frame node using CosseratGeometry."""
    cosserat_in_sofa_frame_node = p_node.addChild(node_name)
    bending_node.addChild(cosserat_in_sofa_frame_node)

    frames_mo = cosserat_in_sofa_frame_node.addObject(
        "MechanicalObject",
        template="Rigid3d",
        name="FramesMO",
        position=geometry.frames,  # Use geometry data
        showIndices=1,
        showObject=1,
        showObjectScale=0.8,
    )

    cosserat_in_sofa_frame_node.addObject("UniformMass", totalMass=beam_mass)

    cosserat_in_sofa_frame_node.addObject(
        "DiscreteCosseratMapping",
        curv_abs_input=geometry.curv_abs_sections,  # Use geometry data
        curv_abs_output=geometry.curv_abs_frames,  # Use geometry data
        name="cosseratMapping",
        input1=bending_node.cosserat_state.getLinkPath(),
        input2=p_node.cosserat_base_mo.getLinkPath(),
        output=frames_mo.getLinkPath(),
        debug=0,
        radius=beam_radius,
    )
    return cosserat_in_sofa_frame_node

def add_mini_header(root_node):
    root_node.addObject("RequiredPlugin", name="Sofa.Component.Mass")
    root_node.addObject("RequiredPlugin", name="Sofa.Component.SolidMechanics.Spring")
    root_node.addObject("RequiredPlugin", name="Sofa.Component.StateContainer")
    root_node.addObject("RequiredPlugin", name="Sofa.Component.Visual")
    root_node.addObject("RequiredPlugin", name="Cosserat")

    # Configure scene
    root_node.addObject(
        "VisualStyle",
        displayFlags="showBehaviorModels showCollisionModels showMechanicalMappings",
    )


def createScene(root_node):
    """Create a basic Cosserat beam scene using the new CosseratGeometry class."""
    # Load required plugins
    add_mini_header(root_node)

    # add gravity
    root_node.gravity = [0, 0.0, 0]

    # === NEW APPROACH: Use CosseratGeometry ===
    # Define beam geometry parameters
    beam_geometry_params = BeamGeometryParameters(
        beam_length=30.0,  # Total beam length
        nb_section=3,  # Number of sections for physics
        nb_frames=3,  # Number of frames for visualization
    )

    # Create geometry object - this automatically calculates all the geometry!
    beam_geometry = CosseratGeometry(beam_geometry_params)

    print(f"✨ Created beam with:")
    print(f"   - Length: {beam_geometry.get_beam_length()}")
    print(f"   - Sections: {beam_geometry.get_number_of_sections()}")
    print(f"   - Frames: {beam_geometry.get_number_of_frames()}")
    print(f"   - Section lengths: {beam_geometry.section_lengths}")

    # Create rigid base
    base_node = _add_rigid_base(root_node, node_name="rigid_base")

    # Custom bending states for this tutorial (slight bend)
    # Each vector [κ₁, κ₂, κ₃] represents the material curvature at each section.
    # Here, we apply a small constant curvature in the z-direction (κ₃ = 0.1)
    # to give the beam a slight initial bend. This helps visualize the
    # beam's orientation and response to forces.
    custom_bending_states = [
        [0.0, 0.0, 0.1],  # Section 1: slight bend in y and z
        [0.0, 0.0, 0.1],  # Section 2: slight bend in y and z
        [0.0, 0.0, 0.1],  # Section 3: slight bend in y and z
    ]

    # Create cosserat state using the geometry object
    bending_node = _add_cosserat_state(root_node, beam_geometry,
                                       custom_bending_states=custom_bending_states)

    # Create cosserat frame using the geometry object
    _add_cosserat_frame(base_node, bending_node, beam_geometry, beam_mass=0.0)

    return root_node