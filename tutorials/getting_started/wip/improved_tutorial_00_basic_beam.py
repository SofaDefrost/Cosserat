# -*- coding: utf-8 -*-
"""
Tutorial 00: Basic Cosserat Beam
===============================

This tutorial introduces the fundamentals of creating a Cosserat beam in SOFA.
We'll build a simple, static beam with a fixed base and a slight bend.

Key concepts covered:
- Creating a rigid base for the beam
- Setting up Cosserat coordinates (sections) with the BeamHookeLawForceField
- Mapping from Cosserat coordinates to SOFA frames
- Understanding the difference between sections and frames

Note: This is a static scene with no gravity or external forces.
"""

import os
import sys

# Add the python package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "python"))

from cosserat import (BeamGeometryParameters, BeamPhysicsBaseParameters,
                      CosseratGeometry)

# Global parameters
stiffness_param: float = 1.0e10  # Stiffness for the base spring
beam_radius: float = 1.0         # Radius of the beam visualization


def add_required_plugins(root_node):
    """Add all required SOFA plugins to the scene."""
    root_node.addObject("RequiredPlugin", name="Sofa.Component.Mass")
    root_node.addObject("RequiredPlugin", name="Sofa.Component.SolidMechanics.Spring")
    root_node.addObject("RequiredPlugin", name="Sofa.Component.StateContainer")
    root_node.addObject("RequiredPlugin", name="Sofa.Component.Visual")
    root_node.addObject("RequiredPlugin", name="Cosserat")  # Our special plugin!

    # Configure visualization
    root_node.addObject(
        "VisualStyle",
        displayFlags="showBehaviorModels showCollisionModels showMechanicalMappings",
    )


def create_rigid_base(parent_node, node_name="rigid_base", positions=None):
    """
    Create a rigid base node that anchors the beam.
    
    Parameters:
        parent_node: The SOFA node to attach this to
        node_name: Name for the node
        positions: Initial position and orientation [x,y,z,qx,qy,qz,qw], defaults to origin
    
    Returns:
        The created rigid base node
    """
    if positions is None:
        positions = [0, 0, 0, 0, 0, 0, 1]  # Default at origin with identity quaternion
        
    rigid_base_node = parent_node.addChild(node_name)
    
    # Add a mechanical object to represent the rigid base
    rigid_base_node.addObject(
        "MechanicalObject",
        template="Rigid3d",            # Rigid body template (position + orientation)
        name="cosserat_base_mo",       # Name we'll reference later
        position=positions,            # Initial position and orientation
        showObject=True,               # Display in the viewer
        showObjectScale="0.1",         # Size of the display
    )
    
    # Add a spring to fix the base in place
    rigid_base_node.addObject(
        "RestShapeSpringsForceField",
        name="spring",
        stiffness=stiffness_param,         # Linear stiffness
        angularStiffness=stiffness_param,  # Angular stiffness
        external_points="0",               # External point index
        mstate="@cosserat_base_mo",        # Which mechanical state to use
        points="0",                        # Which point to fix
        template="Rigid3d",                # Template must match the mechanical object
    )
    
    return rigid_base_node


def create_cosserat_state(parent_node, geometry, node_name="cosserat_coordinates", 
                         custom_bending_states=None):
    """
    Create the Cosserat coordinate node that holds the reduced coordinates (strains).
    
    Parameters:
        parent_node: The SOFA node to attach this to
        geometry: CosseratGeometry object with precalculated geometry
        node_name: Name for the node
        custom_bending_states: Optional custom strain values
        
    Returns:
        The created Cosserat coordinates node
    """
    cosserat_coordinate_node = parent_node.addChild(node_name)

    # Use geometry data or custom bending states
    bending_states = (
        custom_bending_states if custom_bending_states else geometry.bendingState
    )

    # Add a mechanical object to hold the strain variables
    cosserat_coordinate_node.addObject(
        "MechanicalObject",
        template="Vec3d",           # 3D vector template
        name="cosserat_state",      # Name we'll reference in mapping
        position=bending_states,    # Initial strain values
    )
    
    # Add a force field that implements Hooke's Law for the beam
    cosserat_coordinate_node.addObject(
        "BeamHookeLawForceField",
        crossSectionShape="circular",        # Cross-section shape
        length=geometry.section_lengths,     # Length of each section
        radius=2.0,                          # Physical radius
        youngModulus=1.0e3,                  # Material stiffness
        poissonRatio=0.4,                    # Material property
    )
    
    return cosserat_coordinate_node


def create_cosserat_frame(
    parent_node, bending_node, geometry, node_name="cosserat_in_Sofa_frame_node", 
    beam_mass=0.0
):
    """
    Create the node that maps from Cosserat coordinates to SOFA frames.
    
    Parameters:
        parent_node: First parent node (for rigid base)
        bending_node: Second parent node (for Cosserat coordinates)
        geometry: CosseratGeometry object with precalculated geometry
        node_name: Name for the node
        beam_mass: Mass to distribute across the beam
        
    Returns:
        The created frames node
    """
    # Create a node that will be a child of both the rigid base and Cosserat coordinates
    cosserat_in_sofa_frame_node = parent_node.addChild(node_name)
    bending_node.addChild(cosserat_in_sofa_frame_node)

    # Add a mechanical object to represent the frames along the beam
    frames_mo = cosserat_in_sofa_frame_node.addObject(
        "MechanicalObject",
        template="Rigid3d",              # Rigid body template
        name="FramesMO",                 # Name for referencing
        position=geometry.frames,        # Precalculated frame positions
        showIndices=1,                   # Show indices in the viewer
        showObject=1,                    # Show visual representation
        showObjectScale=0.8,             # Size of visual objects
    )

    # Add mass if specified
    cosserat_in_sofa_frame_node.addObject("UniformMass", totalMass=beam_mass)

    # Add the mapping from Cosserat coordinates to frames
    cosserat_in_sofa_frame_node.addObject(
        "DiscreteCosseratMapping",
        curv_abs_input=geometry.curv_abs_sections,  # Curvilinear abscissa for sections
        curv_abs_output=geometry.curv_abs_frames,   # Curvilinear abscissa for frames
        name="cosseratMapping",
        input1=bending_node.cosserat_state.getLinkPath(),  # Link to strain variables
        input2=parent_node.cosserat_base_mo.getLinkPath(), # Link to rigid base
        output=frames_mo.getLinkPath(),                    # Link to output frames
        debug=0,                                           # Debug level
        radius=beam_radius,                                # Visualization radius
    )
    
    return cosserat_in_sofa_frame_node


def createScene(root_node):
    """
    Create a basic Cosserat beam scene using the CosseratGeometry class.
    
    This scene demonstrates:
    1. How to create a basic beam with the Cosserat plugin
    2. The relationship between sections (physics) and frames (visualization)
    3. How to set up the beam with a slight bend
    """
    # Load required plugins
    add_required_plugins(root_node)

    # No gravity in this simple example
    root_node.gravity = [0, 0.0, 0]

    # Define beam geometry parameters
    beam_geometry_params = BeamGeometryParameters(
        beam_length=30.0,  # Total beam length
        nb_section=3,      # Number of sections for physics
        nb_frames=3,       # Number of frames for visualization (matches sections)
    )

    # Create geometry object - this automatically calculates all needed values
    beam_geometry = CosseratGeometry(beam_geometry_params)

    # Print info about the beam
    print(f"✨ Created beam with:")
    print(f"   - Length: {beam_geometry.get_beam_length()}")
    print(f"   - Sections: {beam_geometry.get_number_of_sections()}")
    print(f"   - Frames: {beam_geometry.get_number_of_frames()}")
    print(f"   - Section lengths: {beam_geometry.section_lengths}")

    # Create rigid base
    base_node = create_rigid_base(root_node)

    # Define custom bending states for this tutorial (slight bend in z-direction)
    custom_bending_states = [
        [0.0, 0.0, 0.1],  # Section 1: bend around z-axis
        [0.0, 0.0, 0.1],  # Section 2: bend around z-axis
        [0.0, 0.0, 0.1],  # Section 3: bend around z-axis
    ]

    # Create Cosserat state using the geometry object and custom bending
    bending_node = create_cosserat_state(
        root_node, beam_geometry, custom_bending_states=custom_bending_states
    )

    # Create Cosserat frame mapping
    create_cosserat_frame(base_node, bending_node, beam_geometry, beam_mass=0.0)

    return root_node
