# -*- coding: utf-8 -*-
"""
Debug script to isolate the segmentation fault
"""

import os
import sys

# Add the python package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "python"))

from cosserat import (BeamGeometryParameters, BeamPhysicsBaseParameters,
                      CosseratGeometry)

# Global parameters
stiffness_param: float = 1.0e10
beam_radius: float = 1.0

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

def _add_rigid_base(p_node, node_name="rigid_base", positions=None):
    """Create a rigid base node for the beam."""
    if positions is None:
        positions = [0, 0, 0, 0, 0, 0, 1]
    
    print(f"DEBUG: Creating rigid base node: {node_name}")
    rigid_base_node = p_node.addChild(node_name)
    
    if rigid_base_node is None:
        print(f"ERROR: Failed to create rigid base node: {node_name}")
        return None
    
    print(f"DEBUG: Adding MechanicalObject to {node_name}")
    rigid_base_node.addObject(
        "MechanicalObject",
        template="Rigid3d",
        name="cosserat_base_mo",
        position=positions,
        showObject=True,
        showObjectScale="0.1",
    )
    
    print(f"DEBUG: Adding RestShapeSpringsForceField to {node_name}")
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
    print(f"DEBUG: Successfully created rigid base node: {node_name}")
    return rigid_base_node

def _add_cosserat_state_v2_debug(p_node, geometry: CosseratGeometry, node_name="cosserat_coordinates", custom_bending_states=None):
    """Create the cosserat coordinate node using CosseratGeometry - DEBUG VERSION."""
    print(f"DEBUG: Creating cosserat state node: {node_name}")
    
    if p_node is None:
        print("ERROR: Parent node is None")
        return None
    
    if geometry is None:
        print("ERROR: Geometry is None")
        return None
        
    cosserat_coordinate_node = p_node.addChild(node_name)
    
    if cosserat_coordinate_node is None:
        print(f"ERROR: Failed to create cosserat coordinate node: {node_name}")
        return None

    # Use geometry data or custom bending states
    bending_states = (
        custom_bending_states if custom_bending_states else geometry.bendingState
    )
    
    print(f"DEBUG: Bending states: {bending_states}")
    print(f"DEBUG: Section lengths: {geometry.section_lengths}")

    print(f"DEBUG: Adding MechanicalObject to {node_name}")
    cosserat_coordinate_node.addObject(
        "MechanicalObject",
        template="Vec3d",
        name="cosserat_state",
        position=bending_states,
    )
    
    print(f"DEBUG: Adding HookeSeratPCSForceField to {node_name}")
    try:
        cosserat_coordinate_node.addObject(
            "HookeSeratPCSForceField",
            crossSectionShape="circular",
            length=geometry.section_lengths,  # Use geometry data
            radius=2.0,
            youngModulus=1.0e3,
            poissonRatio=0.4,
        )
        print(f"DEBUG: Successfully added HookeSeratPCSForceField")
    except Exception as e:
        print(f"ERROR: Failed to add HookeSeratPCSForceField: {e}")
        return None
    
    print(f"DEBUG: Successfully created cosserat state node: {node_name}")
    return cosserat_coordinate_node

def _add_cosserat_frame_v2_debug(
    p_node, bending_node, geometry: CosseratGeometry, node_name="cosserat_in_Sofa_frame_node", beam_mass=0.0
):
    """Create the cosserat frame node using CosseratGeometry - DEBUG VERSION."""
    print(f"DEBUG: Creating Cosserat frame node: {node_name}")
    
    # Check inputs
    if p_node is None:
        print("ERROR: p_node is None")
        return None
    if bending_node is None:
        print("ERROR: bending_node is None")
        return None
    if geometry is None:
        print("ERROR: geometry is None")
        return None
    
    print(f"DEBUG: Geometry frames: {geometry.frames}")
    print(f"DEBUG: Geometry curv_abs_sections: {geometry.curv_abs_sections}")
    print(f"DEBUG: Geometry curv_abs_frames: {geometry.curv_abs_frames}")
    
    cosserat_in_sofa_frame_node = p_node.addChild(node_name)
    if cosserat_in_sofa_frame_node is None:
        print(f"ERROR: Failed to create frame node: {node_name}")
        return None
    
    print(f"DEBUG: Adding bending_node as child")
    bending_node.addChild(cosserat_in_sofa_frame_node)

    print(f"DEBUG: Adding MechanicalObject FramesMO")
    frames_mo = cosserat_in_sofa_frame_node.addObject(
        "MechanicalObject",
        template="Rigid3d",
        name="FramesMO",
        position=geometry.frames,  # Use geometry data
        showIndices=1,
        showObject=1,
        showObjectScale=0.8,
    )
    
    if frames_mo is None:
        print("ERROR: Failed to create FramesMO")
        return None

    print(f"DEBUG: Adding UniformMass")
    cosserat_in_sofa_frame_node.addObject("UniformMass", totalMass=beam_mass)

    # Check link paths before creating mapping
    print(f"DEBUG: Checking link paths")
    if hasattr(bending_node, 'cosserat_state'):
        bending_link = bending_node.cosserat_state.getLinkPath()
        print(f"DEBUG: Bending node link: {bending_link}")
    else:
        print("ERROR: bending_node does not have cosserat_state")
        return None
        
    if hasattr(p_node, 'cosserat_base_mo'):
        base_link = p_node.cosserat_base_mo.getLinkPath()
        print(f"DEBUG: Base node link: {base_link}")
    else:
        print("ERROR: p_node does not have cosserat_base_mo")
        return None
    
    frames_link = frames_mo.getLinkPath()
    print(f"DEBUG: Frames link: {frames_link}")

    print(f"DEBUG: Creating HookeSeratDiscretMapping")
    try:
        mapping = cosserat_in_sofa_frame_node.addObject(
            "HookeSeratDiscretMapping",
            curv_abs_input=geometry.curv_abs_sections,  # Use geometry data
            curv_abs_output=geometry.curv_abs_frames,  # Use geometry data
            name="cosseratMapping",
            input1=bending_node.cosserat_state.getLinkPath(),
            input2=p_node.cosserat_base_mo.getLinkPath(),
            output=frames_mo.getLinkPath(),
            debug=1,
            radius=beam_radius,
        )
        print(f"DEBUG: Successfully created HookeSeratDiscretMapping")
    except Exception as e:
        print(f"ERROR: Failed to create HookeSeratDiscretMapping: {e}")
        return None
    
    print(f"DEBUG: Successfully created cosserat frame node: {node_name}")
    return cosserat_in_sofa_frame_node

def createScene(root_node):
    """Debug version of createScene"""
    print("DEBUG: Starting createScene")
    
    # Load required plugins
    add_mini_header(root_node)

    # add gravity
    root_node.gravity = [0, 0.0, 0]

    print("DEBUG: Creating beam geometry")
    # Define beam geometry parameters
    beam_geometry_params = BeamGeometryParameters(
        beam_length=30.0,  # Total beam length
        nb_section=3,  # Number of sections for physics
        nb_frames=3,  # Number of frames for visualization
    )

    # Create geometry object - this automatically calculates all the geometry!
    beam_geometry = CosseratGeometry(beam_geometry_params)
    
    if beam_geometry is None:
        print("ERROR: Failed to create beam geometry")
        return root_node

    print(f"DEBUG: Created beam with:")
    print(f"   - Length: {beam_geometry.get_beam_length()}")
    print(f"   - Sections: {beam_geometry.get_number_of_sections()}")
    print(f"   - Frames: {beam_geometry.get_number_of_frames()}")
    print(f"   - Section lengths: {beam_geometry.section_lengths}")

    # Custom bending states for this tutorial (slight bend)
    custom_bending_states = [
        [0.0, 0.0, 0.1],  # Section 1: slight bend in y and z
        [0.0, 0.0, 0.1],  # Section 2: slight bend in y and z
        [0.0, 0.0, 0.1],  # Section 3: slight bend in y and z
    ]

    print("DEBUG: Creating rigid base")
    # Create beam with old Component
    base_node_2 = _add_rigid_base(root_node, node_name="rigid_base_2")
    
    if base_node_2 is None:
        print("ERROR: Failed to create base node")
        return root_node

    print("DEBUG: Creating cosserat state")
    # Create cosserat states using old method
    bending_node_2 = _add_cosserat_state_v2_debug(root_node, beam_geometry, node_name="c_state_2", custom_bending_states=custom_bending_states)

    if bending_node_2 is None:
        print("ERROR: Failed to create bending node")
        return root_node

    print("DEBUG: Creating cosserat frame - THIS IS LINE 246 EQUIVALENT")
    # This is the equivalent of line 246 in the original script
    frame_node = _add_cosserat_frame_v2_debug(base_node_2, bending_node_2, beam_geometry, node_name="frame_2", beam_mass=1.0)
    
    if frame_node is None:
        print("ERROR: Failed to create frame node")
        return root_node

    print("DEBUG: Scene created successfully")
    return root_node
