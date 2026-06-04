# -*- coding: utf-8 -*-
"""
Basic functions to simulate scenes using Frames2StrainCosseratMapping
"""

import os
import sys

# Add the python package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "python"))

from cosserat import (BeamGeometryParameters, BeamPhysicsBaseParameters,
                      CosseratGeometry)

# Global parameters
stiffness_param: float = 1.0e10
beam_radius: float = 0.5


def add_mini_header(root_node):
    root_node.addObject("RequiredPlugin", pluginName="Sofa.Component.Mass")
    root_node.addObject("RequiredPlugin", pluginName="Sofa.Component.SolidMechanics.Spring")
    root_node.addObject("RequiredPlugin", pluginName="Sofa.Component.StateContainer")
    root_node.addObject("RequiredPlugin", pluginName="Sofa.Component.Visual")
    root_node.addObject("RequiredPlugin", pluginName="Cosserat")
    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.LinearSolver.Direct') # Needed to use components [SparseLDLSolver]  
    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.ODESolver.Backward') # Needed to use components [EulerImplicitSolver]  
    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.LinearSolver.Iterative') # Needed to use components [CGLinearSolver]  
    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.MechanicalLoad') # Needed to use components [ConstantForceField]  
    root_node.addObject('RequiredPlugin', pluginName='Sofa.Component.Constraint.Projective') # Needed to use components [FixedProjectiveConstraint]

    # Configure scene
    root_node.addObject(
        "VisualStyle",
        displayFlags="showBehaviorModels showCollisionModels showMechanicalMappings",
    )
    root_node.addObject("DefaultAnimationLoop")


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


def _add_cosserat_frame(
    p_node, geometry: CosseratGeometry, beam_mass=0.0, constraint=True, node_name="frame_node"
):
    """Create the cosserat frame node using CosseratGeometry."""
    frame_node = p_node.addChild(node_name)
    frames_mo = frame_node.addObject(
        "MechanicalObject",
        template="Rigid3d",
        name="FramesMO",
        position=geometry.frames,  # Use geometry data
        showIndices=1,
        showObject=1,
        showObjectScale=0.8,
    )

    if constraint:
        frame_node.addObject("FixedProjectiveConstraint", indices="0")
    
    frame_node.addObject("UniformMass", totalMass=beam_mass)

    return frame_node


def _add_cosserat_state(p_node, f_node, geometry: CosseratGeometry, custom_bending_states=None, node_name="bending_node", radius=0.5, youngModulus=1e3, poissonRatio=0.4):
    """Create the cosserat coordinate node using CosseratGeometry."""
    bending_node = f_node.addChild(node_name)

    bending_states = []
    for i in range(geometry.params.nb_section):
        bending_states.append([0, 0, 0, 1, 0, 0]) # strain à l'équilibre 

    if custom_bending_states is not None:
        bending_states = custom_bending_states

    bending_node.addObject(
        "MechanicalObject",
        template="Vec6d",
        name="cosserat_state",
        position=bending_states,
    )
    bending_node.addObject(
        "BeamHookeLawForceField",
        crossSectionShape="circular",
        length=geometry.section_lengths,  # Use geometry data
        radius=radius,
        youngModulus=youngModulus,
        poissonRatio=poissonRatio,
    )


    bending_node.addObject(
    "Frames2StrainCosseratMapping", 
    curv_abs_input=geometry.curv_abs_sections, 
    curv_abs_output=geometry.curv_abs_frames, 
    name="cosseratMapping", 
    input1=f_node.getLinkPath(), 
    input2=p_node.cosserat_base_mo.getLinkPath(), 
    output=bending_node.cosserat_state.getLinkPath(), 
    debug=0,
    radius=beam_radius, 
    color=[0., 1., 0., 0.5], #green   
    )

    return bending_node
