"""
The frame positions are obtained using the Strain2RigidCosseratMapping with strain of section i equal to [0, 0, i*0.1, 0, 0, 0].
Remark that in this mapping (same for the DiscreteCosseratMapping) the default elongation is 0 
but for the force field (BeamHookeLawForceField or HookeSeratPCSForceField) the default elongation is fixed to 1. 
"""

######## Pas de ForceField car on veut trouver l'erreur de la fonction apply ########

import os
import sys

# Add the python package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "python"))

from cosserat import BeamGeometryParameters, CosseratGeometry

from basic_functions import (_add_cosserat_frame, _add_cosserat_state,
                                    _add_rigid_base, add_mini_header)

stiffness_param: float = 1e10
v_damping_param: float =  8e-1  # Damping parameter for dynamics
beam_mass: float = 0
nb_section: int = 16
beam_length: int = 5


def createScene(root):
    """Create a Cosserat beam scene with forces and dynamics."""
    # Configure scene with time integration
    add_mini_header(root)
    root.addObject("RequiredPlugin", pluginName="SofaValidation")

    # Add gravity
    root.gravity = [0, 0, 0]  # Add gravity!

    beam_geometry_params = BeamGeometryParameters(
        beam_length=beam_length,
        nb_section=nb_section,
        nb_frames=nb_section,
    )    
        
    beam_geometry = CosseratGeometry(beam_geometry_params)

    # On génère la liste des indices que l'on veut monitorer (0, 1, 2... 7)
    indices_str = "   ".join([str(i) for i in range(nb_section)])

    ## base node
    rigid_base = _add_rigid_base(root, node_name="rigid_base")
    
    ## frame node

    frame_positions = [[0., 0., 0., 0., 0., 0., 1.],
                    [0.3125, 0., 0., 0., 0., 0., 1.],
                    [0.624886, 0.00732288, 0., 0., 0., 0.0234354, 0.999725],
                    [0.935899, 0.0365661, 0., 0., 0., 0.0702546, 0.997529],
                    [1.24122, 0.101942, 0., 0., 0., 0.140162, 0.990129],
                    [1.53158, 0.216235, 0., 0., 0., 0.232235, 0.97266],
                    [1.79136, 0.388646, 0., 0., 0., 0.344365, 0.938836],
                    [1.99838, 0.621359, 0., 0., 0., 0.472555, 0.881301],
                    [2.1259, 0.905121, 0., 0., 0., 0.61015, 0.792286],
                    [2.14788, 1.21501, 0., 0., 0., 0.747141, 0.664666],
                    [2.04806, 1.5087, 0., 0., 0., 0.869746, 0.4935],
                    [1.83185, 1.73036, 0., 0., 0., 0.960575, 0.278022],
                    [1.53712, 1.82336, 0., 0., 0., 0.999714, 0.023919],
                    [1.23709, 1.75203, 0., 0., 0., 0.967073, -0.254498],
                    [1.02703, 1.5272, 0., 0., 0., -0.846182, 0.532893],
                    [0.990702, 1.22244, 0., 0., 0., -0.629302, 0.777161],
                    [1.15357, 0.963264, 0., 0., 0., -0.323185, 0.946336]
                    ]
    
    beam_geometry.frames = frame_positions
    frame_node = _add_cosserat_frame(root, beam_geometry, beam_mass=beam_mass)

    ## bending node
    custom_bending_states = [[0, 0, 0, 1, 0, 0] for _ in range(nb_section)]
    # strain_node = _add_cosserat_state(rigid_base, frame_node, beam_geometry, custom_bending_states=custom_bending_states)

    strain_node = frame_node.addChild("strain_node")
    strain_node.addObject(
        "MechanicalObject",
        template="Vec6d",
        name="cosserat_state",
        position=custom_bending_states,
    )

    strain_node.addObject(
    "Frames2StrainCosseratMapping", 
    curv_abs_input=beam_geometry.curv_abs_sections, 
    curv_abs_output=beam_geometry.curv_abs_frames, 
    name="cosseratMapping", 
    input1=frame_node.getLinkPath(), 
    input2=rigid_base.getLinkPath(), 
    output=strain_node.cosserat_state.getLinkPath(), 
    debug=0,
    radius=0.5, 
    color=[0., 1., 0., 0.5], #green   
    )


    strain_node.addObject("Monitor", name="Monitor_Frames2Strain", template="Vec6d", 
                           listening=True, indices=indices_str, showPositions=True, 
                           ExportPositions=True, ExportVelocities=False, 
                           ExportForces=False, fileName="monitor_frames2strain")

    #### Beam 2 (Red beam using Strain2RigidCosseratMapping) ####
    
    #Using same beam parameters
    beam_geometry2 = CosseratGeometry(beam_geometry_params)

    base_node2 = root.addChild("base_node2")
    base_node2.addObject("MechanicalObject", template="Rigid3d", name="base_mo2", 
                        position=[0., 0., 0., 0., 0., 0., 1.], showIndices="1", showObject="1", showObjectScale=0.1)
    
    base_node2.addObject("RestShapeSpringsForceField", name="spring2", stiffness=stiffness_param, angularStiffness=stiffness_param,
                         external_points="0", mstate="@base_mo2", points="0", template="Rigid3d")
        
    custom_bending_states = [[0, 0, i*0.15, 0, 0, 0] for i in range(nb_section)]
    
    bending_node = root.addChild("bending_node")
    bending_node.addObject("MechanicalObject", template="Vec6d", position=custom_bending_states, name="bending_state")
    

    bending_node.addObject("Monitor", name="Monitor_Strain2Rigid", template="Vec6d", 
                           listening=True, indices=indices_str, showPositions=True, 
                           ExportPositions=True, ExportVelocities=False, 
                           ExportForces=False, fileName="monitor_strain2rigid")

    frame_node2 = base_node2.addChild("frame_node")
    bending_node.addChild(frame_node2)

    frames_mo2 = frame_node2.addObject(
        "MechanicalObject",
        template="Rigid3d",
        name="FramesMO2",
        position=beam_geometry2.frames,  # Use geometry data
        showIndices=1,
        showObject=1,
        showObjectScale=0.8,
    )

    frame_node2.addObject("UniformMass", totalMass=beam_mass)

    frame_node2.addObject(
        "Strain2RigidCosseratMapping",
        curv_abs_input=beam_geometry2.curv_abs_sections,
        curv_abs_output=beam_geometry2.curv_abs_frames,
        name="cosseratMapping2",
        input1=bending_node.bending_state.getLinkPath(),
        input2=base_node2.base_mo2.getLinkPath(),
        output=frames_mo2.getLinkPath(),
        debug=0,
        radius=0.5,
        color=[1.0, 0.0, 0.0, 0.5], #red
    )

    return root
