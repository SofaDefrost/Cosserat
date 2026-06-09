"""
The frame positions are obtained using the Strain2FramesCosseratMapping with strain of section i equal to [0, 0, i*0.1, 0, 0, 0].
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
nb_section: int = 8
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
    
    # ## frame node
# Frame  : 0 = SE3(R=SO3(quat=[1.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000]), t=[0.0000000000000000, 0.0000000000000000, 0.0000000000000000])
# Frame  : 1 = SE3(R=SO3(quat=[1.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000]), t=[0.6250000000000000, 0.0000000000000000, 0.0000000000000000])
# Frame  : 2 = SE3(R=SO3(quat=[0.9995117584851364, 0.0000000000000000, 0.0000000000000000, 0.0312449139853261]), t=[1.2495931784238019, 0.0195248929990088, 0.0000000000000000])
# Frame  : 3 = SE3(R=SO3(quat=[0.9956086864580018, 0.0000000000000000, 0.0000000000000000, 0.0936127312355129]), t=[1.8693130730232506, 0.0973958809932285, 0.0000000000000000])
# Frame  : 4 = SE3(R=SO3(quat=[0.9824733131012553, 0.0000000000000000, 0.0000000000000000, 0.1864032967622699]), t=[2.4688771807691765, 0.2706148516230317, 0.0000000000000000])
# Frame  : 5 = SE3(R=SO3(quat=[0.9515679480481722, 0.0000000000000000, 0.0000000000000000, 0.3074385145803810]), t=[3.0159390404052129, 0.5694761076407728, 0.0000000000000000])
# Frame  : 6 = SE3(R=SO3(quat=[0.8921336993669943, 0.0000000000000000, 0.0000000000000000, 0.4517714714916838]), t=[3.4579067110456743, 1.0077921964662537, 0.0000000000000000])
# Frame  : 7 = SE3(R=SO3(quat=[0.7922858596771785, 0.0000000000000000, 0.0000000000000000, 0.6101500770757915]), t=[3.7258157917714865, 1.5684110434723633, 0.0000000000000000])
# Frame  : 8 = SE3(R=SO3(quat=[0.6409968581633251, 0.0000000000000000, 0.0000000000000000, 0.7675435022360272]), t=[3.7503292063111386, 2.1879536470985115, 0.0000000000000000])

    frame_positions = [[0, 0, 0, 0, 0, 0, 1],
                       [0.625, 0, 0, 0, 0, 0, 1],
                       [1.2495931784238019, 0.0195248929990088, 0, 0, 0, 0.0312449139853261, 0.9995117584851364],
                       [1.8693130730232506, 0.0973958809932285, 0, 0, 0, 0.0936127312355129, 0.9956086864580018],
                       [2.4688771807691765, 0.2706148516230317, 0, 0, 0, 0.1864032967622699, 0.9824733131012553],
                       [3.0159390404052129, 0.5694761076407728, 0, 0, 0, 0.3074385145803810, 0.9515679480481722],
                       [3.4579067110456743, 1.0077921964662537, 0, 0, 0, 0.4517714714916838, 0.8921336993669943],
                       [3.7258157917714865, 1.5684110434723633, 0, 0, 0, 0.6101500770757915, 0.7922858596771785],
                       [3.7503292063111386, 2.1879536470985115, 0, 0, 0, 0.7675435022360272, 0.6409968581633251]]
                   
    
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

    #### Beam 2 (Red beam using Strain2FramesCosseratMapping) ####
    
    #Using same beam parameters
    beam_geometry2 = CosseratGeometry(beam_geometry_params)

    base_node2 = root.addChild("base_node2")
    base_node2.addObject("MechanicalObject", template="Rigid3d", name="base_mo2", 
                        position=[0., 0., 0., 0., 0., 0., 1.], showIndices="1", showObject="1", showObjectScale=0.1)
    
    base_node2.addObject("RestShapeSpringsForceField", name="spring2", stiffness=stiffness_param, angularStiffness=stiffness_param,
                         external_points="0", mstate="@base_mo2", points="0", template="Rigid3d")
        
    custom_bending_states = [[0, 0, i*0.1, 1, 0, 0] for i in range(nb_section)]
    
    bending_node = root.addChild("bending_node")
    bending_node.addObject("MechanicalObject", template="Vec6d", position=custom_bending_states, name="bending_state")
    

    # bending_node.addObject("Monitor", name="Monitor_Strain2Rigid", template="Vec6d", 
    #                        listening=True, indices=indices_str, showPositions=True, 
    #                        ExportPositions=True, ExportVelocities=False, 
    #                        ExportForces=False, fileName="monitor_strain2rigid")

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
        "Strain2FramesCosseratMapping",
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
