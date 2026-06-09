"""
The frame positions are obtained using the Strain2FramesCosseratMapping with strain of section i equal to [0, 0, i*0.1, 0, 0, 0].
Remark that in this mapping (same for the DiscreteCosseratMapping) the default elongation is 0 
but for the force field (BeamHookeLawForceField or HookeSeratPCSForceField) the default elongation is fixed to 1. 
"""


import os
import sys

# Add the python package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "python"))

from cosserat import BeamGeometryParameters, CosseratGeometry

from basic_functions import (_add_cosserat_frame, _add_cosserat_state,
                                    _add_rigid_base, add_mini_header)

stiffness_param: float = 1e10
v_damping_param: float =  8e-1  # Damping parameter for dynamics
beam_mass: float = 1
nb_section: int = 8
beam_length: int = 5


def createScene(root):
    """Create a Cosserat beam scene with forces and dynamics."""
    # Configure scene with time integration
    add_mini_header(root)
    root.addObject("RequiredPlugin", pluginName="SofaValidation")

    # Add gravity
    root.gravity = [0, -9.81, 0]  # Add gravity!

    # Configure time integration and solver

    ## solver node
    solver = root.addChild("solver_node")
    solver.addObject("EulerImplicitSolver", firstOrder="0", rayleighMass=0.1, rayleighStiffness=0.1, vdamping=v_damping_param)
    solver.addObject("CGLinearSolver", iterations=1000, tolerance=1e-10, threshold=1e-10)

    beam_geometry_params = BeamGeometryParameters(
        beam_length=beam_length,
        nb_section=nb_section,
        nb_frames=nb_section,
    )    
        
    beam_geometry = CosseratGeometry(beam_geometry_params)

    # On génère la liste des indices que l'on veut monitorer (0, 1, 2... 7)
    indices_str = "   ".join([str(i) for i in range(nb_section)])

    ## base node
    rigid_base = _add_rigid_base(solver, node_name="rigid_base")
    
    ## frame node
    frame_positions = [[0., 0., 0., 0., 0., 0., 1.],
                    [0.625, 0., 0., 0., 0., 0., 1.],
                    [1.24959, 0.0195249, 0., 0., 0., 0.0312449, 0.999512],
                    [1.86931, 0.0973959, 0., 0., 0., 0.0936127, 0.995609],
                    [2.46888, 0.270615, 0., 0., 0., 0.186403, 0.982473],
                    [3.01594, 0.569476, 0., 0., 0., 0.307439, 0.951568],
                    [3.45791, 1.00779, 0., 0., 0., 0.451771, 0.892134],
                    [3.72582, 1.56841, 0., 0., 0., 0.61015, 0.792286],
                    [3.75033, 2.18795, 0., 0., 0., 0.767544, 0.640997]]
    
    beam_geometry.frames = frame_positions
    frame_node = _add_cosserat_frame(solver, beam_geometry, beam_mass=beam_mass)

    ## bending node
    custom_bending_states = [[0, 0, i*0.1, 1, 0, 0] for i in range(nb_section)]
    strain_node = _add_cosserat_state(rigid_base, frame_node, beam_geometry, custom_bending_states=custom_bending_states)

    strain_node.addObject("Monitor", name="Monitor_Frames2Strain", template="Vec6d", 
                           listening=True, indices=indices_str, showPositions=True, 
                           ExportPositions=True, ExportVelocities=False, 
                           ExportForces=False, fileName="monitor_frames2strain")

    #### Beam 2 (Red beam using Strain2FramesCosseratMapping) ####
    
    #Using same beam parameters

    base_node2 = solver.addChild("base_node2")
    base_node2.addObject("MechanicalObject", template="Rigid3d", name="base_mo2", 
                        position=[0., 0., 0., 0., 0., 0., 1.], showIndices="1", showObject="1", showObjectScale=0.1)
    
    base_node2.addObject("RestShapeSpringsForceField", name="spring2", stiffness=stiffness_param, angularStiffness=stiffness_param,
                         external_points="0", mstate="@base_mo2", points="0", template="Rigid3d")
        
    custom_bending_states = [[0, 0, i*0.1, 0, 0, 0] for i in range(nb_section)]
    
    bending_node = solver.addChild("bending_node")
    bending_node.addObject("MechanicalObject", template="Vec6d", position=custom_bending_states, name="bending_state")
    bending_node.addObject(
        "BeamHookeLawForceField",
        crossSectionShape="circular",
        length=beam_geometry.section_lengths,  # Use geometry data
        radius=0.5,
        youngModulus=1.0e3,
        poissonRatio=0.4,
    )    

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
        position=beam_geometry.frames,  # Use geometry data
        showIndices=1,
        showObject=1,
        showObjectScale=0.8,
    )

    frame_node2.addObject("UniformMass", totalMass=beam_mass)

    frame_node2.addObject(
        "Strain2FramesCosseratMapping",
        curv_abs_input=beam_geometry.curv_abs_sections,
        curv_abs_output=beam_geometry.curv_abs_frames,
        name="cosseratMapping2",
        input1=bending_node.bending_state.getLinkPath(),
        input2=base_node2.base_mo2.getLinkPath(),
        output=frames_mo2.getLinkPath(),
        debug=0,
        radius=0.5,
        color=[1.0, 0.0, 0.0, 0.5], #red
    )


    return root

# (in) Rigid base: 0 0 0 0 0 0 1
# (in) Strain [0]: 0 0 0 0 0 0
# (in) Strain [1]: 0 0 0.1 0 0 0
# (in) Strain [2]: 0 0 0.2 0 0 0
# (in) Strain [3]: 0 0 0.3 0 0 0
# (in) Strain [4]: 0 0 0.4 0 0 0
# (in) Strain [5]: 0 0 0.5 0 0 0
# (in) Strain [6]: 0 0 0.6 0 0 0
# (in) Strain [7]: 0 0 0.7 0 0 0
# (out) Frame [0] :0 0 0 0 0 0 1
# (out) Frame [1] :0.625 0 0 0 0 0 1
# (out) Frame [2] :1.24959 0.0195249 0 0 0 0.0312449 0.999512
# (out) Frame [3] :1.86931 0.0973959 0 0 0 0.0936127 0.995609
# (out) Frame [4] :2.46888 0.270615 0 0 0 0.186403 0.982473
# (out) Frame [5] :3.01594 0.569476 0 0 0 0.307439 0.951568
# (out) Frame [6] :3.45791 1.00779 0 0 0 0.451771 0.892134
# (out) Frame [7] :3.72582 1.56841 0 0 0 0.61015 0.792286
# (out) Frame [8] :3.75033 2.18795 0 0 0 0.767544 0.640997
