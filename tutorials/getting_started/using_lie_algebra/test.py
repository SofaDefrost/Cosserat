import os
import sys

# Add the python package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "python"))

from cosserat import BeamGeometryParameters, CosseratGeometry

from introduction_and_setup import (_add_cosserat_frame_v2, _add_cosserat_state,
                                    _add_rigid_base, add_mini_header)

stiffness_param: float = 1e10
beam_mass: float = 0
nb_section: int = 5
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


    ## base node
    rigid_base = _add_rigid_base(root, node_name="rigid_base")
    
    ## bending node
    bending_state = [[0, 0, 0.2] for _ in range(nb_section)]
    bending_node = _add_cosserat_state(root, beam_geometry, custom_bending_states=bending_state)

    ## frame node
    frame_node = _add_cosserat_frame_v2(rigid_base, bending_node, beam_geometry, beam_mass=beam_mass)


    return root