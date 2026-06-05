import os
import sys

# Add the python package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "python"))

from cosserat import BeamGeometryParameters, CosseratGeometry

from basic_functions import (_add_cosserat_frame, _add_cosserat_state,
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
    
    ## frame node
# Frame  : 0 = SE3(R=SO3(quat=[1.0000000000000000, 0.0000000000000000, 0.0000000000000000, 0.0000000000000000]), t=[0.0000000000000000, 0.0000000000000000, 0.0000000000000000])
# Frame  : 1 = SE3(R=SO3(quat=[0.9950041652780258, 0.0000000000000000, 0.0000000000000000, 0.0998334166468282]), t=[0.9933466539753060, 0.0996671107937919, 0.0000000000000000])
# Frame  : 2 = SE3(R=SO3(quat=[0.9800665778412416, 0.0000000000000000, 0.0000000000000000, 0.1986693307950612]), t=[1.9470917115432522, 0.3946950299855747, 0.0000000000000000])
# Frame  : 3 = SE3(R=SO3(quat=[0.9553364891256060, 0.0000000000000000, 0.0000000000000000, 0.2955202066613396]), t=[2.8232123669751763, 0.8733219254516087, 0.0000000000000000])
# Frame  : 4 = SE3(R=SO3(quat=[0.9210609940028851, 0.0000000000000000, 0.0000000000000000, 0.3894183423086505]), t=[3.5867804544976130, 1.5164664532641732, 0.0000000000000000])
# Frame  : 5 = SE3(R=SO3(quat=[0.8775825618903728, 0.0000000000000000, 0.0000000000000000, 0.4794255386042030]), t=[4.2073549240394819, 2.2984884706593016, 0.0000000000000000])


    frames_pos = [[0, 0, 0, 0, 0, 0, 1],
                  [0.9933466539753060, 0.0996671107937919, 0, 0, 0, 0.0998334166468282, 0.9950041652780258],
                  [1.9470917115432522, 0.3946950299855747, 0, 0, 0, 0.1986693307950612, 0.9800665778412416],
                  [2.8232123669751763, 0.8733219254516087, 0, 0, 0, 0.2955202066613396, 0.9553364891256060],
                  [3.5867804544976130, 1.5164664532641732, 0, 0, 0, 0.3894183423086505, 0.9210609940028851],
                  [4.2073549240394819, 2.2984884706593016, 0, 0, 0, 0.4794255386042030, 0.8775825618903728]]
    
    beam_geometry.frames = frames_pos
    frame_node = _add_cosserat_frame(rigid_base, beam_geometry, beam_mass=beam_mass)

    bending_node = _add_cosserat_state(rigid_base, frame_node, beam_geometry, youngModulus=0., poissonRatio=0.)
    return root