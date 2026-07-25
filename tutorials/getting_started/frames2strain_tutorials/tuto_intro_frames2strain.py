"""
Initially curved beam:

The frame positions are obtained using the Strain2FramesCosseratMapping with strain of section i equal to [0, 0, i*0.1, 0, 0, 0].
Remark that in this mapping (same for the DiscreteCosseratMapping) the default elongation is 0 
but for the force field (BeamHookeLawForceField or HookeSeratPCSForceField) the default elongation is fixed to 1. 
"""
beam_mass: float = 0
nb_section: int = 8
beam_length: float = 10

import os
import sys

# Add the python package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "python"))


from cosserat import BeamGeometryParameters, CosseratGeometry

from basic_functions import (add_mini_header, _add_cosserat_frame, _add_cosserat_state, _add_rigid_base)
def createScene(root):

    ##required plugins
    add_mini_header(root)

    root.gravity = [0, 0, 0]

    ## solver node
    solver = root.addChild("solver_node")
    solver.addObject("EulerImplicitSolver", firstOrder="0", rayleighMass=0.1, rayleighStiffness=0.1)
    solver.addObject("CGLinearSolver", iterations=100, tolerance=1e-7, threshold=1e-7)

    beam_geometry_params = BeamGeometryParameters(
        beam_length=beam_length,  # Total beam length
        nb_section=nb_section,
        nb_frames=nb_section,
    )    
        
    beam_geometry = CosseratGeometry(beam_geometry_params)
   

    ## base node
    rigid_base = _add_rigid_base(solver, node_name="rigid_base")
    
    ## frame node
    # (out) Frame [0] :0 0 0 0 0 0 1
    # (out) Frame [1] :1.25 0 0 0 0 0 1
    # (out) Frame [2] :2.49675 0.0780233 0 0 0 0.0624593 0.998048
    # (out) Frame [3] :3.70474 0.386474 0 0 0 0.186403 0.982473
    # (out) Frame [4] :4.75596 1.0492 0 0 0 0.366273 0.930508
    # (out) Frame [5] :5.42432 2.09012 0 0 0 0.585097 0.810963
    # (out) Frame [6] :5.43452 3.31983 0 0 0 0.806081 0.591805
    # (out) Frame [7] :4.66758 4.26979 0 0 0 0.966827 0.255434
    # (out) Frame [8] :3.46086 4.36543 0 0 0 0.983986 -0.178246

    frame_positions = [[0, 0, 0, 0, 0, 0, 1], 
                       [1.25, 0, 0, 0, 0, 0, 1],
                       [2.49675, 0.0780233, 0, 0, 0, 0.0624593, 0.998048],
                       [3.70474, 0.386474, 0, 0, 0, 0.186403, 0.982473],
                       [4.75596, 1.0492, 0, 0, 0, 0.366273, 0.930508],
                       [5.42432, 2.09012, 0, 0, 0, 0.585097, 0.810963],
                       [5.43452, 3.31983, 0, 0, 0, 0.806081, 0.591805],
                       [4.66758, 4.26979, 0, 0, 0, 0.966827, 0.255434],
                       [3.46086, 4.36543, 0, 0, 0, 0.983986, -0.178246]]
    
    beam_geometry.frames = frame_positions
    frame_node = _add_cosserat_frame(solver, beam_geometry, beam_mass=beam_mass)

    ## bending node
    strain_node = _add_cosserat_state(rigid_base, frame_node, beam_geometry)


    return root
