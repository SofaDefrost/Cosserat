# -*- coding: utf-8 -*-
"""
Tutorial 02: Cosserat Beam with Forces
=====================================

This tutorial builds on Tutorial 01 by adding:
- Gravity forces
- Applied forces at the beam tip
- Mass distribution
- Solver configuration for dynamic simulation

"""

import os
import sys

# Add the python package to the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "python"))

from cosserat import BeamGeometryParameters, CosseratGeometry

from basic_functions import (_add_cosserat_frame, _add_cosserat_state,
                                    _add_rigid_base, add_mini_header)

v_damping_param: float =  8e-1  # Damping parameter for dynamics
beam_mass: float = 5
nb_section: int = 12
beam_length: int = 10


def createScene(root):
    """Create a Cosserat beam scene with forces and dynamics."""
    # Configure scene with time integration
    add_mini_header(root)

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


    ## base node
    rigid_base = _add_rigid_base(solver, node_name="rigid_base")
    
    ## frame node
    frame_node = _add_cosserat_frame(solver, beam_geometry, beam_mass=beam_mass)

    ## bending node
    strain_node = _add_cosserat_state(rigid_base, frame_node, beam_geometry)


    return root

