"""Cosserat Beam Python Module

This module provides the core classes and functions for Cosserat beam simulations.
"""

from .beam import CosseratBase
from .geometry import (CosseratGeometry, calculate_beam_parameters,
                       calculate_frame_parameters, generate_edge_list)
from .header import addHeader, addSolverNode, addVisual
from .params import (BeamGeometryParameters, BeamPhysicsBaseParameters,
                     BeamPhysicsParametersNoInertia, Parameters)
from .utils import addEdgeCollision, addPointsCollision, create_rigid_node

__all__ = [
    'CosseratBase',
    'CosseratGeometry',
    'calculate_beam_parameters',
    'calculate_frame_parameters',
    'generate_edge_list',
    'BeamGeometryParameters',
    'BeamPhysicsBaseParameters',
    'Parameters',
    'addEdgeCollision',
    'addPointsCollision',
    'create_rigid_node',
    'addHeader',
    'addVisual',
    'addSolverNode'
]

