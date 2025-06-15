"""Cosserat Beam Python Module

This module provides the core classes and functions for Cosserat beam simulations.
"""

from .beam import CosseratBase
from .geometry import CosseratGeometry, calculate_beam_parameters, calculate_frame_parameters, generate_edge_list
from .params import BeamGeometryParameters, BeamPhysicsParameters, Parameters
from .utils import addEdgeCollision, addPointsCollision, create_rigid_node
from .header import addHeader, addVisual, addSolverNode

__all__ = [
    'CosseratBase',
    'CosseratGeometry',
    'calculate_beam_parameters',
    'calculate_frame_parameters', 
    'generate_edge_list',
    'BeamGeometryParameters',
    'BeamPhysicsParameters',
    'Parameters',
    'addEdgeCollision',
    'addPointsCollision',
    'create_rigid_node',
    'addHeader',
    'addVisual',
    'addSolverNode'
]

