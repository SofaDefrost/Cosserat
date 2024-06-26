from dataclasses import dataclass, field
from typing import List, Literal, Union

@dataclass
class BeamGeometryParameters:
    """Cosserat Beam Geometry parameters"""
    beam_length: float = 1.0  # beam length in m
    nb_sections: int = 5  # number of sections along the beam length
    nb_frames: int = 30  # number of frames along the beam
    build_collision_model: bool = False
    init_pos: List[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])
    show_frames_object: bool = True
    show_rigid_base_object: bool = True

@dataclass
class BeamPhysicsParameters:
    """Cosserat Beam Physics parameters"""
    young_modulus: float = 1.205e8
    poisson_ratio: float = 0.3
    use_inertia: bool = False
    GI: float = 1.5708
    GA: float = 3.1416e4
    EI: float = 0.7854
    EA: float = 3.1416e4
    beam_mass: float = 1.0
    beam_radius: float = 0.01  # beam radius in m
    beam_length: float = 1.0  # beam length in m along the X axis
    beam_shape: Literal["circular", "rectangular"] = "circular"
    length_Y: float = 0.1  # length of the beam in the Y direction (for rectangular shape)
    length_Z: float = 0.1  # length of the beam in the Z direction (for rectangular shape)

@dataclass
class SimulationParameters:
    """Simulation parameters"""
    rayleigh_stiffness: float = 0.2
    rayleigh_mass: float = 0.1
    first_order: bool = False

@dataclass
class VisualParameters:
    """Visual parameters"""
    show_object: bool = True
    show_object_scale: float = 1.0
    show_object_color: List[float] = field(default_factory=lambda: [1.0, 0.0, 0.0, 1.0])

@dataclass
class ContactParameters:
    """Contact parameters"""
    response_params: str = "mu=0.8"
    response: str = "FrictionContactConstraint"
    alarm_distance: float = 0.05
    contact_distance: float = 0.01
    is_multithreading: bool = False
    tolerance: float = 1.0e-8
    max_iterations: int = 100
    epsilon: float = 1.0e-6

@dataclass
class CosseratBeamParameters:
    """Parameters for the Cosserat Beam"""
    beam_physics: BeamPhysicsParameters = field(default_factory=BeamPhysicsParameters)
    simulation: SimulationParameters = field(default_factory=SimulationParameters)
    contact: ContactParameters = field(default_factory=ContactParameters)
    beam_geometry: BeamGeometryParameters = field(default_factory=BeamGeometryParameters)
    visual: VisualParameters = field(default_factory=VisualParameters)

def create_cosserat_object(params: CosseratBeamParameters = None) -> CosseratBeamParameters:
    """
    Create and return a Cosserat beam object with default or custom parameters.
    
    Args:
        params (CosseratBeamParameters, optional): Custom parameters for the Cosserat beam.
    
    Returns:
        CosseratBeamParameters: An object containing all parameters for the Cosserat beam.
    """
    return params or CosseratBeamParameters()

# Example usage
cosserat_beam = create_cosserat_object()