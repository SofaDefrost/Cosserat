# @todo use this dataclass to create the cosserat object

from dataclasses import dataclass, field
from typing import List, Literal


# @TODO: Improve this function, remove hard coding.
@dataclass
class BeamGeometryParameters:
    """Cosserat Beam Geometry parameters"""

    beam_length: float = 1.0  # beam length in m
    nb_section: int = 5  # number of sections, sections along the beam length
    nb_frames: int = 30  # number of frames along the beam
    build_collision_model: int = 0

    def validate(self):
        assert self.beam_length > 0, "Beam length must be positive"
        assert self.nb_section > 0, "Number of sections must be positive"
        assert self.nb_frames > 0, "Number of frames must be positive"
        assert self.nb_frames >= self.nb_section, "Number of frames must be positive"

@dataclass
class BeamPhysicsBaseParameters:
    """Base class for Cosserat Beam Physics parameters"""

    young_modulus: float = 1.205e8
    poisson_ratio: float = 0.3
    beam_mass: float = 1.0
    beam_radius: float = 0.01  # default radius of 0.02 / 2.0
    beam_length: float = 1.0  # default length along the X axis
    beam_shape: Literal['circular', 'rectangular'] = 'circular'
    length_Y: float = 0.1  # length in Y direction for rectangular beam
    length_Z: float = 0.1  # length in Z direction for rectangular beam
    useInertia: bool = False

    def validate(self):
        assert self.young_modulus > 0, "Young's modulus must be positive"
        assert self.poisson_ratio > 0 and self.poisson_ratio < 0.5, "Poisson's ratio must be between 0 and 0.5"
        assert self.beam_mass > 0, "Beam mass must be positive"
        assert self.beam_radius > 0, "Beam radius must be positive"
        assert self.beam_length > 0, "Beam length must be positive"


@dataclass
class BeamPhysicsParametersNoInertia(BeamPhysicsBaseParameters):
    """Parameters for a Cosserat Beam without inertia"""
    pass


@dataclass
class BeamPhysicsParametersWithInertia(BeamPhysicsBaseParameters):
    """Parameters for a Cosserat Beam with inertia"""

    GI: float = 1.5708
    GA: float = 3.1416e4
    EI: float = 0.7854
    EA: float = 3.1416e4

    def validate(self):
        super().validate()
        assert self.GI > 0, "GI must be positive"
        assert self.GA > 0, "GA must be positive"
        assert self.EI > 0, "EI must be positive"
        assert self.EA > 0, "EA must be positive"


@dataclass
class SimulationParameters:
    """Simulation parameters"""

    rayleigh_stiffness: float = 0.2
    rayleigh_mass: float = 0.1
    firstOrder: bool = False

    def validate(self):
        assert self.rayleigh_stiffness >= 0, "Rayleigh stiffness must be non-negative"
        assert self.rayleigh_mass >= 0, "Rayleigh mass must be non-negative"


@dataclass
class VisualParameters:
    """Visual parameters"""

    showObject: int = 1
    show_object_scale: float = 1.0
    show_object_color: List[float] = field(default_factory=lambda: [1.0, 0.0, 0.0, 1.0])

    def validate(self):
        assert len(self.show_object_color) == 4, "Color must have four components (RGBA)"
        assert all(0.0 <= x <= 1.0 for x in self.show_object_color), "Color components must be in range [0, 1]"


@dataclass
class ContactParameters:
    """Contact parameters"""

    responseParams: str = "mu=0.0"
    response: str = "FrictionContactConstraint"
    alarmDistance: float = 0.05
    contactDistance: float = 0.01
    isMultithreading: bool = False
    tolerance: float = 1.0e-8
    maxIterations: int = 100
    epsilon: float = 1.0e-6


@dataclass
class Parameters:
    """Parameters for the Cosserat Beam"""

    beam_physics_params: BeamPhysicsBaseParameters = field(default_factory=BeamPhysicsParametersNoInertia)
    simu_params: SimulationParameters = field(
        default_factory=SimulationParameters
    )  # SimulationParameters()
    contact_params: ContactParameters = field(
        default_factory=ContactParameters
    )  # ContactParameters()
    beam_geo_params: BeamGeometryParameters = field(
        default_factory=BeamGeometryParameters
    )
    visual_params: VisualParameters = field(default_factory=VisualParameters)

    def validate(self):
        self.beam_physics_params.validate()
        self.simu_params.validate()
        self.contact_params.validate()
        self.beam_geo_params.validate()
        self.visual_params.validate()