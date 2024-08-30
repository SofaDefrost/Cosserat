# @todo use this dataclass to create the cosserat object

from dataclasses import dataclass, field
from typing import List


# @TODO: Improve this function, remove hard coding.
@dataclass
class BeamGeometryParameters:
    """Cosserat Beam Geometry parameters"""

    beamLength: float = 1.0  # beam length in m
    nbSection: int = 5  # number of sections, sections along the beam length
    nbFrames: int = 30  # number of frames along the beam
    buildCollisionModel: int = 0
    init_pos: List[float] = field(
        default_factory=lambda: [0.0, 0.0, 0.0]
    )  # The beam rigid base position as a list [x, y, z]

    """Parameters for the visualisation of the beam"""
    show_frames_object: bool = True
    show_frames_indices: bool = False
    showRigidBaseObject: int = 1


#Todo: two objects from the same base class to define different instead of useInertia
@dataclass
class BeamPhysicsParameters:
    """Cosserat Beam Physics parameters"""

    """First set of parameters"""
    youngModulus: float = 1.205e8
    poissonRatio: float = 0.3

    """Second set of parameters"""
    useInertia: bool = False
    GI: float = 1.5708
    GA: float = 3.1416e4
    EI: float = 0.7854
    EA: float = 3.1416e4

    """Common parameters used in both cases"""
    beamMass: float = 1.0
    beamRadius: float = 0.02 / 2.0  # beam radius in m
    beamLength: float = 1.0  # beam length in m along the X axis
    # Todo : add complex beam shape
    beamShape: str = "circular"  # beam shape, circular or rectangular
    """"If using rectangular beam shape"""
    length_Y: float = 0.1  # length of the beam in the Y direction
    length_Z: float = 0.1  # length of the beam in the Z direction


@dataclass
class SimulationParameters:
    """Simulation parameters"""

    rayleighStiffness: float = 0.2
    rayleighMass: float = 0.1
    firstOrder: bool = False


@dataclass
class VisualParameters:
    """Visual parameters"""

    showObject: int = 1
    showObjectScale: float = 1.0
    showObjectColor: List[float] = field(default_factory=lambda: [1.0, 0.0, 0.0, 1.0])


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

    beamPhysicsParams: BeamPhysicsParameters = field(
        default_factory=BeamPhysicsParameters
    )
    simuParams: SimulationParameters = field(
        default_factory=SimulationParameters
    )  # SimulationParameters()
    contactParams: ContactParameters = field(
        default_factory=ContactParameters
    )  # ContactParameters()
    beamGeoParams: BeamGeometryParameters = field(
        default_factory=BeamGeometryParameters
    )
    visualParams: VisualParameters = field(default_factory=VisualParameters)
