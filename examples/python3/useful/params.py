# @todo use this dataclass to create the cosserat object

from dataclasses import dataclass
@dataclass
class BeamPhysicsParameters:
    """one can only use one of the following two parameters"""
    """First set of parameters"""
    youngModulus: float = 1.205e11
    poissonRatio: float = 0.499

    """Second set of parameters"""
    useInertia: bool = False
    GI: float = 1.5708
    GA: float = 3.1416e4
    EI: float = 0.7854
    EA: float = 3.1416e4

    """Common parameters, this parameters are used in both cases"""
    beamMass: float = 1.
    beamLength: float = 1.  # beam length in m
    beamRadius: float = 0.02 / 2.  # beam radius in m

@dataclass
class BeamGeometryParameters:
    """cosserat beam Geometry parameters"""
    init_pos: list[float]  # = [0., 0., 0.], The beam rigid base position
    beamLength: float = 1. # beam length in m
    nbSection: int = 5 # number of sections, here sections are not cross-sections but sections along the beam length
    nbFrames: int = 30
    buildCollisionModel: int = 0

@dataclass
class SimulationParameters:
    """Simulation parameters"""
    rayleighStiffness: float = 0.2
    rayleighMass: float = 0.1

@dataclass
class ContactParameters:
    """Contact parameters"""
    responseParams: str = 'mu=0.8'
    response: str = 'FrictionContactConstraint'
    alarmDistance: float = 0.05
    contactDistance: float = 0.01
    isMultithreading: bool = False
    tolerance: float = 1.e-8
    maxIterations: int = 100
    epsilon: float = 1.e-6