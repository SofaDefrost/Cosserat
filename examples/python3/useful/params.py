# @todo use this dataclass to create the cosserat object

from dataclasses import dataclass
@dataclass
class CableParameters:
    youngModulus: float = 1.205e11
    poissonRatio: float = 0.499
    radius: float = 0.004
    density: float = 7.850e3
    maxDisplacement: float = 0.3

    GI: float = 1.5708
    GA: float = 3.1416e4
    EI: float = 0.7854
    EA: float = 3.1416e4
    L: float = 1.  # beam length in m
    Rb: float = 0.02 / 2.  # beam radius in m


@dataclass
class CableGeometryParameters:
    init_pos: list[float]  # = [0., 0., 0.]
    tot_length: float = 1.
    nbSectionS: int = 5
    nbFramesF: int = 30
    buildCollisionModel: int = 0
    beamMass: float = 1.