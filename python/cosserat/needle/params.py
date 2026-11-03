"""Scenario parameters for the needle-insertion example.

Ported from the legacy ``examples/python3/cosserat/needle/params.py`` into the
maintained ``cosserat`` package. Units are cm, kg, s.
"""

from dataclasses import dataclass, field


@dataclass
class GeometryParams:
    radius: float = 0.1
    nbSections: int = 16
    nbFrames: int = 15
    totalLength: float = 15.0


@dataclass
class PhysicsParams:
    youngModulus: float = 1.20e9
    poissonRatio: float = 0.4
    mass: float = 0.3
    rayleighStiffness: float = 0.1


@dataclass
class FemParams:
    youngModulus: float = 100
    poissonRatio: float = 0.48
    mass: float = 0.3
    rayleigh: float = 0.1
    minVol: str = "16. -8. -5."
    maxVol: str = "40. 8. 5."
    mesh: str = "6 6 6"
    box: str = "15 -10 -10 41 -6 10"


@dataclass
class ContactParams:
    contactDistance: float = 0.01
    alarmDistance: float = 0.1
    dataMu: str = "mu=0.1"
    angleCone: float = 0.01
    coneFactor: float = 0.0


@dataclass
class NeedleParameters:
    Geometry: GeometryParams = field(default_factory=GeometryParams)
    Physics: PhysicsParams = field(default_factory=PhysicsParams)
    Fem: FemParams = field(default_factory=FemParams)
    contact: ContactParams = field(default_factory=ContactParams)


@dataclass
class ConstraintsParams:
    constraintDistance: float = 1.3  # distance between two constraint points
    entryForce: float = 0.3  # the required force to penetrate the volume
