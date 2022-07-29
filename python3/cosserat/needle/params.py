
from dataclasses import dataclass
import string

# Units are cm, kg, s
@dataclass
class Geometry:
    radius: float = 0.1
    nbSections: int = 16
    nbFrames: int = 15
    totalLength: float = 15.


@dataclass
class Physics:
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
    minVol: string = "16. -8. -5."
    maxVol: string = "40. 8. 5."
    mesh: string = "6 6 6"
    box: string = "15 -10 -10 41 -6 10"

@dataclass
class ContactParams:
    contactDistance: float = 0.01
    alarmDistance: float = 0.1
    dataMu: string = "mu=0.01"
    angleCone: float = "0.01"
    coneFactor: float = "0"

@dataclass
class NeedleParameters:
    Geometry: Geometry = Geometry()
    Physics: Physics = Physics()
    FemParams: FemParams = FemParams()
    contact: ContactParams = ContactParams()
