# -*- coding: utf-8 -*-
"""
    Cosserat class in SofaPython3.
"""

__authors__ = "Yinoussa"
__contact__ = "adagolodjo@protonmail.com"
__version__ = "1.0.0"
__copyright__ = "(c) 2021,Inria"
__date__ = "July, 20 2023"

from useful.header import addHeader
from cosserat.cosseratObject import Cosserat
from dataclasses import dataclass
from useful.params import Parameters, BeamGeometryParameters


# [@info] ================ Unit: N, m, Kg, Pa  ================

# @todo use this dataclass to create the cosserat object
@dataclass
class CableParameters:
    youngModulus: float = 1.205e11
    poissonRatio: float = 0.499
    radius: float = 0.004
    density: float = 7.850e3
    maxDisplacement: float = 0.3

    GI: float =1.5708
    GA: float =3.1416e4
    EI: float =0.7854
    EA: float =3.1416e4
    L:  float = 1.  # beam length in m
    Rb: float = 0.02/2. # beam radius in m

@dataclass
class CableGeometryParameters:
    init_pos: list[float]  # = [0., 0., 0.]
    tot_length: float = 1.
    nbSectionS: int = 5
    nbFramesF: int = 30
    buildCollisionModel: int = 0
    beamMass: float = 1.





YM = 1.0e8
PR = 0.38
rayleighStiffness = 1.e-3  # Nope
firstOrder = 1
EI = 1.e2
coeff = 1



# beam parameters
length = 1.  # in m
Rb = 0.02/2.  # beam radius in m

geoParams = CableGeometryParameters(init_pos=[0., 0., 0.], tot_length=length,
                                    nbSectionS=5, nbFramesF=30, buildCollisionModel=0, beamMass=1.)
# inertialParams = CableParameters(youngModulus=YM, poissonRatio=PR, radius=Rb, density=7.850e3, maxDisplacement=0.3, nbSections=20)
nonLinearConfig = {'init_pos': [0., 0., 0.], 'tot_length': length, 'nbSectionS': 5,
                   'nbFramesF': 30, 'buildCollisionModel': 0, 'beamMass': 1.}
inertialParams = {'GI': 1.5708, 'GA': 3.1416e4,
                  'EI': 0.7854, 'EA': 3.1416e4, 'L': length,  'Rb': Rb}

isCollisionModel = 0

Params = Parameters(beamGeoParams=BeamGeometryParameters(init_pos=[0, 0, 0], beamLength=length,
                                    nbSection=5, nbFrames=30, buildCollisionModel=0))

def createScene(rootNode):
    addHeader(rootNode)

    node = rootNode.addChild('rootNode')
    node.gravity = [0., -9.81, 0.]
    node.addChild(
        Cosserat(parent=node, cosseratGeometry=nonLinearConfig, inertialParams=inertialParams, radius=Rb,
                 useCollisionModel=isCollisionModel, name="cosserat", youngModulus=YM, poissonRatio=PR,
                 rayleighStiffness=rayleighStiffness))