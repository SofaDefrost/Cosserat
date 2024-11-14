# -*- coding: utf-8 -*-
"""
    Scene robot ISIR.
"""
__authors__ = "Yinoussa"
__contact__ = "adagolodjo@protonmail.com"
__version__ = "1.0.0"
__copyright__ = "(c) 2021,Inria"
__date__ = "July, 20 2023"

from useful.header import addHeader, addSolverNode
from useful.params import BeamPhysicsParameters, BeamGeometryParameters, SimulationParameters
from useful.params import Parameters
from cosserat.CosseratBase import CosseratBase
from cable import PullingCable
from utils import createRigidDisk, create_cable_points, FingerController
from splib3.loaders import loadPointListFromFile

import os

"""  @info
Comme tu l'aura compris notre prototype est composé d'une colonne centrale (fils de fer tressés avec une gaine)
 et de 13 disques en iglidur.
 
 -La colonne mesure 65.5 cm pour un diamètre de 6.2mm

# -Les disques ont un diamètres de 5.2cm, une épaisseur de 2mm et une masse de 6.35g
#
# -Les disques sur la colonne sont espacés d'environ 4.7cm
#
# -Concernant les trous des disques, il y en a 3 sortes :
#
# les trous intérieurs sont au nombre de 3, ont un diamètre de 2mm et sont placés à environ 6mm du centre
# les trous intermédiaires sont au nombre de 3, ont un diamètre de 2.5mm et sont placés à environ 2cm du centre
# les trous extérieurs sont au nombre de 12, ont un diamètre de 1mm et sont placés à environ 2.3cm du centre
# -Enfin concernant le routage des câbles :
#
# Le manipulateur est constitué de 3 sections, chaque section est composé de 3 fils (disposés à 120° les uns des autres)
# Les câbles sont routés de manière parallèle
# La première section (celle du haut) ainsi que la deuxième (intermédiaire) sont constituées de 4 disques et la
# troisième (celle du bas) de 5 disques
# Les câbles ne sont déployés que sur leur propre section. Par exemple, les câbles qui fonctionnent sur la section 3
# passent par les trous extérieurs sur cette section et par les trous intérieurs sur les sections 1 et 2.
# [@info] ================ Unit: N, cm, g, Pa  ================
"""

# @todo
geoParams = BeamGeometryParameters(init_pos=[0., 0., 0.], beamLength=65.5, show_frames_object=1,
                                   nbSection=14, nbFrames=14, buildCollisionModel=0)
physicsParams = BeamPhysicsParameters(beamMass=5., youngModulus=1.0e8, poissonRatio=0.38, beamRadius=6.2e-1 / 2.,
                                      beamLength=65.5)
simuParams = SimulationParameters()
Params = Parameters(beamGeoParams=geoParams, beamPhysicsParams=physicsParams, simuParams=simuParams)

rotation = [0.0, 0.0, 0.0]
translation = [0.0, 0.0, 0.0]
pullPointLocation = [0.0, 0.0, 0.0]
youngModulus = 18000
valueType = 'position'


def createScene(rootNode):
    addHeader(rootNode, isConstrained=True)
    rootNode.gravity = [9.81, 0., 0.]

    stemNode = addSolverNode(rootNode, name="stemNode", isConstrained=True)

    cosserat = stemNode.addChild(CosseratBase(parent=stemNode, params=Params))
    createRigidDisk(cosserat.cosseratFrame)

    # create the rigid disk node
    cableBaseState, cable3DState = create_cable_points(geoParams=geoParams, num_segments=14)

    cable1Node = cosserat.cosseratFrame.addChild("cable1Node")
    cable1Node.addObject('MechanicalObject', name="cable1Pos", template='Rigid3d', position=cableBaseState,
                         showObjectScale=0.8, showObject=True)
    cable1Node.addObject('RigidMapping', name="cable1Map", rigidIndexPerPoint=[i for i in range(15)],
                         globalToLocalCoords=1
                         )
    cable1_mechaNode = cable1Node.addChild("cable_mecha1Node")
    cable_Mo = cable1_mechaNode.addObject('MechanicalObject', name="cable1_mechaPos", template='Vec3d',
                                          position=cable3DState,
                                          showIndices=True)
    # cable1_mechaNode.addObject('RigidMapping', name="cable1Map", globalToLocalCoords=1)
    cable1_mechaNode.addObject('SkinningMapping', nbRef='1')

    cable = PullingCable(cable1_mechaNode, "PullingCable",
                         pullPointLocation=pullPointLocation,
                         rotation=rotation,
                         translation=translation,
                         cableGeometry=cable3DState[0],
                         valueType=valueType,
                         input=cable_Mo.getLinkPath()
                         )
    cable1_mechaNode.addObject(FingerController(cable))

    return rootNode
