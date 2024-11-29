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
import os

path = f'{os.path.dirname(os.path.abspath(__file__))}/mesh/'

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
                                    nbSection=5, nbFrames=26, buildCollisionModel=0)
physicsParams = BeamPhysicsParameters(beamMass=5., youngModulus= 1.0e8, poissonRatio=0.38, beamRadius=6.2e-1/2.,
                                      beamLength=65.5)
simuParams = SimulationParameters(rayleighStiffness=1.e-3)
Params = Parameters(beamGeoParams=geoParams, beamPhysicsParams=physicsParams, simuParams=simuParams)


def createIntermediateNode(parent, rigidCentral=None, baseName="rigidState",rigidIndex=2):
    """ The intermediateRigid construction """
    if rigidCentral is None:
        rigidCentral = [0, 0., 0, 0, 0, 0, 1]
    q0 = [0, 0, 0, 1]
    q1 = [0, 0, 0, 1]
    q2 = [0, 0, 0, 1]


    interChildPos = [[0, 1.7, 1., q0[0], q0[1], q0[2], q0[3]],
                     [0, -1.7, 1., q1[0], q1[1], q1[2], q1[3]],
                     [0., 0., -2., q2[0], q2[1], q2[2], q2[3]]]
    intermediateRigid = parent.addChild('intermediateRigid')
    intermediateRigid.addObject('MechanicalObject', name=baseName, template='Rigid3d',
                                                      position=rigidCentral, showObject=True, showObjectScale=0.4)
    # @TODO add the mass of the disk
    intermediateRigid.addObject('RigidMapping', name="interRigidMap", index=rigidIndex)

    """Create intermediate kids"""
    interRigidChild = intermediateRigid.addChild('interRigidChild')
    interRigidChild.addObject('MechanicalObject', name="interRigidChildMo", template='Rigid3d',
                              position=interChildPos, showObject=True, showObjectScale=0.4)
    interRigidChild.addObject('RigidMapping', name="interRigidMap")

    """ Add disk to the scene, actually this disk is only use for the visualisation """
    loadDisk(intermediateRigid)

    return intermediateRigid

def loadDisk(parentNode):
    #### MAPPING of the disks (here some torus) ####    ###########
    diskMapping = parentNode.addChild("diskMapping")
    diskMapping.addObject("MeshOBJLoader", name="diskLoader", filename=f'{path}disqueBlender2.obj')
    diskMapping.addObject("OglModel", name="Visual", src="@diskLoader", color="1.0 0.5 0.25 1.0")
    diskMapping.addObject('RigidMapping', name="diskMap")


def createRigidDisk(parentNode): # @todo add a parameter for the number of disk
    """ Create the rigid disk nodes """
    for i in range(1, 14):
        """ Create the rigid disk nodes """
        createIntermediateNode(parentNode, rigidCentral=[0, 0., 0, 0, 0, 0, 1],
                               baseName=f"rigidState{i}", rigidIndex=i*2)


def createScene(rootNode):
    addHeader(rootNode)
    rootNode.gravity = [0., -9.81, 0.]

    stemNode = addSolverNode(rootNode, name="stemNode")
    cosserat = stemNode.addChild(CosseratBase(parent=stemNode, params=Params))

    # create the rigid disk node
    createRigidDisk(cosserat.cosseratFrame)
    return






