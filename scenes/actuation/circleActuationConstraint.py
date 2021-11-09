# -*- coding: utf-8 -*-

from CosseratActuation import *
import Sofa
import SofaPython
from math import sin,cos, sqrt, pi
from stlib.physics.collision import CollisionMesh
import os
from splib.numerics import sin, cos, to_radians
from stlib.physics.deformable import ElasticMaterialObject
from splib.objectmodel import SofaPrefab, SofaObject
from stlib.physics.mixedmaterial import Rigidify
from stlib.components import addOrientedBoxRoi
from splib.numerics import vec3
from splib.numerics.quat import Quat
path = os.path.dirname(os.path.abspath(__file__))+'/mesh/'

# Gauss Quadrature
# Suurce : https://en.wikipedia.org/wiki/Gaussian_quadrature
# C = [1.0/sqrt(3.0),0.57735]

curv_abs_input = [0, 25, 50, 75]
indices = [0, 1, 2]
curv_abs_output = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75]


###############
# #Rate of angular Deformation  (2 sections)
###############
pos1 = [0.0, 0.0, 0.0]; pos2 = [0.0,0.0,0.0]; pos3 = [0.0,0.0,0.0] ; pos = [pos1, pos2, pos3]
distance1 = [0.0, 0.5, 0.0]; distance2 = [0.0,0.1,0.0]; distance3 = [0.0,0.5,0.0]; _distance = [distance1, distance2, distance3]
ddistance1 = [0.0, 0.0, 0.0]; ddistance2 = [0.0,0.0,0.0]; ddistance3 = [0.0,0.0,0.0]; _ddistance = [ddistance1, ddistance2, ddistance3]
R_b = 1.0
L = 75.0

_tension = 0.0

class DataComputationClass(CosseratActuation):
    """docstring for CosseratActuation.DataComputationClass"""

    def __init__(self, nodeA):
        # print("========================= DataComputationClass ======================= ")
        CosseratActuation.__init__(self)
        self.node = nodeA
        # self.K = [10, 1, 6]

    def initGraph(self, node):
        self.tension = 500
        self.node = node
        self.BeamHookeLawForce = self.node.getObject('BeamHookeLawForce')
        self.rateAngularDeformMO = self.node.getObject('rateAngularDeformMO')
        # self.cableConstraint = self.node.getObject('cableConstraint')
        self.K = self.rateAngularDeformMO.findData('position').value
        self.curv_abs_input = curv_abs_input
        self.X  = self.computeX()

        listIntegral = []
        for i in range(0, len(indices)):
            listIntegral.append([0.0,0.0,0.0])
        self.BeamHookeLawForce.findData('integral').value = listIntegral
        # self.cableConstraint.findData('integral').value = listIntegral

        # self.distance = [] # distance
        # self.d_distance = [] # derivative of the distance

        # self.vecDistance1 = [] # distance at s1 = L_{i-1} + C1(L_{i} - L_{i-1})
        # self.vecDistance2 = [] # distance at s2 = L_{i-1} + C2(L_{i} - L_{i-1})
        # self.vecDDistance1 = [] # derivative of the distance at s1
        # self.vecDDistance2 = [] # derivative of the distance at s2
        # print("================================================++++++> ")
        # CONSTANT parameters ( dy, dz, _dy, _dz)
        self.vec_dy  = [R_b/2.0]; self.vec_dz  = [0.0]
        self.vec_ddy = [0.0];     self.vec_ddz = [0.0]

        #self.vec_dy  = [R_b/2.0]; self.vec_dz  = [0.0]
        #self.vec_ddy = [-R_b/L];     self.vec_ddz = [0.0]

        # self.vec_dy  = [0.0, 0.0, -R_b/2.0];            self.vec_dz  = [R_b/2.0, -R_b/2.0, 0.0]
        # self.vec_ddy = [R_b/(2.0*L), R_b/(2.0*L), 0.0];     self.vec_ddz = [-R_b/(2.0*L), R_b/(2.0*L), 0.0]

        kwarg = [self.vec_dy, self.vec_dz, self.vec_ddy, self.vec_ddz]
        self.computeMultiDistanceVectors(*kwarg)

        # kwargs = [self.vec_dy, self.vec_dz, self.vec_ddy, self.vec_ddz, self.K]
        # self.muti_ActuationIntegral(*kwargs)


        ############################## HELICAL PARAMETERS ####################################""
        # self.BeamHookeLawForce.findData('distance0').value = self.distance[0]
        # self.BeamHookeLawForce.findData('distance1').value = self.distance[1]
        # self.BeamHookeLawForce.findData('ddistance0').value = self.d_distance[0]
        # self.BeamHookeLawForce.findData('ddistance1').value = self.d_distance[1]

    def onBeginAnimationStep(self, dt):
        self.tension = self.tension + 50.0
        self.K = self.rateAngularDeformMO.findData('position').value

        integral = self.muti_ActuationIntegral(self.vec_dy, self.vec_dz, self.vec_ddy, self.vec_ddz, self.K)
        # print ("=++++++++=======+++> 0) muti_ActuationIntegral : ", integral )
        listIntegral = []
        for i in range(0,len(integral)):
            listIntegral.append(integral[i])

        self.BeamHookeLawForce.findData('integral').value = listIntegral
        # self.cableConstraint.findData('integral').value = listIntegral
        # print ("======+++> listIntegral : ", listIntegral)

        if(self.tension < 800000.0):
            # print ("===> Tension : ", self.tension)
            self.BeamHookeLawForce.findData('tension').value = self.tension

def createScene(rootNode):

    rootNode.createObject('RequiredPlugin', pluginName='SoftRobots SofaPython SofaSparseSolver CosseratPlugin BeamAdapter', printLog='0')
    rootNode.createObject('VisualStyle', displayFlags='hideVisualModels showBehaviorModels showCollisionModels hideBoundingCollisionModels showForceFields hideInteractionForceFields hideWireframe')
    rootNode.createObject('FreeMotionAnimationLoop')
    rootNode.createObject('GenericConstraintSolver', tolerance="1e-20", maxIterations=1000, printLog=0)

    rootNode.gravity = "0 0 0"
    rootNode.dt = "0.01"
    rootNode.createObject('EulerImplicitSolver', firstOrder="0", rayleighStiffness="1.0", rayleighMass='0.10')
    rootNode.createObject('SparseLDLSolver', name='solver')
    # rootNode.createObject('GenericConstraintCorrection')

    ###############hresho
    ## RigidBase
    ###############
    rigidBaseNode= rootNode.createChild('rigidBase')
    RigidBaseMO = rigidBaseNode.createObject('MechanicalObject', template='Rigid3d', name="RigidBaseMO", position="0 0 0  0 0 0. 1", showObject='1', showObjectScale='0.1', velocity='0 0 0.0 0.0 0 0' )
    rigidBaseNode.createObject('RestShapeSpringsForceField', name='spring', stiffness="5000000", angularStiffness="500000000", external_points="0", mstate="@RigidBaseMO", points="0", template="Rigid3d"  )

    rateAngularDeformNode = rootNode.createChild('rateAngularDeform')
    rateAngularDeformMO = rateAngularDeformNode.createObject('MechanicalObject', template='Vec3d', name='rateAngularDeformMO', position=pos, velocity='0 0 0 0 0 0 0 0 0') # (2 series of 3 angles for 2 sections. we suppose that the lenght is 10 for each)
    # BeamHookeLawForce = rateAngularDeformNode.createObject('BeamHookeLawForceField', name="BeamHookeLawForce", crossSectionShape='circular', length='25 25 25', radius='1.0', youngModulus='5e6')
    BeamHookeLawForce = rateAngularDeformNode.createObject('CosseratInternalActuation', name="BeamHookeLawForce",  crossSectionShape='circular', length='25 25 25', radius='1.0',
    youngModulus='3.8e4',distance0=_distance, distance1=_distance, tension=_tension)

    # cable_position=["25.0 0.0  0.0 " + " 50.0 0.0 0.0 " + "75   0.0  0.0" ]
    # cable = rateAngularDeformNode.createChild('cable')
    # cable.createObject('MechanicalObject', name="cablePos", position=cable_position, template="", showObject="1", showIndices="1")
    cableConstraint = rateAngularDeformNode.createObject('CosseratActuatorConstraint', name="cableConstraint", indices=indices, value="50", integral="@BeamHookeLawForce.integral")
    # cable.createObject('SkinningMapping', nbRef='1',  mapForces='false', mapMasses='false')

    DataComputationClass(rateAngularDeformNode)

    ##############
    ## Frames
    ##############
    # the node of the frame needs to inherit from rigidBaseMO and rateAngularDeform
    frames = ["0.0 0 0  0 0 0 1  5 0 0 0 0 0 1   10 0 0  0 0 0 1   15 0 0 0 0 0 1  20 0 0  0 0 0 1  25 0 0  0 0 0 1  30 0 0  0 0 0 1  35 0 0  0 0 0 1 " +
                " 40 0 0  0 0 0 1   45 0 0  0 0 0 1  50 0 0  0 0 0 1  55 0 0  0 0 0 1   60 0 0  0 0 0 1    65 0 0  0 0 0 1  70 0 0  0 0 0 1  75 0 0  0 0 0 1" ]

    mappedFrameNode = rigidBaseNode.createChild('MappedFrames')
    rateAngularDeformNode.addChild(mappedFrameNode)
    framesMO = mappedFrameNode.createObject('MechanicalObject', template='Rigid3d', name="FramesMO", position=frames, showObject='1', showObjectScale='1' )

    inputMO = rateAngularDeformMO.getLinkPath() # + " " + RigidBaseMO.getLinkPath()
    #inputMO = rateAngularDeformMO.getLinkPath()
    inputMO_rigid = RigidBaseMO.getLinkPath()
    outputMO = framesMO.getLinkPath()
    # TODO:
    mappedFrameNode.createObject('DiscretCosseratMapping', curv_abs_input=curv_abs_input, curv_abs_output=curv_abs_output, input1=inputMO, input2=inputMO_rigid,output=outputMO, debug='0' )

    ############################
    ## Cable inside the beam  ##
    ############################
    # cable_position=[
    #         "25.0 0.0  0.0 " + " 50.0 0.0 0.0 " + "75   0.0  0.0" ]
    # cable = mappedFrameNode.createChild('cable')
    # cable.createObject('MechanicalObject', name="cablePos", position=cable_position, template="", showObject="1", showIndices="1")
    # cableConstraint = cable.createObject('CosseratActuatorConstraint', name="cableConstraint", indices=indices, value="0.5", force="10", integral="@../../../rateAngularDeform/BeamHookeLawForce.integral")
    # cable.createObject('SkinningMapping', nbRef='1',  mapForces='false', mapMasses='false')

    #### CylinderGridTop
    CylinderCollision = mappedFrameNode.createChild('CylinderCollision')
    # CylinderCollision.createObject('MeshSTLLoader', filename=path+'trunk.stl', name='loader', rotation='0 90 0', scale='0.155')
    CylinderCollision.createObject('CylinderGridTopology', name="loader", nx="8", ny="8", nz="20", length="75", radius="1", axis="1 0 0" )
    CylinderCollision.createObject('Mesh', src='@loader')
    CylinderCollision.createObject('MechanicalObject', template='Vec3d')
    CylinderCollision.createObject('Triangle')
    CylinderCollision.createObject('SkinningMapping', nbRef='2')
    mappedFrameNode.createObject('GenericConstraintCorrection')
    # rigidBaseNode.createObject('LinearSolverConstraintCorrection')
