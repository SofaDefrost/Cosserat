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

#Gauss Quadrature
#Suurce : https://en.wikipedia.org/wiki/Gaussian_quadrature
# C = [1.0/sqrt(3.0),0.57735]

curv_abs_input = [0, 25, 50, 75]
curv_abs_output = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75]

###############
## Rate of angular Deformation  (2 sections)
###############
pos1 = [0.0,0.0,0.0]; pos2 = [0.0,0.0,0.0]; pos3 = [0.0,0.0,0.0] ; pos = [pos1, pos2, pos3]
distance1 = [0.0,0.5,0.0]; distance2 = [0.0,0.1,0.0]; distance3 = [0.0,0.5,0.0]; _distance = [distance1, distance2, distance3]
ddistance1 = [0.0,0.0,0.0]; ddistance2 = [0.0,0.0,0.0]; ddistance3 = [0.0,0.0,0.0]; _ddistance = [ddistance1, ddistance2, ddistance3]
R_b = 1.0
L = 75.0

_tension = 0.0

class DataComputationClass(CosseratActuation):
    """docstring for CosseratActuation.DataComputationClass"""

    def __init__(self):
        # print("The python init================================================++++++> ")
        super(DataComputationClass,Sofa.PythonScriptController.CosseratActuation).__init__()
        # self.K = [10, 1, 6]

    def initGraph(self, node):
        self.tension = 500
        self.node = node;
        self.BeamHookeLawForce = self.node.getObject('BeamHookeLawForce')
        self.rateAngularDeformMO = self.node.getObject('rateAngularDeformMO')
        self.K = self.rateAngularDeformMO.findData('position').value
        self.curv_abs_input = curv_abs_input
        self.X  = self.computeX()
        self.distance = [] # distance
        self.d_distance = [] # derivative of the distance
        # print("================================================++++++> ")

        # #CONSTANT parameters ( dy, dz, _dy, _dz)
        # self.computeDX(R_b/2.0, 0.0, 0.0, 0.0)

        # OBLIQUE parameters ( dy, dz, _dy, _dz)
        # self.computeDX(R_b/2.0, 0.0, -R_b/L, 0.0)

        # Act1 parameters ( dy, dz, _dy, _dz)
        # self.computeDX(0.0, R_b/2.0, R_b/(2*L), (-R_b)/(2*L)) # Act1 parameters ( dy, dz, _dy, _dz)

        # Act2 parameters ( dy, dz, _dy, _dz)
        # self.computeDX(0.0, R_b/2.0, -R_b/(2.0*L),  R_b/(2.0*L)) # Act2 parameters ( dy, dz, _dy, _dz)

        # Act3 parameters ( dy, dz, _dy, _dz)
        # self.computeDX(0.0, R_b/2.0, R_b/(2*L), (-R_b)/(2*L)) # Act3 parameters ( dy, dz, _dy, _dz)

        # CONSTANT parameters ( dy, dz, _dy, _dz)
        # self.vec_dy  = [R_b/2.0]; self.vec_dz  = [0.0]
        # self.vec_ddy = [0.0];     self.vec_ddz = [0.0]

        # self.vec_dy  = [R_b/2.0]; self.vec_dz  = [0.0]
        # self.vec_ddy = [-R_b/L];     self.vec_ddz = [0.0]

        self.vec_dy  = [0.0, 0.0, -R_b/2.0];            self.vec_dz  = [R_b/2.0, -R_b/2.0, 0.0]
        self.vec_ddy = [R_b/(2.0*L), R_b/(2.0*L), 0.0];     self.vec_ddz = [-R_b/(2.0*L), R_b/(2.0*L), 0.0]

        kwargs = [self.vec_dy, self.vec_dz, self.vec_ddy, self.vec_ddz, self.K]
        self.muti_ActuationIntegral(*kwargs)


        ############################## HELICAL PARAMETERS ####################################""
        # d = R_b/2.0
        # alpha = 1.529
        # p = 2.0 * pi * d * tan(alpha)
        # self.computeHelicalParameters(d, p)
        # print ("############################## HELICAL PARAMETERS ####################################")
        # print ("=======+++++++> self.distance : ",self.distance)
        # print ("=======+++++++> self.distance : ",self.d_distance)
        #
        self.BeamHookeLawForce.findData('distance0').value = self.distance[0]
        self.BeamHookeLawForce.findData('distance1').value = self.distance[1]
        self.BeamHookeLawForce.findData('ddistance0').value = self.d_distance[0]
        self.BeamHookeLawForce.findData('ddistance1').value = self.d_distance[1]


    def onBeginAnimationStep(self, dt):
        self.tension = self.tension + 8000.0 * dt;
        self.K = self.rateAngularDeformMO.findData('position').value

        integral = self.muti_ActuationIntegral(self.vec_dy, self.vec_dz, self.vec_ddy, self.vec_ddz, self.K)
        # print ("=++++++++=======+++> 0) muti_ActuationIntegral : ", integral )
        listIntegral = []
        for i in range(0,len(integral)):
            listIntegral.append(integral[i])

        # print("+++++++++++++++++++>>>> listIntegral: ",listIntegral)
        self.BeamHookeLawForce.findData('integral').value = listIntegral
        # print ("=++++++++=======+++> 1) muti_ActuationIntegral : ", integral )

        if(self.tension < 20501.0):
            self.BeamHookeLawForce.findData('tension').value = self.tension
            # print("===============>  Tension : ", self.tension)

        # self.BeamHookeLawForce.findData('tension').value = self.tension
        # print("===============>  Tension : ", self.tension)

def createScene(rootNode):

                rootNode.createObject('RequiredPlugin', pluginName='SoftRobots SofaPython SofaSparseSolver CosseratPlugin BeamAdapter', printLog='0')
                rootNode.createObject('VisualStyle', displayFlags='hideVisualModels showBehaviorModels showCollisionModels hideBoundingCollisionModels showForceFields hideInteractionForceFields hideWireframe')
                rootNode.createObject('FreeMotionAnimationLoop')
                rootNode.createObject('GenericConstraintSolver', tolerance="1e-10", printLog='0')

                rootNode.gravity = "0 0 0"
                rootNode.dt="0.01"
                rootNode.createObject('EulerImplicitSolver', firstOrder="0", rayleighStiffness="1.0", rayleighMass='0.10')
                rootNode.createObject('SparseLDLSolver', name='solver')
                rootNode.createObject('GenericConstraintCorrection')

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
                youngModulus='5.93e4',distance0=_distance, distance1=_distance, ddistance=_ddistance, tension=_tension)
                rateAngularDeformNode.createObject('PythonScriptController', classname="DataComputationClass")

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

                #### CylinderGridTop
                CylinderCollision = mappedFrameNode.createChild('CylinderCollision')
                # CylinderCollision.createObject('MeshSTLLoader', filename=path+'trunk.stl', name='loader', rotation='0 90 0', scale='0.155')
                CylinderCollision.createObject('CylinderGridTopology', name="loader", nx="8", ny="8", nz="20", length="75", radius="1", axis="1 0 0" )
                CylinderCollision.createObject('Mesh', src='@loader')
                CylinderCollision.createObject('MechanicalObject', template='Vec3d')
                CylinderCollision.createObject('Triangle')
                CylinderCollision.createObject('SkinningMapping', nbRef='2')
