# -*- coding: utf-8 -*-

import Sofa
import SofaPython
from math import sin,cos, sqrt, pi
import os
path = os.path.dirname(os.path.abspath(__file__))+'/mesh/'

_tension = 0.0
class TensionComputing(Sofa.PythonScriptController):
    def initGraph(self, node):
        self.tension = 500
        self.node = node;
        self.BeamHookeLawForce = self.node.getObject('BeamHookeLawForce')

    def onBeginAnimationStep(self, dt):
        self.tension = self.tension + 8000.0 * dt;

        if(self.tension < 14501.0):
            self.BeamHookeLawForce.findData('tension').value = self.tension
            print("===============>  Tension : ", self.tension)
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


                ###############
                ## Rate of angular Deformation  (2 sections)
                ###############
                pos1 = [0.0,0.0,0.0]; pos2 = [0.0,0.0,0.0]; pos3 = [0.0,0.0,0.0] ; pos = [pos1, pos2, pos3]
                distance1 = [0.0,0.5,0.0]; distance2 = [0.0,0.5,0.0]; distance3 = [0.0,0.5,0.0]; _distance = [distance1, distance2, distance3]
                ddistance1 = [0.0,0.0,0.0]; ddistance2 = [0.0,0.0,0.0]; ddistance3 = [0.0,0.0,0.0]; _ddistance = [ddistance1, ddistance2, ddistance3]

                rateAngularDeformNode = rootNode.createChild('rateAngularDeform')
                rateAngularDeformMO = rateAngularDeformNode.createObject('MechanicalObject', template='Vec3d', name='rateAngularDeformMO', position=pos, velocity='0 0 0 0 0 0 0 0 0') # (2 series of 3 angles for 2 sections. we suppose that the lenght is 10 for each)
                # BeamHookeLawForce = rateAngularDeformNode.createObject('BeamHookeLawForceField', name="BeamHookeLawForce", crossSectionShape='circular', length='25 25 25', radius='1.0', youngModulus='5e6')
                BeamHookeLawForce = rateAngularDeformNode.createObject('CosseratInternalActuation', name="BeamHookeLawForce",  crossSectionShape='circular', length='25 25 25', radius='1.0',
                youngModulus='5.93e4',distance=_distance, ddistance=_ddistance, tension=_tension)
                rateAngularDeformNode.createObject('PythonScriptController', classname="TensionComputing")

                ##############
                ## Frames
                ##############
                # the node of the frame needs to inherit from rigidBaseMO and rateAngularDeform
                frames = ["0.0 0 0  0 0 0 1  5 0 0 0 0 0 1   10 0 0  0 0 0 1   15 0 0 0 0 0 1  20 0 0  0 0 0 1  25 0 0  0 0 0 1  30 0 0  0 0 0 1  35 0 0  0 0 0 1 " +
                            " 40 0 0  0 0 0 1   45 0 0  0 0 0 1  50 0 0  0 0 0 1  55 0 0  0 0 0 1   60 0 0  0 0 0 1    65 0 0  0 0 0 1  70 0 0  0 0 0 1  75 0 0  0 0 0 1" ]

                mappedFrameNode = rigidBaseNode.createChild('MappedFrames')
                rateAngularDeformNode.addChild(mappedFrameNode)
                framesMO = mappedFrameNode.createObject('MechanicalObject', template='Rigid3d', name="FramesMO", position=frames, showObject='1', showObjectScale='1' )

                inputMO = rateAngularDeformMO.getLinkPath()
                inputMO_rigid = RigidBaseMO.getLinkPath()
                outputMO = framesMO.getLinkPath()
                mappedFrameNode.createObject('DiscretCosseratMapping', curv_abs_input='0 25 50 75', curv_abs_output='0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75', input1=inputMO, input2=inputMO_rigid,output=outputMO, debug='0' )

                #### CylinderGridTop
                CylinderCollision = mappedFrameNode.createChild('CylinderCollision')
                # CylinderCollision.createObject('MeshSTLLoader', filename=path+'trunk.stl', name='loader', rotation='0 90 0', scale='0.155')
                CylinderCollision.createObject('CylinderGridTopology', name="loader", nx="8", ny="8", nz="20", length="75", radius="1", axis="1 0 0" )
                CylinderCollision.createObject('Mesh', src='@loader')
                CylinderCollision.createObject('MechanicalObject', template='Vec3d')
                CylinderCollision.createObject('Triangle')
                CylinderCollision.createObject('SkinningMapping', nbRef='2')
