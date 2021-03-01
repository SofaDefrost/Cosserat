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
        self.tension = self.tension + 8000 * dt;
        self.BeamHookeLawForce.findData('tension').value = self.tension


def createScene(rootNode):

                rootNode.createObject('RequiredPlugin', pluginName='SoftRobots SofaPython SofaSparseSolver CosseratPlugin BeamAdapter', printLog='0')
                rootNode.createObject('VisualStyle', displayFlags='hideVisualModels showBehaviorModels showCollisionModels hideBoundingCollisionModels showForceFields hideInteractionForceFields hideWireframe')
                rootNode.createObject('FreeMotionAnimationLoop')
                rootNode.createObject('GenericConstraintSolver', tolerance="1e-10", printLog='0')

                rootNode.gravity = "0 0 0"
                rootNode.dt="0.01"
                rootNode.createObject('EulerImplicitSolver', firstOrder="0", rayleighStiffness="1.0", rayleighMass='0.10')
                #rootNode.createObject('StaticSolver')
                #rootNode.createObject('CGLinearSolver', name='solver', tolerance='1e-20',  threshold='1e-20', verbose='1')
                rootNode.createObject('SparseLDLSolver', name='solver')

                rootNode.createObject('GenericConstraintCorrection')
                ###############hresho
                ## RigidBase
                ###############
                rigidBaseNode= rootNode.createChild('rigidBase')
                #rigidBaseNode.createObject('GenericConstraintCorrection')

                RigidBaseMO = rigidBaseNode.createObject('MechanicalObject', template='Rigid3d', name="RigidBaseMO", position="0 0 0  0 0 0. 1", showObject='1', showObjectScale='0.1', velocity='0 0 0.0 0.0 0 0' )
                rigidBaseNode.createObject('RestShapeSpringsForceField', name='spring', stiffness="5000000", angularStiffness="500000000", external_points="0", mstate="@RigidBaseMO", points="0", template="Rigid3d"  )

                ###############
                ## Rate of angular Deformation  (2 sections)
                ###############
                pos1 = [0.0,0.0,0.0]
                pos2 = [0.0,0.0,0.0]
                pos3 = [0.0,0.0,0.0]
                pos = [pos1, pos2, pos3]

                distance1 = [0.0,0.2,0.0]
                distance2 = [0.0,0.2,0.0]
                distance3 = [0.0,0.2,0.0]
                _distance = [distance1, distance2, distance3]

                ddistance1 = [0.0,0.0,0.0]
                ddistance2 = [0.0,0.0,0.0]
                ddistance3 = [0.0,0.0,0.0]
                _ddistance = [ddistance1, ddistance2, ddistance3]

                rateAngularDeformNode = rootNode.createChild('rateAngularDeform')
                rateAngularDeformMO = rateAngularDeformNode.createObject('MechanicalObject', template='Vec3d', name='rateAngularDeformMO', position=pos, velocity='0 0 0 0 0 0 0 0 0') # (2 series of 3 angles for 2 sections. we suppose that the lenght is 10 for each)
                # BeamHookeLawForce = rateAngularDeformNode.createObject('CosseratInternalActuation', crossSectionShape='circular', length='10 10 10', radius='0.5', youngModulus='5e6')
                BeamHookeLawForce = rateAngularDeformNode.createObject('CosseratInternalActuation', name="BeamHookeLawForce",  crossSectionShape='circular', length='10 10 10', radius='0.5',
                youngModulus='1e6',distance=_distance, ddistance=_ddistance, tension=_tension)
                rateAngularDeformNode.createObject('PythonScriptController', classname="TensionComputing")

                ##############
                ## Frames
                ##############
                # the node of the frame needs to inherit from rigidBaseMO and rateAngularDeform
                mappedFrameNode = rigidBaseNode.createChild('MappedFrames')
                rateAngularDeformNode.addChild(mappedFrameNode)
                framesMO = mappedFrameNode.createObject('MechanicalObject', template='Rigid3d', name="FramesMO", position="0.5 0 0  0 0 0 1  5 0 0 0 0 0 1   10 0 0  0 0 0 1   15 0 0 0 0 0 1  20 0 0  0 0 0 1 25 0 0  0 0 0 1 30 0 0  0 0 0 1", showObject='1', showObjectScale='1' )

                # The mapping has two inputs: RigidBaseMO and rateAngularDeformMO
                # one output: FramesMO

                inputMO = rateAngularDeformMO.getLinkPath() # + " " + RigidBaseMO.getLinkPath()
                #inputMO = rateAngularDeformMO.getLinkPath()
                inputMO_rigid = RigidBaseMO.getLinkPath()
                outputMO = framesMO.getLinkPath()
                # TODO:
                mappedFrameNode.createObject('DiscretCosseratMapping', curv_abs_input='0 10 20 30', curv_abs_output='0.5 5 10 15 20 25 30', input1=inputMO, input2=inputMO_rigid,output=outputMO, debug='0' )

                #### CylinderGridTop
                CylinderCollision = mappedFrameNode.createChild('CylinderCollision')

                # CylinderCollision.createObject('MeshSTLLoader', filename=path+'trunk.stl', name='loader', rotation='0 90 0', scale='0.155')
                CylinderCollision.createObject('CylinderGridTopology', name="loader", nx="8", ny="8", nz="20", length="30", radius="0.5", axis="1 0 0" )
                CylinderCollision.createObject('Mesh', src='@loader')
                CylinderCollision.createObject('MechanicalObject', template='Vec3d')
                CylinderCollision.createObject('Triangle')
                CylinderCollision.createObject('SkinningMapping', nbRef='2')

                # rootNode.createObject('BilateralInteractionConstraint', template='Rigid3d', object2='@rigidBase/MappedFrames/FramesMO', object1='@targetPos/target', first_point='0', second_point='6')


                return rootNode
