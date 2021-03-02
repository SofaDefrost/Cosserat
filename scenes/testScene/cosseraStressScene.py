# -*- coding: utf-8 -*-

import Sofa
import SofaPython
from math import sin,cos, sqrt, pi
import os
path = os.path.dirname(os.path.abspath(__file__))+'/mesh/'


class Animation(Sofa.PythonScriptController):

    def __init__(self, rigidBaseNode, rateAngularDeformNode):
        self.rigidBaseNode = rigidBaseNode
        self.rateAngularDeformNode = rateAngularDeformNode

        self.rate = 1;
        self.angularRate=0.1;
        return;


    def initGraph(self, nodeRigid):
        self.rigidBaseMO = self.rigidBaseNode.getObject('RigidBaseMO')
        self.rateAngularDeformMO = self.rateAngularDeformNode.getObject('rateAngularDeformMO')


def createScene(rootNode):

        rootNode.createObject('RequiredPlugin', pluginName='SoftRobots SofaPython SofaSparseSolver CosseratPlugin BeamAdapter', printLog='0')
        rootNode.createObject('VisualStyle', displayFlags='hideVisualModels showBehaviorModels showCollisionModels hideBoundingCollisionModels showForceFields hideInteractionForceFields hideWireframe')

        rootNode.createObject('FreeMotionAnimationLoop')
        rootNode.createObject('GenericConstraintSolver', tolerance="1e-10", printLog='0')
        rootNode.gravity = "0 -9810 0"
        rootNode.dt="0.01"
        
        ###############hresho
        ## Solver
        ###############
        rootNode.createObject('EulerImplicitSolver', firstOrder="0", rayleighStiffness="1.0", rayleighMass='0.10')
        rootNode.createObject('SparseLDLSolver', name='solver')
        rootNode.createObject('GenericConstraintCorrection')

        ###############hresho
        ## RigidBase
        ###############
        rigidBaseNode= rootNode.createChild('rigidBase')
        RigidBaseMO = rigidBaseNode.createObject('MechanicalObject', template='Rigid3d', name="RigidBaseMO", position="0 0 0  0 0 0. 1", showObject='1', showObjectScale='0.1', velocity='0 0 0.0 0.0 0 0' )
        rigidBaseNode.createObject('RestShapeSpringsForceField', name='spring', stiffness="50000000", angularStiffness="5000000", external_points="0", mstate="@RigidBaseMO", points="0", template="Rigid3d"  )

        ###############
        ## Rate of angular Deformation  (2 sections)
        ###############
        #pos = pi
        array1 = [0.0,0.0,0.0]
        array2 = [0.0,0.0,0.0]
        array3 = [0.0,0.0,0.0]
        pos = [array1, array2, array3]

        rateAngularDeformNode = rootNode.createChild('rateAngularDeform')
        rateAngularDeformMO = rateAngularDeformNode.createObject('MechanicalObject', template='Vec3d', name='rateAngularDeformMO', position=pos, velocity='0 0 0.0 0 0 0  0 0 0') # (2 series of 3 angles for 2 sections. we suppose that the lenght is 10 for each)
        forceFied = rateAngularDeformNode.createObject('ConstantForceField', indices="0 1 2", force="0 0 -5.2e5")
        BeamHookeLawForce = rateAngularDeformNode.createObject('BeamHookeLawForceField', crossSectionShape='circular', length='10 10 10', radius='0.5', youngModulus='5e6')

        ##############
        ## Frames
        ##############
        # the node of the frame needs to inherit from rigidBaseMO and rateAngularDeform
        mappedFrameNode = rigidBaseNode.createChild('MappedFrames')
        rateAngularDeformNode.addChild(mappedFrameNode)
        framesMO = mappedFrameNode.createObject('MechanicalObject', template='Rigid3d', name="FramesMO", position="0.5 0 0  0 0 0 1  5 0 0 0 0 0 1   10 0 0  0 0 0 1   15 0 0 0 0 0 1  20 0 0  0 0 0 1 25 0 0  0 0 0 1 30 0 0  0 0 0 1", showObject='1', showObjectScale='1' )

        inputMO = rateAngularDeformMO.getLinkPath() 
        inputMO_rigid = RigidBaseMO.getLinkPath()
        outputMO = framesMO.getLinkPath()
        mappedFrameNode.createObject('DiscretCosseratMapping', curv_abs_input='0 10 20 30', curv_abs_output='0.5 5 10 15 20 25 30', input1=inputMO, input2=inputMO_rigid,output=outputMO, debug='0' )

        return rootNode
