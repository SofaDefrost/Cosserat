# -*- coding: utf-8 -*-
"""
    Create the scene with the 
    Units: mm, kg, s.
"""

__authors__ = ("Yinoussa:Younes")
__contact__ = ("adagolodjo@protonmail.com")
__version__ = "1.0.0"
__copyright__ = "(c) 2021, Inria"
__date__ = "Juin 3 2021"

import Sofa
import os


#constants
GRAVITY = 9.810
TOT_MASS = 0.1
YOUNG_MODULUS=7e9
DENSITY = 0.02


def createScene(rootNode):

    rootNode.addObject('RequiredPlugin', name='SoftRobots')
    rootNode.addObject('RequiredPlugin', name='BeamAdapter')
    rootNode.addObject('RequiredPlugin', name='SofaPython3')
    rootNode.addObject('RequiredPlugin', name='SofaSparseSolver')
    rootNode.addObject('RequiredPlugin', name='SofaOpenglVisual')
    rootNode.addObject('RequiredPlugin', name='SofaConstraint')
    rootNode.addObject('RequiredPlugin', name='SofaLoader')
    rootNode.addObject('RequiredPlugin', name='SofaImplicitOdeSolver')
    rootNode.addObject('RequiredPlugin', name='SofaMeshCollision')
    rootNode.addObject('RequiredPlugin', name='SofaRigid')
    rootNode.addObject('RequiredPlugin', name='CosseratPlugin')
    rootNode.addObject('RequiredPlugin', name='SofaDeformable')
    rootNode.addObject('RequiredPlugin', name='SofaGeneralLinearSolver')
    rootNode.addObject('RequiredPlugin', name='SofaGeneralRigid')

    rootNode.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels showCollisionModels hideBoundingCollisionModels showForceFields hideInteractionForceFields hideWireframe')
    rootNode.findData('dt').value=0.01
    rootNode.findData('gravity').value= [0., 0., -GRAVITY]

    rootNode.addObject('FreeMotionAnimationLoop')
    rootNode.addObject('GenericConstraintSolver', tolerance=1e-6, maxIterations=1000)
    

    #Define: the total lenght of the beam
    tot_length = 6.67e-3

    #Define: the number of section, the total lenght and the lenght of each beam.
    nbSectionS = 1
    lengthS = tot_length / nbSectionS

    #Define: the number of frame and the lenght beetween each frame.
    nbFramesF = 15
    lengthF = tot_length /nbFramesF

    points = []
    position = []
    lines = []

    position3D = []
    for i in range(nbFramesF):
        sol = i * lengthF
        points.append(i)
        position.append([sol, 0, 0, 0, 0, 0, 1])
        position3D.append([sol, 0, 0])
        if i!= nbFramesF-1:
            lines+= [i, i+1]
    print("=============> position : ", position)
    
    ############################################################################
    ############################################################################
    ###########                       COSSERAT                       ###########
    ############################################################################
    ############################################################################

    #################################
    ##           RigidBase         ##
    #################################
    cableNode = rootNode.addChild('cableNode')
    cableNode.addObject('EulerImplicitSolver', rayleighStiffness="0.2", rayleighMass='1.')
    cableNode.addObject('SparseLDLSolver', name='solver')
    cableNode.addObject('GenericConstraintCorrection')

    rigidBaseNode= cableNode.addChild('rigidBase')
    RigidBaseMO = rigidBaseNode.addObject('MechanicalObject', template='Rigid3d',
                                             name="RigidBaseMO", position= [0.,0.,0., 0,0,0,1], showObject=1,
                                             showObjectScale=2.)
    rigidBaseNode.addObject('RestShapeSpringsForceField', name='spring', stiffness="5.e8", angularStiffness="5.e8",external_points="0", mstate="@RigidBaseMO", points="0", template="Rigid3d")
    
    #################################
    ## Rate of angular Deformation ##
    #################################

    #Define: the lenght of each beam in a list, the positions of each beam
    
    positionS = []
    longeurS = []
    sum = 0.
    curv_abs_inputS = []
    curv_abs_inputS.append(0.0)
    for i in range(nbSectionS):
        positionS.append([0, 0, 0])
        longeurS.append((((i+1)*lengthS) - i*lengthS))
        sum += longeurS[i]
        curv_abs_inputS.append(sum)

    curv_abs_inputS[nbSectionS] = tot_length
    
    
    """ Define: angular rate which is the torsion(x) and bending (y, z) of each section """
    val = 3.14/(2*tot_length)
    print ('===< val :', val)
    rateAngularDeformNode = cableNode.addChild('rateAngularDeform')
    rateAngularDeformMO = rateAngularDeformNode.addObject('MechanicalObject',
                                    template='Vec3d', name='rateAngularDeformMO',
                                    position=positionS, showIndices=0, rest_position=[0.0, 0., 0])
    BeamHookeLawForce = rateAngularDeformNode.addObject('BeamHookeLawForceField',
                                    crossSectionShape='rectangular', length=longeurS,
                                    lengthZ=1.37e-3, lengthY=0.035e-3, radius=2., youngModulus=2e5)
    #rateAngularDeformNode.addObject('RestShapeSpringsForceField', name='spring', stiffness="5.e8",
                               #external_points="0", mstate="@rateAngularDeformMO", points="0", template="Vec3d")


    #################################
    ##             Frame           ##
    #################################
    #Define: local frames related to each section and parameters curv_abs_outputF 
    framesF = []
    curv_abs_outputF = []
    for i in range(nbFramesF):
        sol = i * lengthF
        framesF.append([sol, 0, 0,  0, 0, 0, 1])
        curv_abs_outputF.append(sol)

    framesF.append([tot_length, 0, 0, 0, 0, 0, 1])
    curv_abs_outputF.append(tot_length)
    print("=============> framesF : ", framesF)

    #The node of the frame needs to inherit from rigidBaseMO and rateAngularDeform
    mappedFrameNode = rigidBaseNode.addChild('MappedFrames')
    rateAngularDeformNode.addChild(mappedFrameNode)

    framesMO = mappedFrameNode.addObject('MechanicalObject', template='Rigid3d',
                                                name="FramesMO", position=framesF,
                                            showObject=1, showObjectScale=0.0002)
    mappedFrameNode.addObject('UniformMass', totalMass="2.2e-10", showAxisSizeFactor='0')
    forceFied = mappedFrameNode.createObject('ConstantForceField', indices="15", forces="0 1e-7 0 0 0 0")

    # The mapping has two inputs: RigidBaseMO and rateAngularDeformMO
    #                 one output: FramesMO
    inputMO = rateAngularDeformMO.getLinkPath()
    inputMO_rigid = RigidBaseMO.getLinkPath()
    outputMO = framesMO.getLinkPath()

    mappedFrameNode.addObject('DiscreteCosseratMapping', curv_abs_input= curv_abs_inputS,
                                 curv_abs_output=curv_abs_outputF, input1=inputMO, input2=inputMO_rigid,
                                 output=outputMO, debug=0)
    cableNode.addObject('MechanicalMatrixMapper', template='Vec3,Rigid3', object1=inputMO,
                        object2=inputMO_rigid, nodeToParse=mappedFrameNode.getLinkPath())

    ############################################################################
    ############################################################################
    ###########                        MAPPING                       ###########
    ############################################################################
    ############################################################################

    cosCollisionPoints = mappedFrameNode.addChild('cosCollisionPoints')
    cosColliMeca = cosCollisionPoints.addObject('MechanicalObject', name="cosColliPoins", template="Vec3d", 
                                                position = position3D)
    cosCollisionPoints.addObject('SkinningMapping', nbRef='2')

    return rootNode
