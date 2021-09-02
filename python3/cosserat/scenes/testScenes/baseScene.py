# -*- coding: utf-8 -*-
"""Basic scene using Cosserat in SofaPython3.

Based on the work done with SofaPython. See POEMapping.pyscn.
"""

__authors__ = ("emenager, yadagolo")
__contact__ = ("etienne.menager@ens-rennes.fr, yinoussa.adagolodjo@inria.fr")
__version__ = "1.0.0"
__copyright__ = "(c) 2020,Inria"
__date__ = "March 1 2020"
import Sofa

def createScene(rootNode):
    """Classical function to create a scene with SofaPython3.

    Parameters:
    ----------
        rootNode: <Sofa object>
            The node of the scene.

    Returns:
    -------
        None.

    """

    rootNode.addObject('RequiredPlugin', pluginName='SoftRobots SofaPython3 SofaSparseSolver CosseratPlugin SofaConstraint SofaDeformable SofaImplicitOdeSolver', printLog='0')
    rootNode.addObject('VisualStyle', displayFlags='hideVisualModels showBehaviorModels showCollisionModels hideBoundingCollisionModels showForceFields hideInteractionForceFields hideWireframe')


    rootNode.addObject('FreeMotionAnimationLoop')
    rootNode.addObject('GenericConstraintSolver', tolerance=1e-10, maxIterations=5000, printLog=0)

    rootNode.gravity = [0, -9810, 0]
    rootNode.dt= 0.01

    rootNode.addObject('EulerImplicitSolver', firstOrder=0, rayleighStiffness=1.0, rayleighMass=0.10)
    rootNode.addObject('SparseLDLSolver', name='solver')
    rootNode.addObject('GenericConstraintCorrection')

    #################################
    ##           RigidBase         ##
    #################################
    rigidBaseNode = rootNode.addChild('rigidBase')
    RigidBaseMO = rigidBaseNode.addObject('MechanicalObject', template='Rigid3d',
                                             name="RigidBaseMO", position=[0, 0, 0, 0, 0, 0, 1], showObject=1,
                                             showObjectScale=2.)
    rigidBaseNode.addObject('RestShapeSpringsForceField', name='spring', stiffness=5000,
                               angularStiffness=5000, external_points=0, mstate="@RigidBaseMO", points=0,
                               template="Rigid3d")
    # rigidBaseNode.addObject('UniformMass', totalMass='0.1' )

    #################################
    ## Rate of angular Deformation ##
    #################################

    # Define: the number of section, the total lenght and the lenght of each beam.
    nbSectionS = 3
    tot_length = 30.0
    lengthS = tot_length / nbSectionS

    # Define: the length of each beam in a list, the positions of eahc beam
    # (flexion, torsion), the abs of each section
    positionS = []
    longeurS = []
    temp = 0.
    curv_abs_inputS = [0.0]
    for i in range(nbSectionS):
        positionS.append([0, 0, 0])
        longeurS.append((((i+1)*lengthS) - i*lengthS))
        temp += longeurS[i]
        curv_abs_inputS.append(temp)
    curv_abs_inputS[nbSectionS] = tot_length

    # Define: sofa elements
    rateAngularDeformNode = rootNode.addChild('rateAngularDeform')
    rateAngularDeformMO = rateAngularDeformNode.addObject('MechanicalObject', template='Vec3d',
                                                          name='rateAngularDeformMO',
                                                          position=positionS, showIndices=0)
    rateAngularDeformNode.addObject('BeamHookeLawForceField',
                                    crossSectionShape='circular', length=longeurS,
                                    radius=0.50, youngModulus=5e6)
    # This is an option, one can add mass to be in complete dynamics simulation
    # rateAngularDeformNode.addObject('UniformMass', totalMass='0.1' )

    #################################
    ##             Frame           ##
    #################################
    # Define: the number of frame and the length between each frame.
    nbFramesF = 6
    lengthF = tot_length / nbFramesF

    # Define: the abs of each frame and the position of each frame.
    framesF = []
    curv_abs_outputF = []
    for i in range(nbFramesF):
        sol = i * lengthF
        framesF.append([sol, 0, 0,  0, 0, 0, 1])
        curv_abs_outputF.append(sol)

    framesF.append([tot_length, 0, 0, 0, 0, 0, 1])
    curv_abs_outputF.append(tot_length)

    # The node of the frame needs to inherit from rigidBaseMO and rateAngularDeform
    mappedFrameNode = rigidBaseNode.addChild('MappedFrames')
    rateAngularDeformNode.addChild(mappedFrameNode)
    framesMO = mappedFrameNode.addObject('MechanicalObject', template='Rigid3d',
                                         name="FramesMO", position=framesF,
                                         showObject=1, showObjectScale=1)

    # The mapping has two inputs: RigidBaseMO and rateAngularDeformMO one outputMO: FramesMO
    inputMO = rateAngularDeformMO.getLinkPath()
    inputMO_rigid = RigidBaseMO.getLinkPath()
    outputMO = framesMO.getLinkPath()

    mappedFrameNode.addObject('DiscreteCosseratMapping', curv_abs_input=curv_abs_inputS,
                              curv_abs_output=curv_abs_outputF, input1=inputMO, input2=inputMO_rigid,
                              output=outputMO, debug=0)
    # Target to reach with the tip of the model (cable or rod)
    targetPosNode = rootNode.addChild('targetPos')
    targetPosNode.addObject('MechanicalObject', template='Rigid3d', name='target', position='5 15 10 0 0 0 1',
                            showObject='1', showObjectScale='0.1')
    targetPosNode.addObject('RestShapeSpringsForceField', name='spring', stiffness="5000000000",
                            angularStiffness="500000000000", external_points="0", mstate="@target", points="0",
                            template="Rigid3d")
    # Constraint to drag the tip to the target, this can also be replace with the RestShapeSpringsForceField
    rootNode.addObject('BilateralInteractionConstraint', template='Rigid3d', object2='@rigidBase/MappedFrames/FramesMO',
                       object1='@targetPos/target', first_point='0', second_point=str(len(curv_abs_outputF) - 1))

    return rootNode
