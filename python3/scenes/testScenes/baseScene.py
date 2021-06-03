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
import os
import numpy as np

TempPath = "/home/stefan/Repos/Code/sofa/plugins/SoftRobots.Inverse/docs/examples/component/constraint/BendLabsEffector/Temp/"
class OrientationSweepController(Sofa.Core.Controller):

    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)

        self.RootNode = kwargs["RootNode"]

        self.LiveAngleDataPath = TempPath + "AngleData.txt"
        self.Data = np.array([])
        self.LastData = self.Data
        
    def onAnimateBeginEvent(self, dt):
        

#        self.RotationAngle = (self.RotationAngle + 0.5)%360
#        
#        
#        
#        Amplitude = 45
#        self.InclinationAngle = Amplitude #* np.abs(np.sin(np.deg2rad(2*self.RotationAngle)))
#        self.RInXY = np.cos(np.deg2rad(90-self.InclinationAngle))
#        self.Height = np.cos(np.deg2rad(self.InclinationAngle))
#        
#        (alpha, beta) = calcAlphaAndBeta(self.RotationAngle, self.Height, self.RInXY)
               
        try:
            self.Data = np.loadtxt(self.LiveAngleDataPath)
        except Exception as e:
            print("Warning: couldn't read pressures from file")
#            self.Data = self.LastData
#            self.RealPressures = self.Data[:6] # in kPa   
            return
        
        if self.Data.shape[0] == 0:
            return
        
        alpha = self.Data[1] # Actually, the sensor (or the arduino lib) doesn't follow our convention of the order of angles (first angle is in xz plane, second in yz plane)
        beta = self.Data[0]
        
        print("setting alpha to: " + str(alpha))
        print("setting beta to: " + str(beta))
        
#        self.BendLabsEffector.alpha.value = alpha
#        self.BendLabsEffector.beta.value = beta
#        
        
        #print("Constraints: " + str(self.FPAMO.constraint.value))
        
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

    rootNode.addObject('RequiredPlugin', pluginName='SoftRobots SofaPython3 SofaSparseSolver CosseratPlugin  SofaConstraint SofaDeformable SofaImplicitOdeSolver', printLog='0')
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
                                             name="RigidBaseMO", position= [0,0,0,0,0,0,1], showObject=1,
                                             showObjectScale=2.)
    rigidBaseNode.addObject('RestShapeSpringsForceField', name='spring', stiffness=5000,
                               angularStiffness=5000, external_points=0, mstate="@RigidBaseMO", points=0,
                               template="Rigid3d")
    #rigidBaseNode.addObject('UniformMass', totalMass='0.1' )

    #################################
    ## Rate of angular Deformation ##
    #################################

    #Define: the number of section, the total lenght and the lenght of each beam.
    nbSectionS = 3
    tot_length = 30.0
    lengthS = tot_length / nbSectionS

    #Define: the longueur of each beam in a list, the positions of eahc beam
    #(flexion, torsion), the abs of each section
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

    longeurS[nbSectionS-1] = longeurS[nbSectionS-1] + 1.
    curv_abs_inputS[nbSectionS] = tot_length

    print("=============> positionS : ", positionS)
    print("=============> longeurS : ", longeurS)
    print("=============> curv_abs_inputS : ", curv_abs_inputS)


    #Define: sofa elements
    rateAngularDeformNode = rootNode.addChild('rateAngularDeform')
    rateAngularDeformMO = rateAngularDeformNode.addObject('MechanicalObject',
                                    template='Vec3d', name='rateAngularDeformMO',
                                    position=positionS, showIndices=0)
    BeamHookeLawForce = rateAngularDeformNode.addObject('BeamHookeLawForceField',
                                    crossSectionShape='circular', length=longeurS,
                                    radius=0.50, youngModulus=5e6)
    #rateAngularDeformNode.addObject('UniformMass', totalMass='0.1' )


    #################################
    ##             Frame           ##
    #################################

    #Define: the number of frame and the lenght beetween each frame.
    nbFramesF = 6
    lengthF = tot_length /nbFramesF


    #Define: the abs of each frame and the position of each frame.
    framesF = []
    curv_abs_outputF = []
    for i in range(nbFramesF):
        sol = i * lengthF
        framesF.append([sol, 0, 0,  0, 0, 0, 1])
        curv_abs_outputF.append(sol)

    framesF.append([tot_length, 0, 0, 0, 0, 0, 1])
    curv_abs_outputF.append(tot_length)


    #The node of the frame needs to inherit from rigidBaseMO and rateAngularDeform
    mappedFrameNode = rigidBaseNode.addChild('MappedFrames')
    rateAngularDeformNode.addChild(mappedFrameNode)

    framesMO = mappedFrameNode.addObject('MechanicalObject', template='Rigid3d',
                                                name="FramesMO", position=framesF,
                                            showObject=1, showObjectScale=1)

    # The mapping has two inputs: RigidBaseMO and rateAngularDeformMO
    #                 one output: FramesMO
    inputMO = rateAngularDeformMO.getLinkPath()
    inputMO_rigid = RigidBaseMO.getLinkPath()
    outputMO = framesMO.getLinkPath()

    mappedFrameNode.addObject('DiscretCosseratMapping', curv_abs_input= curv_abs_inputS,
                                 curv_abs_output=curv_abs_outputF, input1=inputMO, input2=inputMO_rigid,
                                 output=outputMO, debug=0)




    targetPosNode= rootNode.addChild('targetPos')
    targetPosNode.addObject('MechanicalObject', template='Rigid3d', name='target', position='5 15 10 0 0 0 1', showObject='1', showObjectScale='0.1')
    targetPosNode.addObject('RestShapeSpringsForceField', name='spring', stiffness="5000000000", angularStiffness="500000000000", external_points="0", mstate="@target", points="0", template="Rigid3d"  )


    rootNode.addObject('BilateralInteractionConstraint', template='Rigid3d', object2='@rigidBase/MappedFrames/FramesMO', object1='@targetPos/target', first_point='0', second_point=str(len(curv_abs_outputF)-1))
    
    rootNode.addObject(OrientationSweepController(name="OrientationSweeController", RootNode=rootNode))


    return rootNode
