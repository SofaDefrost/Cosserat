# -*- coding: utf-8 -*-
"""
Created on May 5 2022

@author: ckrewcun
"""

# -*- coding: utf-8 -*-

"""
    Python controller to simulate coaxial Cosserat beam chains, typically in the
    context of interventional instruments (e.g.: guidewire/catheter)
"""

import Sofa
import Sofa.constants.Key as Key
import numpy as np
import math
import logging
import warnings
from pyquaternion import Quaternion as Quat
from instrument import Instrument

# Python controller to simulate the navigation of coaxial instruments represented
# by two chains of Cosserat beam elements.
# Required structure of the scene :
# * rootNode
#
#   * Instrument0
#       * rigidBase
#           RigidBaseMO (Rigid)
#           spring (RestShapeSpringsForceField)
#           * MappedFrames
#               FramesMO (Rigid)
#               DiscreteCosseratMapping
#       * rateAngularDeform
#           rateAngularDeformMO (Vec3)
#           beamForceField (BeamHookeLawForceField, BeamPlasticLawForceField)
#           FixedConstraintOnStock (FixedConstraint)
#           * MappedFrames (same node as in the rigidBase node)
#
#   * Instrument1
#       [Same structure and names as in Instrument, except for constraintSprings]
#       * rateAngularDeform
#           * constraintWith0
#               constraintSpringsWith0 (StiffSpringForceField)
#                   NB: object1 = Instrument1, object2 = Instrument0
#
#   * Instrument2
#       [Same structure and names as in Instrument1]
#       * rateAngularDeform
#           * constraintWith0
#               constraintSpringsWith0 (StiffSpringForceField)
#                   NB: object1 = Instrument2, object2 = Instrument0
#           * constraintWith1
#               constraintSpringsWith1 (StiffSpringForceField)
#                  NB: object1 = Instrument2, object2 = Instrument1
#
# Init arguments :
#  - nbInstruments: number of simulated coaxial instruments
#  - instrumentFrameNumberVect : for each instrument, vector containing a number
#    of rigid frames spread uniformely on the instrument total length
#  - incrementDistance : distance of pushing/pulling the instrument at user
#    interaction
#  - incrementDirection : direction (vec3) along which the instruments are
#    navigated
#  - instrumentList : vector containing instances of Instrument objects (as defined
#    in instrument.py), to characterise each instrument properties
#  - curvAbsTolerance : distance threshold used to determine if two close nodes
#    should be merged (and considered as one)
#  - nbIntermediateConstraintFrames : number of intermediate coaxial frames added
#    when coaxial beam segments are detected. A higher number means a finer
#    application of constraints on the coaxial beam segments.
#  - constrainWithSprings : boolean indicating if the constrain enforced between
#    the coaxial beam segments shoudl be enforced by a spring ForceField instead
#    of an actual constraint. NB: this requires the scene to have a
#    RestShapeSpringsForceField component with a direction option (activeDirections)
#  - outputFilename : complete path of a file in which the sequence of inputs
#    made by the user during the simulation will be written. Each input is
#    associated with a time stamp, in order for the simulation to be automatically
#    played from this output file. By default, a file named "navigationScript.txt"
#    will be generated in the working directory if the export option (Ctrl+e) is
#    used
#  - inputFilename : complete path of a file describing a navigation simulation
#    with a list of user inputs and time stamps (indicating when the inputs have
#    to be replayed). Such a file can be created during a manually controlled
#    simulation using the 'outputFilename' parameter of this controller.
#
# /!\ In this controller, as it is the case in BeamAdapter, we assume that the
# instruments are given in the ordre of decreasing diameter, meaning that the
# outmost instrument should have index 0. This assumption is mostly used when
# dynamically recomputing the instruments discretisation: in overlapping segments,
# only the outmost instrument should be visible.
#
class CombinedInstrumentsController(Sofa.Core.Controller):

    # -------------------------------------------------------------------- #
    # -----                      Initialisation                      ----- #
    # -------------------------------------------------------------------- #

    def __init__(self, rootNode, solverNode,
                 nbInstruments,
                 instrumentFrameNumberVect,
                 incrementDistance,
                 incrementAngle,
                 incrementDirection,
                 instrumentList,
                 curvAbsTolerance,
                 nbIntermediateConstraintFrames = 0,
                 constrainWithSprings = False,
                 outputFilename="",
                 inputFilename="",
                 *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)

        ### Checking for the scene graph structure ###

        self.rootNode = rootNode
        self.solverNode = solverNode

        ### Reading the insertion velocity parameters ###

        self.incrementDistance = incrementDistance
        # TO DO : pass the minimal distance for constraints as an input parameter ?
        self.minimalDistanceForConstraint = incrementDistance * 5.0
        # TO DO: check that the number of coaxial frames provided is coherent with beam number and nbIntermediateConstraintFrames
        self.nbIntermediateConstraintFrames = nbIntermediateConstraintFrames
        self.incrementAngle = incrementAngle
        self.incrementDirection = incrementDirection
        self.instrumentFrameNumberVect = instrumentFrameNumberVect
        self.curvAbsTolerance = curvAbsTolerance

        self.constrainWithSprings = constrainWithSprings

        ### Controller settings ###

        self.instrumentList = instrumentList

        # Number of simulated instruments
        self.nbInstruments = nbInstruments

        # Curvilinear Abscissa of the tip of each instrument (modified by up/down arrows)
        self.tipCurvAbsVect = np.zeros(nbInstruments)

        # Rotation angle for each instrument (modified by left/right arrows)
        self.rotationAngleVect = np.zeros(nbInstruments)

        # Index of the currently navigated instrument, and associated rigidBase
        self.currentInstrumentId = 0

        instrumentNodeName = "Instrument0"
        instrumentNode = self.solverNode.getChild(str(instrumentNodeName))
        if instrumentNode is None:
            raise NameError("[CombinedInstrumentsController]: Node \'{}\' not found. Your scene should "
                            "contain a node named \'{}\' among the children of the root node in order "
                            "to use this controller".format(instrumentNodeName, instrumentNodeName))

        controlPointNodeName = "controlPointNode0"
        controlPointNode = self.rootNode.getChild(str(controlPointNodeName))
        if controlPointNode is None:
            raise NameError("[CombinedInstrumentsController]: Node \'{}\' not found. Your scene should "
                            "contain a node named \'{}\' among the children of the root node in order "
                            "to use this controller".format(controlPointNodeName, controlPointNodeName))

        self.currentInstrumentControlPointNode = controlPointNode

        rigidBaseNode = instrumentNode.getChild('rigidBase')
        if rigidBaseNode is None:
            raise NameError("[CombinedInstrumentsController]: Node \'rigidBase\' "
                            "not found. Your scene should contain a node named "
                            "\'rigidBase\' among the children of the \'{}\' "
                            "node, where the base and rigid frames of the "
                            "Cosserat model are defined".format(instrumentNodeName))

        self.currentInstrumentRigidBaseNode = rigidBaseNode

        ### Checking if there is an input script file to control the simulation ###

        self.inputFilename = inputFilename
        if inputFilename != "":
            self.hasControlScript = True
            parsingOutput = self.parseControlScript(inputFilename)
            self.actionTimeList = parsingOutput['actionTimeList']
            self.actionList = parsingOutput['actionList']

            self.nbActions = len(self.actionTimeList)
            self.nextActionId = 0
        else:
            self.hasControlScript = False

        ### Checking if there is an output script file to record user inputs during simulation ###

        self.exportNavigationScript = False
        if outputFilename == "":
            self.outputFilename = "navigationScript.txt"
            warningStr = "[CosseratNavigationController]: No filename was provided to export "
            warningStr += "the simulation controls in the form of a scripted file (with export "
            warningStr += "option Ctrl+e). If the export option is used anyway, a control script file "
            warningStr += "named \"navigationScript.txt\" will be generated in your working directory. "
            warningStr += "If the file already exists, it will be replaced. In such a case, consider "
            warningStr += "saving the existing file before using Ctrl+e"
            logging.info(warningStr)
        else:
            self.outputFilename = outputFilename

        ### Additional settings ###

        # Computing the incremental quaternions for rotation
        # Taking the direction of insertion as rotation axis
        qw = math.cos(math.radians(self.incrementAngle)/2)
        plusQuat = self.incrementDirection * math.sin(math.radians(self.incrementAngle)/2)
        minusQuat = -plusQuat
        self.plusQuat = Quat(np.insert(plusQuat, 0, qw))
        self.minusQuat = Quat(np.insert(minusQuat, 0, qw))

        # constructs a grid of indices to access only position DoFs of the rigid particle
        self.posDoFsIdGrid = np.ix_([0], [0, 1, 2])
        # constructs a grid of indices to access only orientation DoFs of the rigid particle
        self.quatDoFsIdGrid = np.ix_([0], [3, 4, 5, 6])

        self.totalTime = 0.0
        self.dt = self.rootNode.findData('dt').value
        self.nbIterations = 0

        ### Verification of parameter coherence ###

        # Number of coaxial frames

        for instrumentId in range(self.nbInstruments):

            # Retrieving the coaxial frame node
            instrumentNodeName = "Instrument" + str(instrumentId)
            instrumentNode = self.solverNode.getChild(str(instrumentNodeName))
            if instrumentNode is None:
                raise NameError("[CombinedInstrumentsController]: Node \'{}\' not found. Your scene should "
                                "contain a node named \'{}\' among the children of the root node in order "
                                "to use this controller".format(instrumentNodeName, instrumentNodeName))

            cosseratMechanicalNode = instrumentNode.getChild('rateAngularDeform')
            if cosseratMechanicalNode is None:
                raise NameError("[CombinedInstrumentsController]: Node \'rateAngularDeform\' "
                                "not found. Your scene should contain a node named "
                                "\'rateAngularDeform\' among the children of the \'{}\' "
                                "node, gathering the mechanical Cosserat components "
                                "(MechanicalObject, Cosserat forcefield)".format(instrumentNodeName))

            rigidBaseNode = instrumentNode.getChild('rigidBase')
            if rigidBaseNode is None:
                raise NameError("[CombinedInstrumentsController]: Node \'rigidBase\' "
                                "not found. Your scene should contain a node named "
                                "\'rigidBase\' among the children of the \'{}\' "
                                "node, where the base and rigid frames of the "
                                "Cosserat model are defined".format(instrumentNodeName))

            coaxialFramesNode = rigidBaseNode.getChild('coaxialSegmentFrames')
            if coaxialFramesNode is None:
                raise NameError("[CombinedInstrumentsController]: Node \'coaxialSegmentFrames\' "
                                "not found. The \'rigidBase\' node should have a child "
                                "node called \'coaxialSegmentFrames\', in which the Cosserat "
                                "mapping tracking coaxial segments should be defined.")

            nbBeamsInMechanicalObject = len(cosseratMechanicalNode.rateAngularDeformMO.position)
            nbCoaxialFrames = len(coaxialFramesNode.coaxialFramesMO.position)
            if (nbCoaxialFrames != (nbIntermediateConstraintFrames+1) * nbBeamsInMechanicalObject + 1):
                raise ValueError("[CombinedInstrumentsController]: The number of coaxial frames ({}) "
                                 "is not coherent with the number of beams ({}) provided for instrument{}, "
                                 "and this controller nbIntermediateConstraintFrames parameter. For a "
                                 "given number of beams *nbBeams* and a given *nbIntermediateConstraintFrames* "
                                 "parameter, the number of coaxial frames should be : "
                                 "(nbIntermediateConstraintFrames+1)*nbBeams + 1."
                                 "".format(nbCoaxialFrames, nbBeamsInMechanicalObject, instrumentId))

        # TO DO: check if a RestShapeSpringsForceField component is provided
        # if the input parameter is True



    # -------------------------------------------------------------------- #
    # -----                     Animation events                     ----- #
    # -------------------------------------------------------------------- #

    def onAnimateBeginEvent(self, event):  # called at each begin of animation step

        # Before running the dynamic navigation methods, we check if the simulation
        # is scripted. If so, we execute the corresponding user inputs

        timeAtStepEnd = self.totalTime + self.rootNode.findData('dt').value

        if self.hasControlScript:
            # The user actions are read from a script file

            nextActionTime = self.actionTimeList[self.nextActionId]
            if timeAtStepEnd >= nextActionTime:
                # Executing the next action, otherwise nothing to do
                nextAction = self.actionList[self.nextActionId]
                if nextAction == 'upArrow':
                    self.moveForward(self.incrementDistance, self.incrementDirection)
                elif nextAction == 'downArrow':
                    self.moveBackward(self.incrementDistance, self.incrementDirection)
                elif nextAction == 'rightArrow':
                    self.rotateClockwise()
                elif nextAction == 'leftArrow':
                    self.rotateCounterclockwise()
                elif nextAction == 'switchInstrument0':
                    self.switchInstrument(0)
                elif nextAction == 'switchInstrument1':
                    self.switchInstrument(1)
                elif nextAction == 'switchInstrument2':
                    self.switchInstrument(2)
                else:
                    pass # We made sure that this doesn't happen during parsing

                # Checking if this was the last action
                if self.nextActionId == self.nbActions-1:
                    endMessage = "[CombinedInstrumentsController]: all actions of "
                    endMessage += "script file \'{}\' ".format(self.inputFilename)
                    endMessage += "were executed. Pausing the animation"
                    print(endMessage)
                    self.rootNode.findData('animate').value = 0
                    self.hasControlScript = False # In case the simulation is animated again
                else:
                    # Moving to next action
                    self.nextActionId += 1


        ### Dynamic navigation ###

        # Define the vector which contains the curvilinear abscissas of all
        # represented nodes of the beam elements (i.e. similar to curb_abs_input
        # in the Cosserat mapping)
        simulatedNodeCurvAbs = []

        # Define the vector which contains the curvilinear abscissas of all
        # represented rigid frames for the Cosserat components (i.e. similar to
        # curv_abs_output in the Cosserat mapping)
        simulatedFrameCurvAbs = []

        #----- Step 0 : base check on instrument's lengths -----#

        # Check that all instruments haven't been pushed further than their
        # rest length. Otherwise the tip curvilinear abscissa of such an instrument
        # is set back to the instrument's length

        for instrumentId in range (0, self.nbInstruments):
            instrumentLength = self.instrumentList[instrumentId].getTotalLength()
            if self.tipCurvAbsVect[instrumentId] > instrumentLength:
                self.tipCurvAbsVect[instrumentId] = instrumentLength


        #----- Step 1 : instrument's tip curvilinear abscissa and combined length -----#

        # Find:
        #    - which instruments are 'out' and simulated (for which tipCurvAbs > 0)
        #    - tipCurvAbs of the most distal instrument = the length of the combined
        #      instruments

        xBeginVect = []

        combinedInstrumentLength = self.getInstrumentConfiguration(xBeginVect)

        # if the totalLength is 0, we move the first instrument and update the
        # information on the instruments configuration
        if combinedInstrumentLength < self.curvAbsTolerance:
            xBeginVect = []
            self.tipCurvAbsVect[0] = self.incrementDistance
            combinedInstrumentLength = self.getInstrumentConfiguration(xBeginVect)

        # Adding the base point if at least one instrument is out
        if combinedInstrumentLength > 0.:
            simulatedNodeCurvAbs.append(0.)
            simulatedFrameCurvAbs.append(0.)


        #----- Step 2 : computation of the points of interest -----#

        # Computation of the curvilinear abscissas for the points of interest,
        # and Cosserat Rigid frames, which are to be simulated.

        instrumentIdsForNodeVect = []
        instrumentIdsForFrameVect = []

        result = self.computeNodeCurvAbs(xBeginVect, simulatedNodeCurvAbs,
                                         simulatedFrameCurvAbs,
                                         instrumentIdsForNodeVect,
                                         instrumentIdsForFrameVect)

        instrumentIdsForNodeVect = result['instrumentIdsForNodeVect']
        instrumentIdsForFrameVect = result['instrumentIdsForFrameVect']
        decimatedNodeCurvAbs = result['decimatedNodeCurvAbs']
        decimatedFrameCurvAbs = result['decimatedFrameCurvAbs']

        # print("decimatedNodeCurvAbs : {}".format(decimatedNodeCurvAbs))
        # print("instrumentIdsForNodeVect : {}".format(instrumentIdsForNodeVect))


        #----- Step 3 :  -----#

        # Once a discretisation - common to all instruments - is computed, we
        # apply it to the Cosserat components. For more details, see comments
        # on the method definition.
        result = self.updateInstrumentComponents(decimatedNodeCurvAbs, decimatedFrameCurvAbs,
                                                 instrumentIdsForNodeVect, instrumentIdsForFrameVect)

        instrumentLastNodeIds = result['instrumentLastNodeIds']


        #----- Step 4 :  -----#

        # Interpolation of positions, velocities, and rest shapes.
        # Once the scene components are updated accordingly to the new discretisation,
        # we have to make sure that this discretisation is coherent with mechanical
        # quantities at the previous timestep.

        self.interpolateMechanicalQuantities(decimatedNodeCurvAbs, instrumentLastNodeIds)


        self.totalTime = timeAtStepEnd
        self.nbIterations = self.nbIterations + 1



    def onAnimateEndEvent(self, event):  # called at each end of animation step
        pass



    # -------------------------------------------------------------------- #
    # -----                    Auxiliary methods                     ----- #
    # -------------------------------------------------------------------- #

    # This method retrieves the global configuration of the instruments.
    #
    # Parameters:
    #    - xBeginVect: [output] vector containing for each instrument the curvilinear
    # abscissa of the instrument's begin point (proximal). This curvilinear abscissa
    # is expressed relatively to the base point, from which the instruments are
    # pushed/pulled
    #
    # Returns:
    #    - combinedInstrumentLength: furthest instrument tip (distal), among all
    # the insruments, which defines the length of insertion at a given time.
    def getInstrumentConfiguration(self, xBeginVect):

        combinedInstrumentLength = 0.

        # Computation of the curvilinear abscissas of the instrument tip
        # (distal end) and beginning (proximal end).
        # /!\ All curvilinear abscissas are expressed relatively to the base
        # point from which the instruments are pushed/pulled. Therefore, for
        # the tip, the curvilinear abscissa is >0, but for the beginning point
        # it is <0
        for instrumentId in range(0,self.nbInstruments):

            instrumentLength = self.instrumentList[instrumentId].getTotalLength()
            xEnd = self.tipCurvAbsVect[instrumentId]
            xBegin = xEnd - instrumentLength
            xBeginVect.append(xBegin)

            if xEnd > combinedInstrumentLength:
                combinedInstrumentLength = xEnd

        return combinedInstrumentLength


    # During navigation, when instruments are in motion, this method computes
    # the curvilinear abscissas for the points of interest, and Cosserat Rigid
    # frames, which are to be simulated. The points of interest correspond
    # to changes in the instrument shape/configuration. If the instrument is
    # entirely straight, the only two key points are the proximal end (begin
    # point) and the distal end (end point). The key points are stored by
    # means of their curvilinear abscissa. As before, this curvilinear abscissa
    # is expressed relatively to the combined configuration starting point.
    # Points for which the curvilinear abscissa is <0 are ignored.
    # Additionaly, the method completes two structures:
    #    - xBeginVect: vector containing for each instrument the curvilinear
    #      abscissa of the proximal end of the instrument. As the curvilinear
    #      abscissa is computed relatively to the point from which the
    #      instruments are pushed, it is negative for the proximal end
    #    - instrumentIdsForNodeVect: vector containing for each node a list
    #      of the indices of all instruments passing through the node
    #
    # Parameters:
    #    - xBeginVect: [input] vector containing for each instrument the curvilinear
    # abscissa of the instrument's begin point (proximal), expressed relatively
    # to the base point from which the instruments are pushed/pulled
    #    - simulatedNodeCurvAbs: [output] list of curvilinear abscissas corresponding
    # to the new nodes (i.e. the beam element extremities, inputs of the Cosserat
    # mapping)
    #    - simulatedFrameCurvAbs: [output] list of curvilinear abscissas corresponding
    # to the new Cosserat frames (i.e. outputs of the Cosserat mapping)
    #    - instrumentIdsForNodeVect: [output] vector containing for each node
    # the list of instruments which the node belongs to
    #    - instrumentIdsForFrameVect: [output] vector containing for each Cosserat
    # rigid frame the list of instruments which the frame belongs to
    def computeNodeCurvAbs(self, xBeginVect, simulatedNodeCurvAbs, simulatedFrameCurvAbs, instrumentIdsForNodeVect, instrumentIdsForFrameVect):

        maxCurvilinearAbscissa = 0. # Max curvilinear abscissa among the key points

        for instrumentId in range(0,self.nbInstruments):

            currentInstrument = self.instrumentList[instrumentId]

            # NB: in this method, we name as 'keyPoint' values which represent
            # curvilinear abscissas expressed in the instrument local frame,
            # ranging from 0 to the length of the instrument.
            # The values we name as 'curvAbs', however, are curvilinear abscissas
            # expressed relatively to the frame representing the insertion point
            # (common to all instruments). In this context, curvilinear abscissas
            # range between -instrumentLength and instrumentLength, negative values
            # representing the instrument parts which have not been deployed yet.

            proximalEndKeyPoint = 0.
            distalEndKeyPoint = currentInstrument.getTotalLength()

            instrumentKeyPoints = [proximalEndKeyPoint]
            instrumentKeyPoints += currentInstrument.getKeyPoints()
            instrumentKeyPoints.append(distalEndKeyPoint)

            nbBeamDistribution = currentInstrument.getNbBeamDistribution()
            instrumentLength = currentInstrument.getTotalLength()

            nbFrameSegmentsOnInstrument = self.instrumentFrameNumberVect[instrumentId]-1

            instrumentBeginCurvAbs = xBeginVect[instrumentId]

            for keyPointId in range( len(instrumentKeyPoints)-1 ):

                keyPoint = instrumentKeyPoints[keyPointId]
                nextKeyPoint = instrumentKeyPoints[keyPointId+1]

                keyPointCurvAbs = instrumentBeginCurvAbs + keyPoint
                nextKeyPointCurvAbs = instrumentBeginCurvAbs + nextKeyPoint

                # In any case, the key points are added to simulatedNodeCurvAbs
                # as soon as they are deployed

                if (keyPointCurvAbs > 0):
                    simulatedNodeCurvAbs.append(keyPointCurvAbs)
                    simulatedFrameCurvAbs.append(keyPointCurvAbs)

                # Additionaly, if part of the instrument between the two key points
                # is visible, we add a number of beam nodes according to the
                # interval beam density and visible length

                visibleIntervalLength = nextKeyPointCurvAbs - maxCurvilinearAbscissa

                if (visibleIntervalLength > 0):

                    # Compute the number of new nodes to add
                    intervalLength = nextKeyPoint - keyPoint
                    ratio = visibleIntervalLength / intervalLength
                    nbBeamsOnInterval = nbBeamDistribution[keyPointId]
                    nbNewNodes = int(nbBeamsOnInterval * ratio)

                    # Add the new nodes
                    for newNodeId in range(0, nbNewNodes):
                        newNodeCurvAbs = nextKeyPointCurvAbs - (newNodeId+1) * (intervalLength / nbBeamsOnInterval)
                        simulatedNodeCurvAbs.append(newNodeCurvAbs)

                    # Compute the number of new frames to add
                    # Contrarily to beam nodes, the output frames of the Cosserat mapping
                    # are distributed uniformly over the instrument.
                    # As we impose frames on the instrument keyPoints, we have to compute
                    # the number of frame segments on each keyPoint interval. The
                    # computation has to be done individually for each interval to avoid
                    # rounding errors.
                    nbFrameSegmentsOnInterval = int(nbFrameSegmentsOnInstrument * intervalLength / instrumentLength)
                    nbNewFrames = int(nbFrameSegmentsOnInterval * ratio)

                    # Add the new frames
                    for newFrameId in range(0, nbNewFrames):
                        newFrameCurvAbs = nextKeyPointCurvAbs - (newFrameId+1) * (intervalLength / nbFrameSegmentsOnInterval)
                        simulatedFrameCurvAbs.append(newFrameCurvAbs)

                    # Update the max curv. abs.
                    maxCurvilinearAbscissa = nextKeyPointCurvAbs

            # After the end of the for loop above, we just have to process the
            # instrument last key point
            lastKeyPoint = instrumentKeyPoints[len(instrumentKeyPoints)-1]
            lastKeyPointCurvAbs = instrumentBeginCurvAbs + lastKeyPoint
            if (lastKeyPointCurvAbs > 0):
                simulatedNodeCurvAbs.append(lastKeyPointCurvAbs)
                simulatedFrameCurvAbs.append(lastKeyPointCurvAbs)

        # endfor instrumentId in range(0,self.nbInstruments)


        # When all points of interest have been detected, we sort and filter
        # ther curv. abs' list.

        # First: sort the curv. abs. values
        # - Nodes
        sortedNodeCurvAbs = np.sort(simulatedNodeCurvAbs, kind='quicksort') # quicksort, heapsort, mergesort, timsort
        # - Frames
        sortedFrameCurvAbs = np.sort(simulatedFrameCurvAbs, kind='quicksort')

        # Second: remove the duplicated values, according to self.curvAbsTolerance
        # - Nodes
        indicesToRemove = []
        for curvAbsId in range(1, len(sortedNodeCurvAbs)):
            diffWithPrevious = abs(sortedNodeCurvAbs[curvAbsId] - sortedNodeCurvAbs[curvAbsId -1])
            if diffWithPrevious < self.curvAbsTolerance:
                indicesToRemove.append(curvAbsId)
        decimatedNodeCurvAbs = np.delete(sortedNodeCurvAbs, indicesToRemove)

        # - Frames
        indicesToRemove = []
        for curvAbsId in range(1, len(sortedFrameCurvAbs)):
            diffWithPrevious = abs(sortedFrameCurvAbs[curvAbsId] - sortedFrameCurvAbs[curvAbsId -1])
            if diffWithPrevious < self.curvAbsTolerance:
                indicesToRemove.append(curvAbsId)
        decimatedFrameCurvAbs = np.delete(sortedFrameCurvAbs, indicesToRemove)

        # Finally, we complete instrumentIdsForNodeVect and instrumentIdsForFrameVect
        # - Nodes
        for newCurvAbs in decimatedNodeCurvAbs:
            instrumentList = []
            # For the node at newCurvAbs, we test each instrument
            for instrumentId in range(0, self.nbInstruments):
                tipCurvAbs = self.tipCurvAbsVect[instrumentId] # In combined instrument
                beginNodeCurvAbs = xBeginVect[instrumentId] # In combined instrument
                if (beginNodeCurvAbs < newCurvAbs + self.curvAbsTolerance) and (tipCurvAbs > newCurvAbs - self.curvAbsTolerance):
                    # Then this instrument passes throught newCurvAbs
                    instrumentList.append(instrumentId)
            # Once all the instruments were tested, we update instrumentIdsForNodeVect
            instrumentIdsForNodeVect.append(instrumentList)

        # - Frames
        for newCurvAbs in decimatedFrameCurvAbs:
            instrumentList = []
            # For the frame at newCurvAbs, we test each instrument
            for instrumentId in range(0, self.nbInstruments):
                tipCurvAbs = self.tipCurvAbsVect[instrumentId] # In combined instrument
                beginNodeCurvAbs = xBeginVect[instrumentId] # In combined instrument
                if (beginNodeCurvAbs < newCurvAbs + self.curvAbsTolerance) and (tipCurvAbs > newCurvAbs - self.curvAbsTolerance):
                    # Then this instrument passes throught newCurvAbs
                    instrumentList.append(instrumentId)
            # Once all the instruments were tested, we update instrumentIdsForFrameVect
            instrumentIdsForFrameVect.append(instrumentList)

        return {'instrumentIdsForNodeVect': instrumentIdsForNodeVect,
                'decimatedNodeCurvAbs': decimatedNodeCurvAbs,
                'decimatedFrameCurvAbs': decimatedFrameCurvAbs,
                'instrumentIdsForFrameVect': instrumentIdsForFrameVect}


    # This method is meant to apply a new discretisation into Cosserat beams and
    # frames, to a set of instruments. For each instrument, this is done in two
    # steps:
    #  - Updating the beam information, both in the Cosserat mapping component
    # (e.g.: DiscreteCosseratMapping) and the Cosserat forcefield component
    # (e.g.: BeamHookeLawForceField). In the mapping, we change the
    # *curv_abs_input* data field according to the new curvilinear abscissas
    # computed in step 2 (decimatedNodeCurvAbs). We start with the 'last' beam
    # (i.e. the one with the higher index) in order to account for undeployed
    # beams at the proximal end of the instruments. In the force field, we
    # change the *length* data field accordingly.
    #  - Updating the frames information. Based on the new discretisation,
    # new Cosserat frame curvilinear abscissas were also computed in step 2
    # (decimatedFrameCurvAbs). We apply these in the *curv_abs_output* data
    # field of the Cosserat mapping component. Nothing more is to be done
    # regarding the frames, as the associted mechanical object is automatically
    # update by the mapping.
    #
    # Parameters:
    #    - decimatedNodeCurvAbs: [input] set of curvilinear abscissas associated
    # to the Cosserat beam elements (i.e. inputs of the Cosserat mapping)
    #    - decimatedFrameCurvAbs: [input] set of curvilinear abscissas associated
    # to the Cosserat rigid frames (i.e. outputs of the Cosserat mapping)
    #    - instrumentIdsForNodeVect: [input] vector containing for each node
    # the list of instruments which the node belongs to
    #    - instrumentIdsForFrameVect: [input] vector containing for each Cosserat
    # rigid frame the list of instruments which the frame belongs to
    def updateInstrumentComponents(self, decimatedNodeCurvAbs, decimatedFrameCurvAbs, instrumentIdsForNodeVect, instrumentIdsForFrameVect):

        # 'Global' variables for this scope, filled while iterating over the instruments
        nbInstruments = self.nbInstruments
        nbNewNodes = len(decimatedNodeCurvAbs)

        # Precomputation, to analyse the deplyment configuration of the different
        # instruments. The purpose of this precomputation is to fill the
        # instrumentLastNodeIds list, defined below, which contains for each
        # instrument the index of the last beam of the instrument. If the instrument
        # is not deployed yet, the corresponding beam index is 0.
        # The idea to fill the list is the following : we go through the elements
        # of instrumentIdsForNodeVect, and stops whenever the size of the list
        # of instruments passing through the nodes decreases. By turning the lists
        # into sets, we retrieve the indices of the instruments which are no longer
        # in the list. For each of these instruments, we store the corresponding
        # last beam index in instrumentLastNodeIds.
        instrumentLastNodeIds = [0]*nbInstruments
        accumulatedNodeNumber = 0

        if len(instrumentIdsForNodeVect) >= 2: # Requires at least one beam
            # Keeping track of the number of instruments whose distal end we haven't
            # reached yet.
            instrumentIterator = nbInstruments

            while (instrumentIterator > 1 and accumulatedNodeNumber < nbNewNodes-1):
                while (accumulatedNodeNumber < nbNewNodes-1 and len(instrumentIdsForNodeVect[accumulatedNodeNumber+1]) >= instrumentIterator):
                    accumulatedNodeNumber += 1

                if (accumulatedNodeNumber == nbNewNodes-1):
                    # In this case, more than one instrument are ending on the last new node
                    # NB: instrumentIterator can't be equal to 1 here
                    break

                # Retrieving the index of the instruments which distal end we reached
                previousInstrumentList = instrumentIdsForNodeVect[accumulatedNodeNumber]
                lessInstrumentList = instrumentIdsForNodeVect[accumulatedNodeNumber+1]
                instrumentDifferenceSet = set(previousInstrumentList) - set(lessInstrumentList)

                for instrumentId in list(instrumentDifferenceSet):
                    instrumentLastNodeIds[instrumentId] = accumulatedNodeNumber

                # Accordingly decreasing the instrumentIterator
                nbStoppedInstrument = len(instrumentDifferenceSet)
                instrumentIterator -= nbStoppedInstrument

            # When leaving the two loops above, two scenarios are possible :
            #   - either only one instrument remains, meaning that we reached
            # the node from which only this instrument remains
            #   - or more than one instrument remain, meaning that these instruments
            # are coaxial until the last node
            # In both cases, we have to fill instrumentLastNodeIds for all instruments
            # remaining on the last node
            lastNodeInstrumentList = instrumentIdsForNodeVect[nbNewNodes-1]
            for instrumentId in lastNodeInstrumentList:
                instrumentLastNodeIds[instrumentId] = nbNewNodes-1

        # Updating the beam and Cosserat mapping components

        for instrumentId in range(0, nbInstruments):

            # First, retrieving the nodes of the current instrument

            instrumentNodeName = "Instrument" + str(instrumentId)
            instrumentNode = self.solverNode.getChild(str(instrumentNodeName))
            if instrumentNode is None:
                raise NameError("[CombinedInstrumentsController]: Node \'{}\' not found. Your scene should "
                                "contain a node named \'{}\' among the children of the root node in order "
                                "to use this controller".format(instrumentNodeName, instrumentNodeName))

            cosseratMechanicalNode = instrumentNode.getChild('rateAngularDeform')
            if cosseratMechanicalNode is None:
                raise NameError("[CombinedInstrumentsController]: Node \'rateAngularDeform\' "
                                "not found. Your scene should contain a node named "
                                "\'rateAngularDeform\' among the children of the \'{}\' "
                                "node, gathering the mechanical Cosserat components "
                                "(MechanicalObject, Cosserat forcefield)".format(instrumentNodeName))

            rigidBaseNode = instrumentNode.getChild('rigidBase')
            if rigidBaseNode is None:
                raise NameError("[CombinedInstrumentsController]: Node \'rigidBase\' "
                                "not found. Your scene should contain a node named "
                                "\'rigidBase\' among the children of the \'{}\' "
                                "node, where the base and rigid frames of the "
                                "Cosserat model are defined".format(instrumentNodeName))

            mappedFramesNode = rigidBaseNode.getChild('MappedFrames')
            if mappedFramesNode is None:
                raise NameError("[CombinedInstrumentsController]: Node \'MappedFrames\' "
                                "not found. The \'rigidBase\' node should have a child "
                                "node called \'MappedFrames\', in which the Cosserat "
                                "rigid frames and the Cosserat mapping are defined.")

            coaxialFramesNode = rigidBaseNode.getChild('coaxialSegmentFrames')
            if coaxialFramesNode is None:
                raise NameError("[CombinedInstrumentsController]: Node \'coaxialSegmentFrames\' "
                                "not found. The \'rigidBase\' node should have a child "
                                "node called \'coaxialSegmentFrames\', in which the Cosserat "
                                "mapping tracking coaxial segments should be defined.")

            # Retrieving the components
            # TO DO : existance check ? What is the appropriate binding method ?
            beamForceFieldComponent = cosseratMechanicalNode.beamForceField
            fixedConstraintOnStock = cosseratMechanicalNode.FixedConstraintOnStock
            instrumentMapping = mappedFramesNode.DiscreteCosseratMapping
            coaxialMapping = coaxialFramesNode.CoaxialCosseratMapping
            ouputFrameMO = mappedFramesNode.FramesMO

            # Updating the beam information (cf comment of function description)

            # TO DO : replace nbForceFieldBeams by nbInputBeamNodes-1 ?
            nbForceFieldBeams = len(beamForceFieldComponent.length)

            # We keep count of the nodes which are not part of the current
            # instrument, in order to make sure that the beams are added
            # starting from the end
            nbNodesNotOnInstrument = 0

            with beamForceFieldComponent.length.writeable() as forceFieldBeamLengths:
                with instrumentMapping.curv_abs_input.writeable() as curv_abs_input:

                    nbInputBeamNodes = len(curv_abs_input)

                    for curvAbsIterator in range(0, nbNewNodes-1):
                        # TO DO : is it necessary to check that nbNewNodes <= nbInputBeamNodes ?
                        # Probably not, given that the total number of beams is taken into account when
                        # computing the new abscissas
                        # NB: this loop stops one node before the last one, to update
                        # both the curvilinear abscissas in curv_abs_ouput (nbNewNodes)
                        # and the new beam lengths in the force field (nbNewNodes-1).
                        # The last curvilinear abscissa is changed afterwards

                        currentKeyPointCurvAbsId = nbNewNodes-1-curvAbsIterator
                        # Check if the new node belongs to the current instrument,
                        # before applying changes
                        if instrumentId in instrumentIdsForNodeVect[currentKeyPointCurvAbsId]:
                            currentBeamCurvAbsId = nbInputBeamNodes-1-curvAbsIterator+nbNodesNotOnInstrument
                            currentInputBeamId = nbForceFieldBeams-1-curvAbsIterator+nbNodesNotOnInstrument
                            # Modifying curv_abs_input in the Cosserat mapping
                            curv_abs_input[currentBeamCurvAbsId] = decimatedNodeCurvAbs[currentKeyPointCurvAbsId]

                            # Modifying beam lengths in the Cosserat beam ForceField(s)
                            currentBeamLength = decimatedNodeCurvAbs[currentKeyPointCurvAbsId] - decimatedNodeCurvAbs[currentKeyPointCurvAbsId-1]
                            # TO DO : is it necessary to check that nbNewNodes <= nbForceFieldBeams ?
                            forceFieldBeamLengths[currentInputBeamId] = currentBeamLength
                        else:
                            nbNodesNotOnInstrument += 1

                    # Last curv_abs_input, associated with the last beam of the chain
                    curv_abs_input[nbInputBeamNodes-nbNewNodes] = decimatedNodeCurvAbs[0]

                    # Forcing the forceFieldBeamLengths elements which are not
                    # affected to actual beams to 0
                    nbInstrumentNewBeams = nbNewNodes - 1 - nbNodesNotOnInstrument
                    nbUnaffectedBeams = nbForceFieldBeams - nbInstrumentNewBeams
                    forceFieldBeamLengths[0:nbUnaffectedBeams] = [0.] * nbUnaffectedBeams

                    # Same for curv_abs_input
                    # We have to count the first element of curv_abs_input, which is
                    # always 0, as 'unaffected', even though it is technically
                    # part of the newly computed curvilinear abscissas. For this
                    # reason, nbUnaffectedInputCurvAbs is computed using nbInstrumentNewBeams,
                    # instead of nbInstrumentNewNodes (= nbInstrumentNewBeams + 1)
                    nbUnaffectedInputCurvAbs = nbInputBeamNodes - nbInstrumentNewBeams
                    curv_abs_input[0:nbUnaffectedInputCurvAbs] = [0.] * nbUnaffectedInputCurvAbs

                    # Updating the fixed constraint
                    newFixedIndices = list(range(0, nbUnaffectedBeams))
                    fixedConstraintOnStock.indices = newFixedIndices

                    # Updating the second (coaxial) Cosserat mapping input
                    coaxialMapping.curv_abs_input = curv_abs_input

            # Updating the frame information (cf comment of function description)

            nbNewFrames = len(decimatedFrameCurvAbs)
            nbTotalFrames = len(instrumentMapping.curv_abs_output)

            # We keep count of the frames which are not on the current instrument
            # in order to make sure that the new frames are added starting from
            # the end
            nbFramesNotOnInstrument = 0

            with instrumentMapping.curv_abs_output.writeable() as curv_abs_output:

                for frameIterator in range(0, nbNewFrames):
                    # TO DO : is it necessary to check that nbNewFrames <= nbTotalFrames ?
                    # Probably not, given that the total number of frames is taken into account when
                    # computing the new abscissas
                    currentFrameCurvAbsId = nbNewFrames-1-frameIterator
                    currentOutputFrameId = nbTotalFrames-1-frameIterator + nbFramesNotOnInstrument
                    # Check if the new frame belongs to the current instrument,
                    # before applying changes
                    if instrumentId in instrumentIdsForFrameVect[currentFrameCurvAbsId]:
                        curv_abs_output[currentOutputFrameId] = decimatedFrameCurvAbs[currentFrameCurvAbsId]
                    else:
                        nbFramesNotOnInstrument += 1


        # Updating constraints on the coaxial segments

        instrumentIndicesSortedByLength = sorted(range(len(instrumentLastNodeIds)), key=lambda k: instrumentLastNodeIds[k])

        if (nbInstruments > 1 and instrumentLastNodeIds[instrumentIndicesSortedByLength[nbInstruments-2]] > 0):
            # The above condition checks that at least two instruments are deployed.
            # When only one instrument is deployed, all computation on coaxial segments
            # can be skipped.

            # Determining which instruments have not been deployed yet (i.e.: instruments for which
            # instrumentLastNodeIds = 0). These instruments don't have any coaxial beams with
            # other instruments
            shortestDeployedInstrumentRank = 0
            while (instrumentLastNodeIds[instrumentIndicesSortedByLength[shortestDeployedInstrumentRank]] <= 0 and shortestDeployedInstrumentRank < nbInstruments-1):
                shortestDeployedInstrumentRank += 1

            longestDeployedInstrumentRank = nbInstruments-1
            longestDeployedInstrumentId = instrumentIndicesSortedByLength[nbInstruments-1]

            # Loop over all deployed instruments, except the longest
            for instrumentRank in range(shortestDeployedInstrumentRank, longestDeployedInstrumentRank):

                instrumentId = instrumentIndicesSortedByLength[instrumentRank]
                instrumentLastBeamId = instrumentLastNodeIds[instrumentId]

                #--- Retrieving the nodes asscoiated to the current instrument ---#

                instrumentNodeName = "Instrument" + str(instrumentId)
                instrumentNode = self.solverNode.getChild(str(instrumentNodeName))
                if instrumentNode is None:
                    raise NameError("[CombinedInstrumentsController]: Node \'{}\' not found. Your scene should "
                                    "contain a node named \'{}\' among the children of the root node in order "
                                    "to use this controller".format(instrumentNodeName, instrumentNodeName))

                rigidBaseNode = instrumentNode.getChild('rigidBase')
                if rigidBaseNode is None:
                    raise NameError("[CombinedInstrumentsController]: Node \'rigidBase\' "
                                    "not found. Your scene should contain a node named "
                                    "\'rigidBase\' among the children of the \'{}\' "
                                    "node, gathering the Cosserat component RigidBase "
                                    "and frames.".format(instrumentNodeName))

                coaxialFrameNode = rigidBaseNode.getChild('coaxialSegmentFrames')
                if coaxialFrameNode is None:
                    raise NameError("[CombinedInstrumentsController]: Node \'coaxialSegmentFrames\' "
                                    "not found. Your scene should contain a node named "
                                    "\'coaxialSegmentFrames\' among the children of the \'rigidBase\' "
                                    "node, containing the frames used for coaxial constraints.")

                #--- Updating the instrument's coaxial Cosserat mapping (curv_abs_output) ---#

                coaxialBeamCurvAbs = decimatedNodeCurvAbs[0:instrumentLastBeamId+1]
                nbTotalCoaxialFrames = len(coaxialFrameNode.CoaxialCosseratMapping.curv_abs_output)

                # Computing the coaxial frames' curvilinear abscissas
                nbDeployedCoaxialFrames = 0
                coaxialFrameCurvAbs = []
                for beamCurvAbsId in range(len(coaxialBeamCurvAbs)-1):
                    # First we add the proximal end of the beam segment
                    proximalEndCurvAbs = coaxialBeamCurvAbs[beamCurvAbsId]
                    coaxialFrameCurvAbs.append(proximalEndCurvAbs)
                    nbDeployedCoaxialFrames += 1

                    # Then we check wether the beam segment is long enough to add intermediate
                    # coaxial frames. If it is, we compute the additional frames' curvilinear abscissas
                    distalEndCurvAbs = coaxialBeamCurvAbs[beamCurvAbsId+1]
                    intermediateSegmentLength = 0
                    if (self.nbIntermediateConstraintFrames != 0):
                        intermediateSegmentLength = (distalEndCurvAbs - proximalEndCurvAbs) / (self.nbIntermediateConstraintFrames+1)
                    if (intermediateSegmentLength >= self.minimalDistanceForConstraint):
                        for intermediateFrameId in range(self.nbIntermediateConstraintFrames):
                            intermediateCurvAbs = proximalEndCurvAbs + (intermediateFrameId+1) * intermediateSegmentLength
                            coaxialFrameCurvAbs.append(intermediateCurvAbs)
                            nbDeployedCoaxialFrames += 1

                # Finally we add the last beam segment distal end
                coaxialFrameCurvAbs.append(coaxialBeamCurvAbs[len(coaxialBeamCurvAbs)-1])
                nbDeployedCoaxialFrames += 1

                deployedCoaxialFrameIds = list(range(nbTotalCoaxialFrames-nbDeployedCoaxialFrames, nbTotalCoaxialFrames))

                with coaxialFrameNode.CoaxialCosseratMapping.curv_abs_output.writeable() as curv_abs_output:
                    curv_abs_output[nbTotalCoaxialFrames-nbDeployedCoaxialFrames:nbTotalCoaxialFrames] = coaxialFrameCurvAbs
                    curv_abs_output[0:nbTotalCoaxialFrames-nbDeployedCoaxialFrames] = [0.]*(nbTotalCoaxialFrames-nbDeployedCoaxialFrames)

                #--- Second loop over the longer deployed instruments, to enforce constraints on ---#
                #--- the coaxial beam segments                                                   ---#

                for longerInstrumentRank in range(instrumentRank+1, nbInstruments):

                    longerInstrumentId = instrumentIndicesSortedByLength[longerInstrumentRank]
                    longerInstrumentLastBeamId = instrumentLastNodeIds[longerInstrumentId]

                    #--- Retrieving the nodes asscoiated to the second (longer) instrument ---#

                    longerInstrumentNodeName = "Instrument" + str(longerInstrumentId)
                    longerInstrumentNode = self.solverNode.getChild(str(longerInstrumentNodeName))
                    if longerInstrumentNode is None:
                        raise NameError("[CombinedInstrumentsController]: Node \'{}\' not found. Your scene should "
                                        "contain a node named \'{}\' among the children of the root node in order "
                                        "to use this controller".format(longerInstrumentNodeName, longerInstrumentNodeName))

                    rigidBaseNode2 = longerInstrumentNode.getChild('rigidBase')
                    if rigidBaseNode2 is None:
                        raise NameError("[CombinedInstrumentsController]: Node \'rigidBase\' "
                                        "not found. Your scene should contain a node named "
                                        "\'rigidBase\' among the children of the \'{}\' "
                                        "node, gathering the Cosserat component RigidBase "
                                        "and frames.".format(longerInstrumentNodeName))

                    coaxialFrameNode2 = rigidBaseNode2.getChild('coaxialSegmentFrames')
                    if coaxialFrameNode2 is None:
                        raise NameError("[CombinedInstrumentsController]: Node \'coaxialSegmentFrames\' "
                                        "not found. Your scene should contain a node named "
                                        "\'coaxialSegmentFrames\' among the children of the \'rigidBase\' "
                                        "node, containing the frames used for coaxial constraints.")

                    #--- Retrieving the node containing the coupling between the two instruments ---#
                    # NB: this node is both a child of coaxialFrameNode and coaxialFrameNode2

                    if (longerInstrumentId > instrumentId):
                        constraintNodeName = "constraint" + str(longerInstrumentId) + "With" + str(instrumentId)
                    else:
                        constraintNodeName = "constraint" + str(instrumentId) + "With" + str(longerInstrumentId)

                    constraintNode = coaxialFrameNode.getChild(str(constraintNodeName))
                    if constraintNode is None:
                        raise NameError("[CombinedInstrumentsController]: Node \'{}\' "
                                        "not found. The \'coaxialSegmentFrames\' node should have a child "
                                        "node called \'{}\', containing the components "
                                        "which implement its coupling with instrument {}".format(constraintNodeName,
                                        constraintNodeName, longerInstrumentId))

                    # We don't have to update the longer instrument own components : this will be
                    # done when the upper loop reaches it. Here, we just update the coaxial frame
                    # indices in the constraint node corresponding to this pair of instruments.

                    # As the initial instrument is shorter, all its coaxial frames are actually common
                    # with the longer instrument. We can reuse the index range computed above
                    shorterInstrumentCoaxialFrameIds = deployedCoaxialFrameIds

                    # We remove the first index (i.e. the first coaxial beam proximal end frame) if the
                    # corresponding curvilinear abscissa is 0.
                    # NB: this will automatically remove the same first frame for the longer instrument,
                    # because of the computation below
                    if (decimatedNodeCurvAbs[0] < self.curvAbsTolerance):
                        shorterInstrumentCoaxialFrameIds = shorterInstrumentCoaxialFrameIds[1:len(shorterInstrumentCoaxialFrameIds)]

                    # Additionnaly, we remove the first coaxial beam distal end frame, if the
                    # corresponding curvilinear abscissa is too close to the RigidBase. This concretly
                    # means that a coaxial beam segment is only taken into account (i.e. under constraint)
                    # only if its distal end is sufficiently distant from the RigidBase. The purpose of
                    # this check is to avoid overconstraining of a small beam segment.
                    # NB: in this whole section of code, we are under the condition that at least
                    # two instruments are deployed. In this case, there is at least one coaxial beam,
                    # so it is unnecessary to check len(decimatedNodeCurvAbs) before accessing
                    # decimatedNodeCurvAbs[1] or shorterInstrumentCoaxialFrameIds[self.nbIntermediateConstraintFrames+1]
                    if (decimatedNodeCurvAbs[1] < self.minimalDistanceForConstraint):
                        shorterInstrumentCoaxialFrameIds = shorterInstrumentCoaxialFrameIds[self.nbIntermediateConstraintFrames+1:len(shorterInstrumentCoaxialFrameIds)]

                    # For the longer instrument, we have to take into account the additional *coaxial* beams
                    # which are further than the first (shorter) instrument end. The notion of coaxial is
                    # important, because noncoaxial beams won't make a difference in terms of coaxial frame
                    # indices. We therefore have to distinguish two case : if the longer instrument is the
                    # longest instrument, or if it is not.
                    nbAdditionalCoaxialBeams = 0
                    if (longerInstrumentRank == nbInstruments-1):
                        # In this case, the coaxial frame indices are the same as the second longest
                        # deployed instrument.
                        secondLongestInstrumentId = instrumentIndicesSortedByLength[nbInstruments-2]
                        secondLongestInstrumentLastBeamId = instrumentLastNodeIds[secondLongestInstrumentId]
                        nbAdditionalCoaxialBeams = secondLongestInstrumentLastBeamId - instrumentLastBeamId
                    else:
                        # In this case, the coaxial frames of the longer instruments are coincidant
                        nbAdditionalCoaxialBeams = longerInstrumentLastBeamId - instrumentLastBeamId

                    nbTotalCoaxialFramesOfLongerInst = len(coaxialFrameNode2.CoaxialCosseratMapping.curv_abs_output)
                    lastCoaxialFrameIndexWithShorterInstrument = nbTotalCoaxialFramesOfLongerInst - 1 - nbAdditionalCoaxialBeams
                    firstCoaxialFrameIndexWithShorterInstrument = lastCoaxialFrameIndexWithShorterInstrument - len(shorterInstrumentCoaxialFrameIds) + 1
                    longerInstrumentCoaxialFrameIds = list(range(firstCoaxialFrameIndexWithShorterInstrument, lastCoaxialFrameIndexWithShorterInstrument+1))


                    # Once the two sets of indices are computed, we update the RigidDistanceMapping
                    # in the constraint Node accordingly
                    # /!\ We make the assumption that input1 of the mapping is the instrument with
                    # the higher index, so that input2 is the instrument with the lower index. For
                    # instance, if we are considering the constraints between Instrument0 and Instrument1,
                    # then input1 of the RigidDistanceMapping refers to Instrument1, and input2 refers to
                    # Instrument0. This is a convention which have to be followed when describing the
                    # scene structure. It is obviously not dependent on the instruments lengths a time t.
                    # TO DO: Enforce this in a more robust way ?
                    if (longerInstrumentId > instrumentId):
                        constraintNode.coaxialFramesDistanceMapping.first_point = longerInstrumentCoaxialFrameIds
                        constraintNode.coaxialFramesDistanceMapping.second_point = shorterInstrumentCoaxialFrameIds
                    else:
                        constraintNode.coaxialFramesDistanceMapping.first_point = shorterInstrumentCoaxialFrameIds
                        constraintNode.coaxialFramesDistanceMapping.second_point = longerInstrumentCoaxialFrameIds

                    # /!\ TO DO: this call to reinit() is a quick fix to make sure to trigger the update process
                    # of the RigidDistanceMapping component. The reinit() method called here is an empty method
                    # simply checking the componentState in order to trigger a callback. A proper way to
                    # implement this would be to write a binding for the callback method itself.
                    # If the callback is not specifically triggered at this point of the Python execution,
                    # when using a FreeMotionAnimationLoop, (RestShapeSprings)ForceFields applied on the output
                    # of the mapping are used (addForce) before the mapping callback was called. Therefore
                    # leading to a mismatch between the indices in the ForceField (updated directly from the
                    # Python code), and the size of the corresponding MechanicalObject (normally updated by the
                    # callback)
                    constraintNode.coaxialFramesDistanceMapping.reinit()

                    # If the constraint is actually applied through springs (RestShapeSpringsForceField),
                    # We update the corresponding indices in the component
                    if self.constrainWithSprings:
                        nbNewDiffDoFs = len(longerInstrumentCoaxialFrameIds)
                        constraintNode.constraintMappingSpring.points = list(range(nbNewDiffDoFs))

                        # Artificial callback triggering for the MechanicalMatrixMapper components
                        self.solverNode.MMMapperRateAngular0RigidBase0.callbackForMSLinkChanges = not (self.solverNode.MMMapperRateAngular0RigidBase0.callbackForMSLinkChanges)
                        self.solverNode.MMMapperRateAngular0RigidBase1.callbackForMSLinkChanges = not (self.solverNode.MMMapperRateAngular0RigidBase1.callbackForMSLinkChanges)
                        self.solverNode.MMMapperRateAngular0RateAngular1.callbackForMSLinkChanges = not (self.solverNode.MMMapperRateAngular0RateAngular1.callbackForMSLinkChanges)
                        self.solverNode.MMMapperRateAngular1RigidBase0.callbackForMSLinkChanges = not (self.solverNode.MMMapperRateAngular1RigidBase0.callbackForMSLinkChanges)
                        self.solverNode.MMMapperRateAngular1RigidBase1.callbackForMSLinkChanges = not (self.solverNode.MMMapperRateAngular1RigidBase1.callbackForMSLinkChanges)
                        self.solverNode.MMMapperRigidBase0RigidBase1.callbackForMSLinkChanges = not (self.solverNode.MMMapperRigidBase0RigidBase1.callbackForMSLinkChanges)



                # Once all other instruments have been handled, we update the longest deployed
                # instrument. For this instrument, we only have to update the coaxial frame
                # curvilinear abscissas (in the coaxial Cosserat Mapping). Constraints with
                # shorter instruments sharing coaxial beams have already been updated when
                # handling those instruments.

                instrumentNodeName = "Instrument" + str(longestDeployedInstrumentId)
                instrumentNode = self.solverNode.getChild(str(instrumentNodeName))
                if instrumentNode is None:
                    raise NameError("[CombinedInstrumentsController]: Node \'{}\' not found. Your scene should "
                                    "contain a node named \'{}\' among the children of the root node in order "
                                    "to use this controller".format(instrumentNodeName, instrumentNodeName))

                rigidBaseNode = instrumentNode.getChild('rigidBase')
                if rigidBaseNode is None:
                    raise NameError("[CombinedInstrumentsController]: Node \'rigidBase\' "
                                    "not found. Your scene should contain a node named "
                                    "\'rigidBase\' among the children of the \'{}\' "
                                    "node, gathering the Cosserat component RigidBase "
                                    "and frames.".format(instrumentNodeName))

                coaxialFrameNode = rigidBaseNode.getChild('coaxialSegmentFrames')
                if coaxialFrameNode is None:
                    raise NameError("[CombinedInstrumentsController]: Node \'coaxialSegmentFrames\' "
                                    "not found. Your scene should contain a node named "
                                    "\'coaxialSegmentFrames\' among the children of the \'rigidBase\' "
                                    "node, containing the frames used for coaxial constraints.")

                #--- Updating the instrument's coaxial Cosserat mapping (curv_abs_output) ---#

                # For the longest instrument, the coaxial beams are all the beams of the second
                # longest instrument (as the descretisation is common to all instruments on coaxial
                # segments).
                secondLongestInstrumentId = instrumentIndicesSortedByLength[nbInstruments-2]
                secondLongestInstrumentLastBeamId = instrumentLastNodeIds[secondLongestInstrumentId]
                coaxialBeamCurvAbs = decimatedNodeCurvAbs[0:secondLongestInstrumentLastBeamId+1]
                nbTotalCoaxialFrames = len(coaxialFrameNode.CoaxialCosseratMapping.curv_abs_output)

                # Computing the coaxial frames' curvilinear abscissas
                coaxialFrameCurvAbs = []
                nbDeployedCoaxialFrames = 0
                for beamCurvAbsId in range(len(coaxialBeamCurvAbs)-1):
                    # First we add the proximal end of the beam segment
                    proximalEndCurvAbs = coaxialBeamCurvAbs[beamCurvAbsId]
                    coaxialFrameCurvAbs.append(proximalEndCurvAbs)
                    nbDeployedCoaxialFrames += 1

                    # Then we check wether the beam segment is long enough to add intermediate
                    # coaxial frames. If it is, we compute the additional frames' curvilinear abscissas
                    distalEndCurvAbs = coaxialBeamCurvAbs[beamCurvAbsId+1]
                    intermediateSegmentLength = 0
                    if (self.nbIntermediateConstraintFrames != 0):
                        intermediateSegmentLength = (distalEndCurvAbs - proximalEndCurvAbs) / (self.nbIntermediateConstraintFrames+1)
                    if (intermediateSegmentLength >= self.minimalDistanceForConstraint):
                        for intermediateFrameId in range(self.nbIntermediateConstraintFrames):
                            intermediateCurvAbs = proximalEndCurvAbs + (intermediateFrameId+1) * intermediateSegmentLength
                            coaxialFrameCurvAbs.append(intermediateCurvAbs)
                            nbDeployedCoaxialFrames += 1

                # Finally we add the last beam segment distal end
                coaxialFrameCurvAbs.append(coaxialBeamCurvAbs[len(coaxialBeamCurvAbs)-1])
                nbDeployedCoaxialFrames += 1

                # Updating the Cosserat mapping component
                deployedCoaxialFrameIds = list(range(nbTotalCoaxialFrames-nbDeployedCoaxialFrames, nbTotalCoaxialFrames))
                with coaxialFrameNode.CoaxialCosseratMapping.curv_abs_output.writeable() as curv_abs_output:
                    curv_abs_output[nbTotalCoaxialFrames-nbDeployedCoaxialFrames:nbTotalCoaxialFrames] = coaxialFrameCurvAbs
                    curv_abs_output[0:nbTotalCoaxialFrames-nbDeployedCoaxialFrames] = [0.]*(nbTotalCoaxialFrames-nbDeployedCoaxialFrames)

        return {'instrumentLastNodeIds': instrumentLastNodeIds}



    def interpolateMechanicalQuantities(self, decimatedNodeCurvAbs, instrumentLastNodeIds):

        # 'Global' variables for this scope, filled while iterating over the instruments
        nbInstruments = self.nbInstruments
        nbNewNodes = len(decimatedNodeCurvAbs)

        for instrumentId in range(nbInstruments):

            instrument = self.instrumentList[instrumentId]
            instrumentTotalLength = instrument.getTotalLength()

            instrumentLastNodeId = instrumentLastNodeIds[instrumentId]
            nbDeployedBeams = instrumentLastNodeId
            instrumentDeployedLength = decimatedNodeCurvAbs[instrumentLastNodeId]
            instrumentUndeployedLength = instrumentTotalLength - instrumentDeployedLength

            # Retrieving the instrument node components
            instrumentNodeName = "Instrument" + str(instrumentId)
            instrumentNode = self.solverNode.getChild(str(instrumentNodeName))
            if instrumentNode is None:
                raise NameError("[CombinedInstrumentsController]: Node \'{}\' not found. Your scene should "
                                "contain a node named \'{}\' among the children of the root node in order "
                                "to use this controller".format(instrumentNodeName, instrumentNodeName))

            cosseratMechanicalNode = instrumentNode.getChild('rateAngularDeform')
            if cosseratMechanicalNode is None:
                raise NameError("[CombinedInstrumentsController]: Node \'rateAngularDeform\' "
                                "not found. Your scene should contain a node named "
                                "\'rateAngularDeform\' among the children of the \'{}\' "
                                "node, gathering the mechanical Cosserat components "
                                "(MechanicalObject, Cosserat forcefield)".format(instrumentNodeName))

            mappedFramesNode = cosseratMechanicalNode.getChild('MappedFrames')
            if mappedFramesNode is None:
                raise NameError("[CombinedInstrumentsController]: Node \'MappedFrames\' "
                                "not found. The \'rateAngularDeform\' node should have a child "
                                "node called \'MappedFrames\', in which the Cosserat "
                                "rigid frames and the Cosserat mapping are defined.")

            # NB: we retrieve a copy of the instrument Cosserat MO fields, as we are
            # about to update them.
            # previousStepMORestShape = cosseratMechanicalNode.rateAngularDeformMO.rest_position.value.copy()
            previousStepMOPosition = cosseratMechanicalNode.rateAngularDeformMO.position.value.copy()
            previousStepMOVelocity = cosseratMechanicalNode.rateAngularDeformMO.velocity.value.copy()
            nbBeamInAngularRateMO = len(previousStepMOPosition)

            # Retrieving information on the previous discretisation
            prevInstrumentCurvAbs = instrument.getPrevDiscretisation()
            nbDeployedBeamsInPrevDiscretisation = len(prevInstrumentCurvAbs)-1
            instrumentLengthDiff = decimatedNodeCurvAbs[instrumentLastNodeId] - prevInstrumentCurvAbs[len(prevInstrumentCurvAbs)-1]

            for newBeamId in range(nbDeployedBeams):

                # Getting the beam extremities's curvilinear abscissas (expressed
                # in the global instrument configuration, in which 0 corresponds
                # to the 'insertion' point)
                beamBeginCurvAbs = decimatedNodeCurvAbs[newBeamId]
                beamEndCurvAbs = decimatedNodeCurvAbs[newBeamId+1]

                # Expressing the curvilinear abscissas in the instrument local
                # configuration (in which 0 corresponds to the proximal end)
                beamBeginCurvAbsAtRest = beamBeginCurvAbs + instrumentUndeployedLength
                beamEndCurvAbsAtRest = beamEndCurvAbs + instrumentUndeployedLength

                # We use a safety margin to make sure that one of the curvilinear abscissas is not
                # exactly on a key point (i.e. a curvilinear abscissa value for which the rest strain
                # of the instrument changes). This safety margin has to be significantly lower than both the
                # navigation increment distance and the difference between both curvilinear abscissas.
                # We arbitrarily divide the quantities by 10 for comparison.
                # TO DO: use a more rigourous definition ?
                curvAbsSafetyMargin = min(abs(beamBeginCurvAbsAtRest - beamEndCurvAbsAtRest)*1e-1, self.incrementDistance*1e-1)


                ### Interpolating the rest positions ###

                restStrain = instrument.getRestStrain(beamBeginCurvAbsAtRest + curvAbsSafetyMargin,
                                                      beamEndCurvAbsAtRest - curvAbsSafetyMargin)

                # Updating the rest strain in the Cosserat MechanicalObject

                beamIdInMechanicalObject = nbBeamInAngularRateMO - nbDeployedBeams + newBeamId
                with cosseratMechanicalNode.rateAngularDeformMO.rest_position.writeable() as rest_position:
                    rest_position[beamIdInMechanicalObject] = restStrain

                ### Interpolating the position ###

                # NB: we compute the interpolation of position and velocity only if
                # at least one beam of the current instrument was deployed in the
                # last time step
                if (nbDeployedBeamsInPrevDiscretisation > 0):

                    beginCurvAbsBeamIdInPrevDiscretisation = instrument.getBeamIdInPrevDiscretisation(beamBeginCurvAbs+curvAbsSafetyMargin, instrumentLengthDiff)
                    endCurvAbsBeamIdInPrevDiscretisation = instrument.getBeamIdInPrevDiscretisation(beamEndCurvAbs-curvAbsSafetyMargin, instrumentLengthDiff)

                    previousBeamPos = np.array([0., 0., 0.])
                    previousBeamVel = np.array([0., 0., 0.])

                    if (beginCurvAbsBeamIdInPrevDiscretisation == endCurvAbsBeamIdInPrevDiscretisation):
                        # In this case, the new beam was entirely on a unique beam in the last time step.
                        # No interpolation is required: we can directly assign the previous beam's position
                        # and velocity to the new beam
                        beginBeamIdInPreviousMechanicalObject = nbBeamInAngularRateMO - nbDeployedBeamsInPrevDiscretisation + beginCurvAbsBeamIdInPrevDiscretisation
                        previousBeamPos = previousStepMOPosition[beginBeamIdInPreviousMechanicalObject]
                        previousBeamVel = previousStepMOVelocity[beginBeamIdInPreviousMechanicalObject]
                    else:
                        # In this case, the new beam is crossing more than one beam in the last time step.
                        # We have to interpolate the position and velocity of these beams to compute the
                        # position and velocity of the new beam.
                        # As the Cosserat DoFs (angularRate) are expressed in a 'strain' space, we use
                        # a simple linear interpolation.
                        crossedBeamLengths = instrument.getInterpolationBeamLengths(beamBeginCurvAbs + curvAbsSafetyMargin,
                                                                                    beamEndCurvAbs - curvAbsSafetyMargin,
                                                                                    instrumentLengthDiff)

                        totBeamLength = sum(crossedBeamLengths)
                        beginBeamIdInPreviousMechanicalObject = nbBeamInAngularRateMO - nbDeployedBeamsInPrevDiscretisation + beginCurvAbsBeamIdInPrevDiscretisation
                        intermediateBeamId = beginBeamIdInPreviousMechanicalObject
                        for beamLength in crossedBeamLengths:
                            interpolationCoeff = (beamLength / totBeamLength)
                            previousBeamPos += interpolationCoeff * np.array(previousStepMOPosition[intermediateBeamId])
                            previousBeamVel += interpolationCoeff * np.array(previousStepMOVelocity[intermediateBeamId])
                            intermediateBeamId += 1

                    with cosseratMechanicalNode.rateAngularDeformMO.position.writeable() as position:
                        position[beamIdInMechanicalObject] = previousBeamPos
                    with cosseratMechanicalNode.rateAngularDeformMO.velocity.writeable() as velocity:
                        velocity[beamIdInMechanicalObject] = previousBeamVel

            ### Setting the instrument reference discretisation for next timestep ###
            instrumentNewCurvAbs = decimatedNodeCurvAbs[0:instrumentLastNodeId+1]
            instrument.setPrevDiscretisation(instrumentNewCurvAbs)





    # -------------------------------------------------------------------- #
    # -----                   Interaction methods                    ----- #
    # -------------------------------------------------------------------- #

    def onKeypressedEvent(self, event):

        # Record of the user inputs in a new script file
        if event['key'] == 'E':
            if self.exportNavigationScript:
                print("Stops writing actions in {}".format(self.outputFilename))
                self.file.close()
            else:
                print("Starts writing actions in script file {}".format(self.outputFilename))
                self.file=open(self.outputFilename, "w")
            self.exportNavigationScript = not self.exportNavigationScript


        if event['key'] == Key.uparrow:  # Up arrow
            self.moveForward(self.incrementDistance, self.incrementDirection)
            if (self.exportNavigationScript):
                self.file.writelines("{:.2f} upArrow\n".format(self.totalTime))

        if event['key'] == Key.downarrow:  # Down arrow
            self.moveBackward(self.incrementDistance, self.incrementDirection)
            if (self.exportNavigationScript):
                self.file.writelines("{:.2f} downArrow\n".format(self.totalTime))

        if event['key'] == Key.leftarrow:  # Left arrow
            self.rotateCounterclockwise()
            if (self.exportNavigationScript):
                self.file.writelines("{:.2f} leftArrow\n".format(self.totalTime))

        if event['key'] == Key.rightarrow:  # Right arrow
            self.rotateClockwise()
            if (self.exportNavigationScript):
                self.file.writelines("{:.2f} rightArrow\n".format(self.totalTime))

        if event['key'] == '0':
            self.switchInstrument(0)
            if (self.exportNavigationScript):
                self.file.writelines("{:.2f} switchInstrument0\n".format(self.totalTime))

        if event['key'] == '1':
            self.switchInstrument(1)
            if (self.exportNavigationScript):
                self.file.writelines("{:.2f} switchInstrument1\n".format(self.totalTime))

        if event['key'] == '2':
            self.switchInstrument(2)
            if (self.exportNavigationScript):
                self.file.writelines("{:.2f} switchInstrument2\n".format(self.totalTime))

    def moveForward(self, distanceIncrement, direction):
        self.tipCurvAbsVect[self.currentInstrumentId] += distanceIncrement
        # TO DO : check that the RigidBaseMO component exists
        with self.currentInstrumentRigidBaseNode.RigidBaseMO.position.writeable() as rigidBasePos:
            rigidBasePos[self.posDoFsIdGrid] -= direction * distanceIncrement
        # TO DO: check that the tip isn't too far here ?

    def moveBackward(self, distanceIncrement, direction):
        self.tipCurvAbsVect[self.currentInstrumentId] -= distanceIncrement
        # TO DO : check that the RigidBaseMO component exists
        with self.currentInstrumentRigidBaseNode.RigidBaseMO.position.writeable() as rigidBasePos:
            rigidBasePos[self.posDoFsIdGrid] += direction * distanceIncrement
        # TO DO: check that the tip isn't <0 here ?

    def rotateClockwise(self):
        # We apply a rotation of self.incrementAngle degrees around the insertion direction.
        # We have to convert the quaternion part of the control point DoFs to a pyquaternion.Quaternion object,
        # as the order of the quaternion coordinates in pyquaternion (w, x, y, z) is not the same as the Rigid3d
        # objects in Sofa (x, y, z, w)
        with self.currentInstrumentControlPointNode.controlPointMO.position.writeable() as controlPointPos:
            controlPointSofaQuat = controlPointPos[self.quatDoFsIdGrid]  # np matrix of size 1x4
            qx = controlPointSofaQuat[0][0]
            qy = controlPointSofaQuat[0][1]
            qz = controlPointSofaQuat[0][2]
            qw = controlPointSofaQuat[0][3]
            controlPointQuat = Quat(qw, qx, qy, qz)
            newControlPointQuat = self.minusQuat * controlPointQuat
            controlPointPos[self.quatDoFsIdGrid] = np.array([[newControlPointQuat[1], newControlPointQuat[2],
                                                            newControlPointQuat[3], newControlPointQuat[0]]])

    def rotateCounterclockwise(self):
        # Same as above, but with a quaternion corresponding to the opposite angle rotation
        with self.currentInstrumentControlPointNode.controlPointMO.position.writeable() as controlPointPos:
            controlPointSofaQuat = controlPointPos[self.quatDoFsIdGrid]  # np matrix of size 1x4
            qx = controlPointSofaQuat[0][0]
            qy = controlPointSofaQuat[0][1]
            qz = controlPointSofaQuat[0][2]
            qw = controlPointSofaQuat[0][3]
            controlPointQuat = Quat(qw, qx, qy, qz)
            newControlPointQuat = self.plusQuat * controlPointQuat
            controlPointPos[self.quatDoFsIdGrid] = np.array([[newControlPointQuat[1], newControlPointQuat[2],
                                                              newControlPointQuat[3], newControlPointQuat[0]]])

    def switchInstrument(self, newInstrumentId):
        if self.nbInstruments <= newInstrumentId:
            warnings.warn("Instrument number {} doesn't exist. The list of available "
                          "instruments is : {}. Current instrument remains "
                          "{}".format(newInstrumentId, list(range(self.nbInstruments)),
                          self.currentInstrumentId))
        else:
            self.currentInstrumentId = newInstrumentId
            print("Currently controlled: instrument {}".format(newInstrumentId))
            self.changeRefRigidBase(newInstrumentId)
            self.changeControlPoint(newInstrumentId)


    def changeRefRigidBase(self, newInstrumentId):
        instrumentNodeName = "Instrument" + str(newInstrumentId)
        instrumentNode = self.solverNode.getChild(str(instrumentNodeName))
        if instrumentNode is None:
            raise NameError("[CombinedInstrumentsController]: Node \'{}\' not found. Your scene should "
                            "contain a node named \'{}\' among the children of the root node in order "
                            "to use this controller".format(instrumentNodeName, instrumentNodeName))
        rigidBaseNode = instrumentNode.getChild('rigidBase')
        if rigidBaseNode is None:
            raise NameError("[CombinedInstrumentsController]: Node \'rigidBase\' "
                            "not found. Your scene should contain a node named "
                            "\'rigidBase\' among the children of the \'{}\' "
                            "node, where the base and rigid frames of the "
                            "Cosserat model are defined".format(instrumentNodeName))
        self.currentInstrumentRigidBaseNode = rigidBaseNode


    def changeControlPoint(self, newInstrumentId):
        controlPointNodeName = "controlPointNode" + str(newInstrumentId)
        controlPointNode = self.rootNode.getChild(str(controlPointNodeName))
        if controlPointNode is None:
            raise NameError("[CombinedInstrumentsController]: Node \'{}\' not found. Your scene should "
                            "contain a node named \'{}\' among the children of the root node in order "
                            "to use this controller".format(controlPointNodeName, controlPointNodeName))

        self.currentInstrumentControlPointNode = controlPointNode


    def parseControlScript(self, inputFilename):
        actionTimeList = []
        actionList = []

        # Getting file content
        file = open(inputFilename, mode = 'r')
        lines = file.readlines()
        file.close()

        # Parsing
        lineId = 1
        for line in lines:
            currentLine = line.split(' ')

            # Checking line structure
            if len(currentLine) != 2:
                errorMessage = "[CombinedInstrumentsController]: Error in parsing file \'{}: ".format(inputFilename)
                errorMessage += "line {} should contain only two arguments : with the following pattern ".format(lineId)
                errorMessage += "\'time(s) controlAction\'"
                raise ValueError(errorMessage)

            # Reading first argument (timeStamp)
            try:
                timeStamp = float(currentLine[0])
            except ValueError:
                errorMessage = "[CombinedInstrumentsController]: Error in parsing file \'{}: ".format(inputFilename)
                errorMessage += "First argument {} of line {} cannot be converted to float".format(currentLine[0], lineId)
                raise ValueError(errorMessage)
            actionTimeList.append(timeStamp)

            # Reading second argument (action)
            actionStr = currentLine[1].split('\n')[0] # Removing '\n' at the end of currentLine[1]
            if (
                    actionStr == "upArrow" or actionStr == "downArrow" or
                    actionStr == "rightArrow" or actionStr == "leftArrow" or
                    actionStr == "switchInstrument0" or
                    actionStr == "switchInstrument1" or
                    actionStr == "switchInstrument2"
               ):
                actionList.append(actionStr)
            else:
                errorMessage = "[CombinedInstrumentsController]: Error in parsing file \'{}: ".format(inputFilename)
                errorMessage += "control action {} in line {} is not recognised".format(actionStr, lineId)
                raise ValueError(errorMessage)

            # Increment line number
            lineId += 1

        return {'actionTimeList': actionTimeList, 'actionList': actionList}
