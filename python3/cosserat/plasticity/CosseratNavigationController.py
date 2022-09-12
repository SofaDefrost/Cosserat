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

# Python controller to simulate the navigation of coaxial instruments represented
# by two chains of Cosserat beam elements.
# Required structure of the scene :
# * rootNode
#   * controlPointNode
#       controlPointMO (Rigid)
#   * cosseratBeamNode
#       MechanicalMatrixMapper
#       * rigidBaseNode
#           RigidBaseMO (Rigid)
#           * MappedFrames
#               FramesMO (Rigid)
#               DiscreteCosseratMapping
#           * ControlMappedFrame
#               ControlFrameMO (1 Rigid)
#               DiscreteCosseratMapping
#               RestShapeSpringsForceField (linked to controlPointMO)
#       * rateAngularDeform
#           rateAngularDeformMO (Cosserat strains)
#
# Arguments :
#  - nbInstruments: number of simulated coaxial instruments
#  - instrumentBeamDensityVect : vector containing a density of beam elements
#    specific to each instrument
#  - incrementDistance : distance of pushing/pulling the instrument at user
#    interaction
#  - isInstrumentStraightVect : vector of boolean indicating for each instrument
#    if it is entirely straight (true), or if the distal end is nonstraight (false)
#  - curvAbsTolerance : distance threshold used to determine if two close nodes
#    should be merged (and considered as one)
#  - instrumentLengths : vector of double indicating the length of each instrument
#
# /!\ In this controller, as it is the case in BeamAdapter, we assume that the
# instruments are given in the ordre of decreasing diameter, meaning that the
# outmost instrument should have index 0. This assumption mostly used when
# dynamically recomputing the instruments discretisation: in overlapping segments,
# only the outmost instrument should be visible.
#
class CombinedInstrumentsController(Sofa.Core.Controller):

    def __init__(self, rootNode,
                 nbInstruments,
                 instrumentBeamDensityVect,
                 incrementDistance,
                 isInstrumentStraightVect,
                 curvAbsTolerance,
                 instrumentLengths,
                 *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)

        ### Checking for the scene graph structure ###

        self.rootNode = rootNode

        # Cosserat beam node
        beamNode = rootNode.getChild('Instrument0')
        if beamNode is None:
            raise NameError('[CombinedInstrumentsController]: Node \'Instrument0\' not found. Your scene should '
                            'contain a node named \'Instrument0\' among the children of the root node in order '
                            'to use this controller')
        self.rigidBaseNode = beamNode.getChild('rigidBase')
        self.mappedFramesNode = self.rigidBaseNode.getChild('MappedFrames')

        self.cosseratMechanicalNode = beamNode.getChild('rateAngularDeform')

        ### Reading the insertion velocity parameters ###

        self.incrementDistance = incrementDistance
        self.instrumentBeamDensityVect = instrumentBeamDensityVect
        self.curvAbsTolerance = curvAbsTolerance

        ### Controller settings ###

        # Number of simulated instruments
        self.nbInstruments = nbInstruments

        # Total rest length of each instrument
        self.instrumentLengths = instrumentLengths
        if len(instrumentLengths) != nbInstruments:
            raise ValueError('[CombinedInstrumentsController]: size of argument \'instrumentLengths\' '
                             'should be equal to nbInstruments')

        # Curvilinear Abscissa of the tip of each instrument (modified by up/down arrows)
        self.tipCurvAbsVect = np.zeros(nbInstruments)

        # Rotation angle for each instrument (modified by left/right arrows)
        self.rotationAngleVect = np.zeros(nbInstruments)

        # Index of the currently navigated instrument
        self.currentInstrumentId = 0

        ### Additional settings ###

        self.totalTime = 0.0
        self.dt = self.rootNode.findData('dt').value
        self.nbIterations = 0


    def onAnimateBeginEvent(self, event):  # called at each begin of animation step

        # Define the vector which contains the curvilinear abscissas of all
        # represented nodes.
        # TO DO: definition of the nodes
        simulatedCurvAbs = []

        #----- Step 0 : base check on instrument's lengths -----#

        # Check that all instruments haven't been pushed further than their
        # rest length. Otherwise the tip curvilinear abscissa of such an instrument
        # is set back to the instrument's length

        for instrumentId in range (0, self.nbInstruments):
            if self.tipCurvAbsVect[instrumentId] > self.instrumentLengths[instrumentId]:
                self.tipCurvAbsVect[instrumentId] = self.instrumentLengths[instrumentId]


        #----- Step 1 : instrument's tip curvilinear abscissa and combined length -----#

        # Find:
        #    - which instruments are 'out' and simulated (for which tipCurvAbs > 0)
        #    - tipCurvAbs of the most distal instrument = the length of the combined
        #      instruments

        xBeginVect = []

        combinedInstrumentLength = self.getInstrumentConfiguration(xBeginVect)

        # if the totalLength is 0, we move the first instrument and update the
        # information on the instruments configuration
        if combinedInstrumentLength < 0.0001:
            xBeginVect = []
            self.tipCurvAbsVect[0] = 0.0001;
            combinedInstrumentLength = self.getInstrumentConfiguration(xBeginVect)

        # Adding the base point if at least one instrument is out
        if combinedInstrumentLength > 0.:
            simulatedCurvAbs.append(0.)


        #----- Step 2 : computation of the points of interest -----#

        # Computation of the curvilinear abscissas for the points of interest,
        # which are to be simulated. Additionaly, two structures are completed:
        #    - instrumentIdsForNodeVect: vector containing for each node a pluginNameList
        #      of the indices of all instruments passing through the nodes
        #    - xBeginVect: vector containing for each instrument the curvilinear
        #      abscissa of the proximal end of the instrument. As the curvilinear
        #      abscissa is computed relatively to the point from which the
        #      instruments are pushed, it is negative for the proximal end

        instrumentIdsForNodeVect = []

        result = self.computeNodeCurvAbs(simulatedCurvAbs, xBeginVect,
                                         instrumentIdsForNodeVect)

        instrumentIdsForNodeVect = result['instrumentIdsForNodeVect']
        decimatedCurvAbs = result['decimatedCurvAbs']

        print("New decimated curvilinear abscissas : {}".format(decimatedCurvAbs))
        print("New vector for instrument Ids : {}".format(instrumentIdsForNodeVect))


        #----- Step 3 :  -----#

        # Apply changes on the Cosserat components
        # TO DO: from the new curv abs, mimick step3 from the IRC Controller
        #  - apply changes of CA on the controlled instrument (move the Cosserat base, correct beam length, update frames, and apply force on the base according to the motion)
        #  - /!\ Adapt beam length so that the extremity of the beams match between the different isInstrumentStraightVect
        #      => the beam length is no longer determined by the user, but by the configuration
        #      => Alternative = if the beams (sections) are not coincident, is it possible to constrain them properly ?



        self.totalTime = self.totalTime + self.dt
        self.nbIterations = self.nbIterations + 1


    def getInstrumentConfiguration(self, xBeginVect):

        combinedInstrumentLength = 0.

        # Computation of the curvilinear abscissas of the instrument tip
        # (distal end) and beginning (proximal end).
        # /!\ All curvilinear abscissas are expressed relatively to the base
        # point from which the instruments are pushed/pulled. Therefore, for
        # the tip, the curvilinear abscissa is >0, but for the beginning point
        # it is <0
        for instrumentId in range(0,self.nbInstruments):

            xEnd = self.tipCurvAbsVect[instrumentId]
            xBegin = xEnd - self.instrumentLengths[instrumentId]
            xBeginVect.append(xBegin)

            if xEnd > combinedInstrumentLength:
                combinedInstrumentLength = xEnd

        return combinedInstrumentLength


    # Method concretely handling the dynamic navigation of the instruments, and
    # the underlying necessary remeshing.
    # TO DO: complete documentation
    def computeNodeCurvAbs(self, simulatedCurvAbs, xBeginVect,
                           instrumentIdsForNodeVect):

        # Computation of the key points for each instrument. These points correspond
        # to changes in the instrument shape/configuration. If the instrument is
        # entirely straight, the only two key points are the proximal end (begin
        # point) and the distal end (end point). The key points are stored by
        # means of their curvilinear abscissa. As before, this curvilinear abscissa
        # is expressed relatively to the combined configuration starting point.
        # Points for which the curvilinear abscissa is <0 are ignored.
        instrumentKeyPointsVect = [[]]
        maxCurvilinearAbscissa = 0. # Max curvilinear abscissa among the key points

        for instrumentId in range(0,self.nbInstruments):
            # Add first key point = proximal extremity point
            beginNodeCurvAbs = xBeginVect[instrumentId] + 0.0
            # TO DO: are the two check below relevent for the first key point ?
            if (beginNodeCurvAbs > 0):
                simulatedCurvAbs.append(beginNodeCurvAbs)
                # Update the maxCurvilinearAbscissa if beginNodeCurvAbs is higher
                if (beginNodeCurvAbs > maxCurvilinearAbscissa):
                    maxCurvilinearAbscissa = beginNodeCurvAbs

            # TO DO: If the instrument is not entirely straight, add the
            # intermediary key points. For each corresponding interval,
            # also add the points based on the beam density on the interval.

            # Add second key point = distal extremity point
            endNodeCurvAbs = xBeginVect[instrumentId] + self.instrumentLengths[instrumentId]
            if (endNodeCurvAbs > 0):
                # If the distal end of the interval is visible (curv. abs. > 0),
                # it means that a least a part of the interval is out, so the
                # correpsonding curv. abs. has to be added in simulatedCurvAbs.
                #
                # If additionnaly, this curv. abs. is greater than the current
                # maximum abscissa (maxCurvilinearAbscissa), it means that the
                # current interval is visible. In this case, we have to discretise
                # the visible part into several beam elements, based on the
                # desired beam density.

                # Add the end point of the interval
                simulatedCurvAbs.append(endNodeCurvAbs)

                if (endNodeCurvAbs > maxCurvilinearAbscissa):

                    # Compute the number of new nodes to add
                    intervalLength = endNodeCurvAbs - beginNodeCurvAbs
                    visibleIntervalLength = endNodeCurvAbs - maxCurvilinearAbscissa
                    ratio = visibleIntervalLength / intervalLength
                    nbNewNodes = int(self.instrumentBeamDensityVect[instrumentId] * ratio)

                    # Add the new nodes
                    for newNodeId in range(0, nbNewNodes):
                        newNodeCurvAbs = endNodeCurvAbs - ((newNodeId+1) / ratio)
                        simulatedCurvAbs.append(newNodeCurvAbs)

                    # Update the max curv. abs.
                    maxCurvilinearAbscissa = endNodeCurvAbs
        # endfor instrumentId in range(0,self.nbInstruments)

        # When all points of interest have been detected, we sort and filter
        # ther curv. abs' list:
        # First: sort the curv. abs. values
        sortedCurvAbs = np.sort(simulatedCurvAbs, kind='quicksort') # quicksort, heapsort, mergesort, timsort

        # Second: remove the duplicated values, according to self.curvAbsTolerance
        indicesToRemove = []
        for curvAbsId in range(1, len(sortedCurvAbs)):
            diffWithPrevious = abs(sortedCurvAbs[curvAbsId] - sortedCurvAbs[curvAbsId -1])
            if diffWithPrevious < self.curvAbsTolerance:
                indicesToRemove.append(curvAbsId)

        decimatedCurvAbs = np.delete(sortedCurvAbs, indicesToRemove)

        # Finally, complete instrumentIdsForNodeVect
        for newCurvAbs in decimatedCurvAbs:
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

        return {'instrumentIdsForNodeVect': instrumentIdsForNodeVect, 'decimatedCurvAbs': decimatedCurvAbs}




    def onAnimateEndEvent(self, event):  # called at each end of animation step
        pass

    def onKeypressedEvent(self, event):

        if event['key'] == Key.uparrow:  # Up arrow
            self.moveForward(self.incrementDistance)

        if event['key'] == Key.downarrow:  # Down arrow
            self.moveBackward(self.incrementDistance)

        if event['key'] == '0':
            self.currentInstrumentId = 0
            print("Currently controlled: instrument 0")

        if event['key'] == '1':
            if self.nbInstruments <= 1:
                warnings.warn("Instrument number 1 doesn't exist (only one instrument (0) is available).".format(self.nbInstruments))
            else:
                self.currentInstrumentId = 0
                print("Currently controlled: instrument 1")

        if event['key'] == '2':
            if self.nbInstruments <= 2:
                warnings.warn("Instrument number 1 doesn't exist (avalaible instruments are from 0 to {}).".format(self.nbInstruments))
            else:
                self.currentInstrumentId = 2
                print("Currently controlled: instrument 2")

    def moveForward(self, distanceIncrement):
        self.tipCurvAbsVect[self.currentInstrumentId] += distanceIncrement
        # print("New tip curvilinear abscissa for instrument {}: {}".format(self.currentInstrumentId, self.tipCurvAbsVect[self.currentInstrumentId]))
        # TO DO: check that the tip isn't too far here ?

    def moveBackward(self, distanceIncrement):
        self.tipCurvAbsVect[self.currentInstrumentId] -= distanceIncrement
        # TO DO: check that the tip isn't <0 here ?

    def rotateClockwise():
        pass

    def rotateCounterclockwise():
        pass
