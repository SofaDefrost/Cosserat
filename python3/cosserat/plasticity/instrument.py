# -*- coding: utf-8 -*-
"""
Created on December 12 2022

@author: ckrewcun
"""

import Sofa
import os
import numpy as np

class Instrument(object):

    def __init__(self, instrumentNode, totalLength, keyPoints,
                 nbBeamDistribution, restStrain, curvAbsTolerance,
                 *args, **kwargs):
        self.instrumentNode = instrumentNode
        self.totalLength = totalLength
        self.keyPoints = keyPoints
        # Check on nbBeamDistribution
        if len(nbBeamDistribution) != len(keyPoints)+1:
            raise ValueError("[Instrument (class)]: initialisation of object Instrument "
                             "with noncoherent \'keyPoints\' and \'nbBeamDistribution\' parameters. "
                             "The size of \'nbBeamDistribution\' should be one more than the size "
                             "of \'keyPoints\' (as we assume two default keyPoints at the "
                             "instrument extremities).")
        self.nbBeamDistribution = nbBeamDistribution
        self.nbBeams = sum(nbBeamDistribution)
        # Check on restStrain
        if len(restStrain) != len(keyPoints)+1:
            raise ValueError("[Instrument (class)]: initialisation of object Instrument "
                             "with noncoherent \'keyPoints\' and \'restStrain\' parameters. "
                             "The size of \'restStrain\' should be one more than the size "
                             "of \'keyPoints\' (as we assume two default keyPoints at the "
                             "instrument extremities).")
        self.restStrain = restStrain
        self.curvAbsTolerance = curvAbsTolerance
        self.prevCurvAbs = [0.]

    def getTotalLength(self):
        return self.totalLength

    def getKeyPoints(self):
        return self.keyPoints

    def getNbBeamDistribution(self):
        return self.nbBeamDistribution

    def getNbBeams(self):
        return self.nbBeams

    # This method returns the rest DoFs corresponding to a pair of curvilinear abscissas,
    # representing a beam.
    #
    # Parameters:
    #    - beginCurvAbs: curvilinear abscissa corresponding to the beginning of the beam
    #    - endCurvAbs: curvilinear abscissa corresponding to the end of the beam
    #
    # Returns:
    #    - restStrain: Vec3 representing the Cosserat DoFs at rest (torsion,
    # bending in Y, bending in Z), for the beam given as input
    def getRestStrain(self, beginCurvAbs, endCurvAbs):

        nbKeyPoints = len(self.keyPoints)

        if beginCurvAbs > endCurvAbs + self.curvAbsTolerance:
            raise ValueError("[Instrument (class)]: call to method getrestStrain with incorrect "
                             "curvilinear abscissa values. The first argument (beginCurvAbs) "
                             "must be lower than (or equal) the second argument (endCurvAbs).")

        if beginCurvAbs + self.curvAbsTolerance < 0 or endCurvAbs > self.totalLength + self.curvAbsTolerance:
            raise ValueError("[Instrument (class)]: call to method getrestStrain with incorrect "
                             "curvilinear abscissa values. The input values must be between "
                             " 0 and the total Length of the instrument.")

        if nbKeyPoints == 0:
            # Strain at rest is constant
            return self.restStrain[0].copy()
        else:
            keyPointId = 0
            for keyPointCurvAbs in self.keyPoints:
                if (keyPointCurvAbs > beginCurvAbs):
                    break
                keyPointId += 1

            if keyPointId == nbKeyPoints:
                # We didn't encounter a key point curvilinear abscissa higher than
                # the beginning curvilinear abscissa. Consequently, the beam is on
                # the last part of the instrument
                return self.restStrain[len(self.restStrain)-1].copy()
            else:
                # We did encounter a key point
                if endCurvAbs < self.keyPoints[keyPointId]:
                    return self.restStrain[keyPointId].copy()
                else:
                    # We have to interpolate between the current and the next instrument
                    # parts
                    # NB: This case should not happen with regular navigation
                    # TO DO: Implement an interpolation ?
                    return self.restStrain[keyPointId].copy()

    def getBeamIdInPrevDiscretisation(self, newCurvAbs, lengthDiff):

        # print("Previous curv abs : {}".format(self.prevCurvAbs))

        newCurvAbsInPrevConfiguration = newCurvAbs - lengthDiff
        # TO DO: is this correct ? Should we use interpolation with not yet deployed beam instead ?
        if (newCurvAbsInPrevConfiguration < self.curvAbsTolerance):
            newCurvAbsInPrevConfiguration = 0.

        prevDiscrEndBeamCurvAbsId = 0
        for curvAbs in self.prevCurvAbs:
            # print("test {} > {} + {}".format(curvAbs, newCurvAbsInPrevConfiguration, self.curvAbsTolerance))
            if curvAbs > newCurvAbsInPrevConfiguration + self.curvAbsTolerance:
                break
            prevDiscrEndBeamCurvAbsId += 1

        return prevDiscrEndBeamCurvAbsId-1


    def getInterpolationBeamLengths(self, beginCurvAbs, endCurvAbs, lengthDiff):

        if beginCurvAbs > endCurvAbs + self.curvAbsTolerance:
            raise ValueError("[Instrument (class)]: call to method getInterpolationBeamLengths with "
                             "incorrect curvilinear abscissa values. The first argument (beginCurvAbs) "
                             "must be lower than (or equal) the second argument (endCurvAbs).")

        beginCurvAbsInPrevConfiguration = beginCurvAbs - lengthDiff
        endCurvAbsInPrevConfiguration = endCurvAbs - lengthDiff

        # print("beginCurvAbsInPrevConfiguration : {}".format(beginCurvAbsInPrevConfiguration))
        # print("endCurvAbsInPrevConfiguration : {}".format(endCurvAbsInPrevConfiguration))

        # TO DO: is this correct ? Should we use interpolation with not yet deployed beam instead ?
        if (beginCurvAbsInPrevConfiguration < self.curvAbsTolerance):
            newCurvAbsInPrevConfiguration = 0.

        beginBeamId = 0
        endBeamId = 0
        prevCurvAbsId = 0
        while (self.prevCurvAbs[prevCurvAbsId] < beginCurvAbsInPrevConfiguration + self.curvAbsTolerance
               and prevCurvAbsId < len(self.prevCurvAbs)):
            prevCurvAbsId += 1
        beginBeamId = prevCurvAbsId - 1

        # print("beginBeamId : {}".format(beginBeamId))

        while (self.prevCurvAbs[prevCurvAbsId] < endCurvAbsInPrevConfiguration + self.curvAbsTolerance
               and prevCurvAbsId < len(self.prevCurvAbs)):
            prevCurvAbsId += 1
        endBeamId = prevCurvAbsId - 1

        # print("endBeamId : {}".format(endBeamId))

        startBeamLength = self.prevCurvAbs[beginBeamId+1] - beginCurvAbsInPrevConfiguration
        beamLengths = [startBeamLength]
        for intermediateBeamId in range(beginBeamId+1, endBeamId):
            # print("Intermediate beam")
            beamLengths.append(self.prevCurvAbs[intermediateBeamId+1] - self.prevCurvAbs[intermediateBeamId])
        endBeamLength = endCurvAbsInPrevConfiguration - self.prevCurvAbs[endBeamId]
        beamLengths.append(endBeamLength)

        return beamLengths


    def setPrevDiscretisation(self, curvAbs):
        self.prevCurvAbs = curvAbs

    def getPrevDiscretisation(self):
        return self.prevCurvAbs.copy()
