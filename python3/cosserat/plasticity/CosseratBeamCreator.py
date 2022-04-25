# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 2021
@author: ckrewcun
"""

# -*- coding: utf-8 -*-

"""
    Auxiliary methods for Cosserat beam topology generation
"""

def generateRegularSectionsAndFrames(totLength, nbBeams, nbFrames):

    # This method generates the parameters to pass to a Cosserat mode. It basically distributes uniformly a required
    # number of sections (Cosserat beams) and a required number of frames (Rigid3d), over a given length.
    # All coordinated are computed by default along the X axis. From the distribution of sections and frames,
    # different input parameters for the Cosserat components (forcefield, mapping) are returned:
    #   - The length of each section (BeamHookeLawForceField.length)
    #   - The curvilinear coordinate of each section extremities (DiscreteCosseratMapping.curv_abs_input)
    #   - The section DoFs (3 DoFs of strain for torsion and bending)
    #   - The Rigid3d DoFs representing the frames
    #   - The curvilinear coordinate of each frame (DiscreteCosseratMapping.curv_abs_output)
    #   - 3D position of the frames + edge topology connecting them (for collision modelling)

    beamLengths = []
    curvAbsInput = []
    beamStrainDoFs = []

    frameRigidDoFs = []
    frame3DDoFs = []
    frameEdges = []
    curvAbsOutput = []

    # Computing section lengths and sections curvilinear coordinates

    individualBeamLength = totLength / nbBeams
    lengthIncr = 0.0
    curvAbsInput.append(0.0)
    for i in range(0, nbBeams):
        beamLengths.append(individualBeamLength)
        lengthIncr += individualBeamLength
        curvAbsInput.append(lengthIncr)
        beamStrainDoFs.append([0., 0., 0.])
    # NB: there are nbBeams+1 values in curvAbsInput, accounting for all section extremities

    # Computing frames Rigid3d coordinates and curvilinear coordinates

    frameIntervalDim = totLength / (nbFrames-1) # nbFrames-1 because the last frame is at the end of the last frame interval
    lengthIncr = 0.0
    # Adding first frame
    frameRigidDoFs.append([0., 0., 0., 0., 0., 0., 1.])
    frame3DDoFs.append([0., 0., 0.])
    curvAbsOutput.append(0.0)
    # Completing
    for i in range(0, nbFrames-1):
        lengthIncr += frameIntervalDim
        frameRigidDoFs.append([lengthIncr, 0., 0., 0., 0., 0., 1.])
        frame3DDoFs.append([lengthIncr, 0., 0.])
        frameEdges.append(i)
        frameEdges.append(i+1)
        curvAbsOutput.append(lengthIncr)

    return {'sectionLengths': beamLengths, 'curvAbsInput': curvAbsInput, 'sectionDoFs': beamStrainDoFs,
            'frameRigidDoFs': frameRigidDoFs, 'frame3DDoFs': frame3DDoFs, 'frameEdges': frameEdges,
            'curvAbsOutput': curvAbsOutput}


# TO DO : to be verified
def generateRegularBeamsWithSameNbFrames(totLength, nbBeams, nbFramesPerBeam):

    # This method generates the parameters to pass to a Cosserat mode. It basically distributes uniformly a required
    # number of sections (Cosserat beams) and a required number of frames (Rigid3d), over a given length.
    # All coordinated are computed by default along the X axis. From the distribution of sections and frames,
    # different input parameters for the Cosserat components (forcefield, mapping) are returned:
    #   - The length of each section (BeamHookeLawForceField.length)
    #   - The curvilinear coordinate of each section extremities (DiscreteCosseratMapping.curv_abs_input)
    #   - The section DoFs (3 DoFs of strain for torsion and bending)
    #   - The Rigid3d DoFs representing the frames
    #   - The curvilinear coordinate of each frame (DiscreteCosseratMapping.curv_abs_output)
    #   - 3D position of the frames + edge topology connecting them (for collision modelling)

    beamLengths = []
    curvAbsInput = []
    beamStrainDoFs = []

    frame6DDoFs = []
    frame3DDoFs = []
    frameEdgeList = []
    curvAbsOutput = []

    # Computing section lengths and sections curvilinear coordinates

    individualBeamLength = totLength / nbBeams
    lengthIncr = 0.0
    curvAbsInput.append(0.0)
    for i in range(0, nbBeams):
        beamLengths.append(individualBeamLength)
        lengthIncr += individualBeamLength
        curvAbsInput.append(lengthIncr)
        beamStrainDoFs.append([0., 0., 0.])
    # NB: there are nbBeams+1 values in curvAbsInput, accounting for all section extremities

    # Computing frames Rigid3d coordinates and curvilinear coordinates
    frameIntervalDim = individualBeamLength / nbFramesPerBeam

    for beamId in range(nbBeams):
        lengthIncr = 0.0
        beamCurvAbs = curvAbsInput[beamId]
        for frameId in range(nbFramesPerBeam):
            currentFrameCurvAbs = beamCurvAbs+lengthIncr
            frames6DDofs.append([currentFrameCurvAbs, 0., 0., 0., 0., 0., 1.])
            curvOutput.append(currentFrameCurvAbs)
            frames3DDofs.append([currentFrameCurvAbs, 0., 0.])
            lengthIncr+=frameIntervalDim

    # Adding last frame
    curvOutput.append(totLength)
    frames3DDofs.append([totLength, 0., 0.])
    frames6DDofs.append([totLength, 0., 0., 0., 0., 0., 1.])

    for i in range(0, len(frames3DDofs) - 1):
        frameEdgeList.append(i)
        frameEdgeList.append(i + 1)

    return {'beamLengths': beamLengths, 'curvAbsInput': curvAbsInput, 'beamStrainDoFs': beamStrainDoFs,
            'frame6DDoFs': frame6DDoFs, 'frame3DDoFs': frame3DDoFs, 'frameEdgeList': frameEdgeList,
            'curvAbsOutput': curvAbsOutput}
