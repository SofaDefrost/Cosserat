# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 2021

@author: ckrewcun
"""

# -*- coding: utf-8 -*-

"""
    Auxiliary methods for Cosserat beam topology generation
"""

def generateRegularSectionsAndFrames(totLength, nbSections, nbFrames):

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

    sectionLengths = []
    curvAbsInput = []
    sectionDoFs = []

    frameRigidDoFs = []
    frame3DDoFs = []
    frameEdges = []
    curvAbsOutput = []

    # Computing section lengths and sections curvilinear coordinates

    sectionDim = totLength / nbSections
    lengthIncr = 0.0
    curvAbsInput.append(0.0)
    for i in range(0, nbSections):
        sectionLengths.append(sectionDim)
        lengthIncr += sectionDim
        curvAbsInput.append(lengthIncr)
        sectionDoFs.append([0., 0., 0.])
    # NB: there are nbSections+1 values in curvAbsInput, accounting for all section extremities

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

    return {'sectionLengths': sectionLengths, 'curvAbsInput': curvAbsInput, 'sectionDoFs': sectionDoFs,
            'frameRigidDoFs': frameRigidDoFs, 'frame3DDoFs': frame3DDoFs, 'frameEdges': frameEdges,
            'curvAbsOutput': curvAbsOutput}
