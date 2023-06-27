import numpy as np


def addEdgeCollision(parentNode, position3D, edges):
    collisInstrumentCombined = parentNode.addChild('collisInstrumentCombined')
    collisInstrumentCombined.addObject('EdgeSetTopologyContainer', name="collisEdgeSet", position=position3D,
                                       edges=edges)
    collisInstrumentCombined.addObject('EdgeSetTopologyModifier', name="collisEdgeModifier")
    collisInstrumentCombined.addObject('MechanicalObject', name="CollisionDOFs")
    collisInstrumentCombined.addObject('LineCollisionModel', bothSide="1", group='2')
    collisInstrumentCombined.addObject('PointCollisionModel', group='2')
    collisInstrumentCombined.addObject('IdentityMapping', name="mapping")
    return collisInstrumentCombined


# @info: This function is used to build the beam collision node
def addPointsCollision(parentNode, position3D, edges, nodeName):
    collisInstrumentCombined = parentNode.addChild(nodeName)
    collisInstrumentCombined.addObject('EdgeSetTopologyContainer', name="beamContainer", position=position3D,
                                       edges=edges)
    collisInstrumentCombined.addObject('EdgeSetTopologyModifier', name="beamModifier")
    collisInstrumentCombined.addObject('MechanicalObject', name="collisionStats", showObject=False, showIndices=False)
    collisInstrumentCombined.addObject('PointCollisionModel', name="beamColMod", group='2')
    # collisInstrumentCombined.addObject('IdentityMapping', name="beamMapping")
    collisInstrumentCombined.addObject('RigidMapping', name="beamMapping")
    return collisInstrumentCombined


#  @info: This function is used to build the constraint node
def addConstraintPoint(parentNode, beamPath):
    constraintPointsNode = parentNode.addChild('constraintPoints')
    constraintPointsNode.addObject("PointSetTopologyContainer", name="constraintPtsContainer", listening="1")
    constraintPointsNode.addObject("PointSetTopologyModifier", name="constraintPtsModifier", listening="1")
    constraintPointsNode.addObject("MechanicalObject", template="Vec3d", showObject=True, showIndices=True,
                                   name="constraintPointsMo", position=[], showObjectScale=0, listening="1")

    # print(f' ====> The beamTip tip is : {dir(beamPath)}')
    constraintPointsNode.addObject('PointsManager', name="pointsManager", listening="1",
                                   beamPath="/solverNode/needle/rigidBase/cosseratInSofaFrameNode/slidingPoint"
                                            "/slidingPointMO")

    constraintPointsNode.addObject('BarycentricMapping', useRestPosition="false", listening="1")
    return constraintPointsNode


def addSlidingPoints(parenNode, frames3D):
    slidingPoint = parenNode.addChild('slidingPoint')
    slidingPoint.addObject('MechanicalObject', name="slidingPointMO", position=frames3D,
                           showObject=False, showIndices=False)
    slidingPoint.addObject('IdentityMapping')
    return slidingPoint


def getLastConstraintPoint(constraintPointsNode):
    return constraintPointsNode.getObject('constraintPointsMo').position[-1]


def computeDistanceBetweenPoints(constraintPointPos, slidingPointPos):
    if len(constraintPointPos) != 0:
        return np.linalg.norm(constraintPointPos[-1] - slidingPointPos[-1])
    else:
        # print("No constraint points yet")
        return 0

def computePositiveAlongXDistanceBetweenPoints(constraintPointPos, slidingPointPos):
    if len(constraintPointPos) != 0:
        if constraintPointPos[-1][0] > slidingPointPos[-1][0]:
            return np.linalg.norm(constraintPointPos[-1] - slidingPointPos[-1])
        else:
            return 0
    else:
        # print("No constraint points yet")
        return 0


def computeNegativeAlongXDistanceBetweenPoints(constraintPointPos, slidingPointPos):
    if len(constraintPointPos) != 0:
        if constraintPointPos[-1][0] < slidingPointPos[-1][0]:
            return np.linalg.norm(constraintPointPos[-1] - slidingPointPos[-1])
        else:
            return 0
    else:
        # print("No constraint points yet")
        return 0
