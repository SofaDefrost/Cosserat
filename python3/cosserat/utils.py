import numpy as np
import params

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
                                   name="constraintPointsMo", position=[], showObjectScale=0)

    # print(f' ====> The beamTip tip is : {dir(beamPath)}')
    constraintPointsNode.addObject('PointsManager', name="pointsManager", listening="1",
                                   beamPath="/solverNode/needle/rigidBase/cosseratInSofaFrameNode/slidingPoint/slidingPointMO")

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


def computeDistanceBetweenPoints(constraintPointsNode, slidingPoint):
    with constraintPointsNode.position.writeable() as posA:
        if len(posA) != 0:
            print(f' ====> The pointA is : {posA}')
            with slidingPoint.position.writeable() as posB:
                print(f' ====> The pointB is : {posB}')
                return np.linalg.norm(posA[-1] - posB[-1])
        else:
            # print("No constraint points yet")
            return 0
