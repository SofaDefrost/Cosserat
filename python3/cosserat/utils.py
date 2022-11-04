
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
    collisInstrumentCombined.addObject('MechanicalObject', name="collisionStats")
    collisInstrumentCombined.addObject('PointCollisionModel', name="beamColMod", group='2')
    collisInstrumentCombined.addObject('IdentityMapping', name="beamMapping")
    print(f' ====> The beam collision node is : {dir(collisInstrumentCombined)}')
    return collisInstrumentCombined


#  @info: This function is used to build the constraint node
def addConstraintPoint(parentNode):
    constraintPointsNode = parentNode.addChild('constraintPoints')
    constraintPointsNode.addObject("PointSetTopologyContainer", name="constraintPtsContainer",
                                   points=[])
    constraintPointsNode.addObject("PointSetTopologyModifier", name="constraintPtsModifier")
    constraintPointsNode.addObject("MechanicalObject", template="Vec3d", showObject=True, showIndices=True,
                                   name="constraintPointsMo", position=[], showObjectScale=1)
    constraintPointsNode.addObject('BarycentricMapping')
    return constraintPointsNode
