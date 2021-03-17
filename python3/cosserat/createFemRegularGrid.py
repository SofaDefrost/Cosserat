# -*- coding: utf-8 -*-

def createFemCube(parentNode):
    FemNode = parentNode.addChild("FemNode")
    FemNode.addObject('VisualStyle', displayFlags='showBehaviorModels hideCollisionModels hideBoundingCollisionModels '
                                                  'showForceFields hideInteractionForceFields showWireframe')
    gelVolume = FemNode.addChild("gelVolume")
    gelVolume.addObject("RegularGridTopology", name="HexaTop", n="6 6 6", min="40 -6 -10", max="100 30 10")
    gelVolume.addObject("TetrahedronSetTopologyContainer", name="Container", position="@HexaTop.position")
    gelVolume.addObject("TetrahedronSetTopologyModifier", name="Modifier")
    gelVolume.addObject("Hexa2TetraTopologicalMapping", input="@HexaTop", output="@Container", swapping="false")

    GelSurface = FemNode.addChild("GelSurface")
    GelSurface.addObject("TriangleSetTopologyContainer", name="Container", position="@../GelVolume/HexaTop.position")
    # GelSurface.addObject("TriangleSetTopologyModifier", input="@../GelVolume/Container", output="@Container",
    #                      flipNormals="false")

    gelNode = FemNode.addChild("gelNode")
    # gelNode.addObject('VisualStyle', displayFlags='showVisualModels hideBehaviorModels showCollisionModels '
    #                                                  'hideMappings hideForceFields showWireframe '
    #                                                  'showInteractionForceFields hideForceFields')
    gelNode.addObject("EulerImplicitSolver", rayleighMass="0.1", rayleighStiffness="0.1")
    gelNode.addObject('SparseLDLSolver', name='preconditioner')
    gelNode.addObject('TetrahedronSetTopologyContainer', src="@../gelVolume/Container", name='container')
    # gelNode.addObject('TetrahedronSetTopologyModifier')
    gelNode.addObject('MechanicalObject', name='tetras', template='Vec3d')
    gelNode.addObject('TetrahedronFEMForceField', template='Vec3d', name='FEM', method='large',
                      poissonRatio='0.45', youngModulus='500')
    # gelNode.addObject('UniformMass', totalMass='5')
    gelNode.addObject('BoxROI', name='ROI1', box='40 -6 -10 100 -4 10', drawBoxes='true')
    gelNode.addObject('RestShapeSpringsForceField', points='@ROI1.indices', stiffness='1e12')

    surfaceNode = gelNode.addChild("surfaceNode")
    surfaceNode.addObject('TriangleSetTopologyContainer', name="surfContainer", src="@../../GelSurface/Container")
    surfaceNode.addObject('MechanicalObject', name='msSurface')
    surfaceNode.addObject('TriangleCollisionModel', name='surface')
    surfaceNode.addObject('BarycentricMapping')

    gelNode.addObject('LinearSolverConstraintCorrection')

    return FemNode
