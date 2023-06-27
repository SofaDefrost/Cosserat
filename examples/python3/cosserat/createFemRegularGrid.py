# -*- coding: utf-8 -*-
# -*- coding: utf-8 -*-
"""Basic scene using Cosserat in SofaPython3.

Based on the work done with SofaPython. See POEMapping.py
"""

__authors__ = "younesssss"
__contact__ = "adagolodjo@protonmail.com, yinoussa.adagolodjo@inria.fr"
__version__ = "1.0.0"
__copyright__ = "(c) 2021,Inria"
__date__ = "March 16 2021"


def createFemCube(parentNode):
    FemNode = parentNode.addChild("FemNode")
    FemNode.addObject('VisualStyle', displayFlags='showBehaviorModels hideCollisionModels hideBoundingCollisionModels '
                                                  'showForceFields hideInteractionForceFields showWireframe')
    gelVolume = FemNode.addChild("gelVolume")
    gelVolume.addObject("RegularGridTopology", name="HexaTop", n="6 6 6", min="40 -16 -10", max="100 20 10")
    gelVolume.addObject("TetrahedronSetTopologyContainer", name="Container", position="@HexaTop.position")
    gelVolume.addObject("TetrahedronSetTopologyModifier", name="Modifier")
    gelVolume.addObject("Hexa2TetraTopologicalMapping", input="@HexaTop", output="@Container", swapping="false")

    GelSurface = FemNode.addChild("GelSurface")
    GelSurface.addObject("TriangleSetTopologyContainer", name="Container", position="@../GelVolume/HexaTop.position")

    gelNode = FemNode.addChild("gelNode")
    gelNode.addObject("EulerImplicitSolver", rayleighMass="0.1", rayleighStiffness="0.1")
    gelNode.addObject('SparseLDLSolver', name='preconditioner', template="CompressedRowSparseMatrixMat3x3d")
    gelNode.addObject('TetrahedronSetTopologyContainer', src="@../gelVolume/Container", name='container')
    gelNode.addObject('MechanicalObject', name='tetras', template='Vec3d')
    gelNode.addObject('TetrahedronFEMForceField', template='Vec3d', name='FEM', method='large',
                      poissonRatio='0.45', youngModulus='100')
    gelNode.addObject('BoxROI', name='ROI1', box='40 -17 -10 100 -14 10', drawBoxes='true')
    gelNode.addObject('RestShapeSpringsForceField', points='@ROI1.indices', stiffness='1e12')

    surfaceNode = gelNode.addChild("surfaceNode")
    surfaceNode.addObject('TriangleSetTopologyContainer', name="surfContainer", src="@../../GelSurface/Container")
    surfaceNode.addObject('MechanicalObject', name='msSurface')
    surfaceNode.addObject('TriangleCollisionModel', name='surface')
    surfaceNode.addObject('LineCollisionModel', name='line')
    surfaceNode.addObject('BarycentricMapping')

    gelNode.addObject('LinearSolverConstraintCorrection')

    return FemNode


def createFemCubeWithParams(parentNode, geometry):
    FemNode = parentNode.addChild("FemNode")

    gelVolume = FemNode.addChild("gelVolume")
    gelVolume.addObject("RegularGridTopology", name="HexaTop", n=geometry.mesh, min=geometry.minVol,
                        max=geometry.maxVol)
    cont = gelVolume.addObject("TetrahedronSetTopologyContainer", name="TetraContainer", position="@HexaTop.position")
    gelVolume.addObject("TetrahedronSetTopologyModifier", name="Modifier")
    gelVolume.addObject("Hexa2TetraTopologicalMapping", input="@HexaTop", output="@TetraContainer", swapping="false")

    GelSurface = FemNode.addChild("GelSurface")
    GelSurface.addObject("TriangleSetTopologyContainer", name="triangleContainer",
                         position="@../gelVolume/HexaTop.position")
    GelSurface.addObject("TriangleSetTopologyModifier", name="Modifier")
    GelSurface.addObject("Tetra2TriangleTopologicalMapping", input="@../gelVolume/TetraContainer",
                         output="@triangleContainer", flipNormals="false")

    gelNode = FemNode.addChild("gelNode")
    gelNode.addObject("EulerImplicitSolver", rayleighMass=geometry.rayleigh, rayleighStiffness=geometry.rayleigh)
    gelNode.addObject('SparseLDLSolver', name='precond', template="CompressedRowSparseMatrixMat3x3d")
    gelNode.addObject('TetrahedronSetTopologyContainer', src="@../gelVolume/TetraContainer", name='container')
    gelNode.addObject('MechanicalObject', name='tetras', template='Vec3d')
    gelNode.addObject('TetrahedronFEMForceField', template='Vec3d', name='FEM', method='large',
                      poissonRatio=geometry.poissonRatio, youngModulus=geometry.youngModulus)
    gelNode.addObject('BoxROI', name='ROI1', box=geometry.box, drawBoxes='true')
    gelNode.addObject('RestShapeSpringsForceField', points='@ROI1.indices', stiffness='1e12')

    surfaceNode = gelNode.addChild("surfaceNode")
    surfaceNode.addObject('TriangleSetTopologyContainer', name="surfContainer",
                          src="@../../GelSurface/triangleContainer")
    surfaceNode.addObject('MechanicalObject', name='msSurface')
    surfaceNode.addObject('TriangleCollisionModel', name='surface')
    surfaceNode.addObject('LineCollisionModel', name='line')
    surfaceNode.addObject('BarycentricMapping')
    visu = surfaceNode.addChild("visu")

    visu.addObject("OglModel", name="Visual", src="@../surfContainer",  color="0.0 0.1 0.9 0.40" )
    visu.addObject("BarycentricMapping", input="@..", output="@Visual")

    gelNode.addObject('LinearSolverConstraintCorrection')

    return FemNode
