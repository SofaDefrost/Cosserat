"""_summary_ Basic scene using Cosserat in SofaPython3.
    The Finger is modeled using FEM model while de cable is modeled using cosserat theory.
    The link between these two mechanical models is performed by constraint based on Lagrangian Multiplier

Returns:
    _type_: _description_
"""

__authors__ = "Younes"
__contact__ = "adagolodjo@protonmail.com, yinoussa.adagolodjo@inria.fr"
__version__ = "1.0.0"
__copyright__ = "(c) 2020,Inria"
__date__ = "july 2023"

from stlib3.physics.deformable import ElasticMaterialObject
from stlib3.physics.constraints import FixedBox
import os


def addHeader(parentNode, multithreading=False, inverse=False, isConstrained=False, isContact=False, params=None):
    """
    Adds to rootNode the default headers for a simulation with contact. Also adds and returns three nodes:
        - Settings
        - Modelling
        - Simulation

    Args:
        isContact:
        inverse:
        isConstrained:
        parentNode:
        multithreading:

    Usage:
        addHeader(rootNode)

    Returns:
        the three SOFA nodes {settings, modelling, simulation}
    """
    settings = parentNode.addChild('Settings')
    settings.addObject('RequiredPlugin', pluginName=[
        "Cosserat", "Sofa.Component.AnimationLoop",  # Needed to use components FreeMotionAnimationLoop
        "Sofa.Component.Collision.Detection.Algorithm",
        "Sofa.Component.Collision.Detection.Intersection",  # Needed to use components LocalMinDistance
        "Sofa.Component.Collision.Response.Contact",  # Needed to use components RuleBasedContactManager
        "Sofa.Component.Constraint.Lagrangian.Correction",  # Needed to use components GenericConstraintCorrection
        "Sofa.Component.Constraint.Lagrangian.Solver",  # Needed to use components GenericConstraintSolver
        "Sofa.Component.Constraint.Projective",  # Needed to use components FixedConstraint
        "Sofa.Component.IO.Mesh",  # Needed to use components MeshOBJLoader, MeshSTLLoader
        "Sofa.Component.Mass", 'Sofa.Component.LinearSolver.Direct',
        "Sofa.Component.Setting",  # Needed to use components BackgroundSetting
        "Sofa.Component.SolidMechanics.Spring",  # Needed to use components RestShapeSpringsForceField
        "Sofa.Component.Topology.Container.Constant",  # Needed to use components MeshTopology
        "Sofa.Component.Topology.Container.Dynamic", "Sofa.Component.Playback",
        "Sofa.Component.Playback", "Sofa.Component.Visual",  # Needed to use components VisualStyle
        "Sofa.Component.Topology.Container.Grid",  # Needed to use components RegularGridTopology
        "Sofa.Component.Topology.Mapping",  # Needed to use components Edge2QuadTopologicalMapping
        "Sofa.GL.Component.Rendering3D",  # Needed to use components OglGrid, OglModel
        "Sofa.GUI.Component", "Sofa.Component.Collision.Geometry", "Sofa.Component.LinearSolver.Direct",
        "Sofa.Component.Mapping.Linear", "Sofa.Component.MechanicalLoad",
        'Sofa.Component.Engine.Select', 'Sofa.Component.SolidMechanics.FEM.Elastic',
        "Sofa.Component.StateContainer", 'Sofa.Component.ODESolver.Backward',
        'Sofa.Component.SolidMechanics.FEM.HyperElastic'
    ])

    settings.addObject('BackgroundSetting', color=[1, 1, 1, 1])
    # settings.addObject('AttachBodyButtonSetting', stiffness=1e6)

    parentNode.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels showCollisionModels '
                                                     'hideBoundingCollisionModels hideForceFields '
                                                     'showInteractionForceFields hideWireframe showMechanicalMappings')
    if isConstrained:
        parentNode.addObject('FreeMotionAnimationLoop', parallelCollisionDetectionAndFreeMotion=multithreading,
                             parallelODESolving=multithreading)
        if inverse:
            settings.addObject('RequiredPlugin', name="SoftRobots.Inverse")
            parentNode.addObject('QPInverseProblemSolver', name='ConstraintSolver', tolerance=1e-8, maxIterations=100,
                                 multithreading=multithreading, epsilon=1)
        else:
            parentNode.addObject('GenericConstraintSolver', name='ConstraintSolver', tolerance=1e-8, maxIterations=100,
                                 multithreading=multithreading, printLog=1)

    if isContact:
        contactHeader(parentNode, _contact_params=params.contactParams)


# components needed for contact modeling
def contactHeader(parentNode, _contact_params=None):
    parentNode.addObject('CollisionPipeline')
    parentNode.addObject("DefaultVisualManagerLoop")
    parentNode.addObject('BruteForceBroadPhase')
    parentNode.addObject('BVHNarrowPhase')
    if not _contact_params == None:
        parentNode.addObject('RuleBasedContactManager', responseParams=_contact_params.responseParams,
                             response='FrictionContactConstraint')
        parentNode.addObject('LocalMinDistance',  alarmDistance=_contact_params.alarmDistance,
                             contactDistance=_contact_params.contactDistance)
    else :
        parentNode.addObject('RuleBasedContactManager', responseParams='mu=0.1', response='FrictionContactConstraint')
        parentNode.addObject('LocalMinDistance',  alarmDistance=0.05, contactDistance=0.01)


def addVisual(node):
    """
    Adds a visual model to the given node.

    Args:
        node:
    Usage:
        addVisual(node)
    """
    node.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels hideCollisionModels '
                                               'hideBoundingCollisionModels hideForceFields '
                                               'showInteractionForceFields hideWireframe showMechanicalMappings')
    return node


def addSolverNode(node, name='solverNode', template='CompressedRowSparseMatrixd', rayleighMass=0., rayleighStiffness=0.,
                  firstOrder=False,
                  iterative=False, isConstrained=False):
    """
    Adds solvers (EulerImplicitSolver, LDLSolver, GenericConstraintCorrection) to the given node.

    Args:
        name:
        isConstrained:
        node:
        template: for the LDLSolver
        rayleighMass:
        rayleighStiffness:
        firstOrder: for the implicit scheme
        iterative: iterative solver

    Usage:
        addSolversNode(node)
    """
    solverNode = node.addChild(name)
    solverNode.addObject('EulerImplicitSolver', firstOrder=firstOrder, rayleighStiffness=rayleighStiffness,
                         rayleighMass=rayleighMass)
    if iterative:
        solverNode.addObject('CGLinearSolver', name='Solver', template=template)
    else:
        solverNode.addObject('SparseLDLSolver', name='Solver', template=template, printLog=True)
    if isConstrained:
        solverNode.addObject('GenericConstraintCorrection', linearSolver=solverNode.Solver.getLinkPath())

    return solverNode


def addFEMObject(parentNode, path):
    fingerSolver = addSolverNode(parentNode)

    # Load a VTK tetrahedral mesh and expose the resulting topology in the scene .
    loader = fingerSolver.addObject('MeshVTKLoader', name='loader', filename=f'{path}finger.vtk',
                                    translation="-17.5 -12.5 7.5",
                                    rotation="0 180 0")
    fingerSolver.addObject('TetrahedronSetTopologyContainer', position=loader.position.getLinkPath(),
                           tetras=loader.tetras.getLinkPath(), name='container')
    # Create a MechanicalObject component to stores the DoFs of the model
    fingerSolver.addObject('MechanicalObject', template='Vec3', name='dofs')

    # Gives a mass to the model
    fingerSolver.addObject('UniformMass', totalMass='0.075')
    # Add a TetrahedronFEMForceField component which implement an elastic material model
    # solved using the Finite Element Method on
    # tetrahedrons.
    fingerSolver.addObject('TetrahedronFEMForceField', template='Vec3d', name='forceField', method='large',
                           poissonRatio='0.45', youngModulus='500')

    fingerSolver.addObject('BoxROI', name='ROI1', box='-18 -15 -8 2 -3 8', drawBoxes='true')
    fingerSolver.addObject('RestShapeSpringsForceField',
                           points='@ROI1.indices', stiffness='1e12')
    ##########################################
    # Cable points                           #
    ##########################################
    # Mapped points inside the finger volume, these points attached to the FE model
    # are constrained to slide on the cable.

    FEMpos = [" 0.0 0 0 15 0 0 30 0 0 45 0 0 60 0 0 66 0 0 81 0.0 0.0"]
    femPoints = fingerSolver.addChild('femPoints')
    femPoints.addObject('MechanicalObject', name="pointsInFEM", position=FEMpos, showObject="1",
                        showIndices="1")
    femPoints.addObject('BarycentricMapping')


def addMappedPoints(parentNode, name="pointsInFEM", position=None, showObject="1",
                    showIndices="1", femPos="0.0 0 0 15 0 0 30 0 0 45 0 0 60 0 0 66 0 0 81 0.0 0.0"):
    if position is None:
        position = femPos

    femPoints = parentNode.addChild(name)
    femPoints.addObject(
        'MechanicalObject',
        name=f'{name}Mo',
        position=position,
        showObject=showObject,
        showIndices=showIndices,
        template='Vec3d'
    )
    femPoints.addObject('BarycentricMapping')
    return femPoints


def Finger(parentNode=None, name="Finger", rotation=None, translation=None, fixingBox=None, path=None, femPos=None):
    if fixingBox is None:
        fixingBox = [-18, -15, -8, 2, -3, 8]
    if rotation is None:
        rotation = [0., 180., 0]
    if translation is None:
        translation = [-17.5, -12.5, 7.5]
    if path is None:
        path = f'{os.path.dirname(os.path.abspath(__file__))}/'

    # TODO : add physical properties as the finger input
    finger = parentNode.addChild(name)
    e_object = ElasticMaterialObject(finger,
                                     volumeMeshFileName=f'{path}finger.vtk',
                                     poissonRatio=0.48,
                                     youngModulus=500,
                                     totalMass=0.075,
                                     surfaceColor=[0.0, 0.7, 0.7, 0.5],
                                     surfaceMeshFileName=f'{path}finger.stl',
                                     rotation=rotation,
                                     translation=translation)
    FixedBox(e_object,
             doVisualization=True,
             atPositions=fixingBox)
    """Add sliding point to the scene, these points are those on which the cable slides"""
    if femPos is None:
        femPoints = addMappedPoints(e_object)
    else:
        femPoints = addMappedPoints(e_object, femPos=femPos)

    finger.addChild(e_object)
    return finger, femPoints
