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
from useful.params import ContactParameters as DefaultContactParams


def addHeader(parent_node, multithreading=False, inverse=False, is_constrained=False, is_contact=False, contact_params=None):
   
    """
    Adds default headers for a simulation with contact to the parent node.
    Also adds and returns three nodes: Settings, Modelling, Simulation.

    Args:
        parent_node: The parent node to add the headers to.
        multithreading: Enables multithreading (optional, default: False).
        inverse: Enables inverse kinematics (optional, default: False).
        is_constrained: Enables constraints (optional, default: False).
        is_contact: Enables contact simulation (optional, default: False).
        contact_params: Parameters for contact simulation (optional).

    Returns:
        A tuple containing the settings, modelling, and simulation nodes.
    """    

    settings = parent_node.addChild('Settings')
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

    parent_node.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels showCollisionModels '
                                                     'hideBoundingCollisionModels hideForceFields '
                                                     'showInteractionForceFields hideWireframe showMechanicalMappings')

    if is_constrained:
        parent_node.addObject('FreeMotionAnimationLoop', parallelCollisionDetectionAndFreeMotion=multithreading,
                              parallelODESolving=multithreading)
        if inverse:
            settings.addObject('RequiredPlugin', name="SoftRobots.Inverse")
            parent_node.addObject('QPInverseProblemSolver', name='ConstraintSolver', tolerance=1e-8, maxIterations=100,
                                  multithreading=multithreading, epsilon=1)
        else:
            parent_node.addObject('GenericConstraintSolver', name='ConstraintSolver', tolerance=1e-8, maxIterations=100,
                                  multithreading=multithreading)

    if is_contact:
        contactHeader(parent_node, contact_params=contact_params.contact_params)


# components needed for contact modeling

def contactHeader(parent_node, contact_params = DefaultContactParams):
    """
    Adds components for contact simulation to the parent node.

    Args:
        parent_node: The parent node to add the components to.
        contact_params: Optional contact parameters (default: None).
    """
    
    parent_node.addObject('CollisionPipeline')
    parent_node.addObject("DefaultVisualManagerLoop")
    parent_node.addObject('BruteForceBroadPhase')
    parent_node.addObject('BVHNarrowPhase')

    parent_node.addObject('RuleBasedContactManager', responseParams=contact_params.responseParams, response='FrictionContactConstraint')
    parent_node.addObject('LocalMinDistance', alarmDistance=contact_params.alarmDistance, contactDistance=contact_params.contactDistance)


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


def addSolverNode(parent_node, name='solverNode', template='CompressedRowSparseMatrixd', rayleighMass=0., rayleighStiffness=0.,
                  firstOrder=False,
                  iterative=False, isConstrained=False):
    """
    Adds solvers (EulerImplicitSolver, LDLSolver, GenericConstraintCorrection) to the given node.

    Args:
        name:
        isConstrained:
        parent_node:
        template: for the LDLSolver
        rayleighMass:
        rayleighStiffness:
        firstOrder: for the implicit scheme
        iterative: iterative solver

    Usage:
        addSolversNode(node)
    """
    solver_node = parent_node.addChild(name)
    solver_node.addObject('EulerImplicitSolver', firstOrder=firstOrder, rayleighStiffness=rayleighStiffness,
                         rayleighMass=rayleighMass)
    if iterative:
        solver_node.addObject('CGLinearSolver', name='Solver', template=template)
    else:
        solver_node.addObject('SparseLDLSolver', name='Solver', template=template, printLog=True)
    if isConstrained:
        solver_node.addObject('GenericConstraintCorrection', linearSolver=solver_node.Solver.getLinkPath())

    return solver_node


def addFEMObject(parent_node, path):
    finger_solver = addSolverNode(parent_node)

    # Load a VTK tetrahedral mesh and expose the resulting topology in the scene .
    loader = finger_solver.addObject('MeshVTKLoader', name='loader', filename=f'{path}finger.vtk',
                                    translation="-17.5 -12.5 7.5",
                                    rotation="0 180 0")
    finger_solver.addObject('TetrahedronSetTopologyContainer', position=loader.position.getLinkPath(),
                           tetras=loader.tetras.getLinkPath(), name='container')
    # Create a MechanicalObject component to stores the DoFs of the model
    finger_solver.addObject('MechanicalObject', template='Vec3', name='dofs')

    # Gives a mass to the model
    finger_solver.addObject('UniformMass', totalMass='0.075')
    # Add a TetrahedronFEMForceField component which implement an elastic material model
    # solved using the Finite Element Method on
    # tetrahedrons.
    finger_solver.addObject('TetrahedronFEMForceField', template='Vec3d', name='forceField', method='large',
                           poissonRatio='0.45', youngModulus='500')

    finger_solver.addObject('BoxROI', name='ROI1', box='-18 -15 -8 2 -3 8', drawBoxes='true')
    finger_solver.addObject('RestShapeSpringsForceField',
                           points='@ROI1.indices', stiffness='1e12')
    ##########################################
    # Cable points                           #
    ##########################################
    # Mapped points inside the finger volume, these points attached to the FE model
    # are constrained to slide on the cable.

    FEMpos = [" 0.0 0 0 15 0 0 30 0 0 45 0 0 60 0 0 66 0 0 81 0.0 0.0"]
    fem_points = finger_solver.addChild('femPoints')
    fem_points.addObject('MechanicalObject', name="pointsInFEM", position=FEMpos, showObject="1",
                        showIndices="1")
    fem_points.addObject('BarycentricMapping')


def addMappedPoints(parent_node, name="pointsInFEM", position=None, showObject="1",
                    showIndices="1", femPos="0.0 0 0 15 0 0 30 0 0 45 0 0 60 0 0 66 0 0 81 0.0 0.0"):
    if position is None:
        position = femPos

    femPoints = parent_node.addChild(name)
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


def Finger(parent_node=None, name="Finger", rotation=None, translation=None, fixingBox=None, path=None, femPos=None):
    if fixingBox is None:
        fixingBox = [-18, -15, -8, 2, -3, 8]
    if rotation is None:
        rotation = [0., 180., 0]
    if translation is None:
        translation = [-17.5, -12.5, 7.5]
    if path is None:
        path = f'{os.path.dirname(os.path.abspath(__file__))}/'

    # TODO : add physical properties as the finger input
    finger = parent_node.addChild(name)
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


# Test
def createScene(root_node):
    addHeader(root_node)
    addSolverNode(root_node)
