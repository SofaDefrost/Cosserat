
def addHeader(rootnode, multithreading=False, inverse=False, isConstrained=False):
    """
    Adds to rootnode the default headers for a simulation with contact. Also adds and returns three nodes:
        - Settings
        - Modelling
        - Simulation

    Args:
        inverse:
        isConstrained:
        rootnode:
        multithreading:

    Usage:
        addHeader(rootnode)

    Returns:
        the three SOFA nodes {settings, modelling, simulation}
    """
    settings = rootnode.addChild('Settings')
    settings.addObject('RequiredPlugin', pluginName=[
                        "ArticulatedSystemPlugin",
                        "Cosserat",   # Needed to use component RigidDistanceMapping
                        "Sofa.Component.AnimationLoop",  # Needed to use components FreeMotionAnimationLoop
                        "Sofa.Component.Collision.Detection.Algorithm",
                        "Sofa.Component.Collision.Detection.Intersection",  # Needed to use components LocalMinDistance
                        "Sofa.Component.Collision.Response.Contact",  # Needed to use components RuleBasedContactManager
                        "Sofa.Component.Constraint.Lagrangian.Correction",  # Needed to use components GenericConstraintCorrection
                        "Sofa.Component.Constraint.Lagrangian.Solver",  # Needed to use components GenericConstraintSolver
                        "Sofa.Component.Constraint.Projective",  # Needed to use components FixedConstraint
                        "Sofa.Component.IO.Mesh",  # Needed to use components MeshOBJLoader, MeshSTLLoader
                        "Sofa.Component.Mass",  # Needed to use components UniformMass
                        "Sofa.Component.Setting",  # Needed to use components BackgroundSetting
                        "Sofa.Component.SolidMechanics.Spring",  # Needed to use components RestShapeSpringsForceField
                        "Sofa.Component.Topology.Container.Constant",  # Needed to use components MeshTopology
                        "Sofa.Component.Topology.Container.Dynamic",
                        # Needed to use components EdgeSetTopologyContainer, EdgeSetTopologyModifier,
                        # QuadSetTopologyContainer, QuadSetTopologyModifier
                        "Sofa.Component.Topology.Container.Grid",  # Needed to use components RegularGridTopology
                        "Sofa.Component.Topology.Mapping",  # Needed to use components Edge2QuadTopologicalMapping
                        "Sofa.Component.Visual",  # Needed to use components VisualStyle
                        "Sofa.GL.Component.Rendering3D",  # Needed to use components OglGrid, OglModel
                        "Sofa.GUI.Component",  # Needed to use components AttachBodyButtonSetting
    ])
    settings.addObject('BackgroundSetting', color=[1, 1, 1, 1])
    # settings.addObject('AttachBodyButtonSetting', stiffness=1e6)

    if isConstrained:
        rootnode.addObject('VisualStyle')
        rootnode.addObject('CollisionPipeline')
        rootnode.addObject("DefaultVisualManagerLoop")
        rootnode.addObject('RuleBasedContactManager', responseParams='mu=0.8', response='FrictionContactConstraint')
        rootnode.addObject('BruteForceBroadPhase')
        rootnode.addObject('BVHNarrowPhase')
        rootnode.addObject('LocalMinDistance', alarmDistance=0.05, contactDistance=0.01)
        rootnode.addObject('FreeMotionAnimationLoop', parallelCollisionDetectionAndFreeMotion=multithreading, parallelODESolving=multithreading)
        if inverse:
            settings.addObject('RequiredPlugin', name="SoftRobots.Inverse")
            rootnode.addObject('QPInverseProblemSolver', name='ConstraintSolver', tolerance=1e-8, maxIterations=100, multithreading=multithreading,
                               epsilon=1)
        else:
            rootnode.addObject('GenericConstraintSolver', name='ConstraintSolver', tolerance=1e-8, maxIterations=100, multithreading=multithreading)
    rootnode.dt = 0.01

    modelling = rootnode.addChild('Modelling')
    modelling.addChild('Topology')
    simulation = rootnode.addChild('Simulation')
    return settings, modelling, simulation

def addSolvers(node, template='CompressedRowSparseMatrixd', rayleighMass=0., rayleighStiffness=0., firstOrder=False,
               iterative=False):
    """
    Adds solvers (EulerImplicitSolver, LDLSolver, GenericConstraintCorrection) to the given node.

    Args:
        node:
        template: for the LDLSolver
        rayleighMass:
        rayleighStiffness:
        firstOrder: for the implicit scheme
        iterative: iterative solver

    Usage:
        addSolvers(node)
    """
    node.addObject('RequiredPlugin', name='SofaPlugins', pluginName=['Sofa.Component.LinearSolver.Direct',
                                                                     'Sofa.Component.ODESolver.Backward'])
    node.addObject('EulerImplicitSolver', firstOrder=firstOrder, rayleighStiffness=rayleighStiffness, rayleighMass=rayleighMass)
    if iterative:
        node.addObject('CGLinearSolver', name='Solver', template=template)
    else:
        node.addObject('SparseLDLSolver', name='Solver', template=template)
    node.addObject('GenericConstraintCorrection', linearSolver=node.Solver.getLinkPath())

    return node


# Test
def createScene(rootNode):

    addHeader(rootNode)

    node = rootNode.addChild('rootNode')
    solverNode = addSolvers(node)


