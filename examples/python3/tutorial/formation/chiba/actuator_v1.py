import Sofa
from fingerController import FingerController  # for keyboard controlling
from loopstest_function_param3 import looptest  # algorithm for generating the fiber loops

from useful.header import addHeader, addSolverNode


def createScene(root_node):
    addHeader(root_node, isConstrained=True)

    root_node.gravity.value = [-9810, 0, 0]
    # root_node.addObject('FreeMotionAnimationLoop')
    # root_node.addObject('ProjectedGaussSeidelConstraintSolver', tolerance=1e-9, maxIterations=50000)
    root_node.addObject('EulerImplicitSolver', name='odesolver', rayleighStiffness=0.01, rayleighMass=0.01)

    finger = root_node.addChild('finger')
    finger.addObject('SparseLDLSolver')
    finger.addObject('MeshGmshLoader', name='loader', filename='data/mesh/act23mmwithconnector.msh',
                     rotation=[90, 0, -90], translation=[22.5, 0, 0])

    finger.addObject('MeshTopology', src='@loader', name='container')
    finger.addObject('MechanicalObject', name='tetras', template='Vec3', showObject=False, showObjectScale=1)

    youngModulus = 263.8
    poissonRatio = 0.49
    mu_ = youngModulus / (2 * (1 + poissonRatio))
    k0_ = youngModulus / (3 * (1 - 2 * poissonRatio))
    finger.addObject('TetrahedronHyperelasticityFEMForceField', materialName="NeoHookean",
                     ParameterSet=str(mu_) + " " + str(k0_))
    finger.addObject('UniformMass', totalMass=0.0008)
    finger.addObject('BoxROI', name='boxROISubTopo', box=[23, -10, -8, 26, 10, 8], strict=False, drawBoxes=True)
    finger.addObject('BoxROI', name='boxROISubTopo2', box=[0, -10, -4, 28, 10, -2], strict=False, drawBoxes=True)
    finger.addObject('BoxROI', name='boxROI', box=[-2, -10, -20, 2, 10, 20], drawBoxes=True)
    finger.addObject('RestShapeSpringsForceField', points='@boxROI.indices', stiffness=1e12, angularStiffness=1e12)
    finger.addObject('LinearSolverConstraintCorrection')
    # Plastic Part
    modelSubTopo = finger.addChild('modelSubTopo')
    modelSubTopo.addObject('MeshTopology', position='@loader.position', tetrahedra='@boxROISubTopo.tetrahedraInROI',
                           name='container')
    modelSubTopo.addObject('TetrahedronFEMForceField', template='Vec3', name='FEM', method='large', poissonRatio=0.3,
                           youngModulus=4400000)


    # Parameters for the fiber reinforcement 
    num_loops = 8  # number of loops
    loop_distance = 3  # Distance between fiber loops
    loops_data, spring_data, indices1, indices2, lengths = looptest(loop_distance, num_loops)

    # print('-------------------------loops_data-------------------')
    # print(loops_data)
    # print('-------------------------spring_data-------------------')
    # print(spring_data)
    # print('-------------------------indices1-------------------')
    # print(indices1)
    # print('-------------------------indices2-------------------')
    # print(indices2)
    # print('-------------------------lengths-------------------')
    # print(lengths)

    # nbDOFs = 10
    Ks = 1e7  # spring stinffness
    Kd = 5
    position = [[i * loop_distance + 1, 0, -4] for i in range(num_loops)]
    print("Positions: ", position)
    edges = [[i, i + 1] for i in range(num_loops - 1)]
    print("edges: ", edges)

    # Hitching (transversal fiber)
    hitching = finger.addChild("hitching")  # fiber string that crosses the actuator transversally
    hitching.addObject("MechanicalObject", template="Vec3", name="DOF",
                       position=[[i * loop_distance + 1, 0, -4] for i in range(num_loops)],
                       showObject=True, showObjectScale=1)
    hitching.addObject('MeshTopology', name='lines',
                       lines=[[i, i + 1] for i in range(num_loops - 1)])
    hitching.addObject("SpringForceField", template="Vec3d", name="springs", showArrowSize=0.5, drawMode=1,
                       stiffness=Ks, damping=Kd, indices1=[0, 1, 2, 3, 4, 5, 6], indices2=[1, 2, 3, 4, 5, 6, 7],
                       lengths=loop_distance)
    hitching.addObject('BarycentricMapping', name='mapping')
    hitching.addObject('UniformMass', totalMass=0.0008)

    return root_node

    # Mapping

    # hitching.addObject("FixedConstraint", name="FixedConstraint", indices=[0])

    # finger.addObject('MechanicalMatrixMapper', template="Vec3,Vec3",
    #                nodeToParse=hitching.getLinkPath(),  # where to find the forces to map
    #                object1=finger.tetras.getLinkPath())  # in case of multi-mapping, here you can give the second parent

    # Fiber Loops

    for k in range(0, loopnum):
        fiber = finger.addChild("fiber" + str(k))
        fiber.addObject("MechanicalObject", template="Vec3", name="DOF",
                        position=loops[k],
                        showObject=True, showObjectScale=1)
        fiber.addObject('MeshTopology', name='lines', lines=[[i, i + 1] for i in range(loopnum)])

        # fiber.addObject("FixedConstraint", name="FixedConstragetLinkPathint", indices=[0])
        fiber.addObject("StiffSpringForceField", template="Vec3d", name="springs", showArrowSize=0.5, drawMode=1,
                        stiffness=Ks, damping=Kd, indices1=indices1, indices2=indices2, lengths=lengths)

        # Mapping 

        fiber.addObject('BarycentricMapping', name='mapping')
        # finger.addObject("MechanicalMatrixMapper", template="Vec3,Vec3",name="MechanicalMatrixMapper"+str(k), 
        #                 nodeToParse=fiber.getLinkPath(), 
        #                 object1=finger.tetras.getLinkPath(),
        #                 object2=hitching.DOF.getLinkPath()) 

    # Air Chamber(cavity)

    tracker = finger.addChild("tracker")  # trackers  for measuring bending angle
    tracker.addObject("MechanicalObject", template="Vec3", name="DOF",
                      position=[[29, 0, -4], [31, 0, -4]],
                      showObject=True, showObjectScale=1, drawMode=1)
    tracker.addObject("MeshTopology", name="points", points=[[29, 0, 0], [31, 0, 0]])
    tracker.addObject('BarycentricMapping', name='mapping')

    # air Pressure Cavity. Note: 0.1=10kPa
    cavity = finger.addChild('cavity')
    cavity.addObject('MeshSTLLoader', name='cavityLoader', filename='data/mesh/CHAMBER20mm-075mm.stl',
                     rotation=[90, 0, 90])
    cavity.addObject('MeshTopology', src='@cavityLoader', name='cavityMesh')
    cavity.addObject('MechanicalObject', name='cavity')
    cavity.addObject('SurfacePressureConstraint', name='SurfacePressureConstraint', template='Vec3', value=0,
                     triangles='@cavityMesh.triangles', valueType='pressure')
    cavity.addObject('BarycentricMapping', name='mapping')

    rootNode.addObject(FingerController(rootNode))
