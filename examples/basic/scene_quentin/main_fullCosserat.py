# This is a sample Python script.

# Press Maj+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

from backbone import Backbone
from controller import DirectController

def createScene(rootNode):
    rootNode.addObject('RequiredPlugin', name='Cosserat')  # Needed to use components [BeamHookeLawForceField]
    rootNode.addObject('RequiredPlugin',
                   name='Sofa.Component.AnimationLoop')  # Needed to use components [FreeMotionAnimationLoop]
    rootNode.addObject('RequiredPlugin',
                   name='Sofa.Component.Constraint.Lagrangian.Correction')  # Needed to use components [GenericConstraintCorrection]
    rootNode.addObject('RequiredPlugin',
                   name='Sofa.Component.Constraint.Lagrangian.Solver')  # Needed to use components [GenericConstraintSolver]
    rootNode.addObject('RequiredPlugin',
                   name='Sofa.Component.LinearSolver.Direct')  # Needed to use components [SparseLDLSolver]
    rootNode.addObject('RequiredPlugin', name='Sofa.Component.Mass')  # Needed to use components [UniformMass]
    rootNode.addObject('RequiredPlugin',
                   name='Sofa.Component.ODESolver.Backward')  # Needed to use components [EulerImplicitSolver]
    rootNode.addObject('RequiredPlugin',
                   name='Sofa.Component.SolidMechanics.Spring')  # Needed to use components [RestShapeSpringsForceField]
    rootNode.addObject('RequiredPlugin',
                   name='Sofa.Component.StateContainer')  # Needed to use components [MechanicalObject]
    rootNode.addObject('RequiredPlugin', name='Sofa.Component.Visual')  # Needed to use components [VisualStyle]

    rootNode.addObject('FreeMotionAnimationLoop')
    rootNode.addObject('VisualStyle',displayFlags='showVisual showBehavior showMapping')
    rootNode.addObject('GenericConstraintSolver',tolerance=1e-6,maxIterations=1000)

    rootNode.gravity = [0.0,-9.81,0.0]

    modeling = rootNode.addChild('Modeling')

    rod = modeling.addChild('Rod')

    ## Robot parameters
    # Geometry
    L = 200.0
    r_bb = 2.0
    e = 2.0
    w = 2.0
    I = w*(e**3)/12
    J = w*(e**3)*0.312
    A = w*e

    # Material
    E = 3000.0
    nu = 0.33
    G = E/(2*(1+nu))

    Nsd = 5
    Ns = 7
    Nfr = 2
    curvature = []
    s_curv = [0.0]
    seg_length = []
    for k in range(0,Ns*Nsd):
        curvature.append([0.0,0.0,0.0])
        s_curv.append((k+1)*L/(Ns*Nsd))
        seg_length.append(L/(Ns*Nsd))

    frames_cosserat = [[0.0,0.0,0.0,0.0,0.0,0.0,1.0]]
    frames = [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]]

    s_frame = [0.0]
    for k in range(1,Ns*Nfr*Nsd+1):
        frames_cosserat.append([0.0,0.0,0.0,0.0,0.0,0.0,1.0])
        frames.append([k*L/(Ns*Nfr*Nsd),0.0,0.0,0.0,0.0,0.0,1.0])
        s_frame.append(k*L/(Ns*Nfr*Nsd))

    radius = 10.0
    cable_routing = []
    for k in range(0,Nsd):
        cable_routing.append([0.0,radius,0.0])


    # Rigid base
    rigidBase = rod.addChild("RigidBase")
    rigidBase.addObject('MechanicalObject',template='Rigid3',position=[0.0,0.0,0.0,0.0,0.0,0.0,1.0],name = 'rigid_base_mo')
    rigidBase.addObject('FixedProjectiveConstraint',indices = 0)

    # Internal strain states
    internState = rod.addChild("InternalState")
    internState.addObject('MechanicalObject',template='Vec3',position=curvature, name = 'cosserat_state_mo')
    internState.addObject("BeamHookeLawForceField", length = seg_length, crossSectionShape = 'rectangular', lengthY = e, lengthZ = w, youngModulus = E, poissonRatio = nu )


    # Frames along the rod
    frameState = internState.addChild("CosseratFrames")
    frameState.addObject('RegularGridTopology', name='topo', n = [Ns*Nfr*Nsd+1,1,1], min = [0.0,0.0,0.0], max = [L,0.0,0.0])
    mo = frameState.addObject('MechanicalObject',template = 'Rigid3', name = 'cosserat_frame_mo', showObject= True, showObjectScale = 4.0)
    frameState.addObject('UniformMass', totalMass = 0.01)

    # Point force
    forceMag = 0.0
    cff = frameState.addObject('ConstantForceField', indices = Ns*Nfr*(Nsd-3), totalForce = [0.0,-forceMag,0.0,0.0,0.0,0.0])
    frameState.addObject('DiscreteCosseratMapping',curv_abs_input = s_curv, curv_abs_output = s_frame, name = "CosseratMapping", input1 = internState.cosserat_state_mo.getLinkPath(), input2 = rigidBase.rigid_base_mo.getLinkPath(), output = frameState.cosserat_frame_mo.getLinkPath(), debug = 0, radius = r_bb)
    
    # Rigid disks at each constant strain segment extremities
    disks_list = []
    for k in range(1,Nsd+1):
        spacerDiski = frameState.addChild("SpacerDisk"+str(k))
        disks_list.append(spacerDiski)

        spacerDiski.addObject("MechanicalObject", template="Vec3", name="disk_mo",
                           position=cable_routing[k-1], showObject=1, showObjectScale=1)
        spacerDiski.addObject("RigidMapping", input=frameState.cosserat_frame_mo.getLinkPath(),
                           output=disks_list[k - 1].disk_mo.getLinkPath(),index = k*Nfr*Ns)
    
    # Tendon
    cable = frameState.addChild("Cable")
    cable.addObject('MechanicalObject', template='Vec3', name="cable_mo", position = cable_routing)
    cc = cable.addObject('CableConstraint', template='Vec3', name='cable_ctr',
                                hasPullPoint='1',
                                pullPoint = [0.0,radius,0.0],
                                cableInitialLength = L,
                                indices=list(range(0, Nsd)),
                                maxPositiveDisp=100,
                                maxDispVariation=5,
                                minForce=0,
                                value  = 0.0,
                                valueType = 'force')
    indexPairsList = []
    inputList = []
    for k in range(0,Nsd):
        indexPairsList = indexPairsList+[k,0]
        inputList.append(disks_list[k].disk_mo.getLinkPath())

    cable.addObject('SubsetMultiMapping',
                              name="mapping",
                              input=inputList,
                              output=cable.cable_mo.getLinkPath(),
                              indexPairs=indexPairsList)
                              # indexPairs=[0, 0 + i - 1, 1, 0 + i - 1, 2, 0 + i - 1])
                              # [0, 0, 1, 0, 2, 0]

    simulation = rootNode.addChild('Simulation')
    simulation.addObject('EulerImplicitSolver',rayleighMass = 1.0, rayleighStiffness = 1.0, firstOrder=False)
    simulation.addObject('SparseLDLSolver')
    simulation.addObject('GenericConstraintCorrection')

    simulation.addChild(internState)

    rootNode.addObject(DirectController(rootNode = rootNode, cableConstraint = cc, cfMax = 0.0, constantForceField = cff, pfMax = 5.0, risingTime = 10.0, mechanicalObject = mo))

    # Construct the TDCR
    """ bb = Backbone(rotation=[0,0,0],translation=[0,0,0])
    scene.Modelling.addChild(bb)

    scene.addObject(DirectController(scene, [bb.Cables.Cable1,bb.Cables.Cable2,bb.Cables.Cable3]))

    # Add all the elements that need to be solved to the simulation node to ensure the solver find them
    scene.Simulation.addChild(bb.InternalStatesAndForces)
    scene.Simulation.addChild(bb.Frames)
    scene.Simulation.addChild(bb.SpacerDisks.Disk1)
    scene.Simulation.addChild(bb.SpacerDisks.Disk2)
    scene.Simulation.addChild(bb.SpacerDisks.Disk3)
    scene.Simulation.addChild(bb.Cables.Cable1)
    scene.Simulation.addChild(bb.Cables.Cable2)
    scene.Simulation.addChild(bb.Cables.Cable3) """

    # Add a controller for pulling on the tendons






# See PyCharm help at https://www.jetbrains.com/help/pycharm/
