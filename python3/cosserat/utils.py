from robocop.usefulFunction import BuildCosseratGeometry, buildEdges


def addSofaPlugins(node, newPlugins):
    plugins = node.SofaPlugins.pluginName.value
    for plugin in newPlugins.split():
        plugins.append(plugin)
    node.SofaPlugins.pluginName = plugins


def addRigidObject(node, filename, collisionFilename=None, position=None, scale=None,
                   textureFilename='', color=None, density=0.002, name='Object', withSolver=True):
    if color is None:
        color = [1, 1, 1]
    if scale is None:
        scale = [1, 1, 1]
    if position is None:
        position = [0, 0, 0, 0, 0, 0, 1]
    if collisionFilename is None:
        collisionFilename = filename

    _object = node.addChild(name)
    _object.addObject('RequiredPlugin', name='SofaPlugins', pluginName='SofaRigid SofaLoader')
    _object.addObject('MechanicalObject', template='Rigid3', position=position, showObject=False, showObjectScale=5)

    if withSolver:
        _object.addObject('EulerImplicitSolver')
        _object.addObject('CGLinearSolver', tolerance=1e-5, iterations=25)
        _object.addObject('UncoupledConstraintCorrection')

    visu = _object.addChild('Visu')
    visu.addObject('MeshObjLoader', name='loader', filename=filename, scale3d=scale)
    visu.addObject('OglModel', src='@loader', texturename=textureFilename, color=color if textureFilename == '' else '')
    visu.addObject('RigidMapping')

    _object.addObject('GenerateRigidMass', name='mass', density=density, src=visu.loader.getLinkPath())
    _object.mass.init()
    translation = list(_object.mass.centerToOrigin.value)
    _object.addObject('UniformMass', vertexMass="@mass.rigidMass")

    visu.loader.translation = translation

    collision = _object.addChild('Collision')
    collision.addObject('MeshObjLoader', name='loader', filename=collisionFilename, scale3d=scale)
    collision.addObject('MeshTopology', src='@loader')
    collision.addObject('MechanicalObject', translation=translation)
    collision.addObject('TriangleCollisionModel')
    collision.addObject('LineCollisionModel')
    collision.addObject('PointCollisionModel')
    collision.addObject('RigidMapping')

    return _object


def addEdgeCollision(parentNode, position3D, edges):
    CollisInstrumentCombined = parentNode.addChild('CollisInstrumentCombined')
    CollisInstrumentCombined.addObject('EdgeSetTopologyContainer', name="collisEdgeSet", position=position3D,
                                       edges=edges)
    CollisInstrumentCombined.addObject('EdgeSetTopologyModifier', name="collisEdgeModifier")
    CollisInstrumentCombined.addObject('MechanicalObject', name="CollisionDOFs")
    CollisInstrumentCombined.addObject('LineCollisionModel', bothSide="1", group='2')
    CollisInstrumentCombined.addObject('PointCollisionModel', bothSide="1", group='2')
    CollisInstrumentCombined.addObject('IdentityMapping', name="mapping")


def createCosserat(parent, config, config_material, config_simu, name="Cosserat", orientation=None, radius=0.,
                   last_frame=None, rigidBase=None, _showObjectScale=0.02, solver=None, fixedBase=True):
    if last_frame is None:
        last_frame = []
    if orientation is None:
        orientation = [0, 0, 0, 1]
    [x, y, z] = config['init_pos']
    buildCollision = config['collision']

    Young_modulus = config_material['youngModulus']
    poissonRatio = config_material['poissonRatio']
    length_Y = config_material['length_Y']
    length_Z = config_material['length_Z']
    shape = config_material['shape']

    # some parameters
    stiffness = config_simu['stiffness']
    angularStiffness = config_simu['angularStiffness']

    # ################################
    # #           RigidBase         ##
    # ################################

    if parent.hasObject("EulerImplicitSolver") is False:
        base = parent.addChild(name)
        # base.addObject('EulerImplicitSolver',  rayleighStiffness="0.2", rayleighMass='0.1')
        base.addObject('EulerImplicitSolver', rayleighStiffness="1.2", rayleighMass='1.1')
        base.addObject('SparseLDLSolver', name='solver', template="CompressedRowSparseMatrixd")
        base.addObject('GenericConstraintCorrection')
    else:
        # This means you already defined the solver node and give it as parent of this function
        base = parent.addChild(name)

    # base.addObject('EulerImplicitSolver',  rayleighStiffness="0.2", rayleighMass='0.1')

    if rigidBase is None:
        rigidBaseNode = base.addChild('rigidBase')

        RigidBaseMO = rigidBaseNode.addObject('MechanicalObject', template='Rigid3d',
                                              name="RigidBaseMO", position=[x, y, z] + orientation, showObject=1,
                                              showObjectScale=0.2)
        if fixedBase:
            rigidBaseNode.addObject('RestShapeSpringsForceField', name='spring' + str(x), stiffness=stiffness,
                                    angularStiffness=angularStiffness, external_points=0, mstate="@RigidBaseMO",
                                    points=0,
                                    template="Rigid3d")
    else:
        rigidBaseNode = base.addChild(rigidBase)
        RigidBaseMO = rigidBase.RigidBaseMO

    [positionS, curv_abs_inputS, longeurS, framesF, curv_abs_outputF, frames3D] = \
        BuildCosseratGeometry(config)
    # ################################
    # # Rate of angular Deformation ##
    # ################################
    rateAngularDeformNode = base.addChild('rateAngularDeform')
    rateAngularDeformMO = rateAngularDeformNode.addObject('MechanicalObject',
                                                          template='Vec3d', name='rateAngularDeformMO',
                                                          position=positionS, showIndices=0)
    rateAngularDeformNode.addObject('BeamHookeLawForceField', crossSectionShape=shape,
                                    length=longeurS, youngModulus=Young_modulus, poissonRatio=poissonRatio,
                                    lengthY=length_Y, lengthZ=length_Z)
    # ################################
    # # Define Cosserat Frames ##
    # ################################
    # The node of the frame needs to inherit from rigidBaseMO and rateAngularDeform
    mappedFrameNode = rigidBaseNode.addChild('MappedFrames')
    rateAngularDeformNode.addChild(mappedFrameNode)
    framesMO = mappedFrameNode.addObject('MechanicalObject', template='Rigid3d',
                                         name="FramesMO", position=framesF,
                                         showObject=1, showObjectScale=_showObjectScale)
    mappedFrameNode.addObject('UniformMass', totalMass="0.00022", showAxisSizeFactor='0')
    # The mapping has two inputs: RigidBaseMO and rateAngularDeformMO
    #                 one output: FramesMO
    inputMO = rateAngularDeformMO.getLinkPath()
    inputMO_rigid = RigidBaseMO.getLinkPath()
    outputMO = framesMO.getLinkPath()

    mappedFrameNode.addObject('DiscreteCosseratMapping', curv_abs_input=curv_abs_inputS,
                              curv_abs_output=curv_abs_outputF, input1=inputMO, input2=inputMO_rigid,
                              output=outputMO, debug=0, radius=radius)

    """===================================================="""
    """ =+++>           actuator collision        =+++++<=="""
    """===================================================="""
    if buildCollision:
        tab_edges = buildEdges(frames3D)
        addEdgeCollision(mappedFrameNode, frames3D, tab_edges)

    base.addObject('MechanicalMatrixMapper', template='Vec3,Rigid3', object1=inputMO,
                   object2=inputMO_rigid, nodeToParse=mappedFrameNode.getLinkPath())
    return base


def createCosserat_2(parent, config, config_material, name="Cosserat", orientation=[0, 0, 0, 1], radius=0.01,
                     last_frame=[], rigidBase=None, controlePointMo=None):
    [x, y, z, q0, q1, q2, q3] = config['init_pos']
    tot_length = config['tot_length']
    # buildCollision = config['collision']

    Young_modulus = config_material['youngModulus']
    poissonRatio = config_material['poissonRatio']
    length_Y = config_material['length_Y']
    length_Z = config_material['length_Z']
    shape = config_material['shape']

    nbSectionS = config['nbSectionS']
    lengthS = tot_length / nbSectionS

    nbFramesF = config['nbFramesF']
    lengthF = tot_length / nbFramesF

    #################################
    ##           RigidBase         ##
    #################################
    base = parent.addChild(name)
    # base.addObject('EulerImplicitSolver', firstOrder=0, rayleighStiffness=0.2, rayleighMass=0.1)
    # base.addObject('SparseLDLSolver', name='solver')
    # base.addObject('GenericConstraintCorrection')

    if rigidBase is None:
        rigidBaseNode = base.addChild('rigidBase')

        RigidBaseMO = rigidBaseNode.addObject('MechanicalObject', template='Rigid3d',
                                              name="RigidBaseMO", position=[x, y, z, q0, q1, q2, q3], showObject=0,
                                              showObjectScale=2.)
        # rigidBaseNode.addObject('RestShapeSpringsForceField', name='spring' + str(x), stiffness=5000,
        #                         angularStiffness=5000, external_points=0, mstate="@RigidBaseMO", points=0,
        #                         template="Rigid3d")
        if controlePointMo is None:
            rigidBaseNode.addObject('RestShapeSpringsForceField', name='spring', stiffness="1e1",
                                    angularStiffness="1e1",
                                    external_points="0", mstate="@RigidBaseMO", points="0", template="Rigid3d")
        else:
            rigidBaseNode.addObject('RestShapeSpringsForceField', name='spring', stiffness="1.e6",
                                    angularStiffness="1.e6",
                                    external_rest_shape=controlePointMo, external_points="0",
                                    mstate="@RigidBaseMO", points="0", template="Rigid3d")
    else:
        rigidBaseNode = base.addChild(rigidBase)
        RigidBaseMO = rigidBase.RigidBaseMO

    #################################
    ## Rate of angular Deformation ##
    #################################
    # Define: the length of each beam in a list, the positions of each beam
    # (flexion, torsion), the abs of each section
    positionS = []
    longeurS = []
    sum = x
    curv_abs_inputS = [x]
    # curv_abs_inputS.append(x)

    for i in range(nbSectionS):
        positionS.append([0, 0, 0])
        longeurS.append((((i + 1) * lengthS) - i * lengthS))
        sum += longeurS[i]
        curv_abs_inputS.append(sum)

    curv_abs_inputS[nbSectionS] = tot_length + x

    # Define: sofa elements
    rateAngularDeformNode = base.addChild('rateAngularDeform')
    rateAngularDeformMO = rateAngularDeformNode.addObject('MechanicalObject',
                                                          template='Vec3d', name='rateAngularDeformMO',
                                                          position=positionS, showIndices=0)
    rateAngularDeformNode.addObject('BeamHookeLawForceField', crossSectionShape=shape, length=longeurS,
                                    radius=radius, youngModulus=Young_modulus, poissonRatio=poissonRatio)
    # ################################
    # #             Frame           ##
    # ################################
    # Define: the abs of each frame and the position of each frame.
    framesF = []
    frames3D = []
    curv_abs_outputF = []
    for i in range(nbFramesF):
        sol = i * lengthF
        framesF.append([sol + x, y, z, 0, 0, 0, 1])
        frames3D.append([sol + x, y, z])
        curv_abs_outputF.append(sol + x)

    framesF.append([tot_length + x, y, z, 0, 0, 0, 1])
    frames3D.append([tot_length + x, y, z])
    curv_abs_outputF.append(tot_length + x)

    # framesF = [[x, y, z, 0, 0, 0, 1]] + framesF
    # curv_abs_outputF = [x] + curv_abs_outputF

    # The node of the frame needs to inherit from rigidBaseMO and rateAngularDeform
    mappedFrameNode = rigidBaseNode.addChild('MappedFrames')
    rateAngularDeformNode.addChild(mappedFrameNode)
    framesMO = mappedFrameNode.addObject('MechanicalObject', template='Rigid3d',
                                         name="FramesMO", position=framesF,
                                         showObject=1, showObjectScale=3)
    mappedFrameNode.addObject('UniformMass', totalMass="0.00022", showAxisSizeFactor='0')
    # The mapping has two inputs: RigidBaseMO and rateAngularDeformMO
    #                 one output: FramesMO
    inputMO = rateAngularDeformMO.getLinkPath()
    inputMO_rigid = RigidBaseMO.getLinkPath()
    outputMO = framesMO.getLinkPath()

    mappedFrameNode.addObject('DiscreteCosseratMapping', curv_abs_input=curv_abs_inputS,
                              curv_abs_output=curv_abs_outputF, input1=inputMO, input2=inputMO_rigid,
                              output=outputMO, debug=0, radius=0)

    # for i, orient in enumerate(last_frame):
    #     lastFrame = base.addChild('lastFrame_' + str(i))
    #     lastFrameMo = lastFrame.addObject('MechanicalObject', template='Rigid3d',
    #                                       name="RigidBaseMO", position=[x + tot_length, y, z] + orient, showObject=0,
    #                                       showObjectScale=2.)
    #     lastFrame.addObject("RigidRigidMapping", name="mapLastFrame",
    #                         input=framesMO.getLinkPath(),
    #                         output=lastFrameMo.getLinkPath(), index=nbFramesF + 1,
    #                         globalToLocalCoords=True)
    base.addObject('MechanicalMatrixMapper', template='Vec3,Rigid3', object1=inputMO,
                   object2=inputMO_rigid, nodeToParse=mappedFrameNode.getLinkPath())

    # if buildCollision:
    #     tab_edges = buildEdges(frames3D)
    #     CollisInstrumentCombined = mappedFrameNode.addChild('CollisInstrumentCombined')
    #     CollisInstrumentCombined.addObject('EdgeSetTopologyContainer', name="collisEdgeSet", position=frames3D,
    #                                        edges=tab_edges)
    #     CollisInstrumentCombined.addObject('EdgeSetTopologyModifier', name="colliseEdgeModifier")
    #     CollisInstrumentCombined.addObject('MechanicalObject', name="CollisionDOFs")
    #     CollisInstrumentCombined.addObject('LineCollisionModel', bothSide="1", group='2', proximity="0.01")
    #     CollisInstrumentCombined.addObject('PointCollisionModel', bothSide="1", group='2')
    #     CollisInstrumentCombined.addObject('IdentityMapping', name="mapping")
    return base
