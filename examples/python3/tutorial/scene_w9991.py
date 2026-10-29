import Sofa
from useful.params import BeamPhysicsParametersNoInertia, BeamGeometryParameters, SimulationParameters, Parameters

geoParams = BeamGeometryParameters(beam_length=30.0, nb_section=6, nb_frames=24, build_collision_model=0)
physicsParams = BeamPhysicsParametersNoInertia(beam_mass=1.0, young_modulus=1.0e4, poisson_ratio=0.38,
beam_radius=0.2, beam_length=30.0)
simuParams = SimulationParameters()
Params = Parameters(beam_geo_params=geoParams, beam_physics_params=physicsParams, simu_params=simuParams)

stiffness_param = 1.e10

def _add_rigid_base(p_node):
    rigid_base_node = p_node.addChild('rigid_base')
    rigid_base_node.addObject('MechanicalObject', template="Rigid3d", name="cosserat_base_mo",
    position="0 0 0 0 0 0. 1", showObject=1, showObjectScale="0.1")
    rigid_base_node.addObject('RestShapeSpringsForceField', template="Rigid3d", name="spring",
    stiffness=stiffness_param, angularStiffness=stiffness_param,
    mstate="@cosserat_base_mo", external_points="0", points="0")
    return rigid_base_node


def _add_cosserat_state(p_node, bending_states, list_sections_length, _radius=0.2, _innerRadius=0.05):
    cosserat_coordinate_node = p_node.addChild('cosseratCoordinate')
    cosserat_coordinate_node.addObject('MechanicalObject', template="Vec3d", name="cosserat_state",
    position=bending_states)
    cosserat_coordinate_node.addObject('BeamHookeLawForceField', crossSectionShape="circular", innerRadius=_innerRadius,
    length=list_sections_length, radius=_radius,youngModulus=physicsParams.young_modulus, poissonRatio=physicsParams.poisson_ratio)
    return cosserat_coordinate_node

def _add_cosserat_frame(p_node, _bending_node, nb_framesF, _section_curv_abs, _frame_curv_abs, _radius=0.2, _innerRadius=0.05):
    cosserat_in_Sofa_frame_node = p_node.addChild('cosserat_in_Sofa_frame_node')
    _bending_node.addChild(cosserat_in_Sofa_frame_node)
    frames_mo = cosserat_in_Sofa_frame_node.addObject('MechanicalObject', template="Rigid3d", name="FramesMO",
    position=nb_framesF, showIndices=0, showObject=0, showObjectScale=0.8)
    cosserat_in_Sofa_frame_node.addObject('UniformMass', totalMass=physicsParams.beam_mass)
    cosserat_in_Sofa_frame_node.addObject('DiscreteCosseratMapping', name="cosseratMapping",
    curv_abs_input=_section_curv_abs, curv_abs_output=_frame_curv_abs,
    input1=_bending_node.cosserat_state.getLinkPath(), input2=p_node.cosserat_base_mo.getLinkPath(),
    output=frames_mo.getLinkPath(), debug=0, radius=_radius)
    return cosserat_in_Sofa_frame_node

def createScene(root_node):
    root_node.addObject('RequiredPlugin', name="plugins", pluginName=[
    "Cosserat",                                       # Needed to use components [BeamHookeLawForceField,DiscreteCosseratMapping]
    "Sofa.Component.AnimationLoop",                         # Needed to use components [FreeMotionAnimationLoop]
    "Sofa.Component.Collision.Detection.Algorithm",
    "Sofa.Component.Collision.Detection.Intersection",      # Needed to use components [LocalMinDistance]
    "Sofa.Component.Collision.Geometry",
    "Sofa.Component.Collision.Response.Contact",            # Needed to use components [RuleBasedContactManager]
    "Sofa.Component.Constraint.Lagrangian.Correction",      # Needed to use components [GenericConstraintCorrection]
    "Sofa.Component.Constraint.Lagrangian.Solver",          # Needed to use components [BlockGaussSeidelConstraintSolver]
    "Sofa.Component.Constraint.Projective",                 # Needed to use components [FixedConstraint]
    "Sofa.Component.Engine.Select",
    "Sofa.Component.IO.Mesh",                               # Needed to use components MeshOBJLoader, MeshSTLLoader
    "Sofa.Component.ODESolver.Backward",                    # Needed to use components [EulerImplicitSolver]
    "Sofa.Component.LinearSolver.Direct",                   # Needed to use components [SparseLDLSolver]
    "Sofa.Component.Mapping.Linear",
    "Sofa.Component.Mass",                                  # Needed to use components [UniformMass]
    "Sofa.Component.MechanicalLoad",
    "Sofa.Component.MechanicalLoad",
    "Sofa.Component.ODESolver.Backward",
    "Sofa.Component.Playback",
    "Sofa.Component.Setting",                               # Needed to use components [BackgroundSetting]
    "Sofa.Component.SolidMechanics.FEM.Elastic",
    "Sofa.Component.SolidMechanics.FEM.HyperElastic",
    "Sofa.Component.SolidMechanics.Spring",                 # Needed to use components [RestShapeSpringsForceField]
    "Sofa.Component.StateContainer",                        # Needed to use components [MechanicalObject]
    "Sofa.Component.Topology.Container.Constant",           # Needed to use components [MeshTopology]
    "Sofa.Component.Topology.Container.Dynamic",
    "Sofa.Component.Topology.Container.Grid",               # Needed to use components [RegularGridTopology]
    "Sofa.Component.Topology.Mapping",                      # Needed to use components [Edge2QuadTopologicalMapping]
    "Sofa.Component.Visual",                                 # Needed to use components [VisualStyle]
    "Sofa.GL.Component.Rendering3D",                        # Needed to use components [OglGrid, OglModel]
    "Sofa.GUI.Component"])


    root_node.addObject('VisualStyle', displayFlags="showVisualModels showBehaviorModels showMechanicalMappings")

    root_node.gravity = [0.0, -9.18, 0.0]
    root_node.dt = 0.01

    root_node.addObject('FreeMotionAnimationLoop')
    root_node.addObject('EulerImplicitSolver', firstOrder="0", rayleighStiffness="0.0", rayleighMass="0.0")
    root_node.addObject('BlockGaussSeidelConstraintSolver', maxIterations="500", tolerance="1e-20", computeConstraintForces=1, printLog="0")
    root_node.addObject('SparseLDLSolver', name="solver", template="CompressedRowSparseMatrixd")

    base_node = _add_rigid_base(root_node)

    nb_sections = geoParams.nb_section
    beam_length = geoParams.beam_length
    length_s = beam_length/float(nb_sections)
    bending_states = []
    list_sections_length = []
    temp = 0.0
    section_curv_abs = [0.0]

    for i in range(nb_sections):
        bending_states.append([0, 0.0, 0.0])
        list_sections_length.append((((i+1)*length_s)-i*length_s))
        temp += list_sections_length[i]
        section_curv_abs.append(temp)
    bending_states[nb_sections-1] = [0, 0.2, 0.0]
    bending_states[nb_sections-2] = [0, 0.2, 0.0]
    section_curv_abs[nb_sections] = beam_length

    bending_node = _add_cosserat_state(root_node, bending_states, list_sections_length)
    nb_frames = geoParams.nb_frames
    length_f = beam_length/float(nb_frames)
    cosserat_G_frames = []
    frames_curv_abs = []
    cable_position_f = []
    x, y, z = 0, 0, 0

    for i in range(nb_frames):
        sol = i * length_f
        cosserat_G_frames.append([sol + x, y, z, 0, 0, 0, 1])
        cable_position_f.append([sol + x, y, z])
        frames_curv_abs.append(sol + x)

    cosserat_G_frames.append([beam_length + x, y, z, 0, 0, 0, 1])
    cable_position_f.append([beam_length + x, y, z])
    frames_curv_abs.append(beam_length + x)

    _add_cosserat_frame(base_node, bending_node, cosserat_G_frames, section_curv_abs, frames_curv_abs)



    return root_node
