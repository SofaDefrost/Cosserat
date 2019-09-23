import Sofa

from math import *

# from splib.numerics import Quat

import os
path = 'mesh/'

class ConstantForceController(Sofa.PythonScriptController):
    def initGraph(self, node):
        self.time = 0
        self.node = node
        self.constanteForce = self.node.getObject('constanteForce')

    def onBeginAnimationStep(self, dt):
        self.time = self.time + dt
        forces = self.constanteForce.forces
        forces[0][1] = -0.0008 * (sin(self.time))
        self.constanteForce.forces = forces


class ConstantForceController1(Sofa.PythonScriptController):
    def initGraph(self, node):
        self.time = 0
        self.node = node
        self.constanteForce = self.node.getObject('constanteForce')

    def onBeginAnimationStep(self, dt):
        self.time = self.time + dt
        forces = self.constanteForce.forces
        forces[0][1] = 0.0008 * (sin(self.time))
        self.constanteForce.forces = forces


class Animation(Sofa.PythonScriptController):

    def __init__(self, ktNode, beam_0_node, beam_1_node, beam_2_node,beam_3_node):
        self.ktNode = ktNode
        self.beam_0_node = beam_0_node
        self.beam_1_node = beam_1_node
        self.beam_2_node = beam_2_node
        self.beam_3_node = beam_3_node
        return;

    def initGraph(self, nodeRigid):
        self.time = 0
        self.mstate = self.ktNode.getObject('mstate')
        self.mstate_beam0 = self.beam_0_node.getObject('mstate_beam0')
        self.mstate_beam1 = self.beam_1_node.getObject('mstate_beam1')
        self.mstate_beam2 = self.beam_2_node.getObject('mstate_beam2')
        self.mstate_beam3 = self.beam_3_node.getObject('mstate_beam3')

    def onBeginAnimationStep(self, dt):
        self.time = self.time + dt;
        rest_position = self.mstate.rest_position

        pos0 = self.mstate_beam0.findData('position').value;
        rest_position[1] = pos0[2];
        pos1 = self.mstate_beam1.findData('position').value;
        rest_position[2] = pos1[2];
        pos2 = self.mstate_beam2.findData('position').value;
        rest_position[3] = pos2[2];
        pos3 = self.mstate_beam3.findData('position').value;
        rest_position[4] = pos3[2];

        rest_position[0] = [ 0, sin(self.time),  0, 0, 0, 0, 1]
        self.mstate.findData('rest_position').value = rest_position;

        # self.MechanicalState.rest_position = rest_position

pos_str = '0 0 0  20 0 0  40 0 0  60 0 0  80 0 0'
pos_quad_str = '0 0 0 0 0 0 1  20 0 0 0 0 0 1  40 0 0 0 0 0 1   60 0 0 0 0 0 1  80 0 0 0 0 0 1'
line_str = ' 0 1  1 2  2 3  3 4 '


def createScene(rootNode):

    rootNode.createObject(
        'RequiredPlugin', pluginName='SoftRobots SofaPython SofaSparseSolver SofaPreconditioner SofaOpenglVisual CosseratPlugin BeamAdapter')

    rootNode.createObject(
        'VisualStyle', displayFlags='hideVisualModels showBehaviorModels showCollisionModels hideBoundingCollisionModels showForceFields hideInteractionForceFields showWireframe')

    rootNode.findData('dt').value = 0.01
    rootNode.findData('gravity').value = '0 0 0'

    rootNode.createObject('FreeMotionAnimationLoop')
    rootNode.createObject('GenericConstraintSolver',tolerance="1e-6", maxIterations="5000", unbuilt="false")


    beam0 = rootNode.createChild('beam0')
    beam0.createObject(
        'VisualStyle', displayFlags='hideBehaviorModels hideCollisionModels  hideForceFields')
    beam0.createObject('EulerImplicitSolver')
    beam0.createObject('BTDLinearSolver', name="btdsolver")

    beam0.createObject('Mesh', position='0 0 0  10 0 0  20 0 0 ', lines='0 1 1 2')
    beam0.createObject('MechanicalObject', template='Rigid3d', name='mstate_beam0',
                  position='0 0 0 0 0 0 1  10 0 0 0 0 0 1  20 0 0 0 0 0 1', showObject='0')
    beam0.createObject('BeamInterpolation', dofsAndBeamsAligned='1', straight='0',
                  defaultYoungModulus="1000000", radius="0.03", name="interpolation")
    beam0.createObject('ConstantForceField',  name="constanteForce",
                  indices="2", forces="0 -0.0002 0 0 0 0")
    beam0.createObject('AdaptiveBeamForceFieldAndMass',
                  name="BeamForceField", computeMass="1", massDensity="0.0001")

    beam0.createObject('RestShapeSpringsForceField', points="0",
                  template="Rigid3d", stiffness=1e7, angularStiffness=1e7, external_points=0)
    # rigidBaseNode.createObject('RestShapeSpringsForceField', name='spring', stiffness="50000", angularStiffness="50000", external_points="0", mstate="@RigidBaseMO", points="0", template="Rigid3d")

    beam0.createObject('GenericConstraintCorrection', solverName="btdsolver")
    beam0.createObject('PythonScriptController',classname="ConstantForceController")

    beam1 = rootNode.createChild('beam1')
    beam1.createObject(
        'VisualStyle', displayFlags='hideBehaviorModels hideCollisionModels  hideForceFields')
    beam1.createObject('EulerImplicitSolver')
    beam1.createObject('BTDLinearSolver', name="btdsolver")

    beam1.createObject('Mesh', position='20 0 0  30 0 0  40 0 0 ', lines='0 1 1 2')
    beam1.createObject('MechanicalObject', template='Rigid3d', name='mstate_beam1',
                  position='20 0 0 0 0 0 1  30 0 0 0 0 0 1  40 0 0 0 0 0 1', showObject='0')
    beam1.createObject('BeamInterpolation', dofsAndBeamsAligned='1', straight='0',
                  defaultYoungModulus="1000000", radius="0.03", name="interpolation")
    beam1.createObject('ConstantForceField',  name="constanteForce",
                  indices="2", forces="0 -0.0002 0 0 0 0")
    beam1.createObject('AdaptiveBeamForceFieldAndMass',
                  name="BeamForceField", computeMass="1", massDensity="0.0001")

    beam1.createObject('RestShapeSpringsForceField', points="0",
                  template="Rigid3d", stiffness=1e7, angularStiffness=1e7, external_points=0)
    # rigidBaseNode.createObject('RestShapeSpringsForceField', name='spring', stiffness="50000", angularStiffness="50000", external_points="0", mstate="@RigidBaseMO", points="0", template="Rigid3d")

    beam1.createObject('GenericConstraintCorrection', solverName="btdsolver")
    beam1.createObject('PythonScriptController', classname="ConstantForceController1")


    beam2 = rootNode.createChild('beam2')
    beam2.createObject(
        'VisualStyle', displayFlags='hideBehaviorModels hideCollisionModels  hideForceFields')
    beam2.createObject('EulerImplicitSolver')
    beam2.createObject('BTDLinearSolver', name="btdsolver")

    beam2.createObject('Mesh', position='40 0 0  50 0 0  60 0 0 ', lines='0 1 1 2')
    beam2.createObject('MechanicalObject', template='Rigid3d', name='mstate_beam2',
                  position='40 0 0 0 0 0 1  50 0 0 0 0 0 1  60 0 0 0 0 0 1', showObject='0')
    beam2.createObject('BeamInterpolation', dofsAndBeamsAligned='1', straight='0',
                  defaultYoungModulus="1000000", radius="0.03", name="interpolation")
    beam2.createObject('ConstantForceField',  name="constanteForce",
                  indices="2", forces="0 -0.0002 0 0 0 0")
    beam2.createObject('AdaptiveBeamForceFieldAndMass',
                  name="BeamForceField", computeMass="1", massDensity="0.0001")

    beam2.createObject('RestShapeSpringsForceField', points="0",
                  template="Rigid3d", stiffness=1e7, angularStiffness=1e7, external_points=0)
    # rigidBaseNode.createObject('RestShapeSpringsForceField', name='spring', stiffness="50000", angularStiffness="50000", external_points="0", mstate="@RigidBaseMO", points="0", template="Rigid3d")

    beam2.createObject('GenericConstraintCorrection', solverName="btdsolver")
    beam2.createObject('PythonScriptController', classname="ConstantForceController")

    beam3 = rootNode.createChild('beam3')
    beam3.createObject(
        'VisualStyle', displayFlags='hideBehaviorModels hideCollisionModels  hideForceFields')
    beam3.createObject('EulerImplicitSolver')
    beam3.createObject('BTDLinearSolver', name="btdsolver")

    beam3.createObject('Mesh', position='60 0 0  70 0 0  80 0 0 ', lines='0 1 1 2')
    beam3.createObject('MechanicalObject', template='Rigid3d', name='mstate_beam3',
                  position='60 0 0 0 0 0 1  70 0 0 0 0 0 1  80 0 0 0 0 0 1', showObject='0')
    beam3.createObject('BeamInterpolation', dofsAndBeamsAligned='1', straight='0',
                  defaultYoungModulus="1000000", radius="0.03", name="interpolation")
    beam3.createObject('ConstantForceField',  name="constanteForce",
                  indices="2", forces="0 -0.0002 0 0 0 0")
    beam3.createObject('AdaptiveBeamForceFieldAndMass',
                  name="BeamForceField", computeMass="1", massDensity="0.0001")

    beam3.createObject('RestShapeSpringsForceField', points="0",
                  template="Rigid3d", stiffness=1e7, angularStiffness=1e7, external_points=0)
    # rigidBaseNode.createObject('RestShapeSpringsForceField', name='spring', stiffness="50000", angularStiffness="50000", external_points="0", mstate="@RigidBaseMO", points="0", template="Rigid3d")

    beam3.createObject('GenericConstraintCorrection', solverName="btdsolver")
    beam3.createObject('PythonScriptController', classname="ConstantForceController1")

    solverNode = rootNode.createChild('solverNode')

    kt = solverNode.createChild('KTModel')
    kt.createObject('EulerImplicitSolver')
    kt.createObject('BTDLinearSolver', name="btdsolver")

    kt.createObject('Mesh', position=pos_str, lines=line_str)
    kt.createObject('MechanicalObject', template='Rigid3d', name='mstate',
                    position=pos_quad_str, showObject='1', showIndices="1", showIndicesScale="0.2")
    kt.createObject('BeamInterpolation', dofsAndBeamsAligned='1', straight='0',
                    defaultYoungModulus="1000000", radius="0.03", name="interpolation")

    kt.createObject('AdaptiveBeamForceFieldAndMass',
                    name="BeamForceField", computeMass="1", massDensity="0.0001")

    kt.createObject('RestShapeSpringsForceField', points="0", external_rest_shape="@../../Point/mo",
                    template="Rigid3d", stiffness=1e7, angularStiffness=1e7, external_points=0)
    # rigidBaseNode.createObject('RestShapeSpringsForceField', name='spring', stiffness="50000", angularStiffness="50000", external_points="0", mstate="@RigidBaseMO", points="0", template="Rigid3d")

    kt.createObject('LinearSolverConstraintCorrection', solverName="btdsolver")

    Animation(kt, beam0, beam1, beam2, beam3)

    # kt.createObject('PythonScriptController', classname="Animation", beam0, beam1, beam2, beam3)

    # Collision Model
    # collision = kt.createChild('Collis')
    # collision.createObject('Mesh', name='lineMesh', position = pos_str, lines = line_str)
    # collision.createObject('MechanicalObject' , template='Vec3', position= pos_str)
    # collision.createObject('Point', group='2')
    # collision.createObject('Line', group='2')
    # collision.createObject('AdaptiveBeamMapping', name='mapping', mapForces='false', mapMasses='false')

    return rootNode
