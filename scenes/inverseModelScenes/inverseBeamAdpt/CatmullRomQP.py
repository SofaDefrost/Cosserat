# -*- coding: utf-8 -*-
"""
Created on Jun 8 2020

@author: PSC
"""

from CatmullRom3d import make_trajectory

import Sofa
from math import sin, cos
import os
import numpy as np

from tools import Quat, createLinePoints, createLines, transformTableInString, Graph, dijkstra
# from controllers import GraphController, QPPathPlanner, keyboardController

path = os.path.dirname(os.path.abspath(__file__))+'/mesh/'

class CRQPController(Sofa.Core.Controller):
    """ This is a custom controller to perform actions when events are triggered """

    def __init__(self, *args, **kwargs):
        # These are needed (and the normal way to override from a python class)
        Sofa.Core.Controller.__init__(self, *args, **kwargs)

        self.init_done = False

        self.mesh = None
        edges = []
        if kwargs["mesh"]:
            self.mesh = kwargs["mesh"]
            for edge in self.mesh.edges.value:
                edges += [(edge[0], edge[1], 1)]

        skel_graph = Graph()
        for edge in edges:
            skel_graph.add_edge(*edge)

        self.startNode = kwargs["startNode"]
        self.goalNode = kwargs["goalNode"]

        self.path = dijkstra(skel_graph, self.startNode, self.goalNode)
        print("shortest path from {} to {} is {}".format(self.startNode, self.goalNode, self.path))

        self.skelMO = None
        if kwargs["skelMO"]:
            self.skelMO = kwargs["skelMO"]

        self.goalMO = None
        if kwargs["goalMO"]:
            self.goalMO = kwargs["goalMO"]

        self.it = 0
        self.max_it = 1000

    def init_path_pos(self):
        self.path_pos = []
        if self.skelMO:
            # print(self.skelMO.position.value)
            for point in self.path:
                self.path_pos += [self.skelMO.position.value[point][:3]]
        P0 = self.path_pos[-2]
        P1 = self.path_pos[-1]
        last_point = P1 - P0 * np.linalg.norm(P1 - P0) * 2
        x, y, z = make_trajectory([[-10.0, 0.0, 0.0]] + self.path_pos + [last_point], nb_points=self.max_it)
        self.x, self.y, self.z = list(x), list(y), list(z)
        # print('len x ', len(self.x))
        # print('len y ', len(self.y))
        # print('len z ', len(self.z))

        # print(self.path_pos)

    def set_goal_pos(self, x, y, z):
        self.goalMO.rest_position.value[0] = [x, y, z]

    def onAnimateBeginEvent(self, dt):
        if not self.init_done:
            self.init_path_pos()
            self.init_done = True

        if self.it < len(self.x) and self.it < len(self.y) and self.it < len(self.z):
            # print(self.it)
            self.set_goal_pos(self.x[self.it], self.y[self.it], self.z[self.it])
        if self.it < self.max_it - 1:
            self.it += 1

# class RestShapeController(Sofa.Core.Controller):
#     """ This is a custom controller to perform actions when events are triggered """
#
#     def __init__(self, *args, **kwargs):
#         # These are needed (and the normal way to override from a python class)
#         Sofa.Core.Controller.__init__(self, *args, **kwargs)
#
#         self.mo = kwargs["MO"]
#
#     def onAnimateBeginEvent(self, dt):
#         rate = 10.0
#         with self.mo.rest_position.writeable() as rest_pos:
#             # print("rest pos x0 ", rest_pos[0][0])
#             # print("pos x0 ", self.mo.position[0][0])
#             if (abs(rest_pos[0][0]-self.mo.position[0][0]) > rate):
#                 for point in rest_pos:
#                     # print("point before ", point[0])
#                     point[0] += rate
#                     # print("point after ", point[0])

class RestShapeController(Sofa.Core.Controller):
    """ This is a custom controller to perform actions when events are triggered """

    def __init__(self, *args, **kwargs):
        # These are needed (and the normal way to override from a python class)
        Sofa.Core.Controller.__init__(self, *args, **kwargs)

        self.mo = kwargs["MO"]
        self.RSSFFMO = kwargs["RSSFF_MO"]

    def onAnimateBeginEvent(self, dt):
        rate = 10.0
        with self.mo.rest_position.writeable() as rest_pos:
            if (self.mo.position[0][0] - rest_pos[0][0] > rate and rest_pos[0][0] < -rate):
                # print("dist ", self.mo.position[0][0] - rest_pos[0][0])
                # print(rest_pos)
                with self.RSSFFMO.position.writeable() as pos:
                    pos[0] = self.mo.position[0]
                for point in rest_pos:
                    point[0] += rate//2



            elif (rest_pos[0][0] - self.mo.position[0][0] > rate ):
                for point in rest_pos:
                    point[0] -= 3*rate
                with self.RSSFFMO.position.writeable() as pos:
                    pos[0] = self.mo.rest_position[0]


class SpringStiffnessController(Sofa.Core.Controller):
    """ This is a custom controller to perform actions when events are triggered """

    def __init__(self, *args, **kwargs):
        # These are needed (and the normal way to override from a python class)
        Sofa.Core.Controller.__init__(self, *args, **kwargs)

        self.mo = kwargs["MO"]
        self.RSSFFMO = kwargs["RSSFF_MO"]
        self.stiffness = self.RSSFFMO.stiffness.value


    def onAnimateBeginEvent(self, dt):
        self.stiffness *= 0.9992

        # with self.RSSFFMO.stiffness.value.writeable() as stiffness:
        self.RSSFFMO.stiffness.value[0] = self.stiffness




class ActuatorIndexController(Sofa.Core.Controller):
    """ This is a custom controller to perform actions when events are triggered """

    def __init__(self, *args, **kwargs):
        # These are needed (and the normal way to override from a python class)
        Sofa.Core.Controller.__init__(self, *args, **kwargs)

        self.index = kwargs["n_beams"] - 1
        self.trans_actuator = kwargs["trans_actuator"]
        self.rot_actuator = kwargs["rot_actuator"]
        self.actuatorMO = kwargs["actuatorMO"]
        self.constraint = kwargs["constraint"]
        self.constraint_spring = kwargs["constraint_spring"]

        with self.trans_actuator.indices.writeable() as t_index_list:
            t_index_list[0] = self.index
        with self.rot_actuator.indices.writeable() as r_index_list:
            r_index_list[0] = self.index


    def onAnimateBeginEvent(self, dt):
        with self.trans_actuator.indices.writeable() as t_index_list:
            print("t index", t_index_list[0])
            print("full act pos ", self.actuatorMO.position.value)
            print("corresponding pos ", self.actuatorMO.position.value[t_index_list[0]])
            if self.actuatorMO.position[t_index_list[0]][0]>-30:
                if self.index > 0:
                    self.index -= 1
                with self.constraint.indices.writeable() as c_index_list:
                    c_index_list[0] = self.index
                self.constraint_spring.getData("points").value = [i for i in range(self.index//2)]
                # print(self.constraint_spring.getData("points").value)
                # print("changing actuation point for ", self.index)
                t_index_list[0] = self.index


class CostController(Sofa.Core.Controller):
    def __init__(self, *args, **kwargs):
        # These are needed (and the normal way to override from a python class)
        Sofa.Core.Controller.__init__(self, *args, **kwargs)

        self.solver = kwargs["solver"]
        self.goal_pos = kwargs["goal_pos"]
        self.effMO = None
        if not self.effMO:
            self.effMO = kwargs["effMO"]
        self.cost = np.inf
        self.best_cost = np.inf
        self.qp_error = False


    def onAnimateBeginEvent(self, dt):

        qp_error_messages = self.solver.getLoggedMessagesAsString(4)
        if qp_error_messages and qp_error_messages != "":
            self.qp_error = True
            print(">>>>>>>>>>>>>>>>>>>>>>> Caught QP Error ", qp_error_messages)

        if not self.qp_error:
            self.cost = get_distance_to_goal(self.effMO.position[0], self.goal_pos)
            if self.cost < self.best_cost:
                self.best_cost = self.cost


def make_strings(start, end, interval=5, h=3.5, teta=3.14/8):
    n_beams = int(abs(end-start)/interval)
    if start + n_beams * interval != end:
        raise ValueError
    pos_str = ""
    mech_str = ""
    mech_rest_str = ""
    collis_str = ""
    line_str = ""
    for i in range(n_beams+1):
        pos_str += str(start + i*interval)
        pos_str += " 0 0  "
        mech_str += str(start + i*interval)
        mech_str += " 0 0 0 0 0 1  "
        collis_str += str(i*interval)
        collis_str += " 0 0  "

    for i in range(n_beams):
        mech_rest_str += str(start + i*interval)
        mech_rest_str += " 0 0 0 0 0 1  "
        beam = str(i) + ' ' + str(i+1) + '  '
        line_str += beam
    str_last_mech_rest = str(start + (n_beams)*interval - interval/5) + ' ' + str(h) + ' 0 0 0 ' + str(sin(teta)) + ' ' + str(cos(teta))
    mech_rest_str += str_last_mech_rest

    return pos_str, mech_str, mech_rest_str, collis_str, line_str, n_beams


def get_distance_to_goal(kt_endpoint, goal_pos):
    return np.linalg.norm(goal_pos[:3]-kt_endpoint[:3])


configs = [
    {
        "start_node": 3,
        "goal_node": 9,
        "skel": "skel_straight_simplefourche_g2.obj",
        "mesh": "mesh_straight_simplefourche_g2.obj",
        "scale": 10,
        "rotation": [-3.0, 0.0, -100.0]
    },
    {
        "start_node": 7,
        "goal_node": 15,
        "skel": "skel_reversestraight_doublefourche_g4.obj",
        "mesh": "mesh_reversestraight_doublefourche_g4.obj",
        "scale": 10,
        "rotation": [-3.0, 0.0, -100.0]
    }
]

config = configs[0]


def createScene(root, startNode=config["start_node"], goalNode=config["goal_node"], h=3.5, teta=3.14/8):

    root.addObject('RequiredPlugin', pluginName='BeamAdapter')
    root.addObject('RequiredPlugin', pluginName='SoftRobots')
    root.addObject('RequiredPlugin', pluginName='SoftRobots.Inverse')
    root.addObject('RequiredPlugin', pluginName='SofaSparseSolver')

    root.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels showCollisionModels hideMappings '
                                               'showForceFields')

    root.addObject('FreeMotionAnimationLoop')
    root.addObject('CollisionPipeline', verbose=False)
    root.addObject('BruteForceDetection', name="N2")
    root.addObject('CollisionResponse', response="FrictionContact", responseParams="mu=0.")
    root.addObject('LocalMinDistance', name="Proximity", alarmDistance=1, contactDistance=0.1)

    qp_solver = root.addObject('QPInverseProblemSolver', epsilon=0.0, printLog=False, displayTime=0, tolerance=1e-10,
                               maxIterations=10000)

    root.gravity = [0.0, 0.0, 0.0]
    root.dt = 0.01

    #########################################
    # goal                                  #
    #########################################
    goal = root.addChild('goal')
    goal.addObject('VisualStyle', displayFlags="showCollisionModels")
    goal.addObject('EulerImplicit', firstOrder=True)
    goal.addObject('CGLinearSolver', iterations=100, tolerance=1e-5, threshold=1e-4)
    goalMO = goal.addObject('MechanicalObject', name='goalMO', showObject=True, drawMode=1, showObjectScale=1,
                            position=[0.0, 0.0, 0.0])
    goal.addObject('RestShapeSpringsForceField', points=0, angularStiffness=1000, stiffness=1000)
    goal.addObject('UncoupledConstraintCorrection')

    ##########################################
    # Beam Model                             #
    ##########################################
    x_begin = -250
    pos_str, mech_str, mech_rest_str, collis_str, line_str, n_beams = make_strings(x_begin, 5, h=h, teta=teta)

    model = root.addChild('model')
    model.addObject('EulerImplicit', name='odesolver', firstOrder=False, rayleighStiffness=0.1, rayleighMass=0.1)
    model.addObject('SparseLDLSolver', name='ldlsolveur')
    model.addObject('GenericConstraintCorrection', solverName='ldlsolveur')
    model.addObject('Mesh', position=pos_str, lines=line_str)
    beam_mo = model.addObject('MechanicalObject', template='Rigid3d', name='frame1', position=mech_str,
                              rest_position=mech_rest_str)

    model.addObject('BeamInterpolation', dofsAndBeamsAligned=True, straight=False, defaultYoungModulus=100, radius=0.5,
                    name="interpolation")
    model.addObject('AdaptiveBeamForceFieldAndMass',  name="BeamForceField", computeMass=True, massDensity=0.000001)
    constraint = model.addObject('PartialFixedConstraint', indices=0, fixedDirections="0 1 1 0 1 1")
    constraint_spring = model.addObject('RestShapeSpringsForceField', points=0, angularStiffness=10, stiffness=10,
                                        drawSpring=True)

    collision = model.addChild('Collis')

    collision.addObject('Mesh', name='lineMesh', position=collis_str, lines=line_str)
    collision.addObject('MechanicalObject', template='Vec3')
    collision.addObject('Line', group=2)
    collision.addObject('Point', group=2)
    collision.addObject('AdaptiveBeamMapping', name='mapping', mapForces=False, mapMasses=False)

    ##########################################
    # Actuator                               #
    ##########################################

    actuator = model.addChild('actuator')
    actuator.addObject('VisualStyle', displayFlags="showInteractionForceFields")
    actuatorMO = actuator.addObject('MechanicalObject', template='Rigid3d', position='@../frame1.position')
    translation = actuator.addObject('SlidingActuator', template='Rigid3d', indices=n_beams-1, direction="1 0 0 0 0 0",
                        maxNegativeDisp=-2000, maxPositiveDisp=5000, maxDispVariation=0.5, maxForce=100, minForce=-100)
    rotation = actuator.addObject('SlidingActuator', template='Rigid3d', indices=n_beams-1, direction="0 0 0 1 0 0",
                        maxNegativeDisp=-20, maxPositiveDisp=20, maxDispVariation=0.1, maxForce=100, minForce=-100)
    actuator.addObject('AdaptiveBeamMapping', name='mapping', mapForces=False, mapMasses=False)

    ##########################################
    # Effector Model                         #
    ##########################################
    effector = model.addChild('Effector')
    effMO = effector.addObject('MechanicalObject', position=[0.0, 0.0, 0.0])
    effector.addObject('PositionEffector', template='Vec3d', indices=0, effectorGoal="@../../goal/goalMO.position",
                       useDirections='1 1 1')
    effector.addObject('RigidMapping', index=n_beams)

    # ##########################################
    # # Vessel Model                           #
    # ##########################################

    skel = root.addChild("Skeleton")
    skel.addObject('MeshObjLoader', filename=path+config["skel"], flipNormals=True,
                   triangulate=True, name='skelLoader', scale=10, rotation=[-3.0, 0.0, -100.0])
    mesh = skel.addObject('Mesh', position='@skelLoader.position', edges='@skelLoader.edges')
    skelMO = skel.addObject('MechanicalObject', template='Rigid3d', name='frame1', showObject=True, drawMode="1",
                            showObjectScale=2, showColor=[1, 0, 1, 1], position=mesh.position,
                            rest_position=mesh.position)
    skel.addObject('BeamInterpolation', dofsAndBeamsAligned=True, straight=False, defaultYoungModulus=100, radius=0.1,
                   name="interpolation")
    skel.addObject('PartialFixedConstraint', indices=0, fixedDirections="1 1 1 1 1 1")

    ##########################################
    # Obstacle Model                         #
    ##########################################

    obstacle = root.addChild('Obstacle')
    obstacle.addObject('MeshObjLoader', name="loader", filename=path+config["mesh"], scale=config["scale"],
                       rotation=config["rotation"], translation=[0.0, 0.0, 0.0])

    obstacle.addObject('Mesh', src='@loader', name="mesh")
    obstacle.addObject('MechanicalObject', template="Vec3")
    obstacle.addObject('Triangle', group=1)
    obstacle.addObject('Line', group=1)
    obstacle.addObject('Point', group=1)

    root.addObject(CRQPController(name="CRQP", mesh=mesh, skelMO=skelMO, goalMO=goalMO, startNode=startNode,
                                  goalNode=goalNode))

    root.addObject(SpringStiffnessController(name="Spring", MO=beam_mo, RSSFF_MO=constraint_spring))

    skel.init()
    root.addObject(CostController(name="CostController", goal_pos=skel.skelLoader.position.toList()[goalNode],
                                  effMO=effMO, solver=qp_solver))

    # root.addObject(ActuatorIndexController(name="ActuatorIndexController", n_beams=n_beams, trans_actuator=translation, rot_actuator=rotation,
    #             actuatorMO=actuatorMO, constraint=constraint, constraint_spring=constraint_spring))

    return root
