# -*- coding: utf-8 -*-
import Sofa
from stlib.physics.deformable import ElasticMaterialObject
from stlib.physics.constraints import FixedBox
from stlib.scene import Node
from softrobots.actuators import PullingCable
from stlib.physics.collision import CollisionMesh

from splib.loaders import loadPointListFromFile

class FingerController(Sofa.PythonScriptController):
    def __init__(self, node, cable ):
        self.cableconstraintvalue = cable.getObject("CableConstraint").findData('value')
        self.name = "FingerController"

    def onKeyPressed(self,c):
        if (c == "+"):
            self.cableconstraintvalue.value =  self.cableconstraintvalue.value[0][0] + 1.

        elif (c == "-"):
            displacement = self.cableconstraintvalue.value[0][0] - 1.
            if(displacement < 0):
                displacement = 0
            self.cableconstraintvalue.value = displacement

def Finger(parentNode=None, name="Finger",
           rotation=[0.0, 0.0, 0.0], translation=[0.0, 0.0, 0.0],
           fixingBox=[-5.0,0.0,0.0,10.0,15.0,20.0], pullPointLocation=[0.0,0.0,0.0]):

    finger = Node(parentNode, name)
    eobject = ElasticMaterialObject(finger,
                                   volumeMeshFileName="mesh/finger.vtk",
                                   poissonRatio=0.3,
                                   youngModulus=18000,
                                   totalMass=0.5,
                                   surfaceColor=[0.0, 0.8, 0.7],
                                   surfaceMeshFileName="mesh/finger.stl",
                                   rotation=rotation,
                                   translation=translation)

    FixedBox(eobject.node, atPositions=fixingBox, doVisualization=True)

    cable=PullingCable(eobject.node,
                 "PullingCable",
                 pullPointLocation=pullPointLocation,
                 rotation=rotation,
                 translation=translation,
                 cableGeometry=loadPointListFromFile("mesh/cable.json"));

    FingerController(eobject.node, cable)

    CollisionMesh(eobject.node, name="CollisionMesh",
                 surfaceMeshFileName="mesh/finger.stl",
                 rotation=rotation, translation=translation,
                 collisionGroup=[1, 2])

    CollisionMesh(eobject.node, name="CollisionMeshAuto1",
                 surfaceMeshFileName="mesh/fingerCollision_part1.stl",
                 rotation=rotation, translation=translation,
                 collisionGroup=[1])

    CollisionMesh(eobject.node, name="CollisionMeshAuto2",
                 surfaceMeshFileName="mesh/fingerCollision_part2.stl",
                 rotation=rotation, translation=translation,
                 collisionGroup=[2])


    return finger

def createScene(rootNode):
    from stlib.scene import MainHeader, ContactHeader
    MainHeader(rootNode, gravity=[0.0, -981.0, 0.0], plugins=["SoftRobots"])
    ContactHeader(rootNode, alarmDistance=4, contactDistance=3, frictionCoef=0.08)

    Finger(rootNode, translation=[1.0,0.0,0.0])
    return rootNode

