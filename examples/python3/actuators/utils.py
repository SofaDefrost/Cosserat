# -*- coding: utf-8 -*-
from splib3.numerics.quat import Quat
import os
import Sofa
import Sofa.constants.Key as Key
path = f'{os.path.dirname(os.path.abspath(__file__))}/mesh/'


def loadDisk(parentNode, filename):
    # ### MAPPING of the disks (here some torus) ####    ###########
    diskMapping = parentNode.addChild("diskMapping")
    diskMapping.addObject("MeshOBJLoader", name="diskLoader", filename=filename)
    diskMapping.addObject("OglModel", name="Visual", src="@diskLoader", color="1.0 0.5 0.25 1.0")
    diskMapping.addObject('RigidMapping', name="diskMap")


def createIntermediateNode(parent, rigidCentral=None, baseName="rigidState", rigidIndex=0, filename=None):
    """ The intermediateRigid construction """
    if rigidCentral is None:
        rigidCentral = [0, 0., 0, 0, 0, 0, 1]
    q0 = Quat()
    q1 = Quat()
    q2 = Quat()

    interChildPos = [[0, 1.7, 1., q0[0], q0[1], q0[2], q0[3]],
                     [0, -1.7, 1., q1[0], q1[1], q1[2], q1[3]],
                     [0., 0., -2., q2[0], q2[1], q2[2], q2[3]]]
    intermediateRigid = parent.addChild('intermediateRigid')
    intermediateRigid.addObject('MechanicalObject', name=baseName, template='Rigid3d',
                                position=rigidCentral, showObject=True, showObjectScale=0.5)
    # @TODO add the mass of the disk
    intermediateRigid.addObject('RigidMapping', name="interRigidMap", index=rigidIndex)

    # """Create intermediate kids"""
    # interRigidChild = intermediateRigid.addChild('interRigidChild')
    # interRigidChild.addObject('MechanicalObject', name="interRigidChildMo", template='Rigid3d',
    #                           position=interChildPos, showObject=True, showObjectScale=0.5)
    # interRigidChild.addObject('RigidMapping', name="interRigidMap")

    """ Add disk to the scene, actually this disk is only use for the visualisation """
    loadDisk(intermediateRigid, filename=filename)

    return intermediateRigid


def createRigidDisk(parentNode):  # @todo add a parameter for the number of disk
    """ Create the rigid disk nodes """
    for i in range(1, 15):
        """ Create the rigid disk nodes """
        createIntermediateNode(parentNode, rigidCentral=[0, 0., 0, 0, 0, 0, 1],
                               baseName=f"rigidState{i}", rigidIndex=i, filename=f'{path}disqueBlender2.obj')


def create_cable_points(geoParams=None, num_segments=14):
    """
    This function creates cable points based on the given parameters.

    Returns:
        cable_base_state (list): A list of base states for the cable.
        cable_3d_state (list): A list of 3D positions for the cable.
    """

    norm_length = geoParams.beamLength / num_segments

    cable_base_state = []
    cable_3d_state = []

    # for z in [1.7, -1.7, -2]:
    for z in [1.7]:
        base_state = []
        pos_3d = []

        for i in range(num_segments + 1):
            x = i * norm_length
            if z == -2:
                base_state.append([x, 0., z, 0., 0., 0., 1.])
                # pos_3d.append([x, 0., z])
                pos_3d.append([x, 0., z])
            else:
                base_state.append([x, z, 1., 0., 0., 0., 1.])
                # pos_3d.append([x, z, 1.])
                pos_3d.append([x, z, 1])
        cable_base_state.append(base_state)
        cable_3d_state.append(pos_3d)
    return cable_base_state, cable_3d_state


class FingerController(Sofa.Core.Controller):
    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, args, kwargs)
        self.cable = args[0]
        self.name = "FingerController"

    def onKeypressedEvent(self, event):
        displacement = self.cable.CableConstraint.value[0]
        print(f' The displacement is : {displacement}')
        if event["key"] == Key.plus:
            displacement += 0.01

        elif event["key"] == Key.minus:
            displacement -= 0.01
            if displacement < 0:
                displacement = 0
        self.cable.CableConstraint.value = [displacement]