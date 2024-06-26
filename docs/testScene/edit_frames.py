# -*- coding: utf-8 -*-
"""
    Cosserat class in SofaPython3.
"""

__authors__ = "adagolodjo"
__contact__ = "adagolodjo@protonmail.com"
__version__ = "1.0.0"
__copyright__ = "(c) 2021,Inria"
__date__ = "October, 26 2021"

from useful.header import addHeader, addSolverNode
from useful.params import BeamPhysicsParameters, BeamGeometryParameters, SimulationParameters
from useful.params import Parameters
from cosserat.CosseratBase import CosseratBase
from math import sqrt
from splib3.numerics import Quat
import Sofa
from math import pi

geoParams = BeamGeometryParameters(init_pos=[0., 0., 0.], beamLength=30., showFramesObject=1,
                                   nbSection=6, nbFrames=12, buildCollisionModel=0)
physicsParams = BeamPhysicsParameters(beamMass=1., youngModulus=1.0e4, poissonRatio=0.38, beamRadius=0.5,
                                      beamLength=30)
simuParams = SimulationParameters()
Params = Parameters(beamGeoParams=geoParams, beamPhysicsParams=physicsParams, simuParams=simuParams)

force_null = [0., 0., 0., 0., 0., 0.]  # N

class ForceController(Sofa.Core.Controller):
    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        # self.forceNode = kwargs['forceNode']
        self.frames = kwargs['frame_node'].FramesMO
        #self.edite_frames()
        self.edit_pos = True

    def onAnimateEndEvent(self, event):
        if self.edit_pos == True:
            self.edite_frames()
            self.edit_pos = False 
        
        

    def edite_frames(self):
        with self.frames.rest_position.writeable() as position:
            last_frame = len(position)
            print(f'===> The position : {position[last_frame-1]}')
            position[last_frame-1][1] = 2
            
                

    def compute_orthogonal_force(self):
        position = self.frames.position[self.size]  # get the last rigid of the cosserat frame
        orientation = Quat(position[3], position[4], position[5], position[6])  # get the orientation
        # Calculate the direction of the force in order to remain orthogonal to the x_axis
        # of the last frame of the beam.
        with self.forceNode.forces.writeable() as force:
            vec = orientation.rotate([0., self.forceCoeff * 5.e-2, 0.])
            # vec.normalize()
            # print(f' The new vec is : {vec}')
            for count in range(3):
                force[0][count] = vec[count]

    def rotate_force(self):
        if self.forceCoeff <= 100. * pi:
            with self.tip_controller.position.writeable() as position:
                last_frame = self.frames.position[self.size]
                vec = Quat(last_frame[3], last_frame[4], last_frame[5], last_frame[6])  # get the orientation

                vec.rotateFromEuler([self.theta, 0., 0.])  # apply rotation arround x-axis 
                vec.normalize()
                for i, v in enumerate(vec):
                    position[0][i + 3] = v

    def onKeypressedEvent(self, event):
        key = event['key']
        if key == "+":
            self.applyForce = True
        elif key == "-":
            self.applyForce = False


controller_type = 2

def createScene(root_node):
    addHeader(root_node)
    root_node.gravity = [0, 0., 0.]

    solver_node = addSolverNode(root_node, name="solver_node")

    # create cosserat Beam
    cosserat_beam = solver_node.addChild(CosseratBase(parent=solver_node, params=Params))
    cosserat_frames = cosserat_beam.cosseratFrame
    solver_node.addObject(ForceController(frame_node=cosserat_frames))

    

    # solver_node.addObject(ForceController(forceNode=const_force_node, frame_node=cosserat_frames, force_type=controller_type, tip_controller=controller_state))

    return root_node
