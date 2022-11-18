import Sofa
import Sofa.Core
import numpy as np
from cosserat.cosseratObject import Cosserat, cosserat_config


def createScene(root):
    root.gravity = [0, -9.81, 0]

    root.addObject("DefaultAnimationLoop")
    root.addObject("DefaultVisualManagerLoop")
    root.addObject("RequiredPlugin", name="Sofa.Component.Topology.Container.Dynamic")
    root.addObject('VisualStyle', displayFlags='showMechanicalMappings')
    needle = root.addChild(
        Cosserat(parent=root, cosseratGeometry=cosserat_config, name="needle", radius=0.15))
    constraintPointsNode=needle.addChild("constraintPointsNode")

    slidingPoint = needle.addSlidingPoints()
    pathToSlidingMo=slidingPoint.getLinkPath()+str("/slidingPointMO")

    print(f'------> The save path is : {pathToSlidingMo}')

    container = constraintPointsNode.addObject("PointSetTopologyContainer", points=[])
    modifier = constraintPointsNode.addObject("PointSetTopologyModifier")
    state = constraintPointsNode.addObject("MechanicalObject", template="Vec3d", showObject=True, showObjectScale=10)
    pointManager=constraintPointsNode.addObject('PointsManager', name="pointsManager", listening="1",
                                   beamPath="/needle/rigidBase/cosseratInSofaFrameNode/slidingPoint/slidingPointMO")

    root.addObject(PointController(pointManager=pointManager, state=state, baseState=needle))


class PointController(Sofa.Core.Controller):
    def __init__(self, pointManager, state, baseState):
        super().__init__()
        self.pointManager = pointManager
        self.state = state
        self.baseState=baseState


    def onKeypressedEvent(self, event):
        if event["key"] == "M":
            print("Add 10 points")
            print(f'=====> {dir(self.pointManager)}')
            print(f'=====> {dir(self.state)}')
            self.pointManager.addNewPointToState()
            # print("Before setting positions", self.state.position.array())

        elif event["key"] == "D":
            print("Remove a point")
            self.pointManager.removeLastPointfromState()

        elif event["key"] == "B":
            print(f"{len(self.state.position.array())=}")
            with self.state.position.writeable() as state:
                print(f"{len(state)=}")
                # for i in range(len(state)):
                #     state[i] = np.array([i, 0, 0])

            # print("After setting positions", self.state.position.array())import Sofa
