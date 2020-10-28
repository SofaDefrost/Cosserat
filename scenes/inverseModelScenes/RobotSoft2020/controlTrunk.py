from math import cos, sin
import Sofa

class Animation(Sofa.PythonScriptController):
    def __init__(self, targetNode):
        self.target = targetNode        
        self.rate = 1.
        self.time = 0.
        self.radius = 193.
        self.theta = 0.01
        self.dt = 0.01
        self.step = 1
        self.maxTheta = 0.54
        self.minTheta = 0.01
        self.coeff = 0.9
        self.threshold = 0.3
        return

    def initGraph(self, targetNode):
        self.targetMO = self.target.getObject('dofs')
    
    def moveF(self, axes=[0], dt=0.01):
        pos = [0., 0., 0.]
        if self.theta < self.maxTheta:
            pos[2] = self.radius * cos(self.theta)
            for axe in axes:
                pos[axe] = self.coeff * self.radius * sin(self.theta) * (1. - 0.21 * self.theta)
            self.theta += dt
        else:
            self.step += 1
            pos = self.targetMO.findData('position').value

        return pos
    
    def moveBack(self, axes=[0], dt=0.01):
        pos = [0., 0., 0.]
        if self.theta > self.minTheta:
            pos[2] = self.radius * cos(self.theta)
            for axe in axes:
                pos[axe] = self.coeff * self.radius * sin(self.theta) * (1. - 0.21 * self.theta)

            self.theta -= dt
        else:
            self.step += 1
            pos = self.targetMO.findData('position').value
        return pos

    def onBeginAnimationStep(self):
        pos = self.targetMO.findData('position').value
        if self.step == 1:
            axes = [0]
            pos = self.moveF(axes, 0.01)
            self.targetMO.findData('position').value = pos
        elif self.step == 2:
            axes = [0]
            pos = self.moveBack(axes, 0.01)
            self.targetMO.findData('position').value = pos
        elif self.step == 3:
            axes = [1]
            pos = self.moveF(axes, 0.01)
            self.targetMO.findData('position').value = pos
        elif self.step == 4:
            axes = [1]
            pos = self.moveBack(axes, 0.01)
            self.targetMO.findData('position').value = pos
        elif self.step == 5:
            axes = [0, 1]
            self.maxTheta = 0.4
            self.coeff = 0.6
            self.threshold = 1.5
            pos = self.moveF(axes, 0.005)
            self.targetMO.findData('position').value = pos
        elif self.step == 6:
            axes = [0, 1]
            pos = self.moveBack(axes, 0.005)
            self.targetMO.findData('position').value = pos
        elif self.step == 7:
            self.maxTheta = 0.54
            self.minTheta = 0.01
            self.coeff = 0.9
            self.threshold = 0.3
            self.step = 1
        
    def onKeyPressed(self, c):
        if ord(c) == 21:  # up
            pos = self.targetMO.findData('position').value
            pos[0][0] += self.rate
            self.targetMO.findData('position').value = pos

        if ord(c) == 19:  # down
            pos = self.targetMO.findData('position').value
            pos[0][0] -= self.rate
            self.targetMO.findData('position').value = pos

        if c == "+":  # +y
            pos = self.targetMO.findData('position').value
            pos[0][1] += self.rate
            self.targetMO.findData('position').value = pos

        if c == "-":  # -y
            pos = self.targetMO.findData('position').value
            pos[0][1] -= self.rate
            self.targetMO.findData('position').value = pos

        if ord(c) == 18:  # left
            pos = self.targetMO.findData('position').value
            pos[0][2] -= self.rate
            self.targetMO.findData('position').value = pos

        if ord(c) == 20:  # right
            pos = self.targetMO.findData('position').value
            pos[0][2] += self.rate
            self.targetMO.findData('position').value = pos
