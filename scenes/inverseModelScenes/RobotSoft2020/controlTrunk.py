from math import sqrt, pi, cos, sin
from splib.objectmodel import SofaPrefab, SofaObject
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
        self.seuil = 0.3
        
        return

    def initGraph(self, targetNode):
        self.targetMO = self.target.getObject('dofs')
    
    def moveF(self, axes=[0],dt=0.01):
        pos = [0.,0.,0.]
        if (self.theta < self.maxTheta):
            pos[2] = self.radius * cos(self.theta)
            for axe in axes:
                pos[axe] = self.coeff * self.radius * sin(self.theta)* (1. - 0.21*self.theta)
                #if(self.theta > self.seuil) :
                    #pos[axe] = self.coeff * self.radius * sin(self.theta)
                #else:
                    #pos[axe] = self.radius * sin(self.theta)
            self.theta += dt 
#            print ("=================++> ",self.theta)
            
        else :
            self.step +=1 ;
            pos = self.targetMO.findData('position').value
#            print("=================> Passe to knew step")

        return pos
    
    def moveBack(self, axes=[0],dt=0.01):
        pos = [0.,0.,0.]
        if (self.theta > self.minTheta):
            pos[2] = self.radius * cos(self.theta)
            for axe in axes:
                pos[axe] = self.coeff * self.radius * sin(self.theta)* (1. - 0.21*self.theta)
                
#                if(self.theta < self.seuil) :
#                    pos[axe] = self.coeff * self.radius * sin(self.theta)
#                else:
#                    pos[axe] = self.radius * sin(self.theta)
            self.theta -= dt 
#            print ("=================++> ",self.theta)
            
        else :
            self.step +=1 ;
            pos = self.targetMO.findData('position').value
#            print("=================> Passe to knew step")

        return pos
    
        
        
    def onBeginAnimationStep(self, dt):
        
        pos = self.targetMO.findData('position').value
        if self.step == 1:
            axes = [0]
            pos = self.moveF(axes,0.01)
            self.targetMO.findData('position').value = pos
            #if (self.theta < 0.55):
                #z = self.radius * cos(self.theta)
                #if(self.theta > 0.3) :
                    #x = 0.9 * self.radius * sin(self.theta)
                #else:
                    #x = self.radius * sin(self.theta)
                #self.theta += dt 
                #pos[0][0] = x
                #pos[0][2] = z                
                #self.targetMO.findData('position').value = pos
            #elif self.theta > 0.54 :
                #self.step = 2
        elif self.step == 2:
            axes = [0]
            pos = self.moveBack(axes,0.01)
            self.targetMO.findData('position').value = pos
            #if (self.theta > 0.01):
                #z = self.radius * cos(self.theta)
                #if(self.theta < 0.3) :
                    #x = 0.9 * self.radius * sin(self.theta)
                #else:
                    #x = self.radius * sin(self.theta)
                #self.theta -= dt 
                #pos[0][0] = x
                #pos[0][2] = z                
                #self.targetMO.findData('position').value = pos
            #elif self.theta < 0.02 :
                #self.step = 3
                #self.theta = 0.
        elif self.step ==3 :
            axes = [1]
            pos = self.moveF(axes,0.01)
            self.targetMO.findData('position').value = pos
            #if (self.theta < 0.54):
                #z = self.radius * cos(self.theta)
                #if(self.theta < 0.3) :
                    #y = 0.9 * self.radius * sin(self.theta)
                #else:
                    #y = self.radius * sin(self.theta)
                
                #self.theta += dt 
                #pos[0][2] = z
                #pos[0][1] = y                
                #self.targetMO.findData('position').value = pos
            #elif self.theta > 0.53 :
                #self.step = 4
        elif self.step == 4:
            axes = [1]
            pos = self.moveBack(axes,0.01)
            self.targetMO.findData('position').value = pos
            #if (self.theta > 0.01):
                #z = self.radius * cos(self.theta)
                #y = 0
                #if(self.theta < 0.25) :
                    #y = 0.9 * self.radius * sin(self.theta)
                #else:
                    #y = self.radius * sin(self.theta)
                #self.theta -= dt 
                #pos[0][1] = y
                #pos[0][2] = z
                
                #self.targetMO.findData('position').value = pos
            #elif self.theta < 0.02 :
                #self.step = 5
        elif self.step == 5:
#            print("Double variation")
            axes = [0,1]
            self.maxTheta = 0.4
            self.coeff = 0.6
            self.seuil = 1.5
            pos = self.moveF(axes,0.005)
            self.targetMO.findData('position').value = pos
        elif self.step == 6:
#            print("Double variation") 
            axes = [0,1]            
            pos = self.moveBack(axes,0.005)
            self.targetMO.findData('position').value = pos
        elif self.step == 7:
            self.maxTheta = 0.54
            self.minTheta = 0.01
            self.coeff = 0.9
            self.seuil = 0.3
            self.step = 1
        
    def onKeyPressed(self, c):
        if ord(c) == 21:  # up
            pos = self.targetMO.findData('position').value
            pos[0][0] += self.rate
            self.targetMO.findData('position').value = pos
#            print("=======> Position :", pos)
            
        if ord(c) == 19:  # down
            pos = self.targetMO.findData('position').value
            pos[0][0] -= self.rate
            self.targetMO.findData('position').value = pos
#            print("=======> Position :",pos)
             
        if c == "+":  # +y
            pos = self.targetMO.findData('position').value
            pos[0][1] += self.rate
            self.targetMO.findData('position').value = pos
#            print("=======> Position :", pos)
            
        if c == "-":  # -y
            pos = self.targetMO.findData('position').value
            pos[0][1] -= self.rate
            self.targetMO.findData('position').value = pos
#            print("=======> Position :",pos)

        if ord(c) == 18:  # left
            pos = self.targetMO.findData('position').value
            pos[0][2] -= self.rate
            self.targetMO.findData('position').value = pos
#            print("=======> Position :", pos)

        if ord(c) == 20:  # right
            pos = self.targetMO.findData('position').value
            pos[0][2] += self.rate
            self.targetMO.findData('position').value = pos
#            print("=======> Position :", pos)
