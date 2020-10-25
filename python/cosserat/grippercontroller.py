# -*- coding: utf-8 -*-
import Sofa
from splib.numerics import Vec3, Quat
from math import pi

def getTranslated(points, vec,angle):
    r=[]
    
    q = Quat.createFromEuler([(angle[0]*pi)/180., (angle[1]*pi)/180., (angle[2]*pi)/180.], 'ryxz')

    
    v = Vec3(vec[0], vec[1], vec[2])
    sol = v.rotateFromQuat(q)
    
    for v in points:
        r.append( [v[0]+sol[0], v[1]+sol[1], v[2]+sol[2]] )
    return r

class GripperController(Sofa.PythonScriptController):
    def __init__(self, node, fingers, angles):
        self.fingers = fingers
        self.angles = angles
        self.name = "FingerController"

    def onKeyPressed(self,c):
        dir = None
        
        # UP key :
        if ord(c)==19:
            dir = [0.0,1.0,0.0]
        # DOWN key : rear
        elif ord(c)==21:
            dir = [0.0,-1.0,0.0]
        # LEFT key : left
        elif ord(c)==18:
            dir = [-1.0,0.0,0.0]
        elif ord(c)==20:
            dir = [1.0,0.0,0.0]
        elif c=='+':   
            dir = [0.0,0.0,1.0]
        elif c=='-':            
            dir = [0.0,0.0,-1.0]

        if dir != None:
            i = 0
            for finger in self.fingers:
                
                cableNode = finger.getChild("cableNode")
                #trans = finger.trans
#                print("The transformation is :", finger)   
                mecaobject = cableNode.getObject("RigidBaseMO")
                angle = self.angles[i] 
                mecaobject.findData('rest_position').value = getTranslated( mecaobject.rest_position,  dir , angle)
                               
#                print(" The rest Position After: ", mecaobject.rest_position)
                
#                pos = self.rigidBaseMO.findData('position').value
#                pos[0][2] -= self.rate
#                self.rigidBaseMO.findData('position').value = pos

#                cable = m.getChild("PullingCable").getObject("CableConstraint")
#                p = cable.pullPoint[0]
#                cable.findData("pullPoint").value = [p[0]+dir[0], p[1]+dir[1], p[2]+dir[2]]

#16 17 18 19 20 21 48 49 50 51 52 54 63 101 102 103 104 105 106 107 116 128 135 143 150
#16 17 18 19 20 21 48 51 52 54 63 103 104 105 106 107 113 116 128 135 143 150