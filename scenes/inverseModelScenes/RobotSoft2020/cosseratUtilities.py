# -*- coding: utf-8 -*-
from math import cos
from math import sin, sqrt, pi
from splib.numerics import Vec3, Quat

import Sofa
#from stlib.scene import MainHeader, ContactHeader

def createCurvAbsOutput(distance):
    curv_abs_output=[[0.]]
    
    for i in range(0,len(distance)):
        temp = 0
        for k in range(0,i+1):
            temp += distance[k]
        curv_abs_output += [temp]
    
    return curv_abs_output


def compute_BeamLenght(positions):
    listDist=[]        
    for i in range(0,len(positions)-1):
        t1 = positions[i]
        t2 = positions[i+1]
        dist = sqrt((t1[0]-t2[0])**2 + (t1[1]-t2[1])**2 + (t1[2]-t2[2])**2)        
        listDist += [dist] 
        
    return listDist

def createFramesList(positions):
        frames=[]
        for position in positions:
            frames += [[position[0],position[1],position[2],0,0,0,1]]            
        return frames

def extractFEMConstraintPoints(positions,):
    femPoints=[]
    for i in range(0,len(positions)-1,2):
        femPoints += [[positions[i][0], positions[i][1], positions[i][2]]]            
    femPoints += [positions[len(positions)-1]]
    
    return femPoints



