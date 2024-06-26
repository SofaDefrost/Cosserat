#!/usr/bin/env python
# -*- coding: utf-8 -*-
import Sofa.Core
import Sofa.constants.Key as Key
import math
import numpy as np
import pandas as pd
import csv
import sys
import os
path=os.path.dirname(os.path.abspath(__file__))+'/data/'
class FingerController(Sofa.Core.Controller):
    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, args, kwargs)
        self.node = args[0]
        self.fingerNode = self.node.getChild('finger')
        self.pressureConstraint = self.fingerNode.cavity.getObject('SurfacePressureConstraint')#air chamber pressure
        self.trackers = self.fingerNode.tracker.getObject('DOF') # trackers
        #self.timestep=0 #for measuring the time elapsed during each step
        self.MechanicalObject=self.fingerNode.getObject('tetras')
        self.counter=0 #iteration counter
        self.data=[] #pressure time angle data 
        self.totaltime=0.0 #total time elapsed
        self.springs=self.fingerNode.hitching.getObject('springs')
        #self.springs.spring.value.value[0].Ks=1000
        self.Ks=0
        print(self.data)
    

    def onAnimateBeginEvent(self, event): # called at each begin of animation step
        stepsize=2 #how many seconds per step
        steps=27 #number of steps per iteration
        totaliterations=1 #total number of iterations
        self.totaltime=self.node.findData("time").value-self.counter*steps*stepsize
        i=0
        if self.totaltime>=steps*stepsize and self.counter<=(totaliterations-1):
            self.totaltime=0.0
            self.MechanicalObject.reset()
            pressureValue=0.0
            self.pressureConstraint.value = [pressureValue]
            pd.DataFrame(self.data).to_csv(path+'sofaangle'+str(self.counter)+'.csv')
            self.data=[]
            self.counter=self.counter+1
            #self.springs.stiffness.value=self.Ks+delta*self.counter
            self.trackers.reset()
            print(self.data)
            #print(angle)
            print(self.counter)
            print("new iteration")
            i=1
            
        if self.totaltime==0.0:
            print("start")
            self.timestep=0
        #print(self.timestep)
        pressureValue = self.pressureConstraint.value.value[0]
        if self.totaltime<=steps*stepsize and self.counter<=totaliterations:
            if self.totaltime-self.timestep>=stepsize :
                pressureValue+=0.05
                #print("step")
                #print(pressureValue)
                self.timestep=self.totaltime
                if pressureValue>1.3:
                    pressureValue=1.3
                    #print("max pressure")
            #print(self.node.findData("time").value)
            self.pressureConstraint.value = [pressureValue]
            position=self.trackers.findData('position').value
            record=position
            vector_1=[1,0]
            x1=position[0][0]
            z1=position[0][1]
            y1=position[0][2]
            x2=position[1][0]
            z2=position[1][1]
            y2=position[1][2]
            vector_2=[x2-x1,y2-y1]
            unit_vector_1 = vector_1 / np.linalg.norm(vector_1)
            unit_vector_2 = vector_2 / np.linalg.norm(vector_2)
            dot_product = np.dot(unit_vector_1, unit_vector_2)
            angle = np.arccos(dot_product)
            angle=math.degrees(angle)
            if i == 1:
                angle=0
                print(angle, ' angle reset')
            self.data.append([pressureValue,self.totaltime,angle,x1,y1,z1,x2,y2,z2])
            #print(self.totaltime)
        if self.counter>totaliterations:
            pd.DataFrame(self.data).to_csv('sofaangle'+str(self.Ks)+'Ks.csv')
            #print(self.counter)
            #print("FINISHED")
            return




    
        
        

