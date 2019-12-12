# -*- coding: utf-8 -*-

# import Sofa
# import SofaPython
from math import sin,cos, sqrt, pi
import os
path = os.path.dirname(os.path.abspath(__file__))+'/mesh/'


# class TensionComputing(Sofa.PythonScriptController):
#     def initGraph(self, node):
#         self.tension = 500
#         self.node = node;
#         self.BeamHookeLawForce = self.node.getObject('BeamHookeLawForce')
#
#     def onBeginAnimationStep(self, dt):
#         self.tension = self.tension + 8000 * dt;
#         self.BeamHookeLawForce.findData('tension').value = self.tension

class DataComputationClass():
    """docstring for DataComputationClass.Sofa.PythonScriptController"""

    def __init__(self):
        # super(DataComputationClass,Sofa.PythonScriptController).__init__()
        X=[0.1 , 0.3, 0.4]
        dy = 1.0
        dz = 0.0
        _dy = 0.0
        _dz = 0.0

        # dX = []
        # for i in range(0,len(X)):
        #     print ([0.0, dy + _dy*X[i], dz + _dz*X[i]])
        #     dX.append([0.0, dy + _dy*X[i], dz + _dz*X[i]])
        # print(dX)
    def initGraph(self, node):
        self.tension = 500
        self.node = node;
        self.BeamHookeLawForce = self.node.getObject('BeamHookeLawForce')

    def computeDX(self, dy, dz, _dy, _dz, X):
        dX = []
        for i in len(X):
            dX.append([0.0, dy + _dy*X[i], dz + _dz*X[i]])

        return dX

    def computeDDX(self, dy, dz, _dy, _dz, X):
        ddX = []
        for i in len(X):
            dX.append([0.0, _dy, _dz])

        return ddX

DataComputationClass();
