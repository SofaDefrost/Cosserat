# -*- coding: utf-8 -*-

import Sofa
import SofaPython
from math import sin,cos, tan, sqrt, pi
import os
import numpy as np
path = os.path.dirname(os.path.abspath(__file__))+'/mesh/'

#Gauss Quadrature
#Suurce : https://en.wikipedia.org/wiki/Gaussian_quadrature
GaussCoeff = [1.0/sqrt(3.0),0.57735] # Coeff
GaussWeights = [1.0,1.0] #weight
#
# curv_abs_input = [0, 25, 50, 75]
curv_abs_input=[0, 15, 30, 45, 60, 66, 81]
#
# ###############
# ## Rate of angular Deformation  (2 sections)
# ###############
# R_b = 1.0
# L = 75.0
#
# _tension = 0.0

class CosseratActuation(Sofa.PythonScriptController):
    """docstring for DataComputationClass.Sofa.PythonScriptController"""

    def __init__(self):
        print("The python init================================================++++++> ")
        super(DataComputationClass,Sofa.PythonScriptController).__init__()

    def computeX(self):
        X = []
        for j in range(1, len(curv_abs_input)) :
            Li = curv_abs_input[j]
            Li_1 = curv_abs_input[j-1]

            #compute X (the curve absisa) base on the Gauss Quadrature
            # Source:  Federico
            s0 = Li_1  +  GaussCoeff[0] * (Li - Li_1)
            s1 = Li_1  +  GaussCoeff[1] * (Li - Li_1)

            # source https://en.wikipedia.org/wiki/Gaussian_quadrature
            # s0 = ((Li - Li_1)/2.0) * C[0] + (Li + Li_1)/2.0
            # s1 = ((Li - Li_1)/2.0) * C[1] + (Li + Li_1)/2.0
            X.append([s0,s1])
        print ("XXXXX ==> :", X)
        return X

    def computeDX(self, dy, dz, _dy, _dz):
        dX0   = []; dX1   = []
        d_dX0 = []; d_dX1 = []
        for i in range(0,len(self.X)):
            dX0.append([0.0, dy + _dy * self.X[i][0], dz + _dz*self.X[i][0]])
            dX1.append([0.0, dy + _dy * self.X[i][1], dz + _dz*self.X[i][1]])
            d_dX0.append([0.0, _dy, _dz])
            d_dX1.append([0.0, _dy, _dz])
        self.distance   = [dX0,dX1]
        self.d_distance = [d_dX0,d_dX1]

        return [dX0,dX1,d_dX0,d_dX1]


    def computeHelicalParameters(self, d, p):
        dX0  = []; dX1  = []
        d_dX0 = []; d_dX1 = []
        a = (2.0 * pi )/p
        X = self.X
        for i in range(0,len(self.X)):
            dX0.append([0.0, d * cos(a * X[i][0]), d * sin(a * X[i][0]) ])
            dX1.append([0.0, d * cos(a * X[i][1]), d * sin(a * X[i][1]) ])

            d_dX0.append([0.0, -d * a * sin(a * X[i][0]), d * a * cos(a * X[i][0]) ])
            d_dX1.append([0.0, -d * a * sin(a * X[i][1]), d * a * cos(a * X[i][1]) ])

        self.distance   = [dX0,dX1]
        self.d_distance = [d_dX0,d_dX1]

    #
    # def computeD_DX(self, _dy, _dz):
    #     ddX0 = []
    #     ddX1 = []
    #     for i in range(0,len(self.X)):
    #         ddX0.append([0.0, _dy, _dz])
    #         ddX1.append([0.0, _dy, _dz])
    #     print ("_DXDXDXDXDX ==> :", [ddX0,ddX1])
    #     return [ddX0,ddX1]

    def computeIntegrale(self, distance, d_distance,K):
        # distance = self.distance
        # d_distance = self.d_distance
        integral = []

        #### Compute Integral
        for id in range(0,len(K)):
            Li   = curv_abs_input[id+1]
            Li_1 = curv_abs_input[id]

            distance0 = distance[0]
            d_distance0 = d_distance[1]
            Vec0 = np.cross(K[id], distance0[id]) + [1.0, 0.0, 0.0] + d_distance0[id]
            fx0  = np.cross(distance0[id],Vec0)/(np.linalg.norm(Vec0))
            # print ("fx0 : "+str(fx0))

            distance1 = distance[1]
            d_distance1 = d_distance[1]
            Vec1 = np.cross(K[id], distance1[id]) + [1.0, 0.0, 0.0] + d_distance1[id]
            fx1  = np.cross(distance1[id],Vec1)/(np.linalg.norm(Vec1))
            # print ("fx1 : "+str(fx1))
            _int = []
            #Integral base on Gauss Quadrature
            _int.append(np.dot(((Li-Li_1)/2.0), np.dot(GaussWeights[0],fx0) + np.dot(GaussWeights[1],fx1)))
            vecInt = list(_int[0])
            integral.append(vecInt)

        # print ("============<<>>>>> integral by actuation :",vecInt)
        return integral
    def computeMultiDistanceVectors(self, vec_dy, vec_dz, vec_ddy, vec_ddz):
        size_of_actution = len(vec_dy)
        for i in range(size_of_actution) :
            # this doesn't change then can be compute at the initial point
            argu = self.computeDX (vec_dy[i], vec_dz[i], vec_ddy[i], vec_ddz[i])

            # print("====+> The argu : "+str(argu))
            self.distance   = [argu[0], argu[1]]
            self.d_distance = [argu[2], argu[3]]


    def muti_ActuationIntegral(self, vec_dy, vec_dz, vec_ddy, vec_ddz, K):
        #For multi actuation what change for each actuation is d(X)
        #the let compute each d(X)
        integral = []
        size_of_actution = len(vec_dy)
        for i in range(size_of_actution) :
            # this doesn't change then can be compute at the initial point
            argu = self.computeDX (vec_dy[i], vec_dz[i], vec_ddy[i], vec_ddz[i])

            # print("====+> The argu : "+str(argu))
            distance   = [argu[0], argu[1]]
            d_distance = [argu[2], argu[3]]

            # print("====+> The argu 0 : "+str(distance))
            # print("====+> The argu 1 : "+str(d_distance))

            #So for each d(X) let compute the integral
            integral += self.computeIntegrale(distance, d_distance,K)
            # print("====+> The integral : "+str(integral))
        # print ("=++++++++=======+++> muti_ActuationIntegral : "+ str(integral))

        return integral
