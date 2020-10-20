# -*- coding: utf-8 -*-
"""
Created on Jun 5 2020

@author: PSC
"""

import numpy
import pylab as plt
from mpl_toolkits.mplot3d import Axes3D 

def CatmullRomSpline(P0, P1, P2, P3, nPoints=100):
    """
    P0, P1, P2, and P3 should be (x,y) point pairs that define the Catmull-Rom spline.
    nPoints is the number of points to include in this curve segment.
    """
    # Convert the points to numpy so that we can do array multiplication
    P0, P1, P2, P3 = map(numpy.array, [P0, P1, P2, P3])

    # Parametric constant: 0.5 for the centripetal spline, 0.0 for the uniform spline, 1.0 for the chordal spline.
    alpha = 0.5
    # Premultiplied power constant for the following tj() function.
    alpha = alpha/2
    def tj(ti, Pi, Pj):
        xi, yi, zi = Pi
        xj, yj, zj = Pj
        return ((xj-xi)**2 + (yj-yi)**2)**alpha + ti

    # Calculate t0 to t4
    t0 = 0
    t1 = tj(t0, P0, P1)
    t2 = tj(t1, P1, P2)
    t3 = tj(t2, P2, P3)

    # Only calculate points between P1 and P2
    t = numpy.linspace(t1, t2, nPoints)

    # Reshape so that we can multiply by the points P0 to P3
    # and get a point for each value of t.
    t = t.reshape(len(t), 1)
    # print(t)
    A1 = (t1-t)/(t1-t0)*P0 + (t-t0)/(t1-t0)*P1
    A2 = (t2-t)/(t2-t1)*P1 + (t-t1)/(t2-t1)*P2
    A3 = (t3-t)/(t3-t2)*P2 + (t-t2)/(t3-t2)*P3
    # print(A1)
    # print(A2)
    # print(A3)
    B1 = (t2-t)/(t2-t0)*A1 + (t-t0)/(t2-t0)*A2
    B2 = (t3-t)/(t3-t1)*A2 + (t-t1)/(t3-t1)*A3

    C = (t2-t)/(t2-t1)*B1 + (t-t1)/(t2-t1)*B2
    # print("C ", C)
    return C

def CatmullRomChain(P, nb_points=1000):
    """
    Calculate Catmull–Rom for a chain of points and return the combined curve.
    """
    sz = len(P)

    # The curve C will contain an array of (x, y) points.
    C = []
    for i in range(sz-3):
        c = CatmullRomSpline(P[i], P[i+1], P[i+2], P[i+3], nPoints=nb_points//(len(P)-3))
        C.extend(c)

    return C

def make_trajectory(points, nb_points=1000):
    # Calculate the Catmull-Rom splines through the points
    c = CatmullRomChain(points, nb_points=nb_points)
    # Convert the Catmull-Rom curve points into x and y arrays and plot
    x, y, z = zip(*c)
    return x, y, z

# # Define a set of points for curve to go through
# Points = [[0, 1.5, 0], [2, 2, 1], [3, 1, 2], [4, 0.5, 1], [5, 1, 0], [6, 2, 1], [7, 3, 2], [2, 3, 1], [0, 0, 0], [-1, -1, 0]]
#
# x, y, z = make_trajectory(Points, nb_points=200)
#
# xp, yp, zp = [], [], []
# for p in Points:
#     t_xp, t_yp, t_zp = p
#     xp.append(t_xp)
#     yp.append(t_yp)
#     zp.append(t_zp)
#
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.plot(x, y, z, label='parametric curve')
# ax.scatter(xp, yp, zp, c='red', label='control points')
# ax.legend()
#
# plt.show()