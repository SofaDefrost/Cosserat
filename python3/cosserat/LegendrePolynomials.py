#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created legendre polynomials in python3
"""

__authors__ = "younesssss"
__contact__ = "adagolodjo@protonmail.com, yinoussa.adagolodjo@inria.fr"
__version__ = "1.0.0"
__copyright__ = "(c) 2021,Inria"
__date__ = "Nov 18 2021"

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

# matplotlib.use('Agg')


def legendrePoly(n, x):
    if n == 0:
        return 1
    elif n == 1:
        return x
    else:
        return (((2 * n) - 1) * x * legendrePoly(n - 1, x) - (n - 1) * legendrePoly(n - 2, x)) / float(n)


def drawLegendre():
    x = np.linspace(0, 1, 100)
    # print(f'x = {x}')
    for i in range(1, 5):
        # print(f"LX : {legendrePoly(i, x)}")
        plt.plot(x, legendrePoly(i, x), label="Legendre" + str(i))

    plt.legend(loc="best")
    plt.xlabel("X")
    plt.ylabel("Pn")
    plt.show()


def buildMState(sizeAbscissa=6, polyDegree=5):
    x = np.linspace(0.0, 1, sizeAbscissa)
    print(f'x = {x}')
    vectorPoly = []
    for i in range(1, polyDegree):
        poly = legendrePoly(i, x)
        print(f'==> {poly}')
        for p in poly:
            vectorPoly.append([p, 0., 0.])

    print(f'vectorPoly:{vectorPoly}')


# buildMState(sizeAbscissa=3, polyDegree=4)
# for i in range(4):
#     print(f'===> {legendrePoly(i, 0.99)}')

drawLegendre()