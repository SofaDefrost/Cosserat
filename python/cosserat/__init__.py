# -*- coding: utf-8 -*-
"""
The Cosserat plugin for SOFA: Templates Library Documentation
===============================================================

Utility functions and scene templates for the real-time simulation framework `SOFA <https://www.sofa-framework.org/>`_
and the `Cosserat`_ plugin.

The library can be used with scenes written in python and `PSL <https://github.com/sofa-framework/sofa/tree/master/applications/plugins/PSL>`_.

Example:
********

.. sourcecode:: python

    from stlib.scene import MainHeader
    from stlib.physics.rigid import Cube, Floor
    from stlib.physics.deformable import ElasticMaterialObject

    from cosserat.cable import PullingCable

    def createScene(rootNode):
        MainHeader(rootNode)
        DefaultSolver(rootNode)

        Cube(rootNode, translation=[5.0,0.0,0.0])
        Floor(rootNode, translation=[0.0,-1.0,0.0])

       target = ElasticMaterialObject(volumeMeshFileName="mesh/liver.msh",
                                       totalMass=0.5,
                                       attachedTo=node)

        PullingCable(target)
        

Contents of the library
**********************
Template for cosserat .

Content:
********
.. autosummary::

   CosseratCable 

.. autofunction:: cable
.. autofunction:: CosseratFinger
.. autofunction:: addConstraintPoints
"""

from cosseratCable import CosseratFinger, CosseratCable
from grippercontroller import GripperController
