# -*- coding: utf-8 -*-
"""
The Cosserat plugin for SOFA: Templates Library Documentation
===============================================================

Utility functions and scene templates for the real-time simulation framework `SOFA <https://www.sofa-framework.org/>`_
and the `Cosserat`_ plugin.

The library can be used with scenes written in python3

Contents of the library
********************** 

. auto_summary::
    :toc_tree: _auto_summary

    cosserat.cosseratObject

"""

__all__ = ["CosseratBase", "usefulFunctions", "cosseratObject", "nonLinearCosserat", "createFemRegularGrid"]

# The legacy examples also use the maintained package `python/cosserat` (header, params, utils, geometry, ...),
# which has the same name: make its modules importable as `cosserat.<module>` (legacy modules take precedence;
# regular modules of the maintained package win over the namespace-style legacy folders such as `utils/`).
import os as _os

_MAINTAINED_PACKAGE = _os.path.normpath(
    _os.path.join(_os.path.dirname(__file__), "..", "..", "..", "python", "cosserat"))
if _os.path.isdir(_MAINTAINED_PACKAGE) and _MAINTAINED_PACKAGE not in __path__:
    __path__.append(_MAINTAINED_PACKAGE)
