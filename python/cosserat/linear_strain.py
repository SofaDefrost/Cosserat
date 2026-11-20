# -*- coding: utf-8 -*-
"""
Linear Strain Element (LSE) rods: nodal Rigid3 poses + one strain slope per element.

Builds the mechanical graph of the geometrically explicit Cosserat-rod element of
Xun & Tamadazte, "Geometrically explicit Cosserat-rod modeling with piecewise linear strain
for complex rod systems", Computers & Structures 325 (2026) 108181.

Graph (the elastic energy is the sum of two standard Hooke energies):

    rod
    ├── nodes      MechanicalObject<Rigid3d>   nodal poses g_i (mass, contacts, boundary conditions)
    ├── slopes     MechanicalObject<Vec6d>     strain slopes beta_e           (BeamHookeLawForceField, length h^3/12)
    └── strains    MechanicalObject<Vec6d>     mean strains xi_bar_e          (BeamHookeLawForceField, length h)
         └── LinearStrainElementMapping(nodes, slopes -> strains)

The energy of one element, U = h/2 (xi_bar - xi0)' K (xi_bar - xi0) + h^3/24 beta' K beta (Eq. 21), is
reproduced exactly: the first term is the Hooke law on the mapped strain, whose rest value xi0 comes from
the mapping applied to the rest poses (Remark 4 of the paper); the second is the same law on the slope
state with the effective length h^3/12 and a zero rest slope.

Elements are node pairs, so branched rods, closed loops and lattices are built by sharing nodes.
"""

from dataclasses import dataclass
from typing import List, Optional, Sequence

import Sofa

__all__ = ["LinearStrainRod", "add_linear_strain_rod"]


@dataclass
class LinearStrainRod:
    """Handles on the nodes of a rod built by :func:`add_linear_strain_rod`."""

    node: "Sofa.Core.Node"  # holds the nodal poses
    nodes_mo: "Sofa.Core.Object"  # MechanicalObject<Rigid3d>
    slopes_node: "Sofa.Core.Node"
    slopes_mo: "Sofa.Core.Object"  # MechanicalObject<Vec6d>
    strains_node: "Sofa.Core.Node"
    strains_mo: "Sofa.Core.Object"  # MechanicalObject<Vec6d>
    mapping: "Sofa.Core.Object"  # LinearStrainElementMapping
    lengths: List[float]


def add_linear_strain_rod(
    parent: "Sofa.Core.Node",
    node_poses: Sequence[Sequence[float]],
    elements: Sequence[Sequence[int]],
    *,
    young_modulus: float,
    poisson_ratio: float,
    radius: float,
    lengths: Optional[Sequence[float]] = None,
    name: str = "lseRod",
) -> LinearStrainRod:
    """Add a Linear Strain Element rod under ``parent``.

    Args:
        parent: node that holds the ODE and linear solvers (it must be an ancestor of all rod states).
        node_poses: Rigid3 poses ``[x, y, z, qx, qy, qz, qw]`` of the nodes (also the rest configuration).
        elements: node index pairs ``[a, b]`` of each element.
        young_modulus, poisson_ratio, radius: circular cross-section material, giving
            ``K = diag(GJ, EIy, EIz, EA, GA, GA)``.
        lengths: reference length of each element. Default: chord length of the initial poses
            (exact for straight elements; pass arc lengths for initially curved elements).
        name: name of the rod node.
    """
    n_elements = len(elements)
    if n_elements == 0:
        raise ValueError("elements cannot be empty")

    if lengths is None:
        lengths = []
        for a, b in elements:
            pa, pb = node_poses[a], node_poses[b]
            lengths.append(sum((pb[k] - pa[k]) ** 2 for k in range(3)) ** 0.5)
    lengths = [float(h) for h in lengths]
    if len(lengths) != n_elements:
        raise ValueError("lengths must hold one value per element")
    slope_lengths = [h**3 / 12.0 for h in lengths]  # effective length of the slope energy, Eq. 21

    material = dict(
        crossSectionShape="circular",
        youngModulus=young_modulus,
        poissonRatio=poisson_ratio,
        radius=radius,
    )
    zeros = [[0.0] * 6 for _ in range(n_elements)]

    rod = parent.addChild(name)
    nodes_mo = rod.addObject(
        "MechanicalObject", template="Rigid3d", name="nodes", position=[list(p) for p in node_poses]
    )

    slopes_node = rod.addChild("slopes")
    slopes_mo = slopes_node.addObject("MechanicalObject", template="Vec6d", name="slopes", position=zeros)
    slopes_node.addObject(
        "BeamHookeLawForceField", template="Vec6d", name="slopeEnergy", length=slope_lengths, **material
    )

    strains_node = rod.addChild("strains")
    slopes_node.addChild(strains_node)  # second parent: the strains depend on both the poses and the slopes
    strains_mo = strains_node.addObject("MechanicalObject", template="Vec6d", name="strains", position=zeros)
    mapping = strains_node.addObject(
        "LinearStrainElementMapping",
        name="mapping",
        input1=nodes_mo.getLinkPath(),
        input2=slopes_mo.getLinkPath(),
        output=strains_mo.getLinkPath(),
        elements=[list(map(int, e)) for e in elements],
        lengths=lengths,
        applyRestPosition=True,  # the rest strain xi0 is the mapping of the rest poses (paper, Remark 4)
    )
    strains_node.addObject("BeamHookeLawForceField", template="Vec6d", name="elastic", length=lengths, **material)

    return LinearStrainRod(
        node=rod,
        nodes_mo=nodes_mo,
        slopes_node=slopes_node,
        slopes_mo=slopes_mo,
        strains_node=strains_node,
        strains_mo=strains_mo,
        mapping=mapping,
        lengths=lengths,
    )
