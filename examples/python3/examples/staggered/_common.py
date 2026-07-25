# -*- coding: utf-8 -*-
"""
_common.py — Helpers partagés pour les scènes staggered Cosserat
=================================================================
Auteur  : Younes Adagolodjo (DEFROST / INRIA / Polytech Lille)
Date    : 2026-05-31
Branche : painless/python-explicit

Objectif
--------
Éliminer la duplication des blocs de calcul de raideur, masse, inertie qui
étaient répétés dans chaque scène Python staggered. Centraliser aussi
l'import de `Sofa.LieGroups.SO3` avec fallback propre.

Pourquoi ce module est crucial
------------------------------
`PainlessBeamForceField` a des valeurs par défaut INTENTIONNELLEMENT incorrectes
(EA=1e6, GA=1e5, GJ=1e5, EIy=1e4 — voir PainlessBeamForceField.cpp ligne 22-28)
pour DÉTECTER les scènes qui oublient de passer les vraies valeurs. Toute scène
qui n'appelle PAS `circular_stiffness` (ou n'utilise pas `add_painless_beam`)
risque de hériter de ces defaults → instabilité explosive.

Convention de section circulaire pleine
---------------------------------------
    A     = π · r²              (aire)
    I_y   = π · r⁴ / 4          (moment quadratique d'inertie, bending)
    J     = π · r⁴ / 2          (moment polaire, torsion)
    EA    = E · A               (raideur axiale)
    GA    = G · A               (raideur de cisaillement)
    EIy   = E · I_y             (raideur de bending Y)
    EIz   = E · I_y             (raideur de bending Z, identique pour section circulaire)
    GJ    = G · J               (raideur de torsion)

Usage minimal
-------------
    from _common import add_painless_beam, mass_per_node, inertia_per_seg, SO3, HAVE_SO3

    m_node = mass_per_node(L=1.0, N=8, r=0.015, rho=500)
    I_seg  = inertia_per_seg(L=1.0, N=8, r=0.015, rho=500)
    ff = add_painless_beam(beamNode, "@state", E=5e5, G=1.65e5, r=0.015,
                           name="ff", printEvery=100)
"""

import math


# ── SO3 (avec fallback) ───────────────────────────────────────────────────────
try:
    from Sofa.LieGroups import SO3       # noqa: F401  (re-exporté)
    HAVE_SO3 = True
except ImportError:
    SO3 = None
    HAVE_SO3 = False
    print("[_common] WARNING: Sofa.LieGroups not found — SO3 integration disabled.")


# ── Raideur — section circulaire pleine ───────────────────────────────────────

def circular_stiffness(E, G, r):
    """
    Calcule (EA, GA, EIy, EIz, GJ) pour une section circulaire pleine.

    Parameters
    ----------
    E : float
        Module de Young [Pa]
    G : float
        Module de cisaillement [Pa] (≈ E / (2·(1+ν)))
    r : float
        Rayon de la section [m]

    Returns
    -------
    tuple of float
        (EA, GA, EIy, EIz, GJ) en unités SI (N, N, N·m², N·m², N·m²)
    """
    A    = math.pi * r**2
    I_y  = math.pi * r**4 / 4.0
    J    = math.pi * r**4 / 2.0
    return E * A, G * A, E * I_y, E * I_y, G * J


# ── Inertie et masse pour beam discret uniforme ───────────────────────────────

def mass_per_node(L, N, r, rho):
    """
    Masse nodale uniforme : masse totale / (N+1) nœuds.

    Parameters
    ----------
    L : float    Longueur totale du beam [m]
    N : int      Nombre de segments
    r : float    Rayon de section [m]
    rho : float  Densité [kg/m³]

    Returns
    -------
    float        Masse par nœud [kg]
    """
    return math.pi * r**2 * L * rho / (N + 1)


def inertia_per_seg(L, N, r, rho):
    """
    Moment d'inertie scalaire d'un segment cylindrique (axe transverse).

    Formule cylindre plein, axe perpendiculaire à l'axe :
        I = m · (L_seg² / 12 + r² / 4)

    Returns
    -------
    float        Moment d'inertie par segment [kg·m²]
    """
    L_seg = L / N
    mass  = math.pi * r**2 * L_seg * rho
    return mass * (L_seg**2 / 12.0 + r**2 / 4.0)


# ── Helper unique pour ajouter PainlessBeamForceField ─────────────────────────

def add_painless_beam(node, state_path, E, G, r, name="ff", verbose=True, **kwargs):
    """
    Ajoute un PainlessBeamForceField au `node` avec raideurs pré-calculées.

    Garantit que EA/GA/EIy/EIz/GJ sont passés explicitement, évitant ainsi
    les valeurs par défaut C++ délibérément incorrectes.

    Parameters
    ----------
    node : Sofa.Core.Node
        Le noeud SOFA auquel ajouter le ForceField.
    state_path : str
        Chemin vers le CosseratIntrinsicState lié (ex: "@state").
    E, G : float
        Modules de Young et de cisaillement [Pa].
    r : float
        Rayon de la section circulaire [m].
    name : str (default "ff")
        Nom du composant dans le graphe SOFA.
    verbose : bool (default True)
        Affiche les valeurs de raideur à la création.
    **kwargs : dict
        Arguments supplémentaires (ex: printEvery=100, printPerSegment=True).

    Returns
    -------
    Le composant PainlessBeamForceField ajouté.

    Example
    -------
        ff = add_painless_beam(beamNode, "@state", E=5e5, G=1.65e5, r=0.015,
                               printEvery=100)
    """
    EA, GA, EIy, EIz, GJ = circular_stiffness(E, G, r)

    if verbose:
        print(f"  [_common] PainlessBeamForceField '{name}' stiffness:")
        print(f"    EA  = {EA:.4e} N    GA = {GA:.4e} N")
        print(f"    EIy = {EIy:.4e} N·m²   EIz = {EIz:.4e} N·m²   GJ = {GJ:.4e} N·m²")

    return node.addObject(
        "PainlessBeamForceField",
        name=name,
        state=state_path,
        EA=EA, GA=GA, EIy=EIy, EIz=EIz, GJ=GJ,
        **kwargs,
    )


# ── Helper : conversion quaternion → rotation vector ──────────────────────────

import numpy as np  # importé tardivement pour ne pas alourdir l'import quand
                    # seules les fonctions sans numpy sont utilisées.


def quat_to_rotvec(q):
    """
    Convertit un quaternion [qx, qy, qz, qw] en vecteur rotation ω = axis × angle ∈ ℝ³.

    Convention : ω = θ · axis  avec θ ∈ [0, π].
    Identité (qw=1) → ω = [0, 0, 0].

    Utile pour écrire dans `state.orientations.value` qui attend des
    vecteurs rotation (ω = log(R)), pas des quaternions.
    """
    qx, qy, qz, qw = float(q[0]), float(q[1]), float(q[2]), float(q[3])
    qw = np.clip(qw, -1.0, 1.0)
    angle = 2.0 * np.arccos(abs(qw))
    sin_half = np.sin(angle / 2.0)
    if sin_half < 1e-10:
        return np.zeros(3)
    axis = np.array([qx, qy, qz]) / sin_half
    if qw < 0:
        axis = -axis    # garde l'angle dans [0, π]
    return axis * angle


# ── Export ────────────────────────────────────────────────────────────────────
__all__ = [
    "SO3",
    "HAVE_SO3",
    "circular_stiffness",
    "mass_per_node",
    "inertia_per_seg",
    "add_painless_beam",
    "quat_to_rotvec",
]
