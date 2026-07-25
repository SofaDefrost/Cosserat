# -*- coding: utf-8 -*-
"""
Staggered Cosserat — Scène SOFA de validation Euler-Bernoulli
==============================================================
Auteur  : Younes Adagolodjo (DEFROST / INRIA / Polytech Lille)
Contact : adagolodjo@protonmail.com
Branche : painless/python-explicit
Référence : Romanyà-Serrasolsas et al., SIGGRAPH 2025. DOI: 10.1145/3730944

Objectif
--------
Valider la discrétisation staggered en comparant la déflexion de pointe
d'un cantilever sous charge répartie (gravité) avec la solution analytique
d'Euler-Bernoulli (régime petit déplacement) :

    δ_EB = q · L⁴ / (8 · E · I)

où q = ρ · π · r² · |g| est la charge linéique (N/m).

Paramètres choisis pour rester dans le régime EB (δ / L ≈ 1.2 %) :
    L = 0.5 m,  r = 0.01 m,  E = 1.0e8 Pa,  ρ = 200 kg/m³
    →  δ_EB ≈ 6.14 mm  =  1.23 % · L  ✓ (régime linéaire)

Cible de validation
-------------------
    |δ_num − δ_EB| / δ_EB  <  5 %   pour  N = 8
    |δ_num − δ_EB| / δ_EB  <  1 %   pour  N = 32
    (ordre attendu : O(h²) — doubler N divise l'erreur par ≈ 4)

API SOFA utilisée
-----------------
Conformément à MEMOIRE.md §4, on accède aux Data fields via `.value`
(pas via les méthodes pybind11 getPositions/setPositions/getNodalForces…
qui ne sont pas exposées par le binding actuel).

Architecture SOFA
-----------------
      rootNode
        ├── RequiredPlugin (Cosserat)
        ├── DefaultAnimationLoop
        ├── DefaultVisualManagerLoop
        └── staggeredBeam  (child)
              ├── EdgeSetTopologyContainer
              ├── CosseratIntrinsicState     — DOFs : N+1 Vec3 + N SO3 (ω = log R)
              ├── CosseratTopologyBuilder    — init géométrique
              ├── PainlessBeamForceField     — forces élastiques staggered
              └── StaggeredCosseratMapping   — affichage + frames Rigid3d

Contrôleur Python
-----------------
StaggeredValidationController :
  - Intégration Euler explicite sur R³ × SO(3) (onAnimateEndEvent)
  - Amortissement fort (DAMPING ~ 0.97) pour convergence rapide vers équilibre
  - Surveillance convergence : |Δδ| / δ < CONVERGENCE_TOL
  - À convergence : compare δ_num vs δ_EB et affiche l'erreur relative

Usage
-----
  runSofa sofa_staggered_validation.py
  runSofa sofa_staggered_validation.py --argv "N=16"

Touche [v] : forcer l'affichage de la comparaison courante
Touche [r] : reset
Touche [q] : quitter
"""

import math
import os
import sys

import numpy as np
import Sofa

# ── Helpers partagés (_common.py dans le même dossier) ────────────────────────
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    SO3, HAVE_SO3 as _HAVE_SO3,
    add_painless_beam, mass_per_node, inertia_per_seg,
)

# ── Paramètres physiques (régime petit déplacement Euler-Bernoulli) ───────────
BEAM_LENGTH  = 0.5         # m
RADIUS       = 0.01        # m   section circulaire
YOUNG_MOD    = 1.0e8       # Pa  régime linéaire (petit déplacement)
SHEAR_MOD    = YOUNG_MOD / (2.0 * (1.0 + 0.49))  # ν ≈ 0.49 (caoutchouc)
DENSITY      = 200.0       # kg/m³
GRAVITY_Y    = -9.81       # m/s²

# Quantités géométriques
AREA         = math.pi * RADIUS**2
I_SECTION    = math.pi * RADIUS**4 / 4.0
J_SECTION    = math.pi * RADIUS**4 / 2.0
Q_LINE       = DENSITY * AREA * abs(GRAVITY_Y)            # N/m
DELTA_EB     = Q_LINE * BEAM_LENGTH**4 / (8.0 * YOUNG_MOD * I_SECTION)  # m

# ── Paramètres de simulation ──────────────────────────────────────────────────
N_SEGMENTS   = 8           # peut être surchargé via argv (ex: --argv "N=16")
DT           = 5e-5        # s  (Euler explicite → dt petit)
DAMPING_POS  = 0.97        # fort amortissement pour convergence rapide
DAMPING_ANG  = 0.97
CONVERGENCE_TOL  = 5e-4    # |Δδ| / δ < TOL → équilibre atteint
WARMUP_STEPS = 500         # ne pas tester la convergence avant ce step
MAX_STEPS    = 100_000     # limite de sécurité
DEBUG_FREQ   = 1000

# ── Lecture d'un paramètre N depuis argv (ex: --argv "N=16") ──────────────────
for arg in sys.argv[1:]:
    if arg.startswith('N='):
        try:
            N_SEGMENTS = int(arg.split('=')[1])
        except ValueError:
            pass


# ── Helpers : importés depuis _common.py ──────────────────────────────────────
# (mass_per_node, inertia_per_seg, circular_stiffness, add_painless_beam)


# ── Contrôleur d'intégration et de validation ─────────────────────────────────

class StaggeredValidationController(Sofa.Core.Controller):
    """
    Contrôleur Euler explicite sur R³×SO(3) + validation Euler-Bernoulli.

    Event ordering :
      onAnimateEndEvent — PainlessBeamForceField (AnimateBeginEvent) a déjà
      calculé les forces depuis l'état courant ; on les lit puis on intègre.

    À chaque pas :
      1. Lit positions, forces et torques via les Data .value.
      2. Intègre positions (nodes 1..N) et orientations (segments 1..N-1)
         par Euler explicite avec fort amortissement.
      3. Surveille la convergence de la déflexion de pointe.
      4. À convergence, affiche l'erreur relative vs. solution analytique.
    """

    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)

        self.state    = kwargs['state']
        self.ff       = kwargs['ff']
        self.N        = kwargs['N']
        self.m        = kwargs['mass_per_node']
        self.I_seg    = kwargs['inertia_per_seg']
        self.delta_EB = kwargs['delta_EB']

        self._vel_pos = np.zeros((self.N + 1, 3))
        self._vel_ang = np.zeros((self.N,     3))

        h = BEAM_LENGTH / self.N
        self._rest_pos = np.array([[i * h, 0., 0.] for i in range(self.N + 1)])

        self._step       = 0
        self._converged  = False
        self._prev_delta = 0.0

    # ── Touches ───────────────────────────────────────────────────────────────

    def onKeypressedEvent(self, event):
        k = event['key']
        if k == 'v':
            self._print_comparison()
        elif k == 'r':
            self._reset()
        elif k == 'q':
            self.getContext().animate = False

    # ── Réinitialisation ──────────────────────────────────────────────────────

    def _reset(self):
        self.state.positions.value = self._rest_pos.tolist()
        if _HAVE_SO3:
            # Identité SO3 → rotation vector = [0, 0, 0]
            self.state.orientations.value = [[0.0, 0.0, 0.0]] * self.N
        self._vel_pos[:] = 0.
        self._vel_ang[:] = 0.
        self._step       = 0
        self._converged  = False
        self._prev_delta = 0.0
        print("  [reset] Configuration droite restaurée.")

    # ── Affichage de la comparaison ───────────────────────────────────────────

    def _print_comparison(self):
        pos   = np.array(self.state.positions.value).reshape(-1, 3)
        tip_y = float(pos[-1][1])
        delta = abs(tip_y)
        err   = abs(delta - self.delta_EB) / self.delta_EB * 100.0

        print("\n" + "─" * 60)
        print(f"  N segments       : {self.N}")
        print(f"  Étape            : {self._step}")
        print(f"  δ numérique      : {delta * 1000:.4f} mm")
        print(f"  δ Euler-Bernoulli: {self.delta_EB * 1000:.4f} mm")
        print(f"  Erreur relative  : {err:.3f} %")
        print(f"  δ / L            : {delta / BEAM_LENGTH * 100:.2f} %"
              f"   [< 5 % → régime petit déplacement]")
        status = "✓ convergé" if self._converged else "⋯ en cours"
        print(f"  État             : {status}")
        print("─" * 60 + "\n")

    # ── Pas d'animation ───────────────────────────────────────────────────────

    def onAnimateEndEvent(self, event):
        # NOTE: onAnimateEndEvent — PainlessBeamForceField a déjà rempli
        # nodalForces / segmentTorques depuis l'état courant (AnimateBeginEvent).
        # Mettre l'intégrateur sur AnimateBeginEvent crée un lag d'un step → zigzag.
        if self._converged:
            return
        if self._step >= MAX_STEPS:
            print(f"\n[validation] MAX_STEPS={MAX_STEPS} atteint sans convergence.")
            self._print_comparison()
            self._converged = True
            return

        dt = float(self.getContext().dt.value)

        # ── Lecture de l'état et des forces via les Data .value ───────────────
        pos    = np.array(self.state.positions.value).reshape(-1, 3)
        f_el   = np.array(self.ff.nodalForces.value).reshape(-1, 3)
        tau_el = np.array(self.ff.segmentTorques.value).reshape(-1, 3)

        if pos.shape[0] != self.N + 1 or f_el.shape[0] != self.N + 1:
            return

        # ── First-step sanity check ───────────────────────────────────────────
        if self._step == 0:
            print("\n── [validation] FIRST STEP SANITY CHECK ──────────────────────")
            print(f"  pos shape    : {pos.shape}     (expected ({self.N+1}, 3))")
            print(f"  f_el shape   : {f_el.shape}     (expected ({self.N+1}, 3))")
            print(f"  tau_el shape : {tau_el.shape}    (expected ({self.N}, 3))")
            print(f"  max|f_el|    : {np.linalg.norm(f_el, axis=1).max():.4e} N"
                  f"   [≈ 0 au repos]")
            print(f"  dt           : {dt:.2e} s   m={self.m:.3e} kg   I_seg={self.I_seg:.3e} kg·m²")
            print(f"  g_eff_y      : {GRAVITY_Y:.2f} m/s²")
            print(f"  δ_EB cible   : {self.delta_EB * 1000:.4f} mm"
                  f"   ({self.delta_EB / BEAM_LENGTH * 100:.2f} % · L)")
            print(f"  SO3 actif    : {_HAVE_SO3}")
            print("── End first-step check ────────────────────────────────────\n")

        # ── Intégration des positions (nodes 1..N, node 0 clamped) ────────────
        new_pos = pos.copy()
        g_eff   = GRAVITY_Y
        for i in range(1, self.N + 1):
            f_total       = f_el[i].copy()
            f_total[1]   += self.m * g_eff               # gravité en −Y
            acc            = f_total / self.m
            self._vel_pos[i] = self._vel_pos[i] * DAMPING_POS + acc * dt
            new_pos[i]    = pos[i] + self._vel_pos[i] * dt
        self.state.positions.value = new_pos.tolist()

        # ── Intégration SO3 (segments 1..N-1, segment 0 clamped) ──────────────
        if _HAVE_SO3:
            omegas = np.array(self.state.orientations.value).reshape(-1, 3)
            R_list = [SO3.exp(omegas[i]) for i in range(self.N)]
            new_R  = list(R_list)
            for i in range(1, self.N):                   # segment 0 clamped
                alpha             = tau_el[i] / self.I_seg
                self._vel_ang[i]  = self._vel_ang[i] * DAMPING_ANG + alpha * dt
                dR                = SO3.exp(self._vel_ang[i] * dt)
                new_R[i]          = R_list[i] * dR
            self.state.orientations.value = [list(new_R[i].log()) for i in range(self.N)]

        # ── Surveillance de la convergence ────────────────────────────────────
        self._step += 1
        tip_y  = float(new_pos[-1][1])
        delta  = abs(tip_y)
        if self._step > WARMUP_STEPS and delta > 1e-8:
            rel_change = abs(delta - self._prev_delta) / delta
            if rel_change < CONVERGENCE_TOL:
                self._converged = True
                print(f"\n  ✓ ÉQUILIBRE ATTEINT au step {self._step} "
                      f"(t={self._step * dt:.3f} s)")
                self._print_comparison()
        self._prev_delta = delta

        # ── Affichage périodique ──────────────────────────────────────────────
        if self._step % DEBUG_FREQ == 0:
            err = abs(delta - self.delta_EB) / self.delta_EB * 100.0
            print(f"  step={self._step:6d}  δ_tip = {delta * 1000:7.4f} mm"
                  f"  (cible {self.delta_EB * 1000:.4f} mm,"
                  f"  err={err:5.2f} %)")


# ── Construction de la scène SOFA ─────────────────────────────────────────────

def createScene(rootNode):
    N  = N_SEGMENTS
    h  = BEAM_LENGTH / N
    Np = N + 1

    m_node = mass_per_node(BEAM_LENGTH, N, RADIUS, DENSITY)
    I_seg  = inertia_per_seg(BEAM_LENGTH, N, RADIUS, DENSITY)

    # ── Configuration globale ─────────────────────────────────────────────────
    rootNode.gravity.value = [0., 0., 0.]   # gravité gérée par le contrôleur
    rootNode.dt.value      = DT
    rootNode.addObject('DefaultAnimationLoop')
    rootNode.addObject('DefaultVisualManagerLoop')
    rootNode.addObject('VisualStyle',
                       displayFlags='showVisualModels showBehaviorModels')
    rootNode.addObject('InteractiveCamera',
                       position=[0.25, -0.15, 0.8],
                       lookAt=[0.25, -0.02, 0.])

    # ── Plugins ───────────────────────────────────────────────────────────────
    rootNode.addObject('RequiredPlugin', pluginName='Cosserat')

    # ── Noeud poutre ──────────────────────────────────────────────────────────
    beamNode = rootNode.addChild('staggeredBeam')

    beamNode.addObject('EdgeSetTopologyContainer', name='topology')
    beamNode.addObject('EdgeSetTopologyModifier')

    node_pos_str = ' '.join(f'{i*h:.6f} 0 0' for i in range(Np))
    identity_orientations_str = ' '.join(['0 0 0'] * N)
    state = beamNode.addObject('CosseratIntrinsicState',
                               name='state',
                               positions=node_pos_str,
                               orientations=identity_orientations_str,
                               template='Vec3d')

    beamNode.addObject('CosseratTopologyBuilder',
                       name='builder',
                       totalLength=BEAM_LENGTH,
                       nbSections=N,
                       radius=RADIUS,
                       youngModulus=YOUNG_MOD,
                       shearModulus=SHEAR_MOD,
                       intrinsicState='@state',
                       topology='@topology')

    ff = add_painless_beam(beamNode, "@state",
                           E=YOUNG_MOD, G=SHEAR_MOD, r=RADIUS,
                           name="ff", printEvery=0)   # silencieux côté C++

    beamNode.addObject('StaggeredCosseratMapping',
                       name='mapping',
                       state='@state',
                       drawBeam=True,
                       drawFrames=True,
                       drawRadius=RADIUS * 1.5,
                       drawAxisLength=0.015,
                       color='0.15 0.65 0.95 0.9')

    # ── Contrôleur Python ─────────────────────────────────────────────────────
    rootNode.addObject(
        StaggeredValidationController(
            name='validator',
            state=state,
            ff=ff,
            N=N,
            mass_per_node=m_node,
            inertia_per_seg=I_seg,
            delta_EB=DELTA_EB,
        )
    )

    # ── Info de démarrage ─────────────────────────────────────────────────────
    print("\n" + "=" * 64)
    print("  Staggered Cosserat — Validation Euler-Bernoulli")
    print("=" * 64)
    print(f"  L = {BEAM_LENGTH} m   N = {N}   h = {h:.5f} m   r = {RADIUS} m")
    print(f"  E = {YOUNG_MOD:.2e} Pa   G = {SHEAR_MOD:.2e} Pa   ρ = {DENSITY} kg/m³")
    print(f"  I = {I_SECTION:.4e} m⁴   J = {J_SECTION:.4e} m⁴   q = {Q_LINE:.4f} N/m")
    print()
    print(f"  ── Solution analytique Euler-Bernoulli ──")
    print(f"     δ_EB = q · L⁴ / (8 E I) = {DELTA_EB * 1000:.4f} mm"
          f"  ({DELTA_EB / BEAM_LENGTH * 100:.2f} % · L)")
    if DELTA_EB / BEAM_LENGTH > 0.05:
        print(f"     ⚠  δ/L > 5 % — risque sortie du régime petit déplacement")
    else:
        print(f"     ✓  δ/L < 5 % — régime petit déplacement OK pour validation EB")
    print()
    print(f"  dt = {DT:.2e} s   amortissement = {DAMPING_POS}")
    print(f"  Convergence : |Δδ| / δ < {CONVERGENCE_TOL:.0e}"
          f"  (après {WARMUP_STEPS} pas)")
    print(f"  Module LieGroups : {'OK' if _HAVE_SO3 else 'MANQUANT (EA/GA seules)'}")
    print(f"  Touches : [v] valider  [r] reset  [q] quitter")
    print(f"  CLI     : --argv \"N=16\"  pour changer le nombre de segments")
    print("=" * 64 + "\n")

    return rootNode
