# -*- coding: utf-8 -*-
"""
Staggered Cosserat — Scène SOFA de validation analytique
=========================================================
Auteur  : Younes Adagolodjo (DEFROST / INRIA / Polytech Lille)
Contact : adagolodjo@protonmail.com
Branche : painless/base-geometry
Référence : Romanyà-Serrasolsas et al., SIGGRAPH 2025. DOI: 10.1145/3730944

Objectif
--------
Valider la discrétisation staggered en comparant la déflexion de pointe
d'un cantilever sous charge répartie (gravité) avec la solution analytique
d'Euler-Bernoulli :

    δ_EB = q · L⁴ / (8 · E · I)

où q = ρ · π · r² · g est la charge linéique (N/m).

Architecture SOFA utilisée
--------------------------
Tous les composants sont instanciés via la factory SOFA (addObject).
L'intégration SO3 est assurée par un contrôleur Python (Euler explicite
sur le groupe de Lie R³×SO(3)).

      rootNode
        ├── RequiredPlugin (Cosserat)
        ├── DefaultAnimationLoop
        ├── DefaultVisualManagerLoop
        └── staggeredBeam  (child)
              ├── EdgeSetTopologyContainer   — topologie SOFA
              ├── CosseratIntrinsicState     — DOFs : N+1 Vec3 + N SO3
              ├── CosseratTopologyBuilder    — init géométrique + rigidités
              ├── PainlessBeamForceField     — forces élastiques staggered
              └── StaggeredCosseratMapping   — affichage + frames Rigid3d

Contrôleur Python
-----------------
StaggeredValidationController :
  - Intégration explicite Euler sur R³ × SO(3)
  - Amortissement fort (DAMPING ~ 0.97) pour convergence rapide vers équilibre
  - Surveillance convergence : |Δδ| / δ < CONVERGENCE_TOL
  - À convergence : compare δ_num vs δ_EB et affiche l'erreur relative

Usage
-----
  runSofa sofa_staggered_validation.py
  runSofa sofa_staggered_validation.py --argv "N=4,8,16,32"

Touche [v] : forcer l'affichage de la comparaison courante
Touche [r] : reset
Touche [q] : quitter
"""

import Sofa
import numpy as np
import sys

# ── Chargement des bindings pybind11 ─────────────────────────────────────────
#
#  Important : importer explicitement le module Cosserat (pybind11) pour que
#  SOFA Python3 effectue le downcast des proxies CosseratIntrinsicState et
#  PainlessBeamForceField, rendant disponibles getPositions(), getOrientations(),
#  getNodalForces(), etc.
#
try:
    import Cosserat as _CosseratModule   # noqa — enregistre le downcast SOFA
    _HAVE_BINDINGS = True
except ImportError:
    _HAVE_BINDINGS = False
    print("[WARNING] Module pybind11 'Cosserat' non trouvé. "
          "Compiler avec -DCOSSERAT_WITH_PYTHON_BINDINGS=ON.")

try:
    from LieGroups import SO3
    _HAVE_SO3 = True
except ImportError:
    _HAVE_SO3 = False
    print("[WARNING] Module 'LieGroups' non trouvé. "
          "L'intégration SO3 sera désactivée.")

# ── Paramètres physiques (Euler-Bernoulli petit déplacement) ──────────────────

BEAM_LENGTH  = 1.0        # m
RADIUS       = 0.01       # m   section circulaire
YOUNG_MOD    = 1.0e6      # Pa  régime linéaire (petit déplacement)
SHEAR_MOD    = YOUNG_MOD / (2.0 * (1.0 + 0.49))  # ν ≈ 0.49 (caoutchouc)
DENSITY      = 1200.0     # kg/m³
GRAVITY_Y    = -9.81      # m/s²

# Moment quadratique I = π r⁴ / 4
I_SECTION    = np.pi * RADIUS**4 / 4.0
# Charge linéique q = ρ · A · |g|
AREA         = np.pi * RADIUS**2
Q_LINE       = DENSITY * AREA * abs(GRAVITY_Y)   # N/m

# Solution analytique Euler-Bernoulli (charge répartie, encastrement-libre)
DELTA_EB     = Q_LINE * BEAM_LENGTH**4 / (8.0 * YOUNG_MOD * I_SECTION)

# ── Paramètres de simulation ──────────────────────────────────────────────────

N_SEGMENTS   = 8          # peut être surchargé via argv
DT           = 1e-4       # s  (Euler explicite → dt petit)
DAMPING_POS  = 0.97       # fort amortissement pour convergence rapide
DAMPING_ANG  = 0.97
CONVERGENCE_TOL  = 5e-4   # |Δδ| / δ < TOL → équilibre atteint
MAX_STEPS    = 50_000     # limite de sécurité

# ── Lecture d'un paramètre N depuis argv ──────────────────────────────────────
for arg in sys.argv[1:]:
    if arg.startswith('N='):
        try:
            N_SEGMENTS = int(arg.split('=')[1])
        except ValueError:
            pass


# ── Helpers ───────────────────────────────────────────────────────────────────

def _mass_per_node(L, N, r, rho):
    """Masse nodale uniforme (masse totale / (N+1))."""
    return np.pi * r**2 * L * rho / (N + 1)


def _inertia_per_seg(L, N, r, rho):
    """Moment d'inertie scalaire d'un segment cylindrique (axe transverse)."""
    L_seg = L / N
    mass  = np.pi * r**2 * L_seg * rho
    return mass * (L_seg**2 / 12.0 + r**2 / 4.0)


# ── Contrôleur d'intégration et de validation ─────────────────────────────────

class StaggeredValidationController(Sofa.Core.Controller):
    """
    Contrôleur Euler explicite sur R³×SO(3) + validation Euler-Bernoulli.

    À chaque pas d'animation :
      1. Lit les forces élastiques depuis PainlessBeamForceField.
      2. Intègre positions (nodes 1..N) et orientations (segments 1..N-1)
         par Euler explicite avec fort amortissement.
      3. Surveille la convergence de la déflexion de pointe.
      4. À convergence, affiche l'erreur relative vs. solution analytique.
    """

    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)

        self.state   = kwargs['state']
        self.ff      = kwargs['ff']
        self.N       = kwargs['N']
        self.m       = kwargs['mass_per_node']
        self.I_seg   = kwargs['inertia_per_seg']
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
        if _HAVE_BINDINGS:
            try:
                self.state.setPositions(self._rest_pos.copy())
                if _HAVE_SO3:
                    self.state.setOrientations([SO3() for _ in range(self.N)])
            except Exception as e:
                print(f"[reset] Erreur setPositions : {e}")
        self._vel_pos[:] = 0.
        self._vel_ang[:] = 0.
        self._step       = 0
        self._converged  = False
        self._prev_delta = 0.0
        print("  [reset] Configuration droite restaurée.")

    # ── Affichage de la comparaison ───────────────────────────────────────────

    def _print_comparison(self):
        try:
            pos   = self.state.getPositions()   # (N+1, 3)
            tip_y = float(pos[-1][1])
            delta = abs(tip_y)
            err   = abs(delta - self.delta_EB) / self.delta_EB * 100.0
        except Exception as e:
            print(f"[validation] Erreur lecture état : {e}")
            return

        print("\n" + "─" * 60)
        print(f"  N segments       : {self.N}")
        print(f"  Étape            : {self._step}")
        print(f"  δ numérique      : {delta * 1000:.4f} mm")
        print(f"  δ Euler-Bernoulli: {self.delta_EB * 1000:.4f} mm")
        print(f"  Erreur relative  : {err:.3f} %")
        status = "✓ convergé" if self._converged else "⋯ en cours"
        print(f"  État             : {status}")
        print("─" * 60 + "\n")

    # ── Pas d'animation ───────────────────────────────────────────────────────

    def onAnimateBeginEvent(self, event):
        if self._converged:
            return
        if self._step >= MAX_STEPS:
            print(f"[validation] MAX_STEPS={MAX_STEPS} atteint sans convergence.")
            self._print_comparison()
            self._converged = True
            return

        dt = float(self.getContext().dt)

        # ── Lecture de l'état ──────────────────────────────────────────────────
        if not _HAVE_BINDINGS:
            return
        try:
            pos = self.state.getPositions()        # numpy (N+1, 3)
        except Exception:
            return

        # ── Lecture des forces élastiques ─────────────────────────────────────
        try:
            f_el  = self.ff.getNodalForces()       # numpy (N+1, 3)
            tau_el = self.ff.getSegmentTorques()   # numpy (N,   3)  cadre corps
        except Exception:
            return

        if f_el.shape[0] != self.N + 1:
            return

        # ── Intégration des positions (nodes 1..N) ────────────────────────────
        new_pos = pos.copy()
        g_eff   = GRAVITY_Y  # m/s²
        for i in range(1, self.N + 1):
            f_total       = f_el[i].copy()
            f_total[1]   += self.m * g_eff
            acc            = f_total / self.m
            self._vel_pos[i] = self._vel_pos[i] * DAMPING_POS + acc * dt
            new_pos[i]    = pos[i] + self._vel_pos[i] * dt

        try:
            self.state.setPositions(new_pos)
        except Exception:
            return

        # ── Intégration SO3 (segments 1..N-1) ─────────────────────────────────
        if _HAVE_SO3:
            try:
                R_list = self.state.getOrientations()
                new_R  = list(R_list)
                for i in range(1, self.N):
                    alpha          = tau_el[i] / self.I_seg
                    self._vel_ang[i] = self._vel_ang[i] * DAMPING_ANG + alpha * dt
                    new_R[i]       = R_list[i] * SO3.exp(self._vel_ang[i] * dt)
                self.state.setOrientations(new_R)
            except Exception:
                pass

        # ── Surveillance de la convergence ─────────────────────────────────────
        self._step += 1
        tip_y  = float(new_pos[-1][1])
        delta  = abs(tip_y)
        if self._step > 200 and delta > 1e-8:
            rel_change = abs(delta - self._prev_delta) / delta
            if rel_change < CONVERGENCE_TOL:
                self._converged = True
                self._print_comparison()
        self._prev_delta = delta

        # Affichage périodique
        if self._step % 1000 == 0:
            print(f"  step={self._step:5d}  δ_tip = {delta * 1000:.4f} mm  "
                  f"(cible : {self.delta_EB * 1000:.4f} mm)")


# ── Construction de la scène SOFA ─────────────────────────────────────────────

def createScene(rootNode):
    N  = N_SEGMENTS
    h  = BEAM_LENGTH / N
    Np = N + 1

    m_node = _mass_per_node(BEAM_LENGTH, N, RADIUS, DENSITY)
    I_seg  = _inertia_per_seg(BEAM_LENGTH, N, RADIUS, DENSITY)

    # ── Configuration globale ─────────────────────────────────────────────────
    rootNode.gravity.value = [0., 0., 0.]   # gravité gérée manuellement
    rootNode.dt.value      = DT
    rootNode.addObject('DefaultAnimationLoop')
    rootNode.addObject('DefaultVisualManagerLoop')
    rootNode.addObject('VisualStyle',
                       displayFlags='showVisualModels showBehaviorModels')
    rootNode.addObject('InteractiveCamera',
                       position=[0.5, -0.5, 2.0],
                       lookAt=[0.5, 0., 0.])

    # ── Plugins ───────────────────────────────────────────────────────────────
    rootNode.addObject('RequiredPlugin', pluginName='Cosserat')
    rootNode.addObject('RequiredPlugin', pluginName='SofaBaseTopology')

    # ── Noeud poutre ─────────────────────────────────────────────────────────
    beamNode = rootNode.addChild('staggeredBeam')

    beamNode.addObject('EdgeSetTopologyContainer',
                       name='topology',
                       position=' '.join(f'{i*h:.6f} 0 0' for i in range(Np)),
                       edges=' '.join(f'{i} {i+1}' for i in range(N)))
    beamNode.addObject('EdgeSetTopologyModifier')

    state = beamNode.addObject('CosseratIntrinsicState',
                               name='state')

    beamNode.addObject('CosseratTopologyBuilder',
                       name='builder',
                       totalLength=BEAM_LENGTH,
                       nbSections=N,
                       radius=RADIUS,
                       youngModulus=YOUNG_MOD,
                       shearModulus=SHEAR_MOD,
                       l_intrinsicState='@state',
                       l_topology='@topology')

    ff = beamNode.addObject('PainlessBeamForceField',
                            name='ff',
                            l_state='@state')

    beamNode.addObject('StaggeredCosseratMapping',
                       name='mapping',
                       l_state='@state',
                       drawBeam=True,
                       drawFrames=True,
                       drawRadius=RADIUS * 1.5,
                       drawAxisLength=0.02,
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
    print(f"  L = {BEAM_LENGTH} m   N = {N}   h = {h:.5f} m")
    print(f"  E = {YOUNG_MOD:.2e} Pa   G = {SHEAR_MOD:.2e} Pa   r = {RADIUS} m")
    print(f"  I = {I_SECTION:.4e} m⁴   q = {Q_LINE:.4f} N/m")
    print(f"  δ_EB = {DELTA_EB*1000:.4f} mm  (solution analytique)")
    print(f"  dt = {DT:.2e} s   amortissement = {DAMPING_POS}")
    print(f"  Bindings Cosserat : {'OK' if _HAVE_BINDINGS else 'MANQUANTS'}")
    print(f"  Module LieGroups  : {'OK' if _HAVE_SO3 else 'MANQUANT'}")
    print(f"  Touches : [v] valider  [r] reset  [q] quitter")
    print("=" * 64 + "\n")

    return rootNode
