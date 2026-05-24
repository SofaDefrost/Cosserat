# -*- coding: utf-8 -*-
"""
Staggered Cosserat — Cantilever full dynamics (Vec3 + SO3)
===========================================================
Author  : Younes Adagolodjo (DEFROST / INRIA / Polytech Lille)
Contact : adagolodjo@protonmail.com
Date    : 2025
Branch  : painless/base-geometry

Reference
---------
Romanyà-Serrasolsas, Casafranca, Otaduy —
"Painless Differentiable Rotation Dynamics", ACM ToG / SIGGRAPH 2025.
DOI: 10.1145/3730944

Description
-----------
Full staggered Cosserat dynamics using the Python bindings for SO3
(module LieGroups) and CosseratIntrinsicState / PainlessBeamForceField
(module Sofa.Cosserat).

Integration: explicit Euler on the product manifold R³ × SO(3)
---------------------------------------------------------------
Position step (nodes i=1..N, node 0 clamped):
    v_{n+1} = v_n * damping + (f_elastic + f_gravity) / m * dt
    x_{n+1} = x_n + v_{n+1} * dt

SO3 step (segments i=0..N-1, segment 0 clamped by zero-torque BC):
    α_n     = τ_elastic / I_seg            (angular acceleration)
    ω_{n+1} = ω_n * damping + α_n * dt    (in body frame)
    R_{n+1} = R_n * SO3.exp(ω_{n+1} * dt) (right-product update on SO(3))

where τ_elastic is provided by PainlessBeamForceField.d_segmentTorques
in the body frame of each segment.

Controls
--------
  [+] / [-]   change gravity coefficient
  [r]         reset beam to straight configuration
  [p]         print tip position and orientation

Requirements
------------
  import Sofa
  from LieGroups import SO3          # built by python/Binding/LieGroups
  from Sofa.Cosserat import (CosseratIntrinsicState,
                              PainlessBeamForceField)  # python/Binding/Cosserat
"""

import Sofa
import numpy as np

# ── Bindings pybind11 ─────────────────────────────────────────────────────────
#  L'import du module 'Cosserat' enregistre le downcast SOFA Python3 qui permet
#  d'appeler state.getPositions(), ff.getNodalForces(), etc. sur les proxies.
#  Il est chargé automatiquement via RequiredPlugin('Cosserat'), mais un import
#  explicite garantit la disponibilité avant le premier accès Python.
try:
    import Cosserat as _CosseratModule   # noqa — enregistre le downcast
    _HAVE_BINDINGS = True
except ImportError:
    _HAVE_BINDINGS = False
    print("[WARNING] Module pybind11 'Cosserat' non trouvé. "
          "Compiler avec -DCOSSERAT_WITH_PYTHON_BINDINGS=ON.")

try:
    from LieGroups import SO3          # pybind11 module (LieGroups)
    _HAVE_SO3 = True
except ImportError:
    _HAVE_SO3 = False
    print("[WARNING] LieGroups module not found — SO3 integration disabled. "
          "Run with the Cosserat plugin compiled with Python bindings.")

# ── Physical parameters ───────────────────────────────────────────────────────
BEAM_LENGTH  = 1.0       # m
N_SEGMENTS   = 8
RADIUS       = 0.015     # m  cross-section radius
YOUNG_MOD    = 5.0e5     # Pa  (soft, to see large deflections quickly)
SHEAR_MOD    = 1.65e5    # Pa  (≈ E / (2*(1+ν)), ν≈0.5)
DENSITY      = 500.0     # kg/m³
DT           = 5e-4      # s   explicit Euler → keep small!
DAMPING_POS  = 0.995     # velocity damping (positions)
DAMPING_ANG  = 0.990     # angular velocity damping
GRAVITY_Y    = -9.81     # m/s²


def _mass_per_node(L, N, r, rho):
    total = np.pi * r**2 * L * rho
    return total / (N + 1)


def _inertia_per_seg(L, N, r, rho):
    """Scalar moment of inertia for a cylindrical segment (rotation about diameter)."""
    L_seg = L / N
    vol   = np.pi * r**2 * L_seg
    mass  = vol * rho
    # I = m/12 * L² + m/4 * r²  (solid cylinder, transverse axis)
    return mass * (L_seg**2 / 12.0 + r**2 / 4.0)


# ── Integrator controller ─────────────────────────────────────────────────────

class StaggeredFullIntegrator(Sofa.Core.Controller):
    """
    Explicit Euler integrator for the full staggered Cosserat state (R³ + SO(3)).

    Position DOFs (nodes 1..N) and SO3 DOFs (segments 1..N-1) are integrated.
    Node 0 and segment 0 are clamped (Dirichlet BC).
    """

    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)

        self.state      = kwargs['state']
        self.ff         = kwargs['ff']
        self.N          = kwargs['N']
        self.m          = kwargs['mass_per_node']
        self.I          = kwargs['inertia_per_seg']
        self.g_coeff    = 1.0
        self._has_so3   = _HAVE_SO3

        # Velocity buffers
        self._vel_pos = np.zeros((self.N + 1, 3))   # translational  (world)
        self._vel_ang = np.zeros((self.N,     3))   # angular (so3, body frame)

        # Rest positions (for reset)
        h = BEAM_LENGTH / self.N
        self._rest_pos = np.array([[i * h, 0., 0.] for i in range(self.N + 1)])

    # ── Key events ────────────────────────────────────────────────────────────

    def onKeypressedEvent(self, event):
        key = event['key']
        if key == '+':
            self.g_coeff = min(self.g_coeff + 0.5, 20.0)
            print(f"  g_eff = {self.g_coeff * GRAVITY_Y:.2f} m/s²")
        elif key == '-':
            self.g_coeff = max(self.g_coeff - 0.5, 0.0)
            print(f"  g_eff = {self.g_coeff * GRAVITY_Y:.2f} m/s²")
        elif key == 'r':
            self._reset()
        elif key == 'p':
            self._print_state()

    def _reset(self):
        """Reset beam to straight, stress-free configuration."""
        self.state.setPositions(self._rest_pos.copy())
        if self._has_so3:
            self.state.setOrientations([SO3() for _ in range(self.N)])
        self._vel_pos[:] = 0.
        self._vel_ang[:] = 0.
        print("  [reset] Beam returned to straight configuration.")

    def _print_state(self):
        pos = self.state.getPositions()
        tip = pos[-1]
        print(f"  Tip position : [{tip[0]:.4f}, {tip[1]:.4f}, {tip[2]:.4f}] m")
        deflection = np.linalg.norm(tip - self._rest_pos[-1])
        print(f"  Tip deflection : {deflection * 100:.2f} cm")
        if self._has_so3:
            R_list = self.state.getOrientations()
            tip_R  = R_list[-1]
            q = tip_R.quaternion()
            print(f"  Tip segment quat : [{q[0]:.4f}, {q[1]:.4f}, {q[2]:.4f}, {q[3]:.4f}]")

    # ── Animation step ────────────────────────────────────────────────────────

    def onAnimateBeginEvent(self, event):
        dt = float(self.getContext().dt)

        # ── Read current state ─────────────────────────────────────────────────
        try:
            pos = self.state.getPositions()         # (N+1, 3) numpy
        except Exception:
            return

        # ── Read elastic forces from force field ──────────────────────────────
        try:
            f_el  = self.ff.getNodalForces()        # (N+1, 3)
            tau_el = self.ff.getSegmentTorques()    # (N, 3)  body frame
        except Exception:
            return

        if f_el.shape[0] != self.N + 1:
            return

        # ── Integrate node positions (explicit Euler) ─────────────────────────
        g_eff = self.g_coeff * GRAVITY_Y
        new_pos = pos.copy()

        for i in range(1, self.N + 1):             # node 0 clamped
            f_total = f_el[i].copy()
            f_total[1] += self.m * g_eff           # gravity in −Y

            acc = f_total / self.m
            self._vel_pos[i] = self._vel_pos[i] * DAMPING_POS + acc * dt
            new_pos[i] = pos[i] + self._vel_pos[i] * dt

        self.state.setPositions(new_pos)

        # ── Integrate SO3 orientations (explicit Euler on SO(3)) ──────────────
        if self._has_so3:
            R_list = self.state.getOrientations()   # list of SO3

            new_R = list(R_list)
            for i in range(1, self.N):              # segment 0 clamped
                tau_i = tau_el[i]                   # body-frame torque (numpy 3)
                alpha = tau_i / self.I              # angular acceleration (body)
                self._vel_ang[i] = (
                    self._vel_ang[i] * DAMPING_ANG + alpha * dt
                )
                # Right-product update: R_{n+1} = R_n * exp(ω * dt)
                dR = SO3.exp(self._vel_ang[i] * dt)
                new_R[i] = R_list[i] * dR

            self.state.setOrientations(new_R)


# ── Scene builder ─────────────────────────────────────────────────────────────

def createScene(rootNode):
    N  = N_SEGMENTS
    h  = BEAM_LENGTH / N
    Np = N + 1

    m_node = _mass_per_node(BEAM_LENGTH, N, RADIUS, DENSITY)
    I_seg  = _inertia_per_seg(BEAM_LENGTH, N, RADIUS, DENSITY)

    # ── Plugins ───────────────────────────────────────────────────────────────
    rootNode.addObject('RequiredPlugin', name='CosseratPlugin',
                       pluginName='Cosserat')
    rootNode.addObject('RequiredPlugin', name='SofaBaseTopology',
                       pluginName='SofaBaseTopology')

    rootNode.addObject('VisualStyle',
                       displayFlags='showVisualModels showBehaviorModels')
    rootNode.gravity.value = [0., 0., 0.]
    rootNode.dt.value      = DT
    rootNode.addObject('DefaultAnimationLoop')
    rootNode.addObject('DefaultVisualManagerLoop')
    rootNode.addObject('InteractiveCamera',
                       position=[0.5, -0.8, 2.5],
                       lookAt=[0.5, 0., 0.],
                       fieldOfView=45)

    # ── Beam node ─────────────────────────────────────────────────────────────
    beamNode = rootNode.addChild('staggeredBeam')
    beamNode.addObject('EdgeSetTopologyContainer', name='topology')
    beamNode.addObject('EdgeSetTopologyModifier')

    node_pos_str = ' '.join(f'{i*h:.6f} 0 0' for i in range(Np))
    state = beamNode.addObject('CosseratIntrinsicState',
                               name='state',
                               positions=node_pos_str,
                               template='Vec3d')

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
                       drawRadius=RADIUS,
                       drawAxisLength=0.03,
                       color='0.15 0.65 0.95 0.9')

    # ── Python integrator ─────────────────────────────────────────────────────
    rootNode.addObject(
        StaggeredFullIntegrator(
            name='integrator',
            state=state,
            ff=ff,
            N=N,
            mass_per_node=m_node,
            inertia_per_seg=I_seg,
        )
    )

    # ── Startup info ──────────────────────────────────────────────────────────
    print("\n" + "=" * 62)
    print("  Staggered Cosserat — Full dynamics (Vec3 + SO3)")
    print("=" * 62)
    print(f"  L={BEAM_LENGTH} m  N={N}  h={h:.4f} m")
    print(f"  E={YOUNG_MOD:.2e} Pa   G={SHEAR_MOD:.2e} Pa   ρ={DENSITY} kg/m³")
    print(f"  m/node={m_node:.4f} kg   I/seg={I_seg:.2e} kg·m²   dt={DT:.2e} s")
    print(f"  SO3 integration: {'ACTIVE (LieGroups module found)' if _HAVE_SO3 else 'DISABLED (import LieGroups failed)'}")
    print(f"  Keys: [+]/[-] gravity  |  [r] reset  |  [p] print state")
    print("=" * 62 + "\n")

    return rootNode
