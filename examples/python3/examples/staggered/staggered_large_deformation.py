# -*- coding: utf-8 -*-
"""
Staggered Cosserat — Large deformation / prescribed curvature (EI/GJ validation)
==================================================================================
Author  : Younes Adagolodjo (DEFROST / INRIA / Polytech Lille)
Contact : adagolodjo@protonmail.com
Date    : 2025
Branch  : painless/base-geometry

Purpose
-------
This scene initialises the beam in a *geometrically non-trivial* configuration
(a circular arc or a helix) by prescribing the SO3 segment orientations, then
lets the system relax from rest under the elastic forces only (no gravity).

The beam should return to the straight rest state, validating that:
  1. SO3 integration is correct (orientations evolve smoothly on SO(3)).
  2. Angular forces (EI, GJ) are non-zero and in the correct direction.
  3. The energy decreases monotonically (dissipative integration).

Initial configurations (choose via INITIAL_SHAPE)
-------------------------------------------------
  'arc'    — quarter-circle in the XY plane (uniform Y-bending curvature).
  'helix'  — combined Y-bending + torsion (GJ coupling).
  'twist'  — pure torsion along X (only GJ active).

Controls
--------
  [+] / [-]   increase / decrease damping
  [i]         (re)initialise to prescribed shape
  [p]         print current energy and tip state

Running
-------
  runSofa staggered_large_deformation.py
"""

import Sofa
import numpy as np

try:
    from Sofa.LieGroups import SO3
    _HAVE_SO3 = True
except ImportError:
    _HAVE_SO3 = False
    print("[WARNING] LieGroups module not found — SO3 integration disabled.")

# ── Parameters ────────────────────────────────────────────────────────────────
BEAM_LENGTH    = 1.0      # m
N_SEGMENTS     = 10
RADIUS         = 0.01     # m
YOUNG_MOD      = 2.0e5    # Pa
SHEAR_MOD      = 6.7e4    # Pa  (≈ E/3)
DENSITY        = 500.0    # kg/m³
DT             = 2e-4     # s
DAMPING        = 0.97     # per-step velocity damping (≈ Rayleigh α)

INITIAL_SHAPE  = 'arc'   # 'arc' | 'helix' | 'twist'


# ── Geometry helpers ──────────────────────────────────────────────────────────

def build_arc_state(N, L):
    """
    Quarter-circle arc in the XY plane.

    Positions along a circle of radius R = 2L/π so that arc length = L.
    Segment orientations: tangent = X, normal = Y → rotation around Z by θ_i.
    """
    R_arc   = 2.0 * L / np.pi              # radius giving arc-length = L
    theta   = np.pi / 2.0                  # total arc angle (quarter circle)
    h       = L / N

    positions    = np.zeros((N + 1, 3))
    quats        = np.zeros((N, 4))        # [qx, qy, qz, qw]

    for i in range(N + 1):
        s  = i * h
        th = s / R_arc                     # angle traversed
        positions[i] = [R_arc * np.sin(th), R_arc * (1.0 - np.cos(th)), 0.0]

    for i in range(N):
        s_mid = (i + 0.5) * h
        th    = s_mid / R_arc
        # Rotation around Z by angle th: tangent aligned with local arc direction
        angle = th
        axis  = np.array([0., 0., 1.])
        half  = angle / 2.0
        quats[i] = [axis[0]*np.sin(half),
                    axis[1]*np.sin(half),
                    axis[2]*np.sin(half),
                    np.cos(half)]

    return positions, quats


def build_helix_state(N, L):
    """
    Helical arc: combined bending (EI) + torsion (GJ).

    Uses a helix with pitch p = L/4 and radius r_helix = L/(4π).
    """
    turns    = 1.0
    pitch    = L / (turns * 2.0 * np.pi)
    r_helix  = L / (turns * 2.0 * np.pi * 2.0)    # small radius
    h        = L / N

    positions = np.zeros((N + 1, 3))
    quats     = np.zeros((N, 4))

    for i in range(N + 1):
        t = i * h / L * (turns * 2.0 * np.pi)
        positions[i] = [r_helix * np.cos(t) - r_helix,
                        r_helix * np.sin(t),
                        pitch * t / (2.0 * np.pi)]

    for i in range(N):
        t     = (i + 0.5) * h / L * (turns * 2.0 * np.pi)
        # Frenet–Serret: tangent direction
        tang  = np.array([-r_helix * np.sin(t),
                           r_helix * np.cos(t),
                           pitch / (2.0 * np.pi)])
        tang /= np.linalg.norm(tang)
        # Align X-axis of body frame with tangent
        e1    = np.array([1., 0., 0.])
        axis  = np.cross(e1, tang)
        s     = np.linalg.norm(axis)
        c     = np.dot(e1, tang)
        if s < 1e-10:
            quats[i] = [0., 0., 0., 1.]
        else:
            axis /= s
            angle = np.arctan2(s, c)
            half  = angle / 2.0
            quats[i] = [axis[0]*np.sin(half),
                        axis[1]*np.sin(half),
                        axis[2]*np.sin(half),
                        np.cos(half)]

    return positions, quats


def build_twist_state(N, L):
    """
    Pure torsion: positions along X (straight), orientations rotating around X.

    Max twist = π/2 (90°) at the free end → only GJ is stressed.
    """
    h         = L / N
    positions = np.zeros((N + 1, 3))
    quats     = np.zeros((N, 4))

    for i in range(N + 1):
        positions[i] = [i * h, 0., 0.]

    max_angle = np.pi / 2.0
    for i in range(N):
        angle = max_angle * (i + 0.5) / N
        half  = angle / 2.0
        # Rotation around X
        quats[i] = [np.sin(half), 0., 0., np.cos(half)]

    return positions, quats


def get_initial_state(shape, N, L):
    if shape == 'helix':
        return build_helix_state(N, L)
    elif shape == 'twist':
        return build_twist_state(N, L)
    else:
        return build_arc_state(N, L)


# ── Controller ────────────────────────────────────────────────────────────────

class LargeDeformationIntegrator(Sofa.Core.Controller):
    """
    Explicit Euler integrator on R³ × SO(3) with energy monitoring.
    No gravity — only elastic relaxation from the prescribed initial state.
    """

    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        self.state   = kwargs['state']
        self.ff      = kwargs['ff']
        self.N       = kwargs['N']
        self.m       = kwargs['mass_per_node']
        self.I       = kwargs['inertia_per_seg']
        self.damping = DAMPING

        self._vel_pos = np.zeros((self.N + 1, 3))
        self._vel_ang = np.zeros((self.N, 3))
        self._step    = 0

    # ── Keyboard ──────────────────────────────────────────────────────────────

    def onKeypressedEvent(self, event):
        key = event['key']
        if key == '+':
            self.damping = min(self.damping + 0.005, 1.0)
            print(f"  Damping = {self.damping:.3f}")
        elif key == '-':
            self.damping = max(self.damping - 0.005, 0.80)
            print(f"  Damping = {self.damping:.3f}")
        elif key == 'i':
            self._init_state()
        elif key == 'p':
            self._print_energy()

    def _init_state(self):
        pos, quats = get_initial_state(INITIAL_SHAPE, self.N, BEAM_LENGTH)
        self.state.setPositions(pos)
        if _HAVE_SO3:
            R_list = [SO3.from_quat_array(quats[i]) for i in range(self.N)]
            self.state.setOrientations(R_list)
        else:
            self.state.setOrientationsQuat(quats)
        self._vel_pos[:] = 0.
        self._vel_ang[:] = 0.
        self._step = 0
        print(f"  [init] Shape '{INITIAL_SHAPE}' applied.")

    def _print_energy(self):
        try:
            f   = self.ff.getNodalForces()
            tau = self.ff.getSegmentTorques()
            # Approximate kinetic energy from velocities
            KE_pos = 0.5 * self.m * np.sum(self._vel_pos**2)
            KE_ang = 0.5 * self.I * np.sum(self._vel_ang**2)
            print(f"  Step {self._step:6d}  KE_pos={KE_pos:.4e}  KE_ang={KE_ang:.4e}"
                  f"  |f_tip|={np.linalg.norm(f[-1]):.4e}"
                  f"  |τ_mid|={np.linalg.norm(tau[self.N//2]):.4e}")
        except Exception as e:
            print(f"  [energy] {e}")

    # ── Animation step ────────────────────────────────────────────────────────

    def onAnimateBeginEvent(self, event):
        dt = float(self.getContext().dt.value)
        self._step += 1

        try:
            pos = self.state.getPositions()
        except Exception:
            return

        try:
            f_el   = self.ff.getNodalForces()
            tau_el = self.ff.getSegmentTorques()
        except Exception:
            return

        if f_el.shape[0] != self.N + 1:
            return

        # ── Positions (nodes 1..N clamped at 0) ───────────────────────────────
        new_pos = pos.copy()
        for i in range(1, self.N + 1):
            acc = f_el[i] / self.m
            self._vel_pos[i] = self._vel_pos[i] * self.damping + acc * dt
            new_pos[i] = pos[i] + self._vel_pos[i] * dt
        self.state.setPositions(new_pos)

        # ── SO3 orientations (segments 1..N-1, segment 0 clamped) ─────────────
        if _HAVE_SO3:
            R_list = self.state.getOrientations()
            new_R  = list(R_list)
            for i in range(1, self.N):
                alpha  = tau_el[i] / self.I
                self._vel_ang[i] = self._vel_ang[i] * self.damping + alpha * dt
                dR     = SO3.exp(self._vel_ang[i] * dt)
                new_R[i] = R_list[i] * dR
            self.state.setOrientations(new_R)

        # Print progress every 500 steps
        if self._step % 500 == 0:
            self._print_energy()


# ── Scene ─────────────────────────────────────────────────────────────────────

def createScene(rootNode):
    N  = N_SEGMENTS
    h  = BEAM_LENGTH / N

    vol      = np.pi * RADIUS**2 * BEAM_LENGTH
    mass_tot = vol * DENSITY
    m_node   = mass_tot / (N + 1)
    I_seg    = m_node * (h**2 / 12.0 + RADIUS**2 / 4.0)

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
                       position=[0.5, 0.5, 2.5],
                       lookAt=[0.5, 0.2, 0.],
                       fieldOfView=45)

    beamNode = rootNode.addChild('relaxBeam')
    beamNode.addObject('EdgeSetTopologyContainer', name='topology')
    beamNode.addObject('EdgeSetTopologyModifier')

    # Initial positions and identity SO3 elements
    # (will be overwritten by controller at init for non-trivial shapes)
    node_pos_str = ' '.join(f'{i*h:.6f} 0 0' for i in range(N + 1))
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

    # Stiffness parameters from cross-section geometry (explicit — no auto-link to TopologyBuilder)
    _A   = np.pi * RADIUS**2
    _I_y = np.pi * RADIUS**4 / 4.0
    _J   = np.pi * RADIUS**4 / 2.0
    _EA  = YOUNG_MOD * _A
    _GA  = SHEAR_MOD * _A
    _EIy = YOUNG_MOD * _I_y
    _EIz = YOUNG_MOD * _I_y
    _GJ  = SHEAR_MOD * _J
    print(f"\n  [scene] PainlessBeamForceField stiffness:")
    print(f"    EA={_EA:.3e} N   GA={_GA:.3e} N   EIy={_EIy:.3e} N·m²   GJ={_GJ:.3e} N·m²")
    ff = beamNode.addObject('PainlessBeamForceField',
                            name='ff',
                            state='@state',
                            EA=_EA, GA=_GA, EIy=_EIy, EIz=_EIz, GJ=_GJ)

    beamNode.addObject('StaggeredCosseratMapping',
                       name='mapping',
                       state='@state',
                       drawBeam=True,
                       drawFrames=True,
                       drawRadius=RADIUS,
                       drawAxisLength=0.04,
                       color='0.9 0.4 0.1 0.9')

    ctrl = LargeDeformationIntegrator(
        name='ctrl',
        state=state,
        ff=ff,
        N=N,
        mass_per_node=m_node,
        inertia_per_seg=I_seg,
    )
    rootNode.addObject(ctrl)

    # Apply initial deformed shape after scene is built
    # (can't call setOrientations before the state is init'd by SOFA)
    # The controller will apply it on first AnimateBegin.
    # We trigger it via a one-shot flag:
    _orig_begin = ctrl.onAnimateBeginEvent

    def _first_step(event):
        ctrl._init_state()
        ctrl.onAnimateBeginEvent = _orig_begin
        _orig_begin(event)

    ctrl.onAnimateBeginEvent = _first_step

    print("\n" + "=" * 62)
    print(f"  Staggered Cosserat — Large deformation relaxation")
    print("=" * 62)
    print(f"  Shape : '{INITIAL_SHAPE}'")
    print(f"  L={BEAM_LENGTH} m  N={N}  h={h:.4f} m")
    print(f"  E={YOUNG_MOD:.2e} Pa   G={SHEAR_MOD:.2e} Pa")
    print(f"  SO3: {'ACTIVE' if _HAVE_SO3 else 'DISABLED'}")
    print(f"  Keys: [+]/[-] damping | [i] re-init | [p] print energy")
    print("=" * 62 + "\n")

    return rootNode
