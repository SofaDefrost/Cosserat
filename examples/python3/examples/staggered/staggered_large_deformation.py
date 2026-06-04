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

import os
import sys

import numpy as np
import Sofa

# ── Helpers partagés (_common.py dans le même dossier) ────────────────────────
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import (
    SO3, HAVE_SO3 as _HAVE_SO3,
    add_painless_beam, mass_per_node, inertia_per_seg,
    quat_to_rotvec as _quat_to_rotvec,
)

# ── Debug printing ────────────────────────────────────────────────────────────
DEBUG_FREQ = 100  # print summary every N steps  (set to 1 for every step)
DEBUG_FREQ_FULL = 500  # print per-segment table every N steps

# ── Parameters ────────────────────────────────────────────────────────────────
BEAM_LENGTH = 1.0  # m
N_SEGMENTS = 10
RADIUS = 0.01  # m
YOUNG_MOD = 2.0e5  # Pa
SHEAR_MOD = 6.7e4  # Pa  (≈ E/3)
DENSITY = 500.0  # kg/m³
DT = 2e-4  # s
DAMPING = 0.97  # per-step velocity damping (≈ Rayleigh α)

INITIAL_SHAPE = "arc"  # 'arc' | 'helix' | 'twist'


# ── Geometry helpers ──────────────────────────────────────────────────────────


def build_arc_state(N, L):
    """
    Quarter-circle arc in the XY plane.

    Positions along a circle of radius R = 2L/π so that arc length = L.
    Segment orientations: tangent = X, normal = Y → rotation around Z by θ_i.
    """
    R_arc = 2.0 * L / np.pi  # radius giving arc-length = L
    theta = np.pi / 2.0  # total arc angle (quarter circle)
    h = L / N

    positions = np.zeros((N + 1, 3))
    quats = np.zeros((N, 4))  # [qx, qy, qz, qw]

    for i in range(N + 1):
        s = i * h
        th = s / R_arc  # angle traversed
        positions[i] = [R_arc * np.sin(th), R_arc * (1.0 - np.cos(th)), 0.0]

    for i in range(N):
        s_mid = (i + 0.5) * h
        th = s_mid / R_arc
        # Rotation around Z by angle th: tangent aligned with local arc direction
        angle = th
        axis = np.array([0.0, 0.0, 1.0])
        half = angle / 2.0
        quats[i] = [
            axis[0] * np.sin(half),
            axis[1] * np.sin(half),
            axis[2] * np.sin(half),
            np.cos(half),
        ]

    return positions, quats


def build_helix_state(N, L):
    """
    Helical arc: combined bending (EI) + torsion (GJ).

    Uses a helix with pitch p = L/4 and radius r_helix = L/(4π).
    """
    turns = 1.0
    pitch = L / (turns * 2.0 * np.pi)
    r_helix = L / (turns * 2.0 * np.pi * 2.0)  # small radius
    h = L / N

    positions = np.zeros((N + 1, 3))
    quats = np.zeros((N, 4))

    for i in range(N + 1):
        t = i * h / L * (turns * 2.0 * np.pi)
        positions[i] = [
            r_helix * np.cos(t) - r_helix,
            r_helix * np.sin(t),
            pitch * t / (2.0 * np.pi),
        ]

    for i in range(N):
        t = (i + 0.5) * h / L * (turns * 2.0 * np.pi)
        # Frenet–Serret: tangent direction
        tang = np.array(
            [-r_helix * np.sin(t), r_helix * np.cos(t), pitch / (2.0 * np.pi)]
        )
        tang /= np.linalg.norm(tang)
        # Align X-axis of body frame with tangent
        e1 = np.array([1.0, 0.0, 0.0])
        axis = np.cross(e1, tang)
        s = np.linalg.norm(axis)
        c = np.dot(e1, tang)
        if s < 1e-10:
            quats[i] = [0.0, 0.0, 0.0, 1.0]
        else:
            axis /= s
            angle = np.arctan2(s, c)
            half = angle / 2.0
            quats[i] = [
                axis[0] * np.sin(half),
                axis[1] * np.sin(half),
                axis[2] * np.sin(half),
                np.cos(half),
            ]

    return positions, quats


def build_twist_state(N, L):
    """
    Pure torsion: positions along X (straight), orientations rotating around X.

    Max twist = π/2 (90°) at the free end → only GJ is stressed.
    """
    h = L / N
    positions = np.zeros((N + 1, 3))
    quats = np.zeros((N, 4))

    for i in range(N + 1):
        positions[i] = [i * h, 0.0, 0.0]

    max_angle = np.pi / 2.0
    for i in range(N):
        angle = max_angle * (i + 0.5) / N
        half = angle / 2.0
        # Rotation around X
        quats[i] = [np.sin(half), 0.0, 0.0, np.cos(half)]

    return positions, quats


def get_initial_state(shape, N, L):
    if shape == "helix":
        return build_helix_state(N, L)
    elif shape == "twist":
        return build_twist_state(N, L)
    else:
        return build_arc_state(N, L)


# _quat_to_rotvec : importé depuis _common.py (alias de quat_to_rotvec)


# ── Controller ────────────────────────────────────────────────────────────────


class LargeDeformationIntegrator(Sofa.Core.Controller):
    """
    Explicit Euler integrator on R³ × SO(3) with energy monitoring.
    No gravity — only elastic relaxation from the prescribed initial state.

    Data access: uses SOFA Data `.value` fields (not pybind11 get/set methods)
    because the Cosserat pybind11 module may not downcast correctly at runtime.

    Event ordering:
      - onAnimateBeginEvent : applies prescribed initial shape (first step only).
      - onAnimateEndEvent   : integrates (PainlessBeamForceField has already
                              computed forces at AnimateBeginEvent).
    """

    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        self.state = kwargs["state"]
        self.ff = kwargs["ff"]
        self.N = kwargs["N"]
        self.m = kwargs["mass_per_node"]
        self.I = kwargs["inertia_per_seg"]
        self.damping = DAMPING

        self._vel_pos = np.zeros((self.N + 1, 3))
        self._vel_ang = np.zeros((self.N, 3))
        self._step = 0
        self._initialized = False  # True after first step applies deformed shape

        h = BEAM_LENGTH / self.N
        self._rest_pos = np.array([[i * h, 0.0, 0.0] for i in range(self.N + 1)])

    # ── Keyboard ──────────────────────────────────────────────────────────────

    def onKeypressedEvent(self, event):
        key = event["key"]
        if key == "+":
            self.damping = min(self.damping + 0.005, 100.0)
            print(f"  Damping = {self.damping:.3f}")
        elif key == "-":
            self.damping = max(self.damping - 0.005, 0.80)
            print(f"  Damping = {self.damping:.3f}")
        elif key == "i":
            self._apply_initial_shape()
            self._vel_pos[:] = 0.0
            self._vel_ang[:] = 0.0
            self._step = 0
        elif key == "p":
            self._print_state()

    # ── Helpers ───────────────────────────────────────────────────────────────

    def _apply_initial_shape(self):
        """Write prescribed deformed configuration to SOFA Data fields."""
        pos, quats = get_initial_state(INITIAL_SHAPE, self.N, BEAM_LENGTH)
        self.state.positions.value = pos.tolist()

        if _HAVE_SO3:
            # Convert quaternions → rotation vectors → write to orientations
            rot_vecs = [_quat_to_rotvec(quats[i]).tolist() for i in range(self.N)]
            self.state.orientations.value = rot_vecs
        else:
            # Fallback: store rotation vectors directly (identity if SO3 unavailable)
            self.state.orientations.value = [[0.0, 0.0, 0.0]] * self.N

        print(f"  [init] Shape '{INITIAL_SHAPE}' applied via Data .value fields.")
        if _HAVE_SO3:
            # Print first few rotation vectors for sanity check
            omegas = np.array(self.state.orientations.value).reshape(-1, 3)
            norms = np.linalg.norm(omegas, axis=1)
            print(
                f"  [init] |ω_seg| : min={norms.min():.3f}  max={norms.max():.3f} rad"
                f"  [0=identity, π/2=90°, π=180°]"
            )

    def _print_state(self):
        pos = np.array(self.state.positions.value).reshape(-1, 3)
        f_el = np.array(self.ff.nodalForces.value).reshape(-1, 3)
        tau_el = np.array(self.ff.segmentTorques.value).reshape(-1, 3)
        tip = pos[-1]
        tip_def = np.linalg.norm(tip - self._rest_pos[-1])
        KE_pos = 0.5 * self.m * np.sum(np.linalg.norm(self._vel_pos, axis=1) ** 2)
        KE_ang = 0.5 * self.I * np.sum(np.linalg.norm(self._vel_ang, axis=1) ** 2)
        print(f"  Step {self._step:6d}  t={self._step * DT:.4f} s")
        print(
            f"  Tip  : ({tip[0]:.5f}, {tip[1]:.5f}, {tip[2]:.5f}) m"
            f"   |tip-rest|={tip_def * 100:.3f} cm"
        )
        print(
            f"  KE   : pos={KE_pos:.4e} J  ang={KE_ang:.4e} J"
            f"  total={KE_pos + KE_ang:.4e} J"
        )
        print(f"  max|f_el|  : {np.linalg.norm(f_el, axis=1).max():.4e} N")
        print(f"  max|tau_el|: {np.linalg.norm(tau_el, axis=1).max():.4e} N·m")

    # ── Animation steps ───────────────────────────────────────────────────────

    def onAnimateBeginEvent(self, event):
        # Apply prescribed initial shape on the very first AnimateBeginEvent.
        # PainlessBeamForceField may fire at the same event — on the NEXT step,
        # forces are guaranteed to come from the deformed configuration.
        if not self._initialized:
            self._initialized = True
            self._apply_initial_shape()

    def onAnimateEndEvent(self, event):
        # PainlessBeamForceField::handleEvent (AnimateBeginEvent) has already
        # computed forces from the current configuration — integrate here.
        dt = float(self.getContext().dt.value)

        # ── Read via SOFA Data fields ─────────────────────────────────────────
        pos = np.array(self.state.positions.value).reshape(-1, 3)
        f_el = np.array(self.ff.nodalForces.value).reshape(-1, 3)
        tau_el = np.array(self.ff.segmentTorques.value).reshape(-1, 3)

        if pos.shape[0] != self.N + 1 or f_el.shape[0] != self.N + 1:
            return

        # ── First-step sanity check ───────────────────────────────────────────
        if self._step == 0:
            print(
                "\n── [large-deform integrator] FIRST STEP SANITY CHECK ──────────────"
            )
            print(f"  pos shape   : {pos.shape}  (expected ({self.N + 1}, 3))")
            print(f"  f_el shape  : {f_el.shape}  (expected ({self.N + 1}, 3))")
            print(f"  tau_el shape: {tau_el.shape}  (expected ({self.N}, 3))")
            max_f = np.linalg.norm(f_el, axis=1).max()
            max_tau = np.linalg.norm(tau_el, axis=1).max()
            print(f"  max|f_el|   : {max_f:.4e} N   max|tau_el|: {max_tau:.4e} N·m")
            print(f"  dt={dt:.2e} s  m={self.m:.3e} kg  I={self.I:.3e} kg·m²")
            print(f"  SO3: {'ACTIVE' if _HAVE_SO3 else 'DISABLED'}")
            if max_f < 1e-10 and max_tau < 1e-10:
                print(
                    f"  ⚠  All forces zero — initial shape may not have been applied yet."
                )
                print(
                    f"     (Normal on step 1 if AnimateBegin ordering puts init before FF.)"
                )
            print(
                "── End first-step check ─────────────────────────────────────────────\n"
            )

        self._step += 1

        # ── Integrate positions (nodes 1..N, node 0 clamped) ─────────────────
        new_pos = pos.copy()
        for i in range(1, self.N + 1):
            acc = f_el[i] / self.m
            self._vel_pos[i] = self._vel_pos[i] * self.damping + acc * dt
            new_pos[i] = pos[i] + self._vel_pos[i] * dt
        self.state.positions.value = new_pos.tolist()

        # ── Integrate SO3 orientations (segments 1..N-1, segment 0 clamped) ──
        if _HAVE_SO3:
            omegas = np.array(self.state.orientations.value).reshape(-1, 3)
            R_list = [SO3.exp(omegas[i]) for i in range(self.N)]
            new_R = list(R_list)
            for i in range(1, self.N):
                alpha = tau_el[i] / self.I
                self._vel_ang[i] = self._vel_ang[i] * self.damping + alpha * dt
                dR = SO3.exp(self._vel_ang[i] * dt)
                new_R[i] = R_list[i] * dR
            self.state.orientations.value = [
                list(new_R[i].log()) for i in range(self.N)
            ]

        # ── Debug summary ─────────────────────────────────────────────────────
        if self._step % DEBUG_FREQ == 0:
            tip = new_pos[-1]
            tip_def = np.linalg.norm(tip - self._rest_pos[-1])
            f_norms = np.linalg.norm(f_el, axis=1)
            t_norms = np.linalg.norm(tau_el, axis=1)
            v_norms = np.linalg.norm(self._vel_pos, axis=1)
            w_norms = np.linalg.norm(self._vel_ang, axis=1)
            KE_pos = 0.5 * self.m * np.sum(v_norms**2)
            KE_ang = 0.5 * self.I * np.sum(w_norms**2)

            print(
                f"\n── [step {self._step:5d}  t={self._step * dt:.4f} s] ──────────────────────────"
            )
            print(
                f"  Tip     : ({tip[0]:.5f}, {tip[1]:.5f}, {tip[2]:.5f}) m"
                f"   |tip-rest|={tip_def * 100:.3f} cm"
            )
            print(
                f"  max|f_el|  : {f_norms.max():.4e} N  (node {int(f_norms.argmax())})"
            )
            print(
                f"  max|tau_el|: {t_norms.max():.4e} N·m  (seg {int(t_norms.argmax())})"
            )
            print(
                f"  max|vel_pos|: {v_norms.max():.4e} m/s  max|vel_ang|: {w_norms.max():.4e} rad/s"
            )
            print(
                f"  KE  : pos={KE_pos:.4e} J  ang={KE_ang:.4e} J  total={KE_pos + KE_ang:.4e} J"
                f"  [should decrease]"
            )
            print(f"  Damping = {self.damping:.3f}")

            if _HAVE_SO3:
                cur_omegas = np.array(self.state.orientations.value).reshape(-1, 3)
                rot_norms = np.linalg.norm(cur_omegas, axis=1)
                print(
                    f"  max|ω_seg| : {rot_norms.max():.4e} rad"
                    f"   [0=identity, warning if >π/2=1.57]"
                )
                if rot_norms.max() > np.pi / 2:
                    print(f"  ⚠  LARGE ROTATION — segment may have flipped!")
                if rot_norms.max() > np.pi:
                    print(f"  ✗  ROTATION > π — log(R) singularity risk!")

        if self._step % DEBUG_FREQ_FULL == 0:
            print(f"\n  ── Per-segment detail (step {self._step}) ────────────────────")
            print(
                f"  {'seg':>4}  {'|tau_el|':>12}  {'|vel_ang|':>12}  {'|ω_seg|':>10}  {'|f_el|':>12}"
            )
            cur_omegas = (
                np.array(self.state.orientations.value).reshape(-1, 3)
                if _HAVE_SO3
                else np.zeros((self.N, 3))
            )
            for i in range(self.N):
                tau_n = np.linalg.norm(tau_el[i]) if i < tau_el.shape[0] else 0.0
                vel_n = np.linalg.norm(self._vel_ang[i])
                rot_n = np.linalg.norm(cur_omegas[i])
                f_n = np.linalg.norm(f_el[i + 1]) if i + 1 < f_el.shape[0] else 0.0
                print(
                    f"  {i:>4}  {tau_n:>12.4e}  {vel_n:>12.4e}  {rot_n:>10.4e}  {f_n:>12.4e}"
                )
            print(f"  ── End per-segment ───────────────────────────────────────────\n")


# ── Scene ─────────────────────────────────────────────────────────────────────


def createScene(rootNode):
    N = N_SEGMENTS
    h = BEAM_LENGTH / N

    m_node = mass_per_node(BEAM_LENGTH, N, RADIUS, DENSITY)
    I_seg  = inertia_per_seg(BEAM_LENGTH, N, RADIUS, DENSITY)

    rootNode.addObject("RequiredPlugin", pluginName=["Cosserat"])
    rootNode.addObject(
        "VisualStyle", displayFlags="showVisualModels showBehaviorModels"
    )
    rootNode.gravity.value = [0.0, 0.0, 0.0]
    rootNode.dt.value = DT
    rootNode.addObject("DefaultAnimationLoop")
    rootNode.addObject("DefaultVisualManagerLoop")
    # rootNode.addObject(
    #    "InteractiveCamera",
    #    position=[0.5, 0.5, 2.5],
    #    lookAt=[0.5, 0.2, 0.0],
    #    fieldOfView=45,
    # )

    beamNode = rootNode.addChild("relaxBeam")
    beamNode.addObject("EdgeSetTopologyContainer", name="topology")
    beamNode.addObject("EdgeSetTopologyModifier")

    # Initial positions and identity SO3 elements
    # (will be overwritten by controller at init for non-trivial shapes)
    node_pos_str = " ".join(f"{i * h:.6f} 0 0" for i in range(N + 1))
    identity_orientations_str = " ".join(["0 0 0"] * N)
    state = beamNode.addObject(
        "CosseratIntrinsicState",
        name="state",
        positions=node_pos_str,
        orientations=identity_orientations_str,
        template="Vec3d",
    )

    beamNode.addObject(
        "CosseratTopologyBuilder",
        name="builder",
        totalLength=BEAM_LENGTH,
        nbSections=N,
        radius=RADIUS,
        youngModulus=YOUNG_MOD,
        shearModulus=SHEAR_MOD,
        intrinsicState="@state",
        topology="@topology",
    )

    # ── PainlessBeamForceField via _common.add_painless_beam ──────────────────
    ff = add_painless_beam(beamNode, "@state",
                           E=YOUNG_MOD, G=SHEAR_MOD, r=RADIUS,
                           name="ff")

    beamNode.addObject(
        "StaggeredCosseratMapping",
        name="mapping",
        state="@state",
        drawBeam=True,
        drawFrames=True,
        drawRadius=RADIUS,
        drawAxisLength=0.04,
        color="0.9 0.4 0.1 0.9",
    )

    ctrl = LargeDeformationIntegrator(
        name="ctrl",
        state=state,
        ff=ff,
        N=N,
        mass_per_node=m_node,
        inertia_per_seg=I_seg,
    )
    rootNode.addObject(ctrl)

    # NOTE: the initial deformed shape is applied via ctrl._initialized flag
    # on the first onAnimateBeginEvent (before integration in onAnimateEndEvent).

    print("\n" + "=" * 62)
    print(f"  Staggered Cosserat — Large deformation relaxation")
    print("=" * 62)
    print(f"  Shape : '{INITIAL_SHAPE}'")
    print(f"  L={BEAM_LENGTH} m  N={N}  h={h:.4f} m")
    print(f"  E={YOUNG_MOD:.2e} Pa   G={SHEAR_MOD:.2e} Pa")
    print(f"  SO3: {'ACTIVE' if _HAVE_SO3 else 'DISABLED'}")
    print(
        f"  Debug: summary every {DEBUG_FREQ} steps, full table every {DEBUG_FREQ_FULL} steps."
    )
    print(f"  Keys: [+]/[-] damping | [i] re-init | [p] print state")
    print("=" * 62 + "\n")

    return rootNode
