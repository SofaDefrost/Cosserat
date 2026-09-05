# -*- coding: utf-8 -*-
"""
Staggered Cosserat — Torsion validation
========================================
Author  : Younes Adagolodjo (DEFROST / INRIA / Polytech Lille)
Contact : adagolodjo@protonmail.com
Date    : 2025
Branch  : painless/python-explicit

Purpose
-------
Validates the torsion term GJ of PainlessBeamForceField by applying a
constant external torque at the free end of a cantilever beam and comparing
the equilibrium twist profile to the analytical (Euler-Bernoulli torsion) solution.

Physics
-------
A cantilever beam of length L is clamped at node 0 / segment 0.
A constant external torque  T  (around the beam axis = local X)
is applied at the free-end segment (segment N-1).

Analytical equilibrium:
  phi(s) = T * s / GJ          (twist angle at arc-length s from the clamped end)
  phi_tip = T * L / GJ

Discrete staggered solution (expected):
  phi_i = T * i * h / GJ       (twist of segment i, with phi_0 = 0 clamped)
  where h = L / N

The difference (phi_tip_discrete vs phi_tip_analytical) is O(h) = O(1/N).

Convergence indicator
---------------------
The kinetic energy in SO3 should decrease monotonically (heavy velocity damping).
When  max|omega_dot| < CONVERGENCE_TOL  the scene prints the final comparison
and exits (if running headless) or keeps displaying.

Controls
--------
  [p]         print current twist profile vs analytical solution
  [r]         reset to straight configuration
  [+] / [-]   increase / decrease applied torque T

Running
-------
  runSofa staggered_torsion_validation.py
  runSofa --batch --nbIterations 5000 staggered_torsion_validation.py
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
    circular_stiffness, add_painless_beam, mass_per_node, inertia_per_seg,
)

# ── Physical parameters ───────────────────────────────────────────────────────
BEAM_LENGTH = 1.0      # m
N_SEGMENTS  = 8        # → N+1 = 9 nodes
RADIUS      = 0.015    # m  cross-section radius
YOUNG_MOD   = 5.0e5    # Pa
SHEAR_MOD   = 1.65e5   # Pa  (≈ E / (2*(1+ν)), ν ≈ 0.5)
DENSITY     = 500.0    # kg/m³
DT          = 5e-4     # s

# ── Torsion load ──────────────────────────────────────────────────────────────
T_APPLIED   = 5e-3     # N·m  external torque at free end (around beam axis = X)

# ── Damping ───────────────────────────────────────────────────────────────────
DAMPING_ANG = 0.95     # per-step angular velocity damping (strong → fast convergence)

# ── Debug ─────────────────────────────────────────────────────────────────────
DEBUG_FREQ        = 200   # summary every N steps
CONVERGENCE_TOL   = 1e-6  # max |angular_velocity| to declare equilibrium [rad/s]


# ── Derived stiffness and inertia ─────────────────────────────────────────────

def _section_properties(r):
    A   = math.pi * r**2
    I_y = math.pi * r**4 / 4.0
    J   = math.pi * r**4 / 2.0          # polar moment of area
    return A, I_y, J


# ── Helpers : section_properties, mass_per_node, inertia_per_seg via _common ──

# ── Controller ────────────────────────────────────────────────────────────────


class TorsionController(Sofa.Core.Controller):
    """
    Explicit Euler integrator for SO(3) DOFs under applied end torque.

    Positions stay straight (no gravity, no transverse forces from pure torsion).
    Only SO3 orientations are integrated.

    Event ordering:
      onAnimateEndEvent — PainlessBeamForceField has already computed elastic
                          torques at AnimateBeginEvent.

    Applied torque:
      T_APPLIED * e1  in the body frame of the free-end segment (segment N-1).
      For a rotation around the beam axis X, the body-frame X-axis remains
      aligned with world X, so the torque stays constant regardless of twist.
    """

    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        self.state     = kwargs["state"]
        self.ff        = kwargs["ff"]
        self.N         = kwargs["N"]
        self.I         = kwargs["I_seg"]        # kg·m²
        self.GJ        = kwargs["GJ"]           # N·m²
        self.t_applied = kwargs["t_applied"]    # float, N·m

        self._vel_ang = np.zeros((self.N, 3))   # angular velocities in so(3)
        self._step    = 0
        self._converged = False

        h = BEAM_LENGTH / self.N
        # Analytical equilibrium: phi_i = T * i * h / GJ
        self._phi_analytical = np.array(
            [self.t_applied * i * h / self.GJ for i in range(self.N)]
        )
        self._phi_tip_analytical = self.t_applied * BEAM_LENGTH / self.GJ

    # ── Keyboard ─────────────────────────────────────────────────────────────

    def onKeypressedEvent(self, event):
        key = event["key"]
        if key == "p":
            self._print_twist_profile()
        elif key == "r":
            self._reset()
        elif key == "+":
            self.t_applied = min(self.t_applied * 2.0, 1.0)
            self._update_analytical()
            print(f"  T = {self.t_applied:.3e} N·m   phi_tip_analytical = "
                  f"{math.degrees(self._phi_tip_analytical):.2f}°")
        elif key == "-":
            self.t_applied = max(self.t_applied / 2.0, 1e-6)
            self._update_analytical()
            print(f"  T = {self.t_applied:.3e} N·m   phi_tip_analytical = "
                  f"{math.degrees(self._phi_tip_analytical):.2f}°")

    def _update_analytical(self):
        h = BEAM_LENGTH / self.N
        self._phi_analytical = np.array(
            [self.t_applied * i * h / self.GJ for i in range(self.N)]
        )
        self._phi_tip_analytical = self.t_applied * BEAM_LENGTH / self.GJ
        self._converged = False

    def _reset(self):
        self.state.orientations.value = [[0.0, 0.0, 0.0]] * self.N
        self._vel_ang[:] = 0.0
        self._step = 0
        self._converged = False
        print("  [reset] Beam returned to zero-twist configuration.")

    # ── Diagnostics ───────────────────────────────────────────────────────────

    def _print_twist_profile(self):
        """Print per-segment twist vs analytical solution."""
        if not _HAVE_SO3:
            print("  [torsion] SO3 not available.")
            return

        omegas = np.array(self.state.orientations.value).reshape(-1, 3)
        phi_sim = omegas[:, 0]   # X-component = torsion angle (exact for pure torsion)

        h = BEAM_LENGTH / self.N
        print(f"\n  ── Torsion profile  (step={self._step}  T={self.t_applied:.3e} N·m) ──")
        print(f"  {'seg':>4}  {'s_center[m]':>12}  {'phi_sim[°]':>12}  {'phi_ana[°]':>12}  {'err%':>8}")
        print(f"  {'----':>4}  {'----------':>12}  {'----------':>12}  {'----------':>12}  {'----':>8}")
        for i in range(self.N):
            s   = i * h                        # arc-length of seg i (staggered convention)
            pa  = math.degrees(self._phi_analytical[i])
            ps  = math.degrees(phi_sim[i])
            err = (ps - pa) / pa * 100.0 if abs(pa) > 1e-10 else 0.0
            print(f"  {i:>4}  {s:>12.4f}  {ps:>12.4f}  {pa:>12.4f}  {err:>8.2f}")

        phi_free = math.degrees(phi_sim[-1])
        phi_ana_tip = math.degrees(self._phi_tip_analytical)
        print(f"\n  Tip (seg N-1) twist : {phi_free:.4f}°")
        print(f"  Analytical (T*L/GJ) : {phi_ana_tip:.4f}°")
        print(f"  Discrete / analytical = {phi_free / phi_ana_tip:.4f}"
              f"  (expected ≈ {(self.N-1)/self.N:.4f} = (N-1)/N)")
        print(f"  ── End torsion profile ──────────────────────────────────────\n")

    # ── Animation step ────────────────────────────────────────────────────────

    def onAnimateEndEvent(self, event):
        if not _HAVE_SO3:
            return

        dt = float(self.getContext().dt.value)

        # ── Read elastic torques (from PainlessBeamForceField) ────────────────
        tau_el = np.array(self.ff.segmentTorques.value).reshape(-1, 3)
        if tau_el.shape[0] != self.N:
            print(f"  [torsion] tau_el.shape={tau_el.shape}, expected ({self.N}, 3)")
            return

        # ── First-step sanity check ───────────────────────────────────────────
        if self._step == 0:
            print("\n── [TorsionController] FIRST STEP SANITY CHECK ──────────────────")
            print(f"  tau_el shape : {tau_el.shape}  (expected ({self.N}, 3))")
            print(f"  max|tau_el|  : {np.linalg.norm(tau_el, axis=1).max():.4e} N·m"
                  f"  [≈ 0 at rest with straight beam]")
            print(f"  T_applied    : {self.t_applied:.4e} N·m  at segment {self.N-1}")
            print(f"  GJ           : {self.GJ:.4e} N·m²")
            print(f"  phi_tip_ana  : {math.degrees(self._phi_tip_analytical):.4f}°"
                  f"  ({self._phi_tip_analytical:.4f} rad)")
            print(f"  I_seg        : {self.I:.4e} kg·m²   dt={dt:.2e} s")
            # Stability check for torsion
            omega_torsion = math.sqrt(self.GJ / (self.I * (BEAM_LENGTH / self.N)**2))
            print(f"  ω_torsion    : {omega_torsion:.1f} rad/s   (Nyquist dt < {2/omega_torsion:.2e} s)")
            print(f"  dt / dt_Nyquist = {dt * omega_torsion / 2:.3f}   [must be < 1 for stability]")
            print("── End first-step check ─────────────────────────────────────────\n")

        self._step += 1

        # ── Total torque = elastic + external applied at free end ─────────────
        tau_total = tau_el.copy()
        tau_total[self.N - 1, 0] += self.t_applied   # X-component = torsion

        # ── Integrate SO3 (segments 1..N-1, segment 0 clamped) ───────────────
        omegas = np.array(self.state.orientations.value).reshape(-1, 3)
        R_list = [SO3.exp(omegas[i]) for i in range(self.N)]
        new_R  = list(R_list)

        max_vel_ang = 0.0
        for i in range(1, self.N):          # segment 0 clamped (identity, fixed)
            alpha             = tau_total[i] / self.I
            self._vel_ang[i]  = self._vel_ang[i] * DAMPING_ANG + alpha * dt
            dR                = SO3.exp(self._vel_ang[i] * dt)
            new_R[i]          = R_list[i] * dR
            max_vel_ang       = max(max_vel_ang, np.linalg.norm(self._vel_ang[i]))

        self.state.orientations.value = [list(new_R[i].log()) for i in range(self.N)]

        # ── Convergence detection ─────────────────────────────────────────────
        if not self._converged and max_vel_ang < CONVERGENCE_TOL:
            self._converged = True
            print(f"\n  ✓ EQUILIBRIUM REACHED at step {self._step}"
                  f"  (t={self._step * DT:.3f} s)  max|ω|={max_vel_ang:.2e} rad/s")
            self._print_twist_profile()

        # ── Summary print ─────────────────────────────────────────────────────
        if self._step % DEBUG_FREQ == 0:
            omegas_new = np.array(self.state.orientations.value).reshape(-1, 3)
            phi_segs   = omegas_new[:, 0]   # torsion angles [rad]

            tau_x   = tau_el[:, 0]          # elastic torsion torques only
            max_tau = np.linalg.norm(tau_el, axis=1).max()
            KE_ang  = 0.5 * self.I * np.sum(np.linalg.norm(self._vel_ang, axis=1)**2)

            print(f"\n── [step {self._step:6d}  t={self._step * dt:.4f} s] "
                  f"T={self.t_applied:.3e} N·m ──────────────────────────────")
            print(f"  phi_free_end  : {math.degrees(phi_segs[-1]):>10.4f}°  "
                  f"(sim)    vs  {math.degrees(self._phi_tip_analytical):>8.4f}°  (analytical T*L/GJ)")
            print(f"  ratio sim/ana : {phi_segs[-1] / self._phi_tip_analytical:.4f}"
                  f"   (expected → {(self.N-1)/self.N:.4f} = (N-1)/N at equilibrium)")
            print(f"  max|tau_el|   : {max_tau:.4e} N·m   max|tau_el_x|={abs(tau_x).max():.4e} N·m")
            print(f"  max|vel_ang|  : {max_vel_ang:.4e} rad/s   KE_ang={KE_ang:.4e} J")
            print(f"  Converged     : {'YES ✓' if self._converged else 'no (still relaxing)'}")

            # Per-segment torsion table
            print(f"\n  {'seg':>4}  {'phi_sim[°]':>12}  {'phi_ana[°]':>12}  "
                  f"{'tau_el_x[N·m]':>16}  {'vel_ang_x[r/s]':>16}")
            h = BEAM_LENGTH / self.N
            for i in range(self.N):
                pa  = math.degrees(i * h * self.t_applied / self.GJ)
                ps  = math.degrees(phi_segs[i])
                tx  = tau_x[i]
                vx  = self._vel_ang[i, 0]
                clamped_str = " ← clamped" if i == 0 else ""
                print(f"  {i:>4}  {ps:>12.4f}  {pa:>12.4f}  {tx:>16.4e}  {vx:>16.4e}"
                      f"{clamped_str}")
            print()


# ── Scene builder ─────────────────────────────────────────────────────────────


def createScene(rootNode):
    N  = N_SEGMENTS
    h  = BEAM_LENGTH / N
    Np = N + 1

    # ── Derived quantities (via _common.py) ───────────────────────────────────
    A, I_y, J = _section_properties(RADIUS)
    EA, GA, EIy, EIz, GJ = circular_stiffness(YOUNG_MOD, SHEAR_MOD, RADIUS)

    I_seg  = inertia_per_seg(BEAM_LENGTH, N, RADIUS, DENSITY)
    m_node = mass_per_node(BEAM_LENGTH, N, RADIUS, DENSITY)

    phi_tip_analytical = T_APPLIED * BEAM_LENGTH / GJ

    # ── Plugins ───────────────────────────────────────────────────────────────
    rootNode.addObject("RequiredPlugin", pluginName=["Cosserat"])

    rootNode.addObject(
        "VisualStyle", displayFlags="showVisualModels showBehaviorModels"
    )
    rootNode.gravity.value  = [0.0, 0.0, 0.0]   # no gravity — pure torsion test
    rootNode.dt.value       = DT
    rootNode.addObject("DefaultAnimationLoop")
    rootNode.addObject("DefaultVisualManagerLoop")
    rootNode.addObject(
        "InteractiveCamera",
        position=[0.5, 0.4, 1.8],
        lookAt=[0.5, 0.0, 0.0],
        fieldOfView=45,
    )

    # ── Beam node ─────────────────────────────────────────────────────────────
    beamNode = rootNode.addChild("torsionBeam")
    beamNode.addObject("EdgeSetTopologyContainer", name="topology")
    beamNode.addObject("EdgeSetTopologyModifier")

    node_pos_str  = " ".join(f"{i * h:.6f} 0 0" for i in range(Np))
    identity_str  = " ".join(["0 0 0"] * N)

    state = beamNode.addObject(
        "CosseratIntrinsicState",
        name="state",
        positions=node_pos_str,
        orientations=identity_str,
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

    ff = add_painless_beam(beamNode, "@state",
                           E=YOUNG_MOD, G=SHEAR_MOD, r=RADIUS,
                           name="ff",
                           printEvery=2000)   # C++ summary every 2000 FF calls

    beamNode.addObject(
        "StaggeredCosseratMapping",
        name="mapping",
        state="@state",
        drawBeam=True,
        drawFrames=True,
        drawRadius=RADIUS,
        drawAxisLength=0.05,     # longer arrows to see torsion clearly
        color="0.2 0.7 0.3 0.95",
    )

    # ── Torsion controller ────────────────────────────────────────────────────
    rootNode.addObject(
        TorsionController(
            name="torsionCtrl",
            state=state,
            ff=ff,
            N=N,
            I_seg=I_seg,
            GJ=GJ,
            t_applied=T_APPLIED,
        )
    )

    # ── Startup info ──────────────────────────────────────────────────────────
    print("\n" + "=" * 66)
    print("  Staggered Cosserat — Torsion validation (GJ)")
    print("=" * 66)
    print(f"  L={BEAM_LENGTH} m   N={N}   h={h:.4f} m   r={RADIUS} m")
    print(f"  E={YOUNG_MOD:.2e} Pa   G={SHEAR_MOD:.2e} Pa   ρ={DENSITY} kg/m³")
    print()
    print(f"  Cross-section : A={A:.4e} m²   I_y={I_y:.4e} m⁴   J={J:.4e} m⁴")
    print(f"  Stiffness     : EA={EA:.4e} N   GA={GA:.4e} N")
    print(f"                  EI={EIy:.4e} N·m²   GJ={GJ:.4e} N·m²")
    print(f"  I_seg         : {I_seg:.4e} kg·m²   m/node={m_node:.4e} kg")
    print()
    print(f"  Applied torque : T={T_APPLIED:.4e} N·m  at free end (segment {N-1})")
    print(f"                   around beam axis X (= pure torsion)")
    print()
    print(f"  ── Analytical equilibrium (Euler-Bernoulli torsion) ──")
    print(f"     phi(s)  = T * s / GJ")
    print(f"     phi_tip = T * L / GJ = {T_APPLIED:.4e} / {GJ:.4e}"
          f" = {phi_tip_analytical:.4f} rad = {math.degrees(phi_tip_analytical):.2f}°")
    print()
    print(f"  ── Expected discrete solution ──")
    print(f"     phi_i = T * i * h / GJ     (segment i, with phi_0 = 0 clamped)")
    print(f"     phi_{{N-1}} = (N-1)/N * phi_tip = {(N-1)/N * phi_tip_analytical:.4f} rad"
          f" = {math.degrees((N-1)/N * phi_tip_analytical):.2f}°")
    print(f"     O(h) error = {1/N*100:.1f}% (= 1/N * 100%)")
    print()
    print(f"  ── Stability check ──")
    omega_tor = math.sqrt(GJ / (I_seg * h**2))
    print(f"     ω_torsion ≈ {omega_tor:.1f} rad/s   dt_Nyquist < {2/omega_tor:.2e} s")
    print(f"     Current dt = {DT:.2e} s   (margin x{2/(omega_tor * DT):.0f})")
    print()
    print(f"  Debug: summary every {DEBUG_FREQ} steps")
    print(f"  Convergence: max|ω_ang| < {CONVERGENCE_TOL:.0e} rad/s")
    print(f"  Keys: [p] twist profile | [r] reset | [+]/[-] torque")
    print("=" * 66 + "\n")

    return rootNode
