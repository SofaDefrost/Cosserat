# -*- coding: utf-8 -*-
"""
Staggered Cosserat — Cantilever beam under gravity (dynamic)
=============================================================
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
A cantilever beam of length L=1 m is clamped at node 0 and free at node N.
Gravity acts in the −Y direction.

Integration strategy (explicit Euler — proof of concept)
---------------------------------------------------------
CosseratIntrinsicState does not yet plug into standard SOFA solvers
(the integration of SO3 DOFs is separate from Vec3 DOFs).  This scene
uses a Python controller that performs a hand-rolled explicit Euler step:

    v_{n+1} = v_n + (f_elastic/m + g) * dt    (positions only)
    x_{n+1} = x_n + v_{n+1} * dt

where f_elastic is read from PainlessBeamForceField.d_nodalForces and
g = [0, -9.81, 0] is applied as a body force on each free node.

Segment orientations are *fixed* to SO3::identity for now — this exercises
the linear stretch/shear stiffness (EA, GA) only.  Angular stiffness
(EI, GJ) will be activated once full SO3 integration is implemented.

Controls
--------
  [+] / [-]   increase / decrease the gravity magnitude

Running with runSofa
--------------------
  runSofa staggered_cantilever_gravity.py
"""

import math

import numpy as np
import Sofa

# ── Physical parameters ───────────────────────────────────────────────────────
BEAM_LENGTH = 1.0  # m
N_SEGMENTS = 8  # → N+1 = 9 nodes
RADIUS = 0.02  # m  (cross-section radius)
YOUNG_MOD = 1.0e6  # Pa  (soft beam — large deformation visible)
SHEAR_MOD = 3.3e5  # Pa
DENSITY = 700.0  # kg/m³
DT = 1e-3  # s   (explicit Euler → keep small)
GRAVITY = -9.81  # m/s²  (−Y direction)
DAMPING = 0.98   # velocity damping per step (numerical dissipation)

# ── Debug printing ────────────────────────────────────────────────────────────
DEBUG_FREQ      = 50   # print summary every N steps  (set to 1 for every step)
DEBUG_FREQ_FULL = 200  # print per-node table every N steps


def _compute_mass_per_node(length, N, radius, density):
    """Distribute total beam mass uniformly across N+1 nodes (lumped mass)."""
    volume = math.pi * radius**2 * length
    total_mass = volume * density
    return total_mass / (N + 1)


# ── SOFA Controller ───────────────────────────────────────────────────────────


class StaggeredIntegrator(Sofa.Core.Controller):
    """
    Explicit Euler integrator for staggered Cosserat positions (position DOFs only).

    Reads elastic nodal forces from PainlessBeamForceField, adds gravity,
    and integrates node positions.  Node 0 is clamped (Dirichlet BC).

    NOTE: SO3 segment orientations are fixed to identity — this scene exercises
    the linear stretch/shear stiffness (EA, GA) only.  Angular stiffness (EI, GJ)
    is activated in staggered_cantilever_full.py once SO3 integration is wired.

    Event ordering: uses onAnimateEndEvent so that PainlessBeamForceField
    (AnimateBeginEvent) has already computed forces from x_n before we integrate
    to x_{n+1}.  Using onAnimateBeginEvent would cause a 1-step force lag that
    creates artificial differential deflections → potential zigzag.
    """

    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        self.state      = kwargs["state"]
        self.forcefield = kwargs["forcefield"]
        self.N          = kwargs["N"]
        self.mass       = kwargs["mass_per_node"]  # kg
        self.grav_coeff = 1.0
        self._step      = 0

        # Numpy velocity buffer  (N+1 nodes × 3)
        self._vel = np.zeros((self.N + 1, 3))

        # Rest positions (for tip-deflection computation)
        h = BEAM_LENGTH / self.N
        self._rest_pos = np.array([[i * h, 0.0, 0.0] for i in range(self.N + 1)])

    # ── Key events ────────────────────────────────────────────────────────────

    def onKeypressedEvent(self, event):
        key = event["key"]
        if key == "+":
            self.grav_coeff = min(self.grav_coeff + 0.5, 20.0)
            print(f"  g_eff = {self.grav_coeff * GRAVITY:.2f} m/s²")
        elif key == "-":
            self.grav_coeff = max(self.grav_coeff - 0.5, 0.0)
            print(f"  g_eff = {self.grav_coeff * GRAVITY:.2f} m/s²")
        elif key == "p":
            self._print_state()

    def _print_state(self):
        pos = np.array(self.state.positions.value).reshape(-1, 3)
        tip = pos[-1]
        deflection = np.linalg.norm(tip - self._rest_pos[-1])
        print(f"  Tip : ({tip[0]:.5f}, {tip[1]:.5f}, {tip[2]:.5f}) m"
              f"  |deflection|={deflection*100:.3f} cm")

    # ── Animation step ────────────────────────────────────────────────────────

    def onAnimateEndEvent(self, event):
        # NOTE: AnimateEndEvent fires AFTER PainlessBeamForceField::handleEvent
        # (AnimateBeginEvent), so elastic forces from x_n are already available.
        dt = float(self.getContext().dt.value)

        # ── Read state via SOFA Data fields ──────────────────────────────────
        pos = np.array(self.state.positions.value).reshape(-1, 3)
        if pos.shape[0] != self.N + 1:
            return

        f_el = np.array(self.forcefield.nodalForces.value).reshape(-1, 3)
        if f_el.shape[0] != self.N + 1:
            print(f"  [integrator] ERROR: f_el.shape={f_el.shape},"
                  f" expected ({self.N+1}, 3) — forces not ready yet")
            return

        # ── First-step sanity check ────────────────────────────────────────
        if self._step == 0:
            print("\n── [gravity integrator] FIRST STEP SANITY CHECK ──────────────────")
            print(f"  pos shape   : {pos.shape}  (expected ({self.N+1}, 3))")
            print(f"  f_el shape  : {f_el.shape}  (expected ({self.N+1}, 3))")
            max_f = np.linalg.norm(f_el, axis=1).max()
            print(f"  max|f_el|   : {max_f:.4e} N  [≈ 0 at rest before any deflection]")
            print(f"  dt          : {dt:.4e} s")
            print(f"  m/node      : {self.mass:.4e} kg")
            print(f"  g_eff       : {self.grav_coeff * GRAVITY:.2f} m/s²")
            # Predict first-step acceleration from gravity alone
            a_grav = abs(GRAVITY) * self.grav_coeff
            print(f"  Expected a_grav : {a_grav:.3f} m/s²   Δv={a_grav*dt:.3e} m/s"
                  f"   Δx={0.5*a_grav*dt**2:.3e} m  (per free node, step 1)")
            print("── End first-step check ────────────────────────────────────────────\n")

        self._step += 1

        # ── Explicit Euler integration (nodes 1..N, node 0 clamped) ─────────
        g_y = self.grav_coeff * GRAVITY
        new_pos = pos.copy()
        for i in range(1, self.N + 1):
            f_total    = f_el[i].copy()
            f_total[1] += self.mass * g_y          # gravity in −Y
            acc = f_total / self.mass
            self._vel[i] = self._vel[i] * DAMPING + acc * dt
            new_pos[i]   = pos[i] + self._vel[i] * dt

        self.state.positions.value = new_pos.tolist()

        # ── Debug prints ──────────────────────────────────────────────────────
        if self._step % DEBUG_FREQ == 0:
            tip     = new_pos[-1]
            tip_def = np.linalg.norm(tip - self._rest_pos[-1])
            f_norms = np.linalg.norm(f_el, axis=1)
            v_norms = np.linalg.norm(self._vel, axis=1)

            print(f"\n── [step {self._step:5d}  t={self._step * dt:.4f} s] "
                  f"──────────────────────────────────────────────")
            print(f"  Tip position : ({tip[0]:.5f}, {tip[1]:.5f}, {tip[2]:.5f}) m"
                  f"   |deflection|={tip_def*100:.3f} cm")
            print(f"  max|f_el|    : {f_norms.max():.4e} N  (at node {int(f_norms.argmax())})"
                  f"   f_tip={np.linalg.norm(f_el[-1]):.4e} N")
            print(f"  max|vel|     : {v_norms.max():.4e} m/s  (at node {int(v_norms.argmax())})")
            print(f"  g_eff        : {g_y:.2f} m/s²   damping={DAMPING}")

            # Stability estimate: kinetic energy should not grow unboundedly
            ke = 0.5 * self.mass * np.sum(v_norms**2)
            print(f"  KE approx    : {ke:.4e} J  [should stabilise or grow slowly]")

        if self._step % DEBUG_FREQ_FULL == 0:
            print(f"\n  ── Per-node detail (step {self._step}) ──────────────────────")
            print(f"  {'node':>5}  {'|f_el|':>12}  {'|vel|':>12}  {'pos_y':>10}")
            for i in range(self.N + 1):
                f_n = np.linalg.norm(f_el[i])
                v_n = np.linalg.norm(self._vel[i])
                print(f"  {i:>5}  {f_n:>12.4e}  {v_n:>12.4e}  {new_pos[i][1]:>10.6f}")
            print(f"  ── End per-node ─────────────────────────────────────────────\n")


# ── Scene builder ─────────────────────────────────────────────────────────────


def createScene(rootNode):
    N = N_SEGMENTS
    h = BEAM_LENGTH / N
    Np = N + 1

    mass_per_node = _compute_mass_per_node(BEAM_LENGTH, N, RADIUS, DENSITY)

    # ── Plugins ───────────────────────────────────────────────────────────────
    rootNode.addObject("RequiredPlugin", pluginName=["Cosserat"])

    # ── Visualisation & animation ──────────────────────────────────────────────
    rootNode.addObject(
        "VisualStyle",
        displayFlags="showVisualModels showBehaviorModels showForceFields",
    )
    rootNode.gravity.value = [0.0, 0.0, 0.0]  # gravity handled by controller
    rootNode.dt.value = DT
    rootNode.addObject("DefaultAnimationLoop")
    rootNode.addObject("DefaultVisualManagerLoop")

    # ── Camera ────────────────────────────────────────────────────────────────
    rootNode.addObject(
        "InteractiveCamera",
        position=[0.5, -0.8, 2.0],
        lookAt=[0.5, 0.0, 0.0],
        distance=2.0,
        fieldOfView=45,
    )

    # ── Beam node ─────────────────────────────────────────────────────────────
    beamNode = rootNode.addChild("cantilever")

    # ── Topology ──────────────────────────────────────────────────────────────
    beamNode.addObject("EdgeSetTopologyContainer", name="topology")
    beamNode.addObject("EdgeSetTopologyModifier")

    # ── CosseratIntrinsicState ─────────────────────────────────────────────────
    # Initial configuration: straight beam along X, N identity SO3 rotations
    node_positions_str = " ".join(f"{i * h:.6f} 0 0" for i in range(Np))
    identity_orientations_str = " ".join(["0 0 0"] * N)
    state = beamNode.addObject(
        "CosseratIntrinsicState",
        name="state",
        positions=node_positions_str,
        orientations=identity_orientations_str,
        template="Vec3d",
    )

    # ── CosseratTopologyBuilder ────────────────────────────────────────────────
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

    # ── PainlessBeamForceField ─────────────────────────────────────────────────
    # Stiffness parameters computed explicitly from cross-section geometry.
    # CosseratTopologyBuilder does NOT auto-link its outputs to PainlessBeamForceField,
    # so we compute and pass them here to avoid using wrong constructor defaults.
    _A = math.pi * RADIUS**2
    _I_y = math.pi * RADIUS**4 / 4.0
    _J = math.pi * RADIUS**4 / 2.0
    _EA = YOUNG_MOD * _A
    _GA = SHEAR_MOD * _A
    _EIy = YOUNG_MOD * _I_y
    _EIz = YOUNG_MOD * _I_y
    _GJ = SHEAR_MOD * _J
    print(f"\n  [scene] PainlessBeamForceField stiffness:")
    print(
        f"    EA={_EA:.3e} N   GA={_GA:.3e} N   EIy={_EIy:.3e} N·m²   GJ={_GJ:.3e} N·m²"
    )
    forcefield = beamNode.addObject(
        "PainlessBeamForceField",
        name="ff",
        state="@state",
        EA=_EA,
        GA=_GA,
        EIy=_EIy,
        EIz=_EIz,
        GJ=_GJ,
    )

    # ── StaggeredCosseratMapping ───────────────────────────────────────────────
    beamNode.addObject(
        "StaggeredCosseratMapping",
        name="mapping",
        state="@state",
        drawBeam=True,
        drawFrames=True,
        drawRadius=RADIUS,
        drawAxisLength=0.04,
        color="0.2 0.6 0.95 0.9",
    )

    # ── Explicit Euler integrator ──────────────────────────────────────────────
    rootNode.addObject(
        StaggeredIntegrator(
            name="integrator",
            state=state,
            forcefield=forcefield,
            N=N,
            mass_per_node=mass_per_node,
        )
    )

    # ── Info ───────────────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("  Staggered Cosserat — Cantilever under gravity")
    print("=" * 60)
    print(f"  Length          : {BEAM_LENGTH:.2f} m")
    print(f"  N segments      : {N}  →  {Np} nodes")
    print(f"  h (segment)     : {h:.4f} m")
    print(f"  Young modulus   : {YOUNG_MOD:.2e} Pa")
    print(f"  Shear modulus   : {SHEAR_MOD:.2e} Pa")
    print(f"  Density         : {DENSITY:.1f} kg/m³")
    print(f"  Mass / node     : {mass_per_node:.4f} kg")
    print(f"  dt              : {DT:.2e} s")
    print(f"  Gravity (−Y)    : {GRAVITY:.2f} m/s²")
    print(f"  Note: SO3 orientations are fixed (identity).")
    print(f"        Only EA/GA (linear) forces are active.")
    print(f"  Debug: summary every {DEBUG_FREQ} steps, full table every {DEBUG_FREQ_FULL} steps.")
    print(f"  Keys: [+]/[-] gravity  |  [p] print tip state")
    print("=" * 60 + "\n")

    return rootNode
