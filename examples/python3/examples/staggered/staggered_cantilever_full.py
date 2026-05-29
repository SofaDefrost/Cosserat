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
  from Sofa.LieGroups import SO3      # built by python/Binding/LieGroups
  from Sofa.Cosserat import (CosseratIntrinsicState,
                              PainlessBeamForceField)  # python/Binding/Cosserat
"""

import numpy as np
import Sofa

# ── Bindings pybind11 ─────────────────────────────────────────────────────────
#  L'import du module 'Cosserat' enregistre le downcast SOFA Python3 qui permet
#  d'appeler state.getPositions(), ff.getNodalForces(), etc. sur les proxies.
#  Il est chargé automatiquement via RequiredPlugin('Cosserat'), mais un import
#  explicite garantit la disponibilité avant le premier accès Python.
try:
    from Sofa import Cosserat as _CosseratModule  # noqa — enregistre le downcast

    _HAVE_BINDINGS = True
except ImportError:
    _HAVE_BINDINGS = False
    print(
        "[WARNING] Module pybind11 'Cosserat' non trouvé. "
        "Compiler avec -DCOSSERAT_WITH_PYTHON_BINDINGS=ON."
    )

try:
    from Sofa.LieGroups import SO3  # pybind11 module (Sofa/LieGroups)

    _HAVE_SO3 = True
except ImportError:
    _HAVE_SO3 = False
    print(
        "[WARNING] LieGroups module not found — SO3 integration disabled. "
        "Run with the Cosserat plugin compiled with Python bindings."
    )

# ── Physical parameters ───────────────────────────────────────────────────────
BEAM_LENGTH = 1.0  # m
N_SEGMENTS = 8
RADIUS = 0.015  # m  cross-section radius
YOUNG_MOD = 5.0e5  # Pa  (soft, to see large deflections quickly)
SHEAR_MOD = 1.65e5  # Pa  (≈ E / (2*(1+ν)), ν≈0.5)
DENSITY = 500.0  # kg/m³
DT = 5e-4  # s   explicit Euler → keep small!
DAMPING_POS = 0.995  # velocity damping (positions)
DAMPING_ANG = 0.990  # angular velocity damping
GRAVITY_Y = -9.81  # m/s²


def _mass_per_node(L, N, r, rho):
    total = np.pi * r**2 * L * rho
    return total / (N + 1)


def _inertia_per_seg(L, N, r, rho):
    """Scalar moment of inertia for a cylindrical segment (rotation about diameter)."""
    L_seg = L / N
    vol = np.pi * r**2 * L_seg
    mass = vol * rho
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

        self.state = kwargs["state"]
        self.ff = kwargs["ff"]
        self.N = kwargs["N"]
        self.m = kwargs["mass_per_node"]
        self.I = kwargs["inertia_per_seg"]
        self.g_coeff = 1.0
        self._has_so3 = _HAVE_SO3

        # Velocity buffers
        self._vel_pos = np.zeros((self.N + 1, 3))  # translational  (world)
        self._vel_ang = np.zeros((self.N, 3))  # angular (so3, body frame)

        # Rest positions (for reset)
        h = BEAM_LENGTH / self.N
        self._rest_pos = np.array([[i * h, 0.0, 0.0] for i in range(self.N + 1)])

    # ── Key events ────────────────────────────────────────────────────────────

    def onKeypressedEvent(self, event):
        key = event["key"]
        if key == "+":
            self.g_coeff = min(self.g_coeff + 0.5, 20.0)
            print(f"  g_eff = {self.g_coeff * GRAVITY_Y:.2f} m/s²")
        elif key == "-":
            self.g_coeff = max(self.g_coeff - 0.5, 0.0)
            print(f"  g_eff = {self.g_coeff * GRAVITY_Y:.2f} m/s²")
        elif key == "r":
            self._reset()
        elif key == "p":
            self._print_state()

    def _reset(self):
        """Reset beam to straight, stress-free configuration."""
        self.state.positions.value = self._rest_pos.tolist()
        if self._has_so3:
            # Identity SO3 → rotation vector = [0, 0, 0]
            self.state.orientations.value = [[0.0, 0.0, 0.0]] * self.N
        self._vel_pos[:] = 0.0
        self._vel_ang[:] = 0.0
        print("  [reset] Beam returned to straight configuration.")

    def _print_state(self):
        pos = np.array(self.state.positions.value).reshape(-1, 3)
        tip = pos[-1]
        print(f"  Tip position : [{tip[0]:.4f}, {tip[1]:.4f}, {tip[2]:.4f}] m")
        deflection = np.linalg.norm(tip - self._rest_pos[-1])
        print(f"  Tip deflection : {deflection * 100:.2f} cm")
        if self._has_so3:
            omegas = np.array(self.state.orientations.value).reshape(-1, 3)
            tip_R = SO3.exp(omegas[-1])
            q = tip_R.quaternion()
            print(
                f"  Tip segment quat : [{q[0]:.4f}, {q[1]:.4f}, {q[2]:.4f}, {q[3]:.4f}]"
            )

    # ── Animation step ────────────────────────────────────────────────────────

    def onAnimateEndEvent(self, event):
        # NOTE: we use AnimateEndEvent so that PainlessBeamForceField::handleEvent
        # (AnimateBeginEvent, fires first) has already computed forces from the
        # current positions x_N before we integrate to x_{N+1}.
        # Wrong ordering (both on AnimateBegin) caused a 1-step force lag that
        # created artificial differential deflections → zigzag instability.
        dt = float(self.getContext().dt.value)

        # ── Read current positions via SOFA Data field ─────────────────────────
        pos = np.array(self.state.positions.value).reshape(-1, 3)
        if pos.shape[0] != self.N + 1:
            return

        # ── Read elastic forces via SOFA Data fields ──────────────────────────
        f_el   = np.array(self.ff.nodalForces.value).reshape(-1, 3)
        tau_el = np.array(self.ff.segmentTorques.value).reshape(-1, 3)

        if f_el.shape[0] != self.N + 1:
            return

        # ── Debug print every 200 steps ───────────────────────────────────────
        if not hasattr(self, '_step'):
            self._step = 0
        self._step += 1
        if self._step % 200 == 1:
            tip   = pos[-1]
            f_tip = f_el[-1]
            tau_mid = tau_el[self.N // 2] if tau_el.shape[0] > 0 else np.zeros(3)
            so3_status = "ACTIVE" if self._has_so3 else "DISABLED"
            ang_norm = np.linalg.norm(self._vel_ang)
            print(f"  step={self._step:5d}  tip=({tip[0]:.4f},{tip[1]:.4f},{tip[2]:.4f})"
                  f"  |f_tip|={np.linalg.norm(f_tip):.3e}"
                  f"  |tau_mid|={np.linalg.norm(tau_mid):.3e}"
                  f"  SO3={so3_status}  |vel_ang|={ang_norm:.3e}")

        # ── Integrate node positions (explicit Euler) ─────────────────────────
        g_eff = self.g_coeff * GRAVITY_Y
        new_pos = pos.copy()

        for i in range(1, self.N + 1):  # node 0 clamped
            f_total    = f_el[i].copy()
            f_total[1] += self.m * g_eff  # gravity in −Y
            acc = f_total / self.m
            self._vel_pos[i] = self._vel_pos[i] * DAMPING_POS + acc * dt
            new_pos[i] = pos[i] + self._vel_pos[i] * dt

        # Write back via SOFA Data field
        self.state.positions.value = new_pos.tolist()

        # ── Integrate SO3 orientations (explicit Euler on SO(3)) ──────────────
        if self._has_so3:
            # orientations.value → list of Vec3d = rotation vectors ω = log(R)
            omegas = np.array(self.state.orientations.value).reshape(-1, 3)
            # Reconstruct SO3 from rotation vectors
            R_list = [SO3.exp(omegas[i]) for i in range(self.N)]

            new_R = list(R_list)
            for i in range(1, self.N):  # segment 0 clamped
                tau_i = tau_el[i]
                alpha = tau_i / self.I
                self._vel_ang[i] = self._vel_ang[i] * DAMPING_ANG + alpha * dt
                dR = SO3.exp(self._vel_ang[i] * dt)
                new_R[i] = R_list[i] * dR

            # Write back as rotation vectors (log of each SO3)
            new_omegas = [list(new_R[i].log()) for i in range(self.N)]
            self.state.orientations.value = new_omegas


# ── Scene builder ─────────────────────────────────────────────────────────────


def createScene(rootNode):
    N = N_SEGMENTS
    h = BEAM_LENGTH / N
    Np = N + 1

    m_node = _mass_per_node(BEAM_LENGTH, N, RADIUS, DENSITY)
    I_seg = _inertia_per_seg(BEAM_LENGTH, N, RADIUS, DENSITY)

    # ── Plugins ───────────────────────────────────────────────────────────────
    rootNode.addObject("RequiredPlugin", pluginName=["Cosserat"])

    rootNode.addObject(
        "VisualStyle", displayFlags="showVisualModels showBehaviorModels"
    )
    rootNode.gravity.value = [0.0, 0.0, 0.0]
    rootNode.dt.value = DT
    rootNode.addObject("DefaultAnimationLoop")
    rootNode.addObject("DefaultVisualManagerLoop")
    rootNode.addObject(
        "InteractiveCamera",
        position=[0.5, -0.8, 2.5],
        lookAt=[0.5, 0.0, 0.0],
        fieldOfView=45,
    )

    # ── Beam node ─────────────────────────────────────────────────────────────
    beamNode = rootNode.addChild("staggeredBeam")
    beamNode.addObject("EdgeSetTopologyContainer", name="topology")
    beamNode.addObject("EdgeSetTopologyModifier")

    node_pos_str = " ".join(f"{i * h:.6f} 0 0" for i in range(Np))
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

    ff = beamNode.addObject("PainlessBeamForceField", name="ff", state="@state")

    beamNode.addObject(
        "StaggeredCosseratMapping",
        name="mapping",
        state="@state",
        drawBeam=True,
        drawFrames=True,
        drawRadius=RADIUS,
        drawAxisLength=0.03,
        color="0.15 0.65 0.95 0.9",
    )

    # ── Python integrator ─────────────────────────────────────────────────────
    rootNode.addObject(
        StaggeredFullIntegrator(
            name="integrator",
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
    print(
        f"  SO3 integration: {'ACTIVE (LieGroups module found)' if _HAVE_SO3 else 'DISABLED (import LieGroups failed)'}"
    )
    print(f"  Keys: [+]/[-] gravity  |  [r] reset  |  [p] print state")
    print("=" * 62 + "\n")

    return rootNode
