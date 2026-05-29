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

import Sofa
import math

# ── Physical parameters ───────────────────────────────────────────────────────
BEAM_LENGTH  = 1.0      # m
N_SEGMENTS   = 8        # → N+1 = 9 nodes
RADIUS       = 0.02     # m  (cross-section radius)
YOUNG_MOD    = 1.0e6    # Pa  (soft beam — large deformation visible)
SHEAR_MOD    = 3.3e5    # Pa
DENSITY      = 700.0    # kg/m³
DT           = 1e-3     # s   (explicit Euler → keep small)
GRAVITY      = -9.81    # m/s²  (−Y direction)
DAMPING      = 0.98     # velocity damping per step (numerical dissipation)


def _compute_mass_per_node(length, N, radius, density):
    """Distribute total beam mass uniformly across N+1 nodes (lumped mass)."""
    volume = math.pi * radius**2 * length
    total_mass = volume * density
    return total_mass / (N + 1)


# ── SOFA Controller ───────────────────────────────────────────────────────────

class StaggeredIntegrator(Sofa.Core.Controller):
    """
    Explicit Euler integrator for staggered Cosserat positions.

    Reads elastic nodal forces from PainlessBeamForceField, adds gravity,
    and integrates node positions.  Node 0 is clamped (Dirichlet BC).

    NOTE: SO3 segment orientations are not integrated here.  The beam
    behaves as a linear-spring chain until full SO3 dynamics is wired.
    """

    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        self.state        = kwargs['state']
        self.forcefield   = kwargs['forcefield']
        self.N            = kwargs['N']
        self.mass         = kwargs['mass_per_node']   # kg
        self.grav_coeff   = 1.0

        # Velocity buffer  (N+1 vectors)
        self._vel = [[0.0, 0.0, 0.0] for _ in range(self.N + 1)]

    # ── Key events ────────────────────────────────────────────────────────────

    def onKeypressedEvent(self, event):
        key = event['key']
        if key == '+':
            self.grav_coeff = min(self.grav_coeff + 0.5, 20.0)
            print(f"  Gravity coefficient: {self.grav_coeff:.1f}  "
                  f"→  g_eff = {self.grav_coeff * GRAVITY:.2f} m/s²")
        elif key == '-':
            self.grav_coeff = max(self.grav_coeff - 0.5, 0.0)
            print(f"  Gravity coefficient: {self.grav_coeff:.1f}  "
                  f"→  g_eff = {self.grav_coeff * GRAVITY:.2f} m/s²")

    # ── Animation step ────────────────────────────────────────────────────────

    def onAnimateBeginEvent(self, event):
        dt = float(self.getContext().dt.value)

        # Read current positions
        try:
            pos = list(self.state.positions.value)
        except Exception as e:
            print(f"[StaggeredIntegrator] Cannot read positions: {e}")
            return

        # Read elastic forces produced by PainlessBeamForceField
        try:
            f_elastic = list(self.forcefield.nodalForces.value)
        except Exception as e:
            # ForceField has not yet computed forces — skip this step
            return

        if len(f_elastic) != self.N + 1:
            # Not yet initialised — skip
            return

        # ── Explicit Euler integration ─────────────────────────────────────────
        g_y = self.grav_coeff * GRAVITY   # effective gravity m/s²

        new_pos = [list(p) for p in pos]

        for i in range(1, self.N + 1):   # node 0 is clamped
            # Total force = elastic + gravity (body force)
            fx = f_elastic[i][0]
            fy = f_elastic[i][1] + self.mass * g_y
            fz = f_elastic[i][2]

            # Acceleration
            ax = fx / self.mass
            ay = fy / self.mass
            az = fz / self.mass

            # Velocity update (with damping)
            self._vel[i][0] = self._vel[i][0] * DAMPING + ax * dt
            self._vel[i][1] = self._vel[i][1] * DAMPING + ay * dt
            self._vel[i][2] = self._vel[i][2] * DAMPING + az * dt

            # Position update
            new_pos[i][0] += self._vel[i][0] * dt
            new_pos[i][1] += self._vel[i][1] * dt
            new_pos[i][2] += self._vel[i][2] * dt

        # Write updated positions back
        try:
            flat = []
            for p in new_pos:
                flat.extend(p)
            self.state.positions.value = flat
        except Exception as e:
            print(f"[StaggeredIntegrator] Cannot write positions: {e}")

    def onAnimateEndEvent(self, event):
        # Optionally print tip displacement
        pass


# ── Scene builder ─────────────────────────────────────────────────────────────

def createScene(rootNode):
    N  = N_SEGMENTS
    h  = BEAM_LENGTH / N
    Np = N + 1

    mass_per_node = _compute_mass_per_node(BEAM_LENGTH, N, RADIUS, DENSITY)

    # ── Plugins ───────────────────────────────────────────────────────────────
    rootNode.addObject('RequiredPlugin', name='CosseratPlugin',
                       pluginName='Cosserat')
    rootNode.addObject('RequiredPlugin', name='SofaBaseTopology',
                       pluginName='SofaBaseTopology')

    # ── Visualisation & animation ──────────────────────────────────────────────
    rootNode.addObject('VisualStyle',
                       displayFlags='showVisualModels showBehaviorModels '
                                    'showForceFields')
    rootNode.gravity.value = [0., 0., 0.]   # gravity handled by controller
    rootNode.dt.value      = DT
    rootNode.addObject('DefaultAnimationLoop')
    rootNode.addObject('DefaultVisualManagerLoop')

    # ── Camera ────────────────────────────────────────────────────────────────
    rootNode.addObject('InteractiveCamera',
                       position=[0.5, -0.8, 2.0],
                       lookAt=[0.5, 0.0, 0.0],
                       distance=2.0,
                       fieldOfView=45)

    # ── Beam node ─────────────────────────────────────────────────────────────
    beamNode = rootNode.addChild('cantilever')

    # ── Topology ──────────────────────────────────────────────────────────────
    beamNode.addObject('EdgeSetTopologyContainer', name='topology')
    beamNode.addObject('EdgeSetTopologyModifier')

    # ── CosseratIntrinsicState ─────────────────────────────────────────────────
    # Initial configuration: straight beam along X, N identity SO3 rotations
    node_positions_str = ' '.join(f'{i * h:.6f} 0 0' for i in range(Np))
    identity_orientations_str = ' '.join(['0 0 0'] * N)
    state = beamNode.addObject('CosseratIntrinsicState',
                               name='state',
                               positions=node_positions_str,
                               orientations=identity_orientations_str,
                               template='Vec3d')

    # ── CosseratTopologyBuilder ────────────────────────────────────────────────
    beamNode.addObject('CosseratTopologyBuilder',
                       name='builder',
                       totalLength=BEAM_LENGTH,
                       nbSections=N,
                       radius=RADIUS,
                       youngModulus=YOUNG_MOD,
                       shearModulus=SHEAR_MOD,
                       intrinsicState='@state',
                       topology='@topology')

    # ── PainlessBeamForceField ─────────────────────────────────────────────────
    # Stiffness parameters will be read from CosseratTopologyBuilder outputs
    # via the linked intrinsic state (the builder writes EA, GA, etc. to it).
    forcefield = beamNode.addObject('PainlessBeamForceField',
                                    name='ff',
                                    state='@state')

    # ── StaggeredCosseratMapping ───────────────────────────────────────────────
    beamNode.addObject('StaggeredCosseratMapping',
                       name='mapping',
                       state='@state',
                       drawBeam=True,
                       drawFrames=True,
                       drawRadius=RADIUS,
                       drawAxisLength=0.04,
                       color='0.2 0.6 0.95 0.9')

    # ── Explicit Euler integrator ──────────────────────────────────────────────
    rootNode.addObject(StaggeredIntegrator(name='integrator',
                                           state=state,
                                           forcefield=forcefield,
                                           N=N,
                                           mass_per_node=mass_per_node))

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
    print(f"  Keys: [+] / [-]  to adjust gravity coefficient.")
    print("=" * 60 + "\n")

    return rootNode
