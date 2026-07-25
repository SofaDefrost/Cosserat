# -*- coding: utf-8 -*-
"""
Staggered Cosserat — Geometry validation test
==============================================
Author  : Younes Adagolodjo (DEFROST / INRIA / Polytech Lille)
Contact : adagolodjo@protonmail.com
Date    : 2025
Branch  : painless/base-geometry

Purpose
-------
Visual and numerical validation of StaggeredCosseratMapping.

The scene initialises a *straight* cantilever beam along the X axis with
N=5 segments (6 nodes) and all segment orientations equal to SO3::identity.

Expected output (frames in world space)
- Node frames:  centres at x_i = [i*h, 0, 0], orientation = identity quat (0,0,0,1)
- Segment frames: centres at midpoints [(i+0.5)*h, 0, 0], orientation = identity

A Python controller validates these values on the first AnimateEnd event and
prints a PASS / FAIL report to the console.

Running with runSofa
--------------------
  runSofa staggered_geometry_test.py

Or headless:
  runSofa --batch --nbIterations 1 staggered_geometry_test.py
"""

import Sofa

# ── Scene parameters ──────────────────────────────────────────────────────────
BEAM_LENGTH = 1.0  # m
N_SEGMENTS = 5  # number of segments  →  N+1 = 6 nodes
RADIUS = 0.02  # m (for visualisation)
YOUNG_MOD = 1.0e8  # Pa
SHEAR_MOD = 3.3e7  # Pa  (≈ E / (2*(1+0.5)))
DT = 0.01  # s


# ── Validation helper ─────────────────────────────────────────────────────────


def _almost_equal(a, b, tol=1e-9):
    return abs(a - b) < tol


def _check_node_frames(node_frames, h, N, tol=1e-9):
    """Verify N+1 node frames of a straight beam along X."""
    errors = []
    for i in range(N + 1):
        frame = node_frames[i]
        # Position: [i*h, 0, 0]
        expected_pos = [i * h, 0.0, 0.0]
        for j, (got, exp) in enumerate(zip(frame[:3], expected_pos)):
            if not _almost_equal(got, exp, tol):
                errors.append(
                    f"  node[{i}] pos[{j}]: got {got:.6f}, expected {exp:.6f}"
                )
        # Orientation: identity quaternion (x=0, y=0, z=0, w=1)
        quat = frame[3:]  # [qx, qy, qz, qw]
        id_q = [0.0, 0.0, 0.0, 1.0]
        for j, (got, exp) in enumerate(zip(quat, id_q)):
            if not _almost_equal(got, exp, tol):
                errors.append(
                    f"  node[{i}] quat[{j}]: got {got:.6f}, expected {exp:.6f}"
                )
    return errors


def _check_segment_frames(seg_frames, h, N, tol=1e-9):
    """Verify N segment frames of a straight beam along X."""
    errors = []
    for i in range(N):
        frame = seg_frames[i]
        # Position: midpoint of edge i
        expected_pos = [(i + 0.5) * h, 0.0, 0.0]
        for j, (got, exp) in enumerate(zip(frame[:3], expected_pos)):
            if not _almost_equal(got, exp, tol):
                errors.append(f"  seg[{i}] pos[{j}]: got {got:.6f}, expected {exp:.6f}")
        # Orientation: identity quaternion
        quat = frame[3:]
        id_q = [0.0, 0.0, 0.0, 1.0]
        for j, (got, exp) in enumerate(zip(quat, id_q)):
            if not _almost_equal(got, exp, tol):
                errors.append(
                    f"  seg[{i}] quat[{j}]: got {got:.6f}, expected {exp:.6f}"
                )
    return errors


# ── SOFA Controller ───────────────────────────────────────────────────────────


class GeometryValidator(Sofa.Core.Controller):
    """
    Runs once on the first AnimateEnd event to validate StaggeredCosseratMapping
    output against analytic values for a straight beam.
    """

    def __init__(self, *args, **kwargs):
        Sofa.Core.Controller.__init__(self, *args, **kwargs)
        self.mapping = kwargs["mapping"]
        self.h = kwargs["h"]
        self.N = kwargs["N"]
        self._done = False

    def onAnimateEndEvent(self, event):
        if self._done:
            return
        self._done = True

        print("\n" + "=" * 60)
        print("  StaggeredCosseratMapping — Geometry validation")
        print("=" * 60)

        try:
            node_frames = self.mapping.nodeFrames.value
            seg_frames = self.mapping.segmentFrames.value
        except Exception as e:
            print(f"  [ERROR] Cannot read mapping output: {e}")
            return

        print(f"  N segments  = {self.N}")
        print(f"  Beam length = {self.N * self.h:.3f} m")
        print(f"  h (segment) = {self.h:.4f} m\n")

        print("  Node frames (centre, quat):")
        for i, f in enumerate(node_frames):
            print(
                f"    node[{i}]  pos=({f[0]:.4f}, {f[1]:.4f}, {f[2]:.4f})"
                f"  quat=({f[3]:.4f}, {f[4]:.4f}, {f[5]:.4f}, {f[6]:.4f})"
            )

        print("\n  Segment frames (midpoint, quat):")
        for i, f in enumerate(seg_frames):
            print(
                f"    seg[{i}]   pos=({f[0]:.4f}, {f[1]:.4f}, {f[2]:.4f})"
                f"  quat=({f[3]:.4f}, {f[4]:.4f}, {f[5]:.4f}, {f[6]:.4f})"
            )

        # Numerical validation
        node_errors = _check_node_frames(node_frames, self.h, self.N)
        seg_errors = _check_segment_frames(seg_frames, self.h, self.N)

        all_errors = node_errors + seg_errors
        print()
        if not all_errors:
            print("  ✅  ALL CHECKS PASSED")
        else:
            print(f"  ❌  {len(all_errors)} CHECK(S) FAILED:")
            for e in all_errors:
                print(e)
        print("=" * 60 + "\n")


# ── Scene builder ─────────────────────────────────────────────────────────────


def createScene(rootNode):
    h = BEAM_LENGTH / N_SEGMENTS  # segment rest length
    Np = N_SEGMENTS + 1  # number of nodes

    # ── Plugins ───────────────────────────────────────────────────────────────

    #    "Cosserat",
    #    "Sofa.Component.AnimationLoop",
    #    "Sofa.Component.LinearSolver.Iterative",
    #    "Sofa.Component.LinearSolver.Direct",

    rootNode.addObject("RequiredPlugin", pluginName=["Cosserat"])

    # ── Display ───────────────────────────────────────────────────────────────
    rootNode.addObject("VisualStyle", displayFlags="showVisualModels")
    rootNode.gravity.value = [0.0, 0.0, 0.0]
    rootNode.dt.value = DT
    rootNode.addObject("DefaultAnimationLoop")

    # ── Beam node ─────────────────────────────────────────────────────────────
    beamNode = rootNode.addChild("straightBeam")

    # ── Topology ──────────────────────────────────────────────────────────────
    beamNode.addObject("EdgeSetTopologyContainer", name="topology")
    beamNode.addObject("EdgeSetTopologyModifier")

    # ── CosseratIntrinsicState ─────────────────────────────────────────────────
    # Positions: N+1 nodes along X, orientations: N identity SO3 (default)
    # node_positions_str = []
    # for i in range(Np):
    # node_positions_str.append([i * 0.2, 0, 0])

    node_positions_str = " ".join(f"{i * h:0.6f} 0 0" for i in range(Np))
    # N identity SO3 elements — format "x y z w" per element (our operator>>)
    # Identity quaternion = (x=0, y=0, z=0, w=1)
    identity_orientations_str = " ".join(["0 0 0"] * N_SEGMENTS)

    print(f"Nodes positions  : {node_positions_str}")
    print(f"Segment SO3 count: {N_SEGMENTS}")
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
        nbSections=N_SEGMENTS,
        radius=RADIUS,
        youngModulus=YOUNG_MOD,
        shearModulus=SHEAR_MOD,
        intrinsicState="@state",
        topology="@topology",
    )

    # ── StaggeredCosseratMapping ───────────────────────────────────────────────
    mapping = beamNode.addObject(
        "StaggeredCosseratMapping",
        name="mapping",
        state="@state",
        drawBeam=True,
        drawFrames=True,
        drawRadius=RADIUS,
        color="0.2 0.5 0.9 0.85",
    )

    # ── Geometry validator ─────────────────────────────────────────────────────
    rootNode.addObject(
        GeometryValidator(name="validator", mapping=mapping, h=h, N=N_SEGMENTS)
    )

    return rootNode
