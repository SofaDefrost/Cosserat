#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Benchmark : convergence et performance de la discrétisation Staggered (Painless)
=================================================================================
Auteur  : Younes Adagolodjo (DEFROST / INRIA / Polytech Lille)
Contact : adagolodjo@protonmail.com
Branche : painless/base-geometry
Référence : Romanyà-Serrasolsas et al., SIGGRAPH 2025. DOI: 10.1145/3730944

Usage
-----
  # Depuis le répertoire de build SOFA (pour trouver les .so) :
  python3 benchmark_staggered_convergence.py

  # Avec matplotlib pour les graphiques :
  python3 benchmark_staggered_convergence.py --plot

  # Sauvegarder les résultats CSV :
  python3 benchmark_staggered_convergence.py --csv results.csv

Description
-----------
Ce script évalue la discrétisation staggered sur le problème canonique du
cantilever encastré-libre sous charge répartie (gravité).

Solution analytique de référence (Euler-Bernoulli, régime petit déplacement) :
    δ_EB = q · L⁴ / (8 · E · I)   où  q = ρ · π · r² · |g|

Pour chaque N ∈ [4, 8, 16, 32, 64] :
  1. Initialisation du beam staggered via pybind11
     (CosseratIntrinsicState + PainlessBeamForceField)
  2. Intégration explicite Euler sur R³×SO(3) avec amortissement fort
     jusqu'à convergence (critère : |Δδ_tip| / δ_tip < TOL)
  3. Mesures :
     • δ_tip numérique
     • Erreur relative vs. δ_EB
     • Nombre de pas jusqu'à convergence
     • Temps CPU total et par pas

Résultats attendus (régime petit déplacement) :
  • Convergence en O(h²) pour la discrétisation staggered
  • Doublement de N → erreur divisée par ~4

Dépendances pybind11
---------------------
  import Cosserat   (module Cosserat : CosseratIntrinsicState, PainlessBeamForceField)
  from LieGroups import SO3
"""

import sys
import time
import math
import argparse

# ── Imports pybind11 ──────────────────────────────────────────────────────────
try:
    from Sofa.Cosserat import CosseratIntrinsicState, PainlessBeamForceField
    _HAVE_COSSERAT = True
except ImportError:
    _HAVE_COSSERAT = False
    print("[ERREUR] Module pybind11 'Cosserat' non trouvé.")
    print("  → Compiler avec -DCOSSERAT_WITH_PYTHON_BINDINGS=ON")
    print("  → Lancer depuis le répertoire de build SOFA")

try:
    from Sofa.LieGroups import SO3
    _HAVE_SO3 = True
except ImportError:
    _HAVE_SO3 = False
    print("[ERREUR] Module 'LieGroups' non trouvé.")

try:
    import numpy as np
    _HAVE_NP = True
except ImportError:
    _HAVE_NP = False
    print("[ERREUR] numpy non trouvé : pip install numpy")

# ── Paramètres physiques ──────────────────────────────────────────────────────

BEAM_LENGTH = 1.0       # m
RADIUS      = 0.01      # m
YOUNG_MOD   = 1.0e6     # Pa  (régime petits déplacements)
POISSON     = 0.49
SHEAR_MOD   = YOUNG_MOD / (2.0 * (1.0 + POISSON))
DENSITY     = 1200.0    # kg/m³
GRAVITY     = 9.81      # m/s²  (sens −Y dans l'intégration)

I_SEC  = math.pi * RADIUS**4 / 4.0
A_SEC  = math.pi * RADIUS**2
Q_LINE = DENSITY * A_SEC * GRAVITY          # N/m
DELTA_EB = Q_LINE * BEAM_LENGTH**4 / (8.0 * YOUNG_MOD * I_SEC)

# Paramètres d'intégration
DT          = 5e-5      # pas de temps (Euler explicite, stable si petit)
DAMPING_POS = 0.98      # amortissement fort → convergence rapide
DAMPING_ANG = 0.98
TOL         = 1e-4      # critère convergence |Δδ|/δ
MAX_STEPS   = 200_000   # limite de sécurité

# Niveaux de maillage à tester
N_VALUES = [4, 8, 16, 32, 64]


# ── Calcul des rigidités ──────────────────────────────────────────────────────

def compute_stiffness(E, G, r):
    """Rigidités canoniques pour section circulaire pleine."""
    A   = math.pi * r**2
    Iy  = math.pi * r**4 / 4.0
    Iz  = Iy
    J   = 2.0 * Iy   # J polaire ≈ 2 I pour cercle
    return {
        'EA':  E * A,
        'GA':  G * A,
        'GJ':  G * J,
        'EIy': E * Iy,
        'EIz': E * Iz,
    }


# ── Intégrateur staggered (sans SOFA GUI) ────────────────────────────────────

def run_staggered_beam(N):
    """
    Simule un cantilever staggered (N segments) sous gravité jusqu'à l'équilibre.

    Returns
    -------
    dict avec :
      - delta_tip   : déflexion de pointe |y_N|  [m]
      - error_pct   : erreur relative vs δ_EB     [%]
      - n_steps     : pas avant convergence
      - cpu_s       : temps CPU total              [s]
      - cpu_per_step: temps CPU par pas            [μs]
      - converged   : bool
    """
    if not (_HAVE_COSSERAT and _HAVE_SO3 and _HAVE_NP):
        return None

    h = BEAM_LENGTH / N
    Np1 = N + 1

    # ── Masse / inertie ────────────────────────────────────────────────────
    m_total  = DENSITY * A_SEC * BEAM_LENGTH
    m_node   = m_total / Np1
    L_seg    = BEAM_LENGTH / N
    m_seg    = DENSITY * A_SEC * L_seg
    I_seg    = m_seg * (L_seg**2 / 12.0 + RADIUS**2 / 4.0)

    # ── Initialisation pybind11 ────────────────────────────────────────────
    state = CosseratIntrinsicState()
    ff    = PainlessBeamForceField()

    # Positions initiales (droite, le long de x)
    rest_pos = np.array([[i * h, 0., 0.] for i in range(Np1)])
    state.setPositions(rest_pos.copy())
    state.setOrientations([SO3() for _ in range(N)])

    # Rigidités
    k = compute_stiffness(YOUNG_MOD, SHEAR_MOD, RADIUS)
    ff.setEA(k['EA'])
    ff.setGA(k['GA'])
    ff.setGJ(k['GJ'])
    ff.setEIy(k['EIy'])
    ff.setEIz(k['EIz'])
    ff.setState(state)
    ff.init()

    # ── Buffers de vitesse ─────────────────────────────────────────────────
    vel_pos = np.zeros((Np1, 3))
    vel_ang = np.zeros((N,   3))

    prev_delta = 0.0
    converged  = False
    n_steps    = 0

    t0 = time.perf_counter()

    for step in range(MAX_STEPS):
        # ── Forces élastiques ──────────────────────────────────────────────
        ff.computeForces()
        f_el   = ff.getNodalForces()    # (Np1, 3)
        tau_el = ff.getSegmentTorques() # (N,   3)

        if f_el.shape[0] != Np1:
            break

        # ── Intégration positions (nodes 1..N, node 0 encastré) ───────────
        pos    = state.getPositions()   # (Np1, 3)
        new_pos = pos.copy()

        for i in range(1, Np1):
            f_total     = f_el[i].copy()
            f_total[1] -= m_node * GRAVITY    # gravité en −Y
            acc          = f_total / m_node
            vel_pos[i]   = vel_pos[i] * DAMPING_POS + acc * DT
            new_pos[i]   = pos[i] + vel_pos[i] * DT

        state.setPositions(new_pos)

        # ── Intégration SO3 (segments 1..N-1, segment 0 encastré) ─────────
        R_list  = state.getOrientations()
        new_R   = list(R_list)
        for j in range(1, N):
            alpha      = tau_el[j] / I_seg
            vel_ang[j] = vel_ang[j] * DAMPING_ANG + alpha * DT
            new_R[j]   = R_list[j] * SO3.exp(vel_ang[j] * DT)
        state.setOrientations(new_R)

        # ── Critère de convergence ─────────────────────────────────────────
        n_steps += 1
        tip_y   = float(new_pos[-1][1])
        delta   = abs(tip_y)
        if step > 500 and delta > 1e-9:
            rel_change = abs(delta - prev_delta) / delta
            if rel_change < TOL:
                converged = True
                break
        prev_delta = delta

    cpu_s = time.perf_counter() - t0

    tip_y = float(state.getPositions()[-1][1])
    delta = abs(tip_y)
    error_pct = abs(delta - DELTA_EB) / DELTA_EB * 100.0 if DELTA_EB > 0 else float('nan')

    return {
        'N':            N,
        'h':            h,
        'delta_tip':    delta,
        'delta_EB':     DELTA_EB,
        'error_pct':    error_pct,
        'n_steps':      n_steps,
        'cpu_s':        cpu_s,
        'cpu_per_step': cpu_s / n_steps * 1e6 if n_steps > 0 else 0.0,
        'converged':    converged,
    }


# ── Affichage ─────────────────────────────────────────────────────────────────

def print_header():
    print("\n" + "=" * 80)
    print("  BENCHMARK — Staggered Cosserat : convergence vs. Euler-Bernoulli")
    print("=" * 80)
    print(f"  L = {BEAM_LENGTH} m   r = {RADIUS} m   E = {YOUNG_MOD:.1e} Pa   ρ = {DENSITY} kg/m³")
    print(f"  q = {Q_LINE:.4f} N/m   δ_EB = {DELTA_EB * 1000:.4f} mm")
    print(f"  dt = {DT:.1e} s   damping = {DAMPING_POS}   tol = {TOL}")
    print("-" * 80)
    print(f"  {'N':>5}  {'h [mm]':>8}  {'δ_num [mm]':>11}  {'δ_EB [mm]':>10}  "
          f"{'err [%]':>8}  {'O(h^p)':>7}  {'steps':>7}  {'t_CPU [s]':>10}  {'μs/step':>8}")
    print("-" * 80)


def print_row(res, prev_res=None):
    if res is None:
        print("  [ERREUR — bindings manquants]")
        return

    # Estimation de l'ordre de convergence
    order_str = "   —"
    if prev_res and prev_res['error_pct'] > 0 and res['error_pct'] > 0:
        ratio = prev_res['error_pct'] / res['error_pct']
        h_ratio = prev_res['h'] / res['h']
        if ratio > 0 and h_ratio > 0:
            p = math.log(ratio) / math.log(h_ratio)
            order_str = f"O(h^{p:.2f})"

    print(f"  {res['N']:>5}  {res['h']*1000:>8.3f}  {res['delta_tip']*1000:>11.5f}  "
          f"{res['delta_EB']*1000:>10.5f}  {res['error_pct']:>8.3f}  "
          f"{order_str:>7}  {res['n_steps']:>7}  {res['cpu_s']:>10.3f}  "
          f"{res['cpu_per_step']:>8.2f}")


def print_footer(results):
    print("-" * 80)
    if len(results) >= 2:
        e1 = results[-2]['error_pct']
        e2 = results[-1]['error_pct']
        h1 = results[-2]['h']
        h2 = results[-1]['h']
        if e1 > 0 and e2 > 0:
            p = math.log(e1 / e2) / math.log(h1 / h2)
            print(f"  Ordre de convergence global : p ≈ {p:.2f}  (attendu : 2.0 pour staggered)")
    print("=" * 80 + "\n")


# ── Export CSV ────────────────────────────────────────────────────────────────

def save_csv(results, path):
    import csv
    keys = ['N', 'h', 'delta_tip', 'delta_EB', 'error_pct',
            'n_steps', 'cpu_s', 'cpu_per_step', 'converged']
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        for r in results:
            if r:
                writer.writerow({k: r[k] for k in keys})
    print(f"  → Résultats sauvegardés dans {path}")


# ── Graphiques matplotlib ─────────────────────────────────────────────────────

def plot_results(results):
    try:
        import matplotlib.pyplot as plt
        import matplotlib.ticker as ticker
    except ImportError:
        print("  [INFO] matplotlib non disponible — graphiques désactivés.")
        return

    valid = [r for r in results if r and r['converged']]
    if not valid:
        return

    N_vals   = [r['N'] for r in valid]
    h_vals   = [r['h'] for r in valid]
    err_vals = [r['error_pct'] for r in valid]
    cpu_vals = [r['cpu_per_step'] for r in valid]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle("Staggered Cosserat — Convergence vs Euler-Bernoulli", fontsize=13)

    # ── Graphe 1 : Erreur relative vs h ──────────────────────────────────────
    ax = axes[0]
    ax.loglog(h_vals, err_vals, 'o-', color='steelblue', linewidth=2, label='Staggered')

    # Reference slope O(h²)
    if len(h_vals) >= 2:
        h_ref = np.array([min(h_vals), max(h_vals)])
        e_ref = err_vals[-1] * (h_ref / h_vals[-1])**2
        ax.loglog(h_ref, e_ref, '--', color='gray', alpha=0.7, label=r'$O(h^2)$ ref')

    ax.set_xlabel("h (longueur de segment) [m]", fontsize=11)
    ax.set_ylabel("Erreur relative |δ_num - δ_EB| / δ_EB  [%]", fontsize=11)
    ax.set_title("Convergence en maillage", fontsize=11)
    ax.legend()
    ax.grid(True, which='both', alpha=0.3)
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())

    # ── Graphe 2 : Temps CPU par pas vs N ────────────────────────────────────
    ax = axes[1]
    ax.semilogy(N_vals, cpu_vals, 's-', color='darkorange', linewidth=2)
    ax.set_xlabel("N (nombre de segments)", fontsize=11)
    ax.set_ylabel("Temps CPU par pas [μs]", fontsize=11)
    ax.set_title("Coût de calcul (Euler explicite + SO3)", fontsize=11)
    ax.grid(True, which='both', alpha=0.3)
    ax.xaxis.set_major_locator(ticker.FixedLocator(N_vals))

    plt.tight_layout()
    plt.savefig("benchmark_staggered.png", dpi=150, bbox_inches='tight')
    print("  → Graphiques sauvegardés dans benchmark_staggered.png")
    plt.show()


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Benchmark convergence discrétisation Staggered Cosserat")
    parser.add_argument('--plot', action='store_true',
                        help="Afficher les graphiques matplotlib")
    parser.add_argument('--csv', type=str, default='',
                        help="Chemin de sortie CSV (ex: results.csv)")
    parser.add_argument('--N', type=str, default='',
                        help="Valeurs de N séparées par des virgules (ex: 4,8,16)")
    args = parser.parse_args()

    # Niveaux de maillage
    n_vals = N_VALUES
    if args.N:
        try:
            n_vals = [int(x) for x in args.N.split(',')]
        except ValueError:
            print(f"[WARNING] --N invalide : '{args.N}', utilisation de N_VALUES par défaut")

    if not (_HAVE_COSSERAT and _HAVE_SO3 and _HAVE_NP):
        print("\n[ERREUR] Dépendances manquantes — abandon.")
        sys.exit(1)

    print_header()
    results  = []
    prev_res = None

    for N in n_vals:
        sys.stdout.write(f"  N={N:3d}  en cours... ")
        sys.stdout.flush()
        res = run_staggered_beam(N)
        print_row(res, prev_res)
        results.append(res)
        prev_res = res

    print_footer(results)

    if args.csv:
        save_csv(results, args.csv)

    if args.plot:
        plot_results(results)


if __name__ == '__main__':
    main()
