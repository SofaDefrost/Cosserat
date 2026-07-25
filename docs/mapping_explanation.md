# Cinématique du mapping Cosserat : `apply`, `applyJ`, `applyJT` (forces et contraintes)

> **Fichiers de référence**
> - Mapping géométrique : `src/liegroups/CosseratBodyJacobian.h`
> - Mapping discret SOFA : `src/Cosserat/mapping/DiscreteCosseratMapping.inl`
> - Mapping dynamique : `src/Cosserat/mapping/DiscreteDynamicCosseratMapping.inl`
> - Groupe de Lie : `src/liegroups/SE3.h`

---

## 1. Modélisation géométrique — rappels

### 1.1 Convention de déformation (Cosserat)

Un bras Cosserat discrétisé en **N éléments** est décrit par :

| Symbole | Espace | Signification |
|---------|--------|---------------|
| $\xi_k \in \mathbb{R}^6$ | section $k$ | déformation locale : $\xi_k = [\varphi_k,\, \rho_k]^T$ |
| $\varphi_k = [\varphi_x,\varphi_y,\varphi_z]$ | | torsion ($x$) + flexion ($y,z$) |
| $\rho_k = [\rho_x,\rho_y,\rho_z]$ | | élongation ($x$) + cisaillement ($y,z$) |
| $L_k$ | $\mathbb{R}_{>0}$ | longueur au repos de l'élément $k$ |
| $g_k^{\text{loc}} \in SE(3)$ | section $k$ | pose locale : $g_k^{\text{loc}} = \exp_{\text{Cos}}(\xi_k, L_k)$ |
| $G_k \in SE(3)$ | nœud $k$ | pose globale cumulée depuis la base |

La **convention tête/queue** dans tout le code est :

$$\xi = \underbrace{[\varphi_x,\varphi_y,\varphi_z]}_{\text{head}\langle 3\rangle\ \text{(angulaire)}}\ \big|\ \underbrace{[\rho_x,\rho_y,\rho_z]}_{\text{tail}\langle 3\rangle\ \text{(linéaire)}}$$

### 1.2 Exponentielles utilisées

| Fonction | Formule | Usage |
|----------|---------|-------|
| `expCosserat(ξ, L)` | $\exp\!\bigl(\hat\xi_{\text{Cos}}\bigr)$, ajoute $+1$ à $\rho_x$ | mapping géométrique (`apply`) |
| `computeExp(ξ)` | $\exp(\hat\xi)$ standard $SE(3)$ | intégrateurs, log/exp purs |

> `expCosserat` encode l'**élongation nominale** $+1$ selon l'axe $x$ via `buildXiHat`. Ne pas mélanger les deux.

### 1.3 Matrice tangente-exponentielle $J_{local}$

Pour l'élément $k$, la matrice $J_{local,k} \in \mathbb{R}^{6\times 6}$ est le résultat de `computeTangExpImplementation(L_k, \xi_k)`. Elle relie incrémentalement la variation de déformation à la variation de pose :

$$\delta g_k^{\text{loc}} = g_k^{\text{loc}} \cdot \hat{J_{local,k} \cdot \delta\xi_k}$$

---

## 2. `apply` — mapping géométrique

### 2.1 Objet

Mappe les déformations $(\xi_1,\ldots,\xi_N)$ vers les **poses globales** $(G_0, G_1,\ldots, G_N)$.

### 2.2 Récurrence directe (racine → bout)

$$\boxed{G_0 = G_{\text{base}} \quad \text{(donné)}}$$

$$\boxed{G_k = G_{k-1} \cdot g_k^{\text{loc}}, \qquad g_k^{\text{loc}} = \exp_{\text{Cos}}(\xi_k,\, L_k), \quad k = 1,\ldots,N}$$

### 2.3 Implémentation (pseudo-code)

```cpp
G[0] = G_base;
for k = 1 to N:
    g_loc[k] = SE3::expCosserat(xi[k], L[k]);
    G[k]     = G[k-1].compose(g_loc[k]);
```

### 2.4 Interprétation géométrique

```
base ──[g_1^loc]──○──[g_2^loc]──○── ··· ──[g_N^loc]──○ bout
G_0              G_1            G_2                    G_N
```

$G_k$ est la **pose absolue** du nœud $k$ dans le repère monde. La pose $G_k^{-1} G_j$ (pour $j > k$) donne la pose du nœud $j$ vue depuis le nœud $k$.

---

## 3. `applyJ` — propagation des vitesses

### 3.1 Objet

Mappe les **taux de déformation** $\dot\xi_k$ vers les **torseurs de vitesse corps** $\eta_k$ à chaque nœud. C'est le produit Jacobien $\eta = J \cdot \dot\xi$.

### 3.2 Grandeurs

| Symbole | Espace | Signification |
|---------|--------|---------------|
| $\dot\xi_k \in \mathbb{R}^6$ | section $k$ | vitesse de déformation |
| $\eta_k \in se(3) \cong \mathbb{R}^6$ | nœud $k$ | torseur de vitesse corps (twist) |
| $Ad_{g^{-1}} \in \mathbb{R}^{6\times 6}$ | — | matrice adjointe inverse de $g$ |

### 3.3 Récurrence directe (racine → bout)

$$\boxed{\eta_0 = \bar\eta_0 \quad \text{(torseur imposé à la base, souvent 0)}}$$

$$\boxed{\eta_k = Ad_{g_k^{-1}} \cdot \bigl(\eta_{k-1} + J_{local,k} \cdot \dot\xi_k\bigr), \qquad k = 1,\ldots,N}$$

### 3.4 Développement explicite

En déroulant la récurrence :

$$\eta_k = \sum_{j=1}^{k} \underbrace{Ad_{G_k^{-1} G_{j-1}}}_{B_{j,k}} \cdot J_{local,j} \cdot \dot\xi_j \;+\; Ad_{G_k^{-1}} \cdot \bar\eta_0$$

avec $G_k = G_0 \cdot g_1^{\text{loc}} \cdots g_k^{\text{loc}}$ et $B_{j,k} = Ad_{g_k^{-1}} \cdots Ad_{g_{j+1}^{-1}}$.

Le Jacobien corps complet à la section $k$ est donc :

$$J_{\text{body}}(k) = \bigl[B_{1,k}\,J_{local,1} \mid B_{2,k}\,J_{local,2} \mid \cdots \mid J_{local,k} \mid 0 \mid \cdots \mid 0\bigr] \in \mathbb{R}^{6 \times 6N}$$

(matrice **triangulaire inférieure par blocs**).

### 3.5 Implémentation (pseudo-code)

```cpp
twists[0] = base_twist;
for k = 0 to N-1:               // section k (0-indexé)
    transported = Ad_inv[k] * (twists[k] + J_local[k] * strain_rates[k]);
    twists[k+1] = transported;
```

### 3.6 Interprétation

Chaque section $k$ :

1. **Transporte** le twist entrant $\eta_{k-1}$ dans son propre repère via $Ad_{g_k^{-1}}$.
2. **Ajoute** la contribution propre de sa déformation $J_{local,k}\cdot\dot\xi_k$.

L'opération $Ad_{g^{-1}}$ est un **changement de repère** : elle convertit le twist exprimé au nœud $k-1$ en un twist exprimé au nœud $k$, dans le repère corps de l'élément $k$.

---

## 4. `applyJT` — rétro-propagation des forces

### 4.1 Objet

Mappe les **torseurs extérieurs** $w_k$ appliqués aux nœuds vers les **forces de déformation** $q_k$ (forces généralisées en espace des déformations). C'est l'opération transposée $q = J^T \cdot w$.

### 4.2 Grandeurs

| Symbole | Espace | Signification |
|---------|--------|---------------|
| $w_k \in se(3)^* \cong \mathbb{R}^6$ | nœud $k$ | torseur extérieur (wrench) : $w_k = [\tau_k, f_k]^T$ |
| $q_k \in \mathbb{R}^6$ | section $k$ | force généralisée de déformation |
| $f_k \in se(3)^*$ | nœud $k$ | **co-torseur accumulé** (grandeur intermédiaire) |
| $b_0 \in se(3)^*$ | nœud $0$ | torseur de réaction à la base |

### 4.3 Principe de puissance virtuelle

Le principe doit être satisfait pour tout $\delta\dot\xi_k$ admissible :

$$\underbrace{\sum_{k=0}^{N} w_k^T \cdot \delta\eta_k}_{\text{puissance externe}} = \underbrace{\sum_{k=1}^{N} q_k^T \cdot \delta\dot\xi_k}_{\text{puissance interne des déformations}} + \underbrace{b_0^T \cdot \delta\eta_0}_{\text{réaction à la base}}$$

### 4.4 Dérivation

En substituant la récurrence de `applyJ` dans chaque terme de puissance externe :

$$w_k^T \cdot \delta\eta_k = w_k^T \cdot Ad_{g_k^{-1}}\bigl(\delta\eta_{k-1} + J_{local,k}\,\delta\dot\xi_k\bigr)$$

$$= \underbrace{\bigl(Ad_{g_k^{-T}} w_k\bigr)^T \cdot \delta\eta_{k-1}}_{\text{transport du torseur vers } k-1} + \underbrace{\bigl(J_{local,k}^T\, Ad_{g_k^{-T}} w_k\bigr)^T \cdot \delta\dot\xi_k}_{\text{contribution à } q_k}$$

En accumulant ces transports de $k = N$ jusqu'à $k = 1$, on définit le **co-torseur accumulé** $f_k$ :

$$f_k = \sum_{j=k}^{N} Ad_{G_j^{-1} G_{k-1}}^T \cdot w_j$$

dont la récurrence locale est :

$$f_N = w_N \qquad\qquad \text{(initialisation au bout)}$$

$$f_{k-1} = w_{k-1} + Ad_{g_k^{-T}} \cdot f_k, \qquad k = N, N-1, \ldots, 1$$

La force de déformation à la section $k$ :

$$\boxed{q_k = J_{local,k}^T \cdot Ad_{g_k^{-T}} \cdot f_k}$$

La réaction à la base :

$$b_0 = f_0 = w_0 + Ad_{g_1^{-T}} \cdot f_1$$

### 4.5 Point crucial — causalité des efforts

> $f_k$ est la somme des torseurs extérieurs **avals** (positions $k$ à $N$) transportés en position $k$.
> Le torseur $w_{k-1}$, appliqué **en amont** de la section $k$, n'entre **jamais** dans $q_k$.
> Il entre seulement dans $q_{k-1}, q_{k-2}, \ldots$ via la rétro-propagation.

Physiquement : la déformation $\dot\xi_k$ de la section $k$ crée un mouvement au **bout** de cette section (nœud $k$), donc elle ne fait du travail qu'avec les efforts appliqués **en aval** (nœuds $k, k+1, \ldots, N$).

### 4.6 Implémentation (pseudo-code)

```cpp
acc = w[N];                              // f_N = w_N
for k = N-1 downto 0:                   // section k (0-indexé)
    transported  = Ad_inv[k].T * acc;   // Ad_{g_k^{-T}} · f_{k+1}
    q[k]         = J_local[k].T * transported;   // force de déformation
    acc          = w[k] + transported;  // f_k = w_k + Ad_{g_k^{-T}} · f_{k+1}

base_wrench_out = acc;                  // b_0 = f_0
```

> **Erreur classique** : calculer `q[k] = J.T * acc` **après** avoir fait `acc = w[k] + Ad.T * acc`.
> Cela pollue $q_k$ avec $J_{local,k}^T w_k$ qui n'y a pas sa place.
> La bonne séquence est : transporter → calculer $q_k$ → ajouter $w_k$.

---

---

## 5. `applyJT` — version contraintes (`MatrixDeriv`)

### 5.1 Pourquoi une deuxième signature ?

SOFA dispose de deux paires `applyJ` / `applyJT` :

| Signature | Rôle | Contexte SOFA |
|-----------|------|---------------|
| `applyJ(VecDeriv)` / `applyJT(VecDeriv)` | propagation **dense** de vitesses / forces nodales | dynamique, forces externes |
| `applyJ(MatrixDeriv)` / `applyJT(MatrixDeriv)` | propagation **creuse** de directions de contraintes | solveur de contraintes (LCP, contact) |

La version `MatrixDeriv` opère sur la **matrice Jacobienne de contraintes** $H$, dont chaque ligne est une direction de contrainte scalaire $c$ portée par plusieurs repères de sortie.

```
H_out ∈ ℝ^{n_c × n_frames}    →    applyJT(MatrixDeriv)    →    H_in ∈ ℝ^{n_c × n_strains}
(directions de contraintes                                         (contraintes ramenées
 sur les repères mappés)                                            dans l'espace des déformations)
```

La relation est : $H_{\text{in}} = H_{\text{out}} \cdot J_{\text{body}}$, soit ligne par ligne :

$$\boxed{h_{c,\text{in}} = J_{\text{body}}^T \cdot h_{c,\text{out}}}$$

### 5.2 Grandeurs supplémentaires

| Symbole | Espace | Signification |
|---------|--------|---------------|
| $\delta_c^s \in \mathbb{R}^6$ | repère $s$ | direction de contrainte $c$ au repère $s$ (dans le repère monde) |
| $P_s \in \mathbb{R}^{6\times 6}$ | repère $s$ | **projecteur** de la convention Rigid6 SOFA → SE(3) |
| $B^T \in \mathbb{R}^{3\times 6}$ | espace des déformations | **sélecteur** des DOF de déformation actifs ($= [I_3 \mid 0_3]$) |
| $b_s \in \mathbb{N}$ | — | indice de l'élément contenant le repère $s$ |
| $\ell_c^s \in \mathbb{R}^6$ | corps | direction de contrainte dans le repère local de l'élément $b_s$ |

### 5.3 Structure à deux niveaux

Il y a deux niveaux de sous-échantillonnage dans le mapping Cosserat :

```
Nœuds de déformation (strain)  ──[applyJ]──►  Nœuds de la poutre  ──[apply frames]──►  Repères de sortie
         ξ_k                                        g_k^{loc}                               G_s, s ∈ b_k
```

Un seul élément poutre $k$ peut contenir **plusieurs** repères de sortie $s \in b_k$ (interpolation intra-élément). La matrice `m_framesExponentialSE3Vectors[s]` porte le transport de l'élément vers le repère $s$, tandis que `m_nodesExponentialSE3Vectors[k]` porte le transport de la base vers le nœud $k$.

### 5.4 Algorithme (par ligne de contrainte $c$)

#### Phase 1 — Contribution directe de chaque repère

Pour chaque repère $s$ ayant une entrée non nulle dans la ligne $c$ :

**Étape 1 — Direction de contrainte → repère local de l'élément**

$$\ell_c^s = \underbrace{CoAdj\!\left(g_s^{-1}\right)}_{\texttt{computeCoAdjoint}\atop\texttt{(framesExp[s])}} \cdot \underbrace{P_s^T}_{\texttt{buildProjector}(s)^T} \cdot \delta_c^s$$

$P_s^T$ gère la convention SOFA (translation avant rotation dans `Rigid3`) ↔ SE(3) (angulaire avant linéaire). $CoAdj(g_s^{-1}) = Ad_{g_s^{-T}}$ est le co-adjoint inverse qui ramène la force du monde vers le corps.

**Étape 2 — Projection dans l'espace des déformations de l'élément $b_s$**

$$q_{c,\, b_s} \mathrel{+}= B^T \cdot \underbrace{J_{local,s}^T}_{\texttt{framesTangExp[s]}^T} \cdot \ell_c^s \quad \in \mathbb{R}^3$$

$B^T = [I_3 \mid 0_3]$ sélectionne les trois premières composantes (angulaires dans la convention $[\varphi,\rho]$), qui correspondent aux DOF de déformation actifs (flexion, torsion, élongation réduite).

**Étape 3 — Stockage pour la rétro-propagation**

La direction locale $\ell_c^s$ est accumulée au nœud $b_s$ pour le passage arrière :

$$\tilde{f}_{b_s} \mathrel{+}= \ell_c^s$$

#### Phase 2 — Rétro-propagation co-adjointe vers la racine

Pour $k = k_{\max}, k_{\max}-1, \ldots, 1$ :

$$\tilde{f} \leftarrow \underbrace{CoAdj\!\left(g_k^{-1}\right)}_{\texttt{computeCoAdjoint}\atop\texttt{(nodesExp[k-1])}} \cdot \tilde{f}$$

$$q_{c,\, k-1} \mathrel{+}= B^T \cdot \underbrace{T_k^T}_{\texttt{nodesTangExp[k-1]}^T} \cdot \tilde{f} \quad \in \mathbb{R}^3$$

À chaque pas, le co-adjoint transporte le torseur de force du nœud $k$ vers le nœud $k-1$, et $T_k^T$ le projette dans l'espace des déformations de l'élément $k-1$.

#### Phase 3 — Réaction à la base

$$b_c = P_0 \cdot \tilde{f} \quad \in \mathbb{R}^6$$

$P_0$ = `buildProjector(frame[0])` reconvertit dans la convention Rigid6 SOFA pour le repère de base rigide.

### 5.5 Schéma complet pour une contrainte $c$ à un repère $s \in b_s$

```
δ_c^s (monde)
    │
    │ P_s^T          (Rigid6 SOFA → SE(3))
    ▼
    │ CoAdj(g_s^{-1}) (monde → corps de l'élément b_s)
    ▼
ℓ_c^s ──────────────────────────── B^T · J_local,s^T ──► q_{c, b_s}  (strain space, élément b_s)
    │
    │ CoAdj(g_{b_s}^{-1})          (nœud b_s → nœud b_s-1)
    ▼
    │ B^T · T_{b_s-1}^T  ──────────────────────────────► q_{c, b_s-1}
    │
    │ CoAdj(g_{b_s-1}^{-1})
    ▼
    │ B^T · T_{b_s-2}^T  ──────────────────────────────► q_{c, b_s-2}
    │
    ·
    ·
    │ CoAdj(g_1^{-1})
    ▼
    │ P_0  ─────────────────────────────────────────────► b_c  (base frame, In2)
```

### 5.6 Différence entre `DiscreteCosseratMapping` et `DiscreteDynamicCosseratMapping`

Les deux implémentent le même calcul mathématique, mais diffèrent dans leur stratégie d'agrégation quand plusieurs repères appartiennent au même élément :

| Aspect | `DiscreteCosseratMapping` | `DiscreteDynamicCosseratMapping` |
|--------|--------------------------|----------------------------------|
| **Agrégation** | `vector<tuple<int,Vec6>>` trié par décroissant puis fusionné | tableau `accumulatedCartesianForces[maxBeamIndex+1]`, accumulation directe $O(1)$ |
| **Complexité** | $O(n_s \log n_s)$ (tri) | $O(n_s + N)$ (accumulateur + passage linéaire) |
| **Phases** | Entrelacées dans la boucle `while (i > 0)` | Séparées : phase 1 (accumulation) puis phase 2 (propagation) |
| **Lisibilité** | Plus difficile à suivre (boucle `while` imbriquée) | Plus claire (deux boucles indépendantes) |

> La version `Dynamic` corrige aussi un bug de compression de l'ancienne version (boucle `while (numNode_i == numNode_i1)` fautive dans certains cas limites).

### 5.7 Rôle de `buildProjector` et de $B^T$

**`buildProjector(frame)`** construit la matrice $P \in \mathbb{R}^{6\times6}$ qui change de convention entre :
- SOFA `Rigid3` : $[t_x, t_y, t_z, \omega_x, \omega_y, \omega_z]$ (translation avant rotation)
- SE(3) Cosserat : $[\varphi_x, \varphi_y, \varphi_z, \rho_x, \rho_y, \rho_z]$ (angulaire avant linéaire)

Elle incorpore aussi la rotation du repère pour exprimer la contrainte dans le bon référentiel.

**$B^T = [I_3 \mid 0_3]$** sélectionne uniquement les 3 composantes angulaires du torseur 6D. Cela reflète que chaque élément Cosserat a **3 DOF actifs** dans l'espace des déformations (flexion $y$, flexion $z$, torsion $x$), les déformations linéaires étant soit contraintes, soit représentées par des paramètres séparés selon le type de mapping.

---

## 6. Tableau de dualité complet — les quatre passes

|  | **`apply`** | **`applyJ`** | **`applyJT` (VecDeriv)** | **`applyJT` (MatrixDeriv)** |
|--|-------------|--------------|--------------------------|------------------------------|
| **Sens** | $k=1\to N$ | $k=1\to N$ | $k=N\to 1$ | Phase 1 : $\forall s$ ; Phase 2 : $k=N\to 1$ |
| **Entrée** | $\xi_k$ | $\dot\xi_k$ | $w_k$ (wrenches) | $\delta_c^s$ (directions de contraintes) |
| **Espace entrée** | $\mathbb{R}^6$ / section | $\mathbb{R}^6$ / section | $se(3)^*$ / nœud | $\mathbb{R}^6$ / repère $s$ |
| **Init** | $G_0 = G_{\text{base}}$ | $\eta_0 = \bar\eta_0$ | $f_N = w_N$ | $\tilde f = 0$, accumulation |
| **Récurrence** | $G_k = G_{k-1} \cdot g_k^{\text{loc}}$ | $\eta_k = Ad_{g_k^{-1}}(\eta_{k-1} + J_k\dot\xi_k)$ | $f_{k-1} = w_{k-1} + Ad_{g_k^{-T}} f_k$ | $\tilde f \leftarrow CoAdj(g_k^{-1})\cdot\tilde f$ |
| **Sortie déformation** | — | — | $q_k = J_k^T Ad_{g_k^{-T}} f_k$ | $q_{c,k} = B^T T_k^T \tilde f$ |
| **Sortie directe** | $G_k \in SE(3)$ | $\eta_k \in se(3)$ | — | $q_{c,b_s} \mathrel{+}= B^T J_s^T \ell_c^s$ |
| **Sortie base** | — | — | $b_0 = f_0$ | $b_c = P_0 \cdot \tilde f$ |
| **Couche supplémentaire** | — | — | — | $P_s^T$ (Rigid6↔SE3) et $B^T$ (sélecteur 3/6) |
| **Complexité** | $O(N)$ | $O(N)$ | $O(N)$ | $O(n_c \cdot (n_s + N))$ |

---

## 8. Relations entre les quantités

### 6.1 Lien `apply` → `applyJ`

La matrice $g_k^{\text{loc}} = \exp_{\text{Cos}}(\xi_k, L_k)$ issue de `apply` est **réutilisée** par `applyJ` via son adjoint inverse $Ad_{g_k^{-1}}$, et par `applyJT` via son co-adjoint $Ad_{g_k^{-T}}$.

$$Ad_{g_k^{-T}} = \left(Ad_{g_k^{-1}}\right)^T$$

### 6.2 Propriété de dualité (vérification numérique)

Pour tout $\dot\xi$ et tout $w$, les trois passes doivent satisfaire :

$$\sum_{k=1}^N w_k^T \cdot (J\dot\xi)_k = \sum_{k=1}^N (J^T w)_k^T \cdot \dot\xi_k$$

soit, en notation compacte : $\langle w,\, J\dot\xi \rangle = \langle J^T w,\, \dot\xi \rangle$.

C'est le test `VirtualPowerDualityForwardTranspose` dans `test_BodyJacobianUncertaintyBezier.cpp`.

### 6.3 Puissance virtuelle élémentaire

Pour un seul élément $k$ isolé, la puissance de $w_k$ contre $\eta_k$ s'écrit :

$$P_k = w_k^T \cdot \eta_k = \tau_k^T \varphi_k + f_k^T \rho_k$$

avec $\eta_k = [\varphi_k, \rho_k]^T$ (twist corps) et $w_k = [\tau_k, f_k]^T$ (wrench).
Ce produit est **invariant** par changement de repère via le transport co-adjoint :

$$w'^T \eta' = (Ad^{-T} w)^T (Ad^{-1} \eta) = w^T (Ad^T Ad^{-1}) \eta = w^T \eta$$

---

## 9. Schéma de flux de données

```
Déformations ξ_k
       │
       ▼  apply
Poses locales g_k^loc ──────────────────────────────────────► Poses globales G_k = G_{k-1}·g_k^loc
       │                                                                │
       │  (pré-calcul : Ad_{g_k^{-1}}, J_local_k)                    │ (positions, repères de sortie)
       │                                                                │
       ├──────────────────────────────────────────────────┐           │
       │                                                  │           │
       ▼  applyJ (VecDeriv)                  applyJT ◄───┘           │
 ξ̇_k ──► η_0 → η_1 → ··· → η_N      [w_0,…,w_N] ──► q_k (forces)  │
           propagation avant              rétro-propagation           │
                                              │                        │
                                              ▼                        │
                                         b_0 (base)                    │
                                                                        │
       ┌────────────────────────────────────────────────────────────── ┘
       │  (solveur de contraintes : contact, bilatéral, …)
       │
       ▼  applyJT (MatrixDeriv)   ─── par ligne de contrainte c ───
δ_c^s (monde)
  │ P_s^T              Rigid6 SOFA → SE(3)
  │ CoAdj(g_s^{-1})    monde → corps de l'élément b_s
  ▼
ℓ_c^s ──── B^T · J_local,s^T ──────────────────────────► q_{c, b_s}  (In1 : strain space)
  │
  │ CoAdj(g_{b_s}^{-1}) · … · CoAdj(g_1^{-1})   rétro-propagation co-adjointe
  ▼
  │ B^T · T_k^T  (à chaque nœud k < b_s) ────────────► q_{c, k}     (In1 : strain space)
  │
  │ P_0          (projecteur base)
  ▼
b_c                                                       (In2 : base frame)
```

---

## 10. Notes d'implémentation

### 8.1 Pré-calcul de $Ad_{g_k^{-1}}$

Dans `CosseratBodyJacobian`, `Ad_inv_[k]` est calculé une seule fois lors de `pushSection()` :

```cpp
Ad_inv_.push_back(g_k.computeAdjoint().inverse());
```

Il est réutilisé tel quel par `applyJ` (multiplication directe) et par `applyJT` (transposé).

### 8.2 Relation entre `computeAdjoint()` et `computeCoAdjoint()`

Pour $g = (R, t) \in SE(3)$ :

$$Ad_g = \begin{pmatrix} R & 0 \\ \hat{t}\,R & R \end{pmatrix}, \qquad Ad_g^{-1} = Ad_{g^{-1}} = \begin{pmatrix} R^T & 0 \\ -R^T\hat{t} & R^T \end{pmatrix}$$

$$Ad_g^{-T} = \begin{pmatrix} R & -\hat{t}R \\ 0 & R \end{pmatrix}$$

> `computeCoAdjoint()` retourne $Ad_g^{-T}$ directement (utile pour le transport des wrenches).

### 8.3 Matrice tangente-exponentielle $J_{local}$

$J_{local,k}$ est la matrice $T(L_k, \xi_k)$ issue de `computeTangExpImplementation()`. Elle dépend de la déformation courante et de la longueur de l'élément. Pour une déformation nulle, elle se réduit à $L_k \cdot I_6$.

---

## 11. Références

- **Simo & Vu-Quoc (1988)** — formulation Cosserat originale pour les poutres 3D
- **Murray, Li & Sastry (1994)** — *A Mathematical Introduction to Robotic Manipulation*, ch. 2-3 (adjoint, twists, wrenches)
- **Munthe-Kaas (1998)** — intégrateurs géométriques RKMK sur les groupes de Lie
- Code source : `src/liegroups/CosseratBodyJacobian.h`, `src/liegroups/LieGroupIntegrators.h`

---

## 12. Comparaison des implémentations `applyJT` — co-adjoint complet vs permutation seule

### 12.1 Contexte

Deux approches ont coexisté dans le code pour projeter la force extérieure du repère SOFA vers l'espace des déformations :

| Implémentation | Fichier de référence | Statut |
|----------------|---------------------|--------|
| **Co-adjoint complet** | `DiscreteCosseratMapping.inl` (version active) | ✅ Correcte |
| **Permutation seule** | version commentée / simplifiée | ⚠️ Approximation (voir ci-dessous) |

### 12.2 Version co-adjoint complet (correcte)

```cpp
// Étape 1 : projecteur de convention SOFA Rigid6 → SE(3) + rotation du repère
Mat6x6 P_trans = buildProjector(frame[s]);   // rotation R du repère incluse
P_trans.transpose();

// Étape 2 : co-adjoint ramène le torseur du monde → repère corps de l'élément
Mat6x6 coAdjoint;
computeCoAdjoint(m_framesExponentialSE3Vectors[s], coAdjoint);

// Force locale = Ad_{g_s^{-T}} · P_s^T · w_s
Vec6 local_F = coAdjoint * P_trans * valueConst;

// Projection dans l'espace des déformations
Vec6 f = matB_trans * J_local[s].T * local_F;
```

Le projecteur $P_s$ incorpore :

1. La **permutation** $[\tau,\,f] \leftrightarrow [\varphi,\,\rho]$ (SOFA Rigid3 vs Cosserat).
2. La **rotation** $R_s$ du repère $s$ : les composantes angulaires/linéaires sont exprimées dans le repère local, pas dans le monde.

Le co-adjoint $Ad_{g_s^{-T}}$ (= `computeCoAdjoint`) effectue ensuite le transport de force depuis la position du repère $s$ jusqu'à l'origine de l'élément $b_s$.

### 12.3 Version permutation seule (approximation)

```cpp
// Version simplifiée : seule la permutation des indices est effectuée
Vec6 local_F;
local_F[0] = w[3]; local_F[1] = w[4]; local_F[2] = w[5];  // torque
local_F[3] = w[0]; local_F[4] = w[1]; local_F[5] = w[2];  // force

// Pas de co-adjoint, pas de rotation de repère
Vec6 f = matB_trans * J_local[s].T * local_F;
```

Cette version applique seulement la **permutation des composantes** (ordre SOFA vs ordre Cosserat), sans tenir compte de :

- La **rotation** $R_s$ du repère courant par rapport au repère monde.
- Le **transport** de la force depuis le repère $s$ vers la base de l'élément.

### 12.4 Quand les deux versions coïncident

La version simplifiée est exacte uniquement si :

- Tous les repères sont alignés avec le repère monde : $R_s = I$ pour tout $s$.
- La poutre est droite et le co-adjoint de transport est l'identité.

En pratique, c'est vérifié pour une **poutre droite en configuration de référence** (non déformée, orientée selon $x$).

### 12.5 Erreur pour les poutres courbes

Dès que la poutre se déforme ou qu'un repère pivote, l'erreur de la version simplifiée est :

$$\Delta f_k = \left(Ad_{g_k^{-T}} \cdot P_k^T - \text{permutation}\right) \cdot w_k$$

proportionnelle à $\|R_k - I\|$ (amplitude de la rotation) et à $\|w_k\|$ (norme de la force appliquée).

**Conclusion** : utiliser systématiquement `buildProjector` + `computeCoAdjoint`. La version simplifiée ne doit pas être utilisée pour des poutres déformées ou courbées.

### 12.6 Cas `matB_trans` — Vec3Types vs Vec6Types

La matrice $B^T$ sélectionne les DOF actifs de l'espace des déformations :

| Type de DOF | `matB_trans` | Interprétation |
|------------|-------------|----------------|
| **Vec3Types** | $[I_3 \mid 0_3]$ (3×6) | seules les 3 composantes angulaires $\varphi$ sont actives (flexion + torsion) |
| **Vec6Types** | $I_6$ (6×6) | les 6 composantes $[\varphi, \rho]$ sont actives (angulaire + linéaire) |

> **Bug corrigé** dans `DiscreteCosseratMapping.cpp` (spécialisation Vec6Types) : la boucle était `for(k<3)` au lieu de `for(k<6)`, rendant $B^T = [I_3 \mid 0_3]$ pour une spécialisation Vec6, ce qui annulait la moitié des forces de déformation.

---

## 13. Implémentation `Strain2RigidCosseratMapping` — variante liegroups

### 13.1 Présentation générale

`Strain2RigidCosseratMapping` est une **nouvelle implémentation** du même mapping Cosserat, construite directement sur la bibliothèque `liegroups` (types Eigen natifs). Elle coexiste avec `DiscreteCosseratMapping` et vise une meilleure lisibilité et maintenabilité.

Fichiers :
- `src/Cosserat/mapping/Strain2RigidCosseratMapping.h` — déclaration de classe
- `src/Cosserat/mapping/Strain2RigidCosseratMapping.inl` — implémentation
- `src/Cosserat/mapping/Strain2RigidCosseratMapping.cpp` — instanciations

### 13.2 Structures de données internes

| Structure | Rôle |
|-----------|------|
| `FrameInfo` | Une frame de sortie : transformation locale, adjoint, index de section |
| `SectionInfo` | Une section (élément) : déformation, longueur, transformation $g_k^{\text{loc}}$, adjoint |
| `m_section_properties[i]` | Tableau de `SectionInfo` (taille N+1, index 0 = base) |
| `m_frameProperties[i]` | Tableau de `FrameInfo` (autant de frames de sortie) |
| `m_indices_vectors[i]` | Index de section associé à la frame $i$ |

### 13.3 `apply`

```cpp
// Mise à jour des transformations locales (section par section)
SE3Types _gx = SE3Types::expCosserat(strain, section_length);
m_section_properties[i+1].setTransformation(_gx);

// Composition pour obtenir la pose globale de la frame i
auto current_frame = base_frame;
for j = 0 to related_beam_idx-1:
    current_frame = current_frame * m_section_properties[j].getTransformation();
current_frame = current_frame * m_frameProperties[i].getTransformation();
```

Différence principale par rapport à `DiscreteCosseratMapping` : le résultat est converti en `sofa::Rigid3` via un quaternion Eigen, sans passer par des matrices SOFA internes.

### 13.4 `applyJ`

La récurrence est la même qu'en §3 mais exprimée en types Eigen :

```cpp
// Base : projection SOFA → SE(3) via buildProjectionMatrix
node_velocities[0] = base_projector * base_vel_local;

// Propagation section par section
for i = 1 to N:
    Ad_gim1 = section.inverse().getAdjoint();
    node_velocities[i] = Ad_gim1 * (node_velocities[i-1] + tang_adj * strain_vel_i);

// Propagation intra-section pour chaque frame de sortie
g_inv = frame.getInverseTransformation();
Ad_gm1 = g_inv.computeAdjoint();
eta_frame = Ad_gm1 * (node_velocities[section_idx] + tang_adj * frame_strain_vel);

// Reprojection SE(3) → SOFA
output_vel = frame_projector * eta_frame;
```

`buildProjectionMatrix(R)` est l'équivalent de `buildProjector(frame)` de `DiscreteCosseratMapping`, mais formulé sur les types Eigen natifs.

### 13.5 `applyJT` (VecDeriv)

```cpp
// Phase 1 : projection globale → locale pour chaque frame
AdjointMatrix P_trans = absoluteFrame.buildProjectionMatrix(R);
TangentVector localForce = P_trans.transpose() * frameForce;
localForces.push_back(localForce);

// Phase 2 : rétro-propagation (bout → racine)
for s = sz-1 downto 0:
    coAdjoint = frame.getCoAdjoint();
    currentLocalForce = coAdjoint * localForces[s];
    f = matB_trans * J_local[s].T * currentLocalForce;

    if section_changed:
        coAdjoint = section.getCoAdjoint();
        totalForce = coAdjoint * totalForce;  // transport co-adjoint
        temp_f = matB_trans * T_section.T * totalForce;
        strainForces[section-1] += temp_f;

    totalForce += currentLocalForce;
    strainForces[section-1] += f;

// Base : reprojection co-adjoint → SOFA
M = absoluteFrame0.buildProjectionMatrix(R0);
baseForces[baseIndex] += M * totalForce;
```

### 13.6 `applyJT` (MatrixDeriv)

Même structure à deux phases que `DiscreteCosseratMapping.inl` :

```cpp
// Pour chaque contrainte c et chaque frame s ∈ c :
P_trans = absoluteFrame.buildProjectionMatrix(R);
coAdjoint = frame.getCoAdjoint();
localForce = coAdjoint * P_trans.T * constraintValue;

f = matB_trans * J_local[s].T * localForce;  // → espace des déformations
o1.addCol(sectionIndex-1, f);

NodesInvolved.push_back({sectionIndex, localForce});

// Tri décroissant + compression des nœuds dupliqués
// Puis rétro-propagation co-adjointe :
while i > 0:
    coAdjoint = section.getCoAdjoint();
    CumulativeF = coAdjoint * CumulativeF;
    temp_f = matB_trans * T_section.T * CumulativeF;
    o1.addCol(i-2, temp_f);
    i--;

// Réaction à la base :
M = absoluteFrame0.buildProjectionMatrix(R0);
o2.addCol(baseIndex, M * CumulativeF);
```

### 13.7 Comparaison avec `DiscreteCosseratMapping`

| Aspect | `DiscreteCosseratMapping` | `Strain2RigidCosseratMapping` |
|--------|--------------------------|-------------------------------|
| **Types** | Types SOFA (`Vec6`, `Mat6x6`, `Frame`) | Types Eigen natifs (`TangentVector`, `AdjointMatrix`, `SE3Types`) |
| **Projecteur** | `buildProjector(frame)` → Mat6x6 SOFA | `buildProjectionMatrix(R)` → Eigen |
| **Co-adjoint** | `computeCoAdjoint(expVec, result)` — procédure | `frame.getCoAdjoint()` — accesseur objet |
| **Helpers** | Tableaux plats (`m_framesExponentialSE3Vectors`) | Classes `FrameInfo` / `SectionInfo` |
| **Lisibilité** | Moyenne (mélange types SOFA/calculs) | Meilleure (types homogènes, encapsulation) |
| **Testé** | Oui (intégration continue) | En cours (tests unitaires partiels) |
| **matB_trans** | Bogue corrigé : `k<6` pour Vec6Types | 3×6 fixe (Vec3Types uniquement par défaut) |

### 13.8 Points d'attention

1. **Inversion de frame** : `frame.getInverseTransformation()` calcule $g_s^{-1}$ à chaque appel dans `applyJ`. Il serait préférable de pré-calculer et de mettre en cache l'adjoint inverse dans `FrameInfo`.

2. **Lecture de `framePositions`** : dans `applyJT(VecDeriv)` et `applyJT(MatrixDeriv)`, les positions absolues des frames sont relues depuis SOFA pour construire $P_{\text{trans}}$. Ces matrices de projection pourraient être pré-calculées lors de `apply` et stockées dans `FrameInfo`.

3. **`matB_trans` = 3×6** : la matrice $B^T$ est fixée à $[I_3 \mid 0_3]$, ce qui convient uniquement pour les DOF Vec3Types. Si un type Vec6Types est instancié, il faudra la passer à $I_6$ (même bug que dans `DiscreteCosseratMapping.cpp`).

4. **Boucle de compression des nœuds** (dans `applyJT(MatrixDeriv)`) : le `while` de compression est fragile — voir le commentaire dans le code signalant un ancien bug sur la condition de sortie. Une réécriture avec un algorithme de `group_by` explicite serait plus sûre.

---

## 14. Plan d'amélioration général des fonctions Apply/ApplyJ/ApplyJT

### 14.1 Problèmes transverses identifiés

| # | Problème | Impact | Fichiers concernés |
|---|----------|--------|--------------------|
| P1 | `buildProjector` / `buildProjectionMatrix` : deux API différentes pour la même opération | Fragmentation, risque d'incohérence | `DiscreteCosseratMapping`, `Strain2RigidCosseratMapping` |
| P2 | Pré-calculs dispersés entre `apply`, `applyJ` (tableaux `m_framesExp…`, `m_nodesTangExp…`) | Couplage fort, pas thread-safe | `DiscreteCosseratMapping.inl` |
| P3 | Accès répété aux positions SOFA dans `applyJT` (reconstruction de la frame à chaque appel) | Performances | `Strain2RigidCosseratMapping.inl` |
| P4 | `matB_trans` (Vec6Types) : boucle `k<3` au lieu de `k<6` | Bug fonctionnel (forces coupées) | `DiscreteCosseratMapping.cpp`, potentiellement `Strain2RigidCosseratMapping.inl` |
| P5 | Boucle `while` de compression des nœuds dans `applyJT(MatrixDeriv)` | Fragile, difficile à maintenir | `DiscreteCosseratMapping.inl`, `Strain2RigidCosseratMapping.inl` |
| P6 | Duplication de la logique de rétro-propagation co-adjointe | Maintenance double | idem |
| P7 | Vérification de dualité `⟨w, Jξ̇⟩ = ⟨J^Tw, ξ̇⟩` absente dans `DiscreteCosseratMapping` | Risque de bug silencieux | `DiscreteCosseratMapping.inl` |

### 14.2 Plan d'action proposé

#### Étape 1 — Corrections de bugs immédiates

- [x] **P4** : `matB_trans` → `k<6` dans `DiscreteCosseratMapping.cpp` (✅ fait dans ce diff)
- [ ] **P4** : Vérifier `matB_trans` dans `Strain2RigidCosseratMapping.inl` pour toute instanciation Vec6

#### Étape 2 — Unification de l'API de projection

Créer une fonction libre (ou méthode de `SE3`) :

```cpp
/// Matrice de projection convention SOFA Rigid6 ↔ SE(3) Cosserat.
/// Retourne P telle que v_cosserat = P * v_sofa.
/// @param frame : repère courant (rotation R incluse)
template<typename Scalar>
Eigen::Matrix<Scalar,6,6> buildProjectionMatrix(const SE3<Scalar>& frame);
```

Remplacer tous les appels à `buildProjector(frame)` (SOFA) et `buildProjectionMatrix(R)` (Eigen) par cette fonction unifiée.

#### Étape 3 — Encapsuler les pré-calculs dans des structures de données

Enrichir `FrameInfo` / `SectionInfo` (ou créer leur équivalent dans `DiscreteCosseratMapping`) :

```cpp
struct SectionCache {
    SE3Types  g_local;          // expCosserat(ξ, L)
    Matrix6   Ad_inv;           // Ad_{g^{-1}} — pour applyJ
    Matrix6   coAdjoint;        // Ad_{g^{-T}} — pour applyJT
    Matrix6   J_local;          // computeTangExp — pour applyJ / applyJT
    Matrix6   P;                // buildProjectionMatrix — pour applyJT
};
```

`apply` remplit le cache. `applyJ` et `applyJT` ne font que lire.

#### Étape 4 — Refactoriser la rétro-propagation en algorithme générique

```cpp
/// Rétro-propagation co-adjointe d'un wrench depuis la tip vers la racine.
/// @param sections  : cache des sections (tip → root)
/// @param tip_wrench : wrench initial au bout
/// @param on_section : callback appelé à chaque section avec le wrench transporté
///                     → renvoie la force de déformation locale
template<typename SectionIt, typename Callback>
void backward_coAdjoint(SectionIt begin, SectionIt end,
                        const TangentVector& tip_wrench,
                        Callback on_section);
```

Les deux `applyJT` (VecDeriv et MatrixDeriv) appellent ce même algorithme avec des callbacks différents.

#### Étape 5 — Tests de dualité automatisés

Ajouter dans le test `test_DiscreteCosseratMapping` (ou équivalent) :

```cpp
// Pour tout vecteur aléatoire ξ̇ et w :
// ⟨w, applyJ(ξ̇)⟩ == ⟨applyJT(w), ξ̇⟩
EXPECT_NEAR(w.dot(eta), q.dot(xi_dot), 1e-10);
```

Ce test détecte immédiatement tout bug dans `applyJT` (ordre des opérations, `matB_trans` incorrect, etc.).

#### Étape 6 — Documentation inline

Ajouter dans chaque fonction `applyJ` / `applyJT` un commentaire court renvoyant aux équations de ce document :

```cpp
// Implémente : η_k = Ad_{g_k^{-1}} · (η_{k-1} + J_local_k · ξ̇_k)
// Voir docs/mapping_explanation.md §3
```
