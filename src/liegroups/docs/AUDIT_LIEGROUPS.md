# Audit de `src/liegroups`

Branche `claude/liegroups-audit-b60c6e` : 24 commits au-dessus de `eb334be` (`git log eb334be..HEAD`).
Le message de chaque commit donne la cause, le correctif, l'erreur mesurée avant/après et le test.

Les erreurs ci-dessous sont absolues, mesurées contre des références en `long double` calculées
à partir de la définition matricielle des groupes (séries de matrices, Jacobiennes par la série
Σ ad^k/(k+1)!), jamais avec les formes fermées de la bibliothèque (`Tests/unit/LieGroupReference.h`).
Angles testés : de 1e-9 à π−1e-9.

## Commits

`8983acd` SE3 hat/vee/isApprox · `4d2e2a2` slerp · `83d962f` BCH ordre 5 · `8f21704` dexp/dlog génériques ·
`5725557` coefficients SO3 · `e056e87` SE3 exp · `b51fa51` SE3 log · `9d1c111` SE3 Q/Jacobiennes ·
`989f0ad` SE2 constructeurs · `17f1961` SE2 ad/Ad · `f4717bb` SE2 exp/log · `69e8bfc` repli Jacobiennes ·
`94ab6d0` SE23 exp/log · `852dc5d` SE23 hat/vee/ad/Ad · `afb70e3` SE23 isApprox · `0e584c1` SGal3 ·
`4f18ce4` fichiers `.inl` · `03c7eec` SE2/RealSpace isApprox · `eab71b9` Bundle · `42ef25e` Types ·
`fc9ac68` expCosserat · `db8a429` tests mal spécifiés · `e0e4ec4` SO3Test · `dd3412a` float.

Les trois corrections de `Types.h` (`SO3Log`, `SE3Exp`, `SE3Log`) sont dans un seul commit.

## Constats

Les tests sont dans `Tests/unit/` sauf mention contraire.

### 1. Fonctions génériques de `LieGroupBase` (`LieGroupBase.inl`)
- `dexp_alternative`, `dexpInv_alternative`, `dlog` : mauvais coefficients (ceux de dlog pour dexp, forme SO(3)
  appliquée à SE(3)). Erreur d'ordre 1 pour SE3 (0,6 à 0,8), 1e-16 après correction.
- `SE3::dlog()` renvoyait la co-adjointe : override supprimé, la version de base s'applique. Aucun appelant.
- BCH ordre 5 : termes ne faisant pas partie de la série de Dynkin ; l'ordre 5 ne faisait pas mieux que l'ordre 4
  (3,75e-9, puis 1,2e-11 après correction).
- `slerp` reparamétré par sin(t d)/(…) : `slerp(0,25)` s'écartait de `interpolate(0,25)` de 0,125 rad, < 1e-12 après.
- Tests : `test_LieGroupBase_Generic`.

### 2. Membres qui ne compilaient pas à l'instanciation
Un membre de template n'échoue qu'à l'instanciation : `-fsyntax-only` sur les headers ne les voit pas.
- `SE3` sans `computeHat`/`computeVee`, `computeIsApprox` sans epsilon par défaut : BCH(order), `matrixLog`,
  `distance`, `squaredDistance`, `slerp` de la base inutilisables pour SE3 (`test_LieGroup_Instantiation`).
- `leftJacobian`/`rightJacobian` (et inverses) : hooks `computeLeftJacobian`… inexistants pour SO2, SE2, RealSpace,
  SE23, SGal3, Bundle. Repli sur la série Σ ad^k/(k+1)! (`test_SE2_Audit`).
- `SE2(Matrix3)` et `SE2(TangentVector)` : non compilables.
- `Bundle` ne compilait pour aucune combinaison de groupes : `Dim`/`AlgebraDim`/`ActionDim` confondus
  (SE3 : 4/6/3, SO2 : 2/1/2, SE2 : 3/3/2), `actionDimension()` inexistant, `exp` non statique, `random` inexistant
  (`test_Bundle_Audit`).
- `SGal3` inutilisable : pas de `compose`, hooks privés, `BCH(X,Y)` déclaré mais défini seulement dans un `.inl`
  non inclus, action lisant les indices 6..9 d'un vecteur de taille 7. Réécrite comme groupe de Kelly
  X = [[R, v, p], [0, 1, t], [0, 0, 1]] (5×5), action de dimension 7 comme déclaré (`test_SGal3_Audit`).

### 3. `SE3.inl`, `SE23.inl`, `SGal3.inl` jamais inclus
`fromEulerZYX`/`toEulerZYX` jamais définies, redéclaration de template avec argument par défaut, `BCH` défini mais
non déclaré, `dualMatrix` avec la mauvaise disposition `[ρ;φ]`. Corrigés et inclus en fin de header
(`SE3InlFreeFunctions`, `SE23InlFreeFunctions`, `SGal3InlFreeFunctions` dans `test_LieGroup_Instantiation`).

### 4. Annulations numériques (`test_SE3_Numerics`)

| Fonction | Avant | Après |
|---|---|---|
| Coefficients SO3 (commutation à 1e-4) | 2,2e-13 | < 3e-15 |
| `SE3::computeExp` | 5,3e-13 | < 4e-15 |
| `SE3::computeLog` près de π (forme (1+cos)/sin) | 1,2e-9 | < 2e-14 |
| `computeQ` et Jacobiennes SE3 (Jl, Jr, inverses) | 7,6e-9 (9,8e-9 vers 2e-4) | < 4e-15 |
| `SE3::expCosserat` à \|φ\| = 1e-9 | 4,6e-10 | < 4e-15 |
| `SE3f` (float) | 1,7e-4 | < 6e-7 |

- `computeQ` avait aussi une branche petit angle fausse (coefficient 1/12 au lieu de 1/24, termes O(t²) perdus).
- `expCosserat` passe maintenant par `computeExp` : l'unité d'élongation O(1) contenue dans ξ̂ multipliait des
  coefficients qui s'annulent.
- Float : la série est utilisée jusqu'à 1 rad (0,25 rad en double).

### 5. SE2 (`SE2.h`, `test_SE2_Audit`)
- Log faux pour tout angle : cot(θ) au lieu de cot(θ/2), signe des termes hors diagonale, facteur 2
  (erreur ~1 à 1,9e9).
- `computeAd` n'était pas `ad` ; colonne de translation de `Ad_g` de mauvais signe (3,1e-1, maintenant 0).
- Exp : 9,5e-10, maintenant < 4e-15.

### 6. SE23 (`SE23.h`, `test_SE23_Audit`)
- Exp : vitesse divisée par θ (1,8e9), le log avait l'erreur miroir (4,7) ; maintenant < 1e-14 / < 2e-14.
- `hat`, `vee`, `ad` étaient des zéros ; `Ad` avait le bloc [u]× R en (0,6) au lieu de (6,3).

### 7. `Types.h` (`test_Types_Audit`)
`SO3Log` (4,3e-7 vers 0,014 rad, jusqu'à 5 rad d'erreur près de π), `SE3Exp` (5e-9), `SE3Log` (de 1e-6 à 5 rad
pour tout angle > 1e-9), `logSO3`. Tous < 5e-15 après correction (borne élargie près de π : la matrice est arrondie).

### 8. `isApprox` relatif d'Eigen contre un vecteur nul
SE3, SE2, SE23, SGal3, RealSpace : g·g⁻¹ n'était jamais « égal » à l'identité (translation de 1e-17 contre 0).
Passé en tolérance absolue pondérée par max(1, |t|).

### Vérifiés sans défaut
`SO3::computeLog` (cas qw<0, branche sous epsilon), `Twist.h`, `Wrench.h`, `SO2.h`.

## Utilisateurs dans `src/Cosserat` de ce qui change de comportement
- `SE3::expCosserat` : `Strain2FramesCosseratMapping.inl:279,289`. Valeurs inchangées à ~1e-15 pour une courbure normale.
- `SE3::computeRightJacobianInverse` : `Frames2StrainCosseratMapping.inl`, trois sites. Plus précis, même convention.
- `SE3::computeIsApprox` : `CosseratBeamGeometry.h:207`. Tolérance absolue pour la translation (avant : relative).
- `improvedBaseCosseratMapping.h` appelle `SE3d::hat`, qui compile maintenant ; `SE3d::dexp` est inchangé (Jacobienne droite).
- `Bundle`, `SE23`, `SGal3`, `Types::SE3Exp/SE3Log/SO3Log`, `SE3::dlog` : aucun utilisateur hors `liegroups`.

## Tests existants

Chaque fichier a été compilé en standalone (gtest, stub de `sofa/testing/BaseTest.h`).

- **Corrigés** (`db8a429`, `e0e4ec4`) : ils échouaient sur du code correct.
  - 5 tests de `differentiation/test_analytical_jacobians.cpp` : différences finies mal construites (dérivée de
    `log` du produit entier, en R = identité) ; `SO3_JacobiansConsistency` attendait −I au lieu de +I.
  - `SO3_LogJacobian` : d/dδ log(exp(δ) R) vaut Jl⁻¹(log R), pas l'identité.
  - `LieGroupIntegrators` (3 tests) : oubli de l'élongation unité (l'erreur de pose 0,4 vaut exactement L) ; bornes
    d'ordre décalées de un (erreur locale en h^(p+1)).
  - `ZeroNoiseStraightRod` : « déformation nulle » est une translation de `length`, Ad ≠ I.
  - `UtilsTest.Sinc` : sinc est paire. `SO3Test.IsApprox` : 2 epsilon d'angle reste dans la tolérance sur le quaternion.
- **Enregistrés** dans `Tests/CMakeLists.txt` : les huit nouveaux fichiers de test, plus `SE3Test`, `SGal3Test`, `SO3Test`.
- **Non listés** : `integration/LieGroupIntegrationTest.cpp`, `test_cosserat_bindings.py`, `test_extended_liegroups.py`
  (`benchmarks/` a son propre CMake).
- **Toujours en échec, non diagnostiqués** : `test_parameter_estimation` (1), `test_ilqr_control` (1),
  `test_trajectory_optimization` (4), `TypesTest` (`IsZero`, `IsApprox`, `NormalizeAngle`), `UtilsTest` (normalisation et
  distance d'angle). Soupçon : tests incohérents avec le code, non vérifié.

## Laissé inchangé volontairement

- `SO3::dexp`/`dlog` sont la Jacobienne **gauche**, alors que `LieGroupBase::dexp` et `SE3::dexp` sont **droits**.
  `DexpIsLeftJacobianAndDlogItsInverse` fixe la version gauche et `SE3::computeLeftJacobian` en dépend.
  `SE3::dlog` suit la base (droit).
- Conventions de perturbation mélangées : `composeJacobians` droite-droite, `inverseJacobian` gauche-gauche
  (−Ad_{g⁻¹}), `actionJacobians` perturbe à gauche, `actionJacobian` numérique de la base perturbe à droite.
  Chaque formule est vérifiée par différences finies dans sa propre convention (`db8a429`), mais l'ensemble n'est pas
  cohérent : convention à trancher.
- Dispositions du vecteur tangent : SE3 `[φ;ρ]`, SE2 `[ρ;θ]`, SE23 `[v;w;a]`, SGal3 `[ρ;w;β;τ]`,
  `Types::SE3Exp/SE3Log` `[ρ;φ]`. SE23 et SGal3 suivent leur documentation, pas la convention « angulaire d'abord ».
- `Types.h` : `normalizeAngle(π)` renvoie −π (intervalle [−π, π)), alors que `TypesTest`/`UtilsTest` attendent π ;
  `isZero` utilise 100·epsilon, les tests supposent 1e-12 ; `cosc` a une branche petit angle qui recopie `sinc`
  (faux), mais `TypesTest.Cosc` fixe `cosc(0) = 1`.
- `SE3::distance` est une somme pondérée ; `SE2::interpolate` n'est pas géodésique ; la base lève
  `NumericalInstabilityException` dans `BCH` dès que la norme dépasse π/2 ; `SE3::computeTangent` a une sémantique
  obscure ; `SE3::projector()` duplique la co-adjointe ; `expCosseratSmall`/`expCosseratGeneral`/`buildXiHat` de
  `SE3.h` sont inutilisés mais gardés (cités par `STRAIN_CONVENTION.md`) ; `Wrench::transformBy` utilise Adᵀ (pull-back).

## Non vérifié

- Pas de build SOFA ni de lien de l'exécutable de tests (le dossier `cmake-build-relwithdebinfo` n'a pas été touché).
- Tests nécessitant SOFA : `AdvancedFeatures`, `HookeSerat_*`, `Phase2`, `MultiSectionBeam`, `jacobian_verification`
  (besoin de `Cosserat/config.h`), `LieGroupIntegrationTest` (besoin de `Cosserat/liegroups/Bundle.h`). Ils utilisent
  SE3/SO3, que j'ai modifiés.
- Passe `-fsyntax-only` sur les 27 unités de compilation du plugin (flags de `compile_commands.json` réécrits vers
  le worktree) : 27/27 code 0, aucune erreur ni nouvel avertissement dans `liegroups`. Ne remplace pas une compilation
  complète ni l'exécution des scènes.
- Float : seul SE3 en float a été testé numériquement.
