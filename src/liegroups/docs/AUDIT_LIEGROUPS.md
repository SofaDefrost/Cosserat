# Audit de `src/liegroups`

Branche `claude/liegroups-audit-b60c6e` : 36 commits au-dessus de `eb334be` (`git log eb334be..HEAD`).
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

Deuxième passe sur SO3 / SE3 / `Types.h` (§9) : `3a3b20c` Types sinc2/sinc3/oneMinusCos · `0882c91` tolérance par défaut
d'`isApprox` · `3869966` tirages aléatoires (Haar) · `f9576b5` constructeurs SE3(q,t) et SO3(Matrice) · `76ecb86` BCH ordres 4 et 5 ·
`e62535a` code mort et `m_quat` privé.

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

### 9. Deuxième passe : SO3, SE3 et le toolkit scalaire de `Types.h`
Il n'existe pas de fichiers `SO3Types`/`SE3Types` : j'ai relu `SO3.h/.inl`, `SE3.h/.inl` et `Types.h` (alias, tolérances,
fonctions sinc, tirages). Tests : `unit/test_Types_Numerics.cpp`, `unit/test_SO3_SE3_Audit.cpp`.

| Défaut | Mesure avant → après | Commit |
|---|---|---|
| `Types::sinc2`, `sinc3`, `oneMinusCos` : formes fermées qui s'annulent, bascule série à 1e-4 seulement (même défaut que celui corrigé dans les coefficients de Jacobienne, qui restait ici) | erreur relative sinc2 5,2e-9 à 1,01e-4, sinc3 7,6e-9, `oneMinusCos` 100 % à 1e-9 et 8,3e-8 à 1e-5 → < 5e-15 | `3a3b20c` |
| Tolérance par défaut d'`isApprox` / `isNearIdentity` = epsilon machine (2,2e-16), donc (g·h)·h⁻¹ ≠ g à cause d'un ulp | échecs sur 10000 comparaisons : SO3 126, SE3 4743, SE2 3563, SE23 6367 → 0 (défaut `Types::tolerance()` = 100 ε) | `0882c91` |
| `SO3::computeRandom` / `SE3::computeRandom` : angle uniforme sur un axe aléatoire, pas la mesure de Haar ; `randomUnitVector` normalisait un vecteur du cube | angle moyen 1,5714 → 2,2074 (Haar : π/2 + 2/π), trace moyenne 1 → 0, E\|axe.x\| 0,516 → 0,500 | `3869966` |
| `SE3(Quaternion, t)` repassait par `toRotationMatrix()` : pour un quaternion non unitaire la matrice n'est pas une rotation ; `SO3(Matrice)` ne normalisait pas le quaternion | norme 1,4142 → 1, erreur de matrice 1,0 → < 1e-15 ; norme 1,0039 pour 1,01·R → 1 | `f9576b5` |
| `SO3::BCH` documente les ordres 1 à 5 mais s'arrêtait à 3 | erreur identique pour les ordres 3, 4, 5 (7,7e-5 à \|v\|,\|w\| ≈ 0,25) → ordre 4 : 3,7e-5, ordre 5 : 3,9e-6, pentes h⁵ et h⁶ ; égal au BCH générique à 1e-15 | `76ecb86` |
| `SO3.h` : `exp_`, `expRodrigues` (doublons d'`exp`), `expCosseratExact` et un bloc « TODO » commentés, `m_quat` public ; `SE3.h` : deux générations de `computeExp/computeLog` en commentaire | aucun changement de comportement ; l'invariant de norme unité ne peut plus être cassé de l'extérieur | `e62535a` |

Vérifiés sans défaut (mesurés) : `SO3::computeLog` de 1e-300 à 1e-10 rad (erreur relative 0), `SO3::expCosserat(k, L)` =
`exp(L k)` à 1e-12 près pour |k| de 1e-14 à 5 et L de 1e-3 à 10, `SO3::interpolate` entre quaternions à 1e-9 (7,9e-17),
quaternions antipodaux (`isApprox` vrai, distance 0), `SE3::expCosserat` sans déformation (translation `[L,0,0]`),
aller-retour `exp(log)` en float pour SE3 (7e-7), symétrie et invariance à gauche de `SE3::distance`, Euler ZYX
(l'erreur de 1e-8 près de ±π/2 est celle du conditionnement des angles d'Euler : atan2 donne le même ordre).

### Vérifiés sans défaut
`SO3::computeLog` (cas qw<0, branche sous epsilon), `Twist.h`, `Wrench.h`, `SO2.h`.

## Utilisateurs dans `src/Cosserat` de ce qui change de comportement

(Mise à jour : le plugin appelle maintenant l'API publique — `exp`, `inverse`, `adjoint`, `log`, `act`, `Identity`,
`rightJacobianInverse`, `coAdjoint` — dans `BeamStateEstimator.h`, `BeamShapeInterpolation.h`, `CosseratBeamGeometry.h/.inl`,
`Frames2StrainCosseratMapping.inl`, `Strain2FramesCosseratMapping.inl`, `CosseratILQRController.inl` et
`CosseratBeamGeometryEngine.cpp`.)
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
- Deuxième passe : `SO3::interpolate` borne `t` à [0, 1] alors que `SE3::interpolate` extrapole ; `SO3::BCH` ne lève pas
  d'exception pour de grandes normes alors que la base le fait ; `SE3::computeTangent` renvoie `[R⁻¹(p−t); log R⁻¹]`
  (position d'abord, contraire à `[φ;ρ]`) et n'a aucun appelant ; `Types::cosc` (voir plus haut, aucun appelant) ;
  `SO3::expCosserat` fait `exp(L k)` avec un test de petit angle sur |k| seul (sans effet : le seuil est epsilon).
- Parmi les propositions de « Questions de conception », ne sont pas appliquées : hooks `protected` (casserait `cosserat_shell`),
  suppression de `SO3::expCosserat*`, déplacement des fichiers Cosserat hors de `liegroups` (§3 et §4, décision : ils servent
  à `cosserat_shell` et à la version « painless »).

## Questions de conception

Statut après la demande de suivi (commits `32ee2c9`, `04a9089` et le commit « hooks » qui les suit) :
- **§1 appliqué, sauf la visibilité des hooks.** `compose` est devenu le hook `computeCompose` (la base expose `compose`
  public), les wrappers dupliqués de `SO3`, `SE3` et `Bundle` sont supprimés, `computeCoAdjoint` → `coAdjoint`, `computeQ` →
  `jacobianQ`, `computeBCH` → `bchSecondOrder`, et le plugin passe par l'API publique (`exp`, `inverse`, `adjoint`, `log`,
  `act`, `Identity`, `rightJacobianInverse`, `coAdjoint`). Les hooks restent **publics** : `cosserat_shell` appelle encore
  `SO3::computeHat`, `SO3::computeRightJacobian` et `SE3::computeAdjoint` (équivalents publics : `hat`, `rightJacobian`,
  `adjoint()`), et ce dépôt ne doit pas être modifié d'ici. On passera les hooks en `protected` + `friend` dès qu'il sera migré.
  La convention est écrite dans le commentaire de `LieGroupBase`.
- **§2 appliqué** : doublons matriciels de `Types.h` supprimés. `static` conservé.
- **§3 et §4 non appliqués, par décision** : ces fichiers servent au projet `cosserat_shell` ou à la version « painless »
  de Cosserat en cours d'implémentation.
- **§6 appliqué** : un seul `hat` par groupe (`SO3::buildAntisymmetric`, `dualMatrix`, `Types::skew3/unskew3` supprimés).
  `Twist::hat()` reste (méthode d'un objet `Twist`).
- Reste en l'état : `SE3::computeTangent` (non-hook nommé `compute*`, sémantique obscure, sans appelant) ; `print` et
  `getTypeName` (hooks sans préfixe).

Les propositions sont décrites ci-dessous telles qu'elles ont été faites.

### 1. Pourquoi `computeXxx` et `xxx` coexistent-ils ?

C'est le motif CRTP de `LieGroupBase`, et le préfixe a une raison précise.
- L'API publique (`exp`, `log`, `inverse`, `adjoint`, `act`, `isApprox`, `Identity`, `hat`, `vee`, `ad`, `BCH`, les
  Jacobiennes) est écrite une fois dans la base. Elle appelle un « point de personnalisation » du groupe dérivé :
  `Derived::computeExp`, `computeInverse`, `computeHat`, etc.
- Si le point de personnalisation s'appelait lui aussi `exp`, un dérivé qui l'oublie ferait que `Derived::exp` résout
  vers la base elle-même : récursion infinie. C'est la classe de défaut trouvée au départ (`ad` ↔ `computeAd` dans
  `SO3`). Le préfixe `compute` évite ça : un hook manquant donne une erreur de compilation, pas une boucle.

Ce n'est donc pas « tout est `compute` » : `compute*` désigne les hooks, le reste l'API. Mais l'application actuelle est
incohérente, sur quatre points :
1. `compose` est un hook sans préfixe (la base appelle `derived().compose(other)`) : pas de risque de récursion, car la base
   n'a pas de `compose` public, mais le nom casse la règle. Idem `print` et `getTypeName`.
2. `SO3` redéclare à la main des wrappers publics (`exp`, `log`, `inverse`, `adjoint`, `act`, `isApprox`, `identity`) qui
   masquent ceux de la base sans rien ajouter ; `SE3` en fait autant pour `log` et `adjoint`.
3. Des fonctions nommées `compute*` ne sont pas des hooks : `computeCoAdjoint`, `computeQ`, `computeBCH`, `computeTangent`.
4. Tous les hooks sont `public` dans les dérivés. Le plugin les appelle donc directement : `SE3Types::computeExp` (9 sites),
   `.computeInverse()` (6), `.computeAdjoint()` / `.computeCoAdjoint()`, `computeRightJacobianInverse` (3), `computeIdentity` (3),
   alors que l'API publique existe.

Proposition : hooks `protected`/`private` avec `friend LieGroupBase<...>` ; renommer `compose` en `computeCompose` ;
supprimer les wrappers dupliqués de `SO3`/`SE3` ; faire passer le plugin par `exp`/`inverse`/`adjoint`/`Identity`.
Coût : une dizaine de sites dans le plugin, mécanique ; risque faible une fois les tests de `test_LieGroup_Instantiation`
en place. `computeCoAdjoint` n'a pas d'équivalent dans la base (il faudrait un `coAdjoint()` public).

### 2. Des fonctions `static` : bonne idée ? Faut-il les sortir de `liegroups` ?

Oui pour `static`, non pour les sortir, avec deux réserves.
- Sont statiques les fonctions qui ne dépendent d'aucun élément du groupe : `exp(ξ)` (fabrique), `hat`, `vee`, `ad`, `BCH`,
  les Jacobiennes d'un vecteur tangent, `dexp`, `Identity`, `Random`. C'est mathématiquement exact (elles agissent sur
  l'algèbre ou construisent un élément) et c'est ce que font Sophus, manif et Pinocchio (`exp6`, `Jexp6`).
- Les sortir en fonctions libres ne changerait rien à l'appel (`SE3d::exp(x)` contre `exp<SE3d>(x)`) et casserait la
  dispatch CRTP. Les ranger dans un type `LieAlgebra<G>` / `Tangent<G>` (la base en déclare un, `LieAlgebra`, jamais défini)
  serait plus propre conceptuellement (groupe : `compose`, `inverse`, `act` ; algèbre : `hat`, `ad`, Jacobiennes, `BCH`),
  mais c'est une refonte d'API sans gain fonctionnel : pas recommandé maintenant.
- Réserve 1 : `Types.h` est un fourre-tout (constantes, RNG, interpolation, et des copies matricielles de `SO3Log`,
  `SE3Exp`, `SE3Log`, `interpolateSE3`, `SE3Distance`, `generateSE3Trajectory`) avec une autre disposition du vecteur
  tangent (`[ρ;φ]`). Aucun appelant dans le plugin, seulement des tests : à supprimer ou à déprécier plutôt qu'à maintenir
  en double (je les ai corrigées parce qu'elles étaient publiques et fausses).
- Réserve 2 : les fonctions libres des `.inl` (`interpolate`, `operator/`, `fromEulerZYX`, `dualMatrix`…) sont déjà hors
  des classes, ce qui est bien ; mais `dualMatrix` n'est qu'un alias de `computeHat` et `interpolate` de la méthode membre.

### 3. Les fonctions « Cosserat » ont-elles leur place dans `liegroups` ?

Elles ont une raison d'être, mais pas à cet endroit.
- `SE3::expCosserat(strain, L)` vaut `exp(L (strain + [0; e_x]))`. Le `+ e_x` est une convention de modélisation (poutre à
  élongation unité au repos, voir `STRAIN_CONVENTION.md`), pas de la théorie des groupes. Elle est utilisée par le plugin
  (`Strain2FramesCosseratMapping.inl:279,289`) et par d'autres en-têtes de `liegroups` (`LieGroupIntegrators.h`,
  `CosseratUncertaintyPropagator.h`, `CosseratBodyJacobian.h`).
  Proposition : une fonction libre `expCosserat(strain, L)` côté plugin (par exemple dans `mapping/SofaLieGroupsUtils.h`),
  qui appelle `SE3::exp`. `liegroups` reste une bibliothèque de groupes.
- `SO3::expCosserat`, `expCosseratFromOmega`, `expCosseratExact` : aucun appelant, et `SO3::expCosserat(k, L)` est exactement
  `exp(L k)` (pas d'élongation sur SO(3)). Le même nom a donc deux sens selon la classe. À supprimer.
- `SE3::expCosseratSmall`, `expCosseratGeneral`, `buildXiHat` : inutilisés depuis `fc9ac68` (seul `STRAIN_CONVENTION.md` les cite).
  À supprimer en même temps que la doc.
- Fichiers de `liegroups/` qui ne sont pas de la théorie des groupes : `CosseratBodyJacobian.h`,
  `CosseratUncertaintyPropagator.h`, `control/CosseratILQRController.h`, `calibration/CosseratParameterEstimator.h`,
  `optimization/CosseratTrajectoryOptimizer.h`, `KinematicChain.h`, `RevoluteJoint.h`, `MultiDofJoints.h`, `BezierSE3.h`,
  `LieGroupIntegrators.h` (environ 3 000 des 8 950 lignes), plus `python/articulated_rigids_py.cpp`. Seuls
  `CosseratBodyJacobian.h` et, via `BeamStateEstimator.h`, `CosseratUncertaintyPropagator.h` et `GaussianOnManifold.h` sont inclus
  par le plugin ; `control/`, `calibration/`, `optimization/` ne sont utilisés que par leurs tests (le plugin a son propre
  `Cosserat/controller/CosseratILQRController.h`, qui n'inclut pas celui de `liegroups`).
  Proposition : un répertoire ou module à part (par exemple `src/Cosserat/geometry/` ou `liegroups_cosserat/`), avec
  `liegroups/` limité à `Types`, `LieGroupBase`, `SO2/SO3/SE2/SE3/SE23/SGal3/RealSpace/Bundle`, `Twist`, `Wrench`,
  `GaussianOnManifold`.

### 4. Utilité des fonctions de `liegroups` dans le plugin

Ce que le plugin utilise réellement (compte de sites dans `src/Cosserat`) :
- `SE3`, `SO3` uniquement. Appels statiques : `computeExp` (9), `exp` (4), `computeRightJacobianInverse` (3), `computeIdentity` (3),
  `hat` (2), `expCosserat` (2), `dexp` (1), `ad` (1). Méthodes d'instance : `rotation` (~31), `matrix` (~26), `inverse` (~25),
  `act` (~17), `log` (~15), `computeInverse` (6), `computeAdjoint` (6), `translation` (5), `compose` (4), `adjoint` (4),
  `computeCoAdjoint` (2), `computeLog` (2), `interpolate`, `distance`, `computeAction`, `quaternion` (1 chacun).
- `distance` : `CosseratBeamGeometry.h` (somme pondérée rotation + translation, pas une distance géodésique).
- `CosseratBodyJacobian`, `Twist`/`Wrench` (alias seulement dans `Strain2FramesCosseratMapping.h`), `GaussianOnManifold`,
  `CosseratUncertaintyPropagator` (via `BeamStateEstimator.h`).
- `Types<double>` : seulement les constantes et `epsilon` dans les forcefields (`BaseBeamForceField.h`, `HookeSeratBaseForceField.h`).

Ce que le plugin n'utilise pas du tout : `SO2`, `SE2`, `SE23`, `SGal3`, `RealSpace`, `Bundle`, `KinematicChain`,
`RevoluteJoint`, `MultiDofJoints`, `AutodiffSupport`, `BezierSE3` et `LieGroupIntegrators` (inclusions commentées dans
`Strain2FramesCosseratMapping.h`), `BCH`, `dlog`, `vee`, `lieBracket`, `composeJacobians`, `inverseJacobian`,
`actionJacobians`, `parallelTransport`, `computeTangent`, `projector`, `Types::SE3Exp/SE3Log/SO3Log`. `improvedBaseCosseratMapping.h`
appelle `SE3d::hat/ad/dexp` et inclut `<Cosserat/liegroups/Utils.h>` qui n'existe pas dans l'arbre : à vérifier s'il est
encore compilé (ce serait un en-tête mort).

Lecture : le plugin s'appuie sur ~15 fonctions de `SE3`/`SO3`. Le reste (environ les trois quarts du code) est de la
généralité non utilisée : sa valeur est une bibliothèque réutilisable (c'est ce que veut le projet CosseratShell qui
dépend de `SO3`/`SE3`), mais chaque ligne non testée a produit un défaut dans cet audit. Soit on teste et on garde (c'est
ce qui a été fait ici pour SE2, SE23, SGal3, Bundle), soit on supprime ce que personne ne veut (`Types::SE3*`, `SO3::expCosserat*`).

### 5. `SE3`, `SE2`… peuvent-ils hériter des types SOFA (`Rigid3d`) ?

Techniquement oui, mais c'est déconseillé, pour quatre raisons vérifiées dans les sources de SOFA.
1. `RigidCoord<3,Real>` stocke un `sofa::type::Vec<3>` et un `sofa::type::Quat` (ordre `[x,y,z,w]`) ; `SE3` stocke un
   vecteur et un quaternion Eigen (ordre `[w,x,y,z]`). Hériter obligerait à choisir un seul stockage, donc à perdre Eigen
   (SIMD, `toRotationMatrix`, `slerp`, autodiff via `Scalar`) ou à convertir à chaque appel.
2. Les opérateurs de `RigidCoord` n'ont pas le sens de ceux d'un groupe : `operator*(RigidCoord)` renvoie un **réel**
   (produit scalaire), `operator+(Deriv)` est une mise à jour « translation + exponentielle de la rotation », pas
   `g · exp(ξ)`, `operator-` est une différence de coordonnées. Les mélanger avec `LieGroupBase::operator*` serait une source
   de défauts silencieux.
3. Les dispositions diffèrent : `RigidDeriv<3>` est `[v_centre ; ω]` (linéaire d'abord), `SE3::TangentVector` est `[φ ; ρ]`
   (angulaire d'abord).
4. `liegroups` est aujourd'hui une bibliothèque d'en-têtes sans dépendance à SOFA (c'est ce qui permet de la compiler et de la
   tester seule, comme ici, et de la faire utiliser par un autre projet). Hériter de `RigidCoord` la lierait à SOFA.

Ce qui existe déjà et suffit : `mapping/SofaLieGroupsUtils.h` (`rigidCoordToSE3`, `se3ToRigidCoord`, utilisées sur 14 sites).
Deux améliorations possibles sans héritage : (a) `rigidCoordToSE3` passe par `Eigen::Quaternion -> toRotationMatrix -> SO3(Matrix)`
alors que `SO3(Quaternion)` existe : on évite un aller-retour matrice ; (b) si l'on veut des `SE3` dans des `Data<>` SOFA, il
faut au minimum `operator<<`/`>>` (présents pour `SO3`, pas pour `SE3`) ; en faire un vrai `DataTypes` (comme `Rigid3Types`)
fermerait la porte aux mappings, forcefields et solveurs SOFA qui attendent `Rigid3Types`.

### 6. Différence entre `computeHat` et `hat`

- `computeHat` est le **hook** du groupe : fonction statique qui construit la représentation matricielle d'un vecteur tangent.
- `hat` est l'**API** : statique, dans `LieGroupBase`, avec un `static_assert` qui vérifie que le dérivé définit
  `computeHat`, puis `return Derived::computeHat(v)`. Même résultat ; l'appelant écrit `SE3d::hat(x)`.
- Le type de retour est `AlgebraMatrix = Matrix<Dim, Dim>` (taille de la représentation matricielle du groupe, pas du
  vecteur tangent) : 3×3 pour SO3 et SE2, 4×4 pour SE3, et pour SE23 (9×9) et SGal3 (10×10) le bloc 5×5 de la vraie
  représentation est placé dans le coin supérieur gauche d'une matrice de la taille `Dim`, pour que les commutateurs de la
  base (BCH, `ad` générique) fonctionnent.
- `vee` / `computeVee` sont leur inverse.
- Il existe d'autres « hat » qui font la même chose pour SO(3) : `SO3::buildAntisymmetric`, `Types::skew3` / `unskew3`,
  `Twist::hat()` (4×4 de se(3)), la fonction libre `dualMatrix`. Quatre noms pour [·]× sur ℝ³ : à unifier (garder
  `SO3::hat`).

## Non vérifié

- Pas de build SOFA ni de lien de l'exécutable de tests (le dossier `cmake-build-relwithdebinfo` n'a pas été touché).
- Tests nécessitant SOFA : `AdvancedFeatures`, `HookeSerat_*`, `Phase2`, `MultiSectionBeam`, `jacobian_verification`
  (besoin de `Cosserat/config.h`), `LieGroupIntegrationTest` (besoin de `Cosserat/liegroups/Bundle.h`). Ils utilisent
  SE3/SO3, que j'ai modifiés.
- Passe `-fsyntax-only` sur les 27 unités de compilation du plugin (flags de `compile_commands.json` réécrits vers
  le worktree) : 27/27 code 0, aucune erreur ni nouvel avertissement dans `liegroups`. Ne remplace pas une compilation
  complète ni l'exécution des scènes.
- Float : seul SE3 en float a été testé numériquement.
