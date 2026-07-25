# TESTING — Cosserat plugin

> Comment compiler et lancer la suite de tests C++ et Python du plugin.
> Auteur : 2026-05-31. 
> Branche : `painless/python-explicit`.

---

## 1. Chemins de référence

| Variable | Valeur |
|---|---|
| Source du plugin | `/Users/yadagolo/travail/plugin/plugin.Cosserat` |
| Build SOFA | `/Users/yadagolo/travail/sofa/src/cmake-build-relwithdebinfo` |
| CMake CLion | `/Users/yadagolo/Applications/CLion.app/Contents/bin/cmake/mac/aarch64/bin/cmake` |
| Binaire de tests | `<BUILD>/bin/Cosserat_tests` |

Pour alléger, on note `BUILD = /Users/yadagolo/travail/sofa/src/cmake-build-relwithdebinfo`.

---

## 2. Activer la compilation des tests (une fois)

Les tests sont conditionnés par `COSSERAT_BUILD_TESTS=ON` (défaut : ON si `SOFA_BUILD_TESTS` est ON).
Si tu obtiens « target `Cosserat_tests` not found » :

```bash
cmake -DCOSSERAT_BUILD_TESTS=ON \
      -S /Users/yadagolo/travail/plugin/plugin.Cosserat \
      -B $BUILD
```

---

## 3. Compiler les tests

```bash
# Compile uniquement l'exécutable de tests (le plus rapide)
cmake --build $BUILD --target Cosserat_tests -j8

# Si tu changes aussi le plugin lui-même (src/Cosserat/*) :
cmake --build $BUILD --target Cosserat -j8
cmake --build $BUILD --target Cosserat_tests -j8
```

En cas d'erreurs de link étranges, recompile en série pour avoir un log propre :

```bash
cmake --build $BUILD --target Cosserat_tests -j1
```

---

## 4. Lancer la suite complète

```bash
cd $BUILD
./bin/Cosserat_tests
```

Ou via CTest depuis le build dir :

```bash
ctest -R Cosserat_tests -V
```

---

## 5. Lancer les tests round-trip log/exp (le nouveau fichier)

Le fichier `tests/unit/LogExpRoundTripTest.cpp` contient 15 tests groupés en deux suites
GTest : `SO3LogExpRoundTrip` et `SE3LogExpRoundTrip`.

```bash
# Toute la suite round-trip (SO3 + SE3)
$BUILD/bin/Cosserat_tests --gtest_filter="*LogExpRoundTrip*"

# Uniquement SO3 ou uniquement SE3
$BUILD/bin/Cosserat_tests --gtest_filter="SO3LogExpRoundTrip.*"
$BUILD/bin/Cosserat_tests --gtest_filter="SE3LogExpRoundTrip.*"

# Le test « A = exp(log(A)) » à la matrice
$BUILD/bin/Cosserat_tests --gtest_filter="*MatrixLevelIdentity*"

# Le sentinelle du bug SE3 corrigé le 2026-05-31 (V/V_inv mauvaise convention K)
$BUILD/bin/Cosserat_tests --gtest_filter="*BUG_SENTINEL*"

# Stress aléatoire (100 ξ tirés au sort)
$BUILD/bin/Cosserat_tests --gtest_filter="*RandomStress100"

# Tests sur la branche Taylor (θ < 1e-4)
$BUILD/bin/Cosserat_tests --gtest_filter="*Taylor*"

# Tests autour de la singularité θ → π
$BUILD/bin/Cosserat_tests --gtest_filter="*NearPi*"
```

### Options gtest utiles

| Option | Effet |
|---|---|
| `--gtest_list_tests` | Liste les tests sans les exécuter |
| `--gtest_print_time=1` | Affiche le temps de chaque test |
| `--gtest_repeat=N` | Re-lance N fois (utile pour les tests aléatoires) |
| `--gtest_shuffle` | Mélange l'ordre d'exécution |
| `--gtest_break_on_failure` | S'arrête sur le premier échec |
| `--gtest_color=yes` | Force la couleur |
| `--gtest_output=xml:rep.xml` | Export rapport JUnit XML |

Exemples combinés :

```bash
# Lister tous les tests round-trip enregistrés
$BUILD/bin/Cosserat_tests --gtest_list_tests | grep -A1 LogExpRoundTrip

# Stress aléatoire répété 10× pour chercher un cas pathologique
$BUILD/bin/Cosserat_tests --gtest_filter="*RandomStress100" --gtest_repeat=10

# Suite complète avec timing et rapport XML
$BUILD/bin/Cosserat_tests --gtest_print_time=1 --gtest_output=xml:tests_report.xml
```

---

## 6. Autres suites de tests existantes

| Suite | Filtre gtest |
|---|---|
| Strain → Rigid mapping | `--gtest_filter="Strain2RigidCosseratMappingTest.*"` |
| Discrete Cosserat applyDJT | `--gtest_filter="DiscreteCosseratMappingApplyDJTTest.*"` |
| Intégration liegroups (Twist, Bezier, Kalman, iLQR…) | `--gtest_filter="TwistWrenchTest.*:BezierSE3Test.*:BeamStateEstimatorTest.*"` |
| Tout SE3 / SO3 / LieGroup* | `--gtest_filter="SE3*:SO3*:LieGroup*"` |

```bash
$BUILD/bin/Cosserat_tests --gtest_filter="Strain2RigidCosseratMappingTest.*"
$BUILD/bin/Cosserat_tests --gtest_filter="DiscreteCosseratMappingApplyDJTTest.*"
timeout 60 $BUILD/bin/Cosserat_tests --gtest_filter="SE3*:SO3*:LieGroup*"
```

---

## 7. Debug avec lldb (macOS)

Quand un test plante (segfault, assertion non-gtest) :

```bash
# Run et backtrace sur la première erreur
lldb $BUILD/bin/Cosserat_tests \
  -o "run --gtest_filter=*BUG_SENTINEL*" \
  -o "bt 20" \
  -o "quit"

# Ouvre lldb en interactif, charge le binaire avec filtre
lldb $BUILD/bin/Cosserat_tests
(lldb) run --gtest_filter="*MatrixLevelIdentity*"
(lldb) bt all
(lldb) frame variable
```

Pour gdb (Linux) :

```bash
gdb --args $BUILD/bin/Cosserat_tests --gtest_filter="*BUG_SENTINEL*"
(gdb) run
(gdb) bt 20
```

---

## 8. Tests Python (scènes SOFA)

> Les 5 scènes staggered importent désormais leurs helpers depuis
> `examples/python3/examples/staggered/_common.py` (raideur, masse, inertie,
> SO3, `add_painless_beam`). Toute modification d'une formule de raideur ne se
> fait qu'à un seul endroit (`_common.py:circular_stiffness`).

Les scènes Python ne sont pas dans `Cosserat_tests` — elles se lancent via `runSofa` :

```bash
cd /Users/yadagolo/travail/plugin/plugin.Cosserat

# Validation Euler-Bernoulli (régime petit déplacement)
runSofa -g qt examples/python3/examples/staggered/sofa_staggered_validation.py
runSofa -g qt examples/python3/examples/staggered/sofa_staggered_validation.py --argv "N=16"
runSofa -g qt examples/python3/examples/staggered/sofa_staggered_validation.py --argv "N=32"

# Cantilever full (positions + SO3, grande déformation)
runSofa -g qt examples/python3/examples/staggered/staggered_cantilever_full.py

# Validation torsion (GJ)
runSofa -g qt examples/python3/examples/staggered/staggered_torsion_validation.py

# Cantilever gravité (positions seules)
runSofa -g qt examples/python3/examples/staggered/staggered_cantilever_gravity.py

# Large déformation (arc, hélice, torsion prescrite)
runSofa -g qt examples/python3/examples/staggered/staggered_large_deformation.py

# Test géométrie pure (frames Rigid3d)
runSofa -g qt examples/python3/examples/staggered/staggered_geometry_test.py
```

Mode batch (headless, sans GUI) :

```bash
runSofa --batch --nbIterations 10000 \
        examples/python3/examples/staggered/sofa_staggered_validation.py
```

---

## 9. Tests des bindings Python

Si configurés avec `-DCOSSERAT_WITH_PYTHON_BINDINGS=ON` :

```bash
python3 /Users/yadagolo/travail/plugin/plugin.Cosserat/tests/unit/test_cosserat_bindings.py
```

Ou via CTest :

```bash
ctest -R CosseratPythonBindings -V
```

---

## 10. Output attendu (référence)

### Suite round-trip log/exp (15 tests)

```
[==========] Running 15 tests from 2 test suites.
[----------] Global test environment set-up.
[----------] 5 tests from SO3LogExpRoundTrip
[ RUN      ] SO3LogExpRoundTrip.Identity
[       OK ] SO3LogExpRoundTrip.Identity
[ RUN      ] SO3LogExpRoundTrip.SmallRotationsTaylorBranch
[       OK ] SO3LogExpRoundTrip.SmallRotationsTaylorBranch
[ RUN      ] SO3LogExpRoundTrip.GeneralBranch
[       OK ] SO3LogExpRoundTrip.GeneralBranch
[ RUN      ] SO3LogExpRoundTrip.NearPiSingularity
[       OK ] SO3LogExpRoundTrip.NearPiSingularity
[ RUN      ] SO3LogExpRoundTrip.RandomStress100
[       OK ] SO3LogExpRoundTrip.RandomStress100
[----------] 10 tests from SE3LogExpRoundTrip
[ RUN      ] SE3LogExpRoundTrip.Identity
[       OK ] SE3LogExpRoundTrip.Identity
[ RUN      ] SE3LogExpRoundTrip.PureRotationNoTranslation
[       OK ] SE3LogExpRoundTrip.PureRotationNoTranslation
[ RUN      ] SE3LogExpRoundTrip.PureTranslationNoRotation
[       OK ] SE3LogExpRoundTrip.PureTranslationNoRotation
[ RUN      ] SE3LogExpRoundTrip.MixedTwistTheta05_BUG_SENTINEL
[       OK ] SE3LogExpRoundTrip.MixedTwistTheta05_BUG_SENTINEL
[ RUN      ] SE3LogExpRoundTrip.MixedTwistTheta1
[       OK ] SE3LogExpRoundTrip.MixedTwistTheta1
[ RUN      ] SE3LogExpRoundTrip.MixedTwistTaylorBranch
[       OK ] SE3LogExpRoundTrip.MixedTwistTaylorBranch
[ RUN      ] SE3LogExpRoundTrip.BranchBoundaryAgreement
[       OK ] SE3LogExpRoundTrip.BranchBoundaryAgreement
[ RUN      ] SE3LogExpRoundTrip.NearPiSingularity
[       OK ] SE3LogExpRoundTrip.NearPiSingularity
[ RUN      ] SE3LogExpRoundTrip.RandomStress100
[       OK ] SE3LogExpRoundTrip.RandomStress100
[ RUN      ] SE3LogExpRoundTrip.MatrixLevelIdentity_A_equals_exp_log_A
[       OK ] SE3LogExpRoundTrip.MatrixLevelIdentity_A_equals_exp_log_A
[----------] Global test environment tear-down
[==========] 15 tests from 2 test suites ran.
[  PASSED  ] 15 tests.
```

---

## 11. Pièges fréquents

- ❌ **« target Cosserat_tests not found »** → `cmake -DCOSSERAT_BUILD_TESTS=ON` puis relancer la compilation.
- ❌ **« No such file `bin/Cosserat_tests` »** → la compilation n'a pas abouti, regarder le log. Souvent un fichier source manquant dans `tests/CMakeLists.txt`.
- ❌ **`MixedTwistTheta05_BUG_SENTINEL` échoue** → le fix SE3 du 2026-05-31 n'a pas été appliqué ou pas recompilé. `cmake --build … --target Cosserat -j8` puis re-tester.
- ❌ **`runSofa` ne trouve pas le plugin** → vérifier `SOFA_PLUGIN_PATH` ou que `PluginManager` charge bien le `Cosserat` plugin.
- ❌ **Tests Python `AttributeError` sur `.getPositions()`** → utiliser `.positions.value` (cf. `MEMOIRE.md §4` et `AUDIT_2026-05-31.md`).

---

## 12. Workflow recommandé après modification du code C++

```bash
# 1. Recompiler le plugin + les tests
cmake --build $BUILD --target Cosserat -j8 && \
cmake --build $BUILD --target Cosserat_tests -j8

# 2. Lancer les tests les plus rapides (sanity check) — quelques secondes
$BUILD/bin/Cosserat_tests --gtest_filter="*LogExpRoundTrip*"

# 3. Si OK, lancer la suite complète
$BUILD/bin/Cosserat_tests

# 4. Si OK, lancer une scène staggered pour vérifier la pile globale
runSofa -g qt examples/python3/examples/staggered/staggered_cantilever_full.py
```

---

## 13. Référence rapide (cheat-sheet)

```bash
# Tout en une commande (recompile + suite round-trip)
cmake --build $BUILD --target Cosserat_tests -j8 && \
  $BUILD/bin/Cosserat_tests --gtest_filter="*LogExpRoundTrip*" --gtest_print_time=1

# Le test « A = exp(log(A)) » seul
$BUILD/bin/Cosserat_tests --gtest_filter="*MatrixLevelIdentity*"

# Stress aléatoire de validation
$BUILD/bin/Cosserat_tests --gtest_filter="*RandomStress100" --gtest_repeat=10
```
