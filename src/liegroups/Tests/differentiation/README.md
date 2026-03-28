# Tests de Différentiabilité - Liegroups

Ce dossier contient l'infrastructure complète pour tester et valider la différentiabilité de la librairie `liegroups`.

## Contenu

### Fichiers

1. **`DifferentiationTestUtils.h`** - Utilitaires de test
   - Différences finies (forward et central)
   - Calcul de jacobiens numériques
   - Comparaison de matrices avec erreur relative
   - Support pour matrices statiques et dynamiques

2. **`test_finite_differences.cpp`** - Tests de base
   - Validation de la précision des différences finies
   - Tests exp/log pour SO2, SO3, SE2, SE3
   - Tests d'action de groupes

3. **`test_analytical_jacobians.cpp`** - Tests des jacobiens analytiques
   - Validation des jacobiens de composition
   - Validation des jacobiens d'inverse
   - Validation des jacobiens d'action
   - Tests de cohérence

## Utilisation

### Compiler les Tests

Les tests utilisent Google Test. Pour compiler (depuis la racine du projet) :

```bash
build_clion
```

### Exécuter les Tests

```bash
./build/bin/Cosserat_liegroups_tests
```

Ou pour des tests spécifiques :

```bash
# Tests de différences finies uniquement
./build/bin/test_finite_differences

# Tests de jacobiens analytiques uniquement  
./build/bin/test_analytical_jacobians
```

## Jacobiens Implémentés

### SO3 (Rotations 3D)

| Méthode | Signature | Description |
|---------|-----------|-------------|
| `composeJacobians(S)` | `pair<Matrix3, Matrix3>` | ∂(R\*S)/∂R et ∂(R\*S)/∂S |
| `inverseJacobian()` | `Matrix3` | ∂(R⁻¹)/∂R |
| `actionJacobians(p)` | `pair<Matrix3, Matrix3>` | ∂(R\*p)/∂R et ∂(R\*p)/∂p |
| `actionJacobianRotation(p)` | `Matrix3` | ∂(R\*p)/∂R uniquement |
| `actionJacobianPoint(p)` | `Matrix3` | ∂(R\*p)/∂p uniquement |

### SE3 (Transformations Rigides 3D)

| Méthode | Signature | Description |
|---------|-----------|-------------|
| `composeJacobians(h)` | `pair<Matrix6, Matrix6>` | ∂(g\*h)/∂g et ∂(g\*h)/∂h |
| `inverseJacobian()` | `Matrix6` | ∂(g⁻¹)/∂g |
| `actionJacobians(p)` | `pair<Matrix3x6, Matrix3>` | ∂(g\*p)/∂g et ∂(g\*p)/∂p |
| `actionJacobianGroup(p)` | `Matrix3x6` | ∂(g\*p)/∂g uniquement |
| `actionJacobianPoint(p)` | `Matrix3` | ∂(g\*p)/∂p uniquement |

## Exemples d'Utilisation

### Tester un Jacobien

```cpp
#include "DifferentiationTestUtils.h"
#include "../../SO3.h"

using namespace sofa::component::cosserat::liegroups;
using TestUtils = testing::DifferentiationTestUtils<double>;

// Créer une rotation
SO3d R = SO3d::exp(Vector3(0.5, 0.3, -0.2));

// Obtenir le jacobien analytique
auto J_analytical = R.inverseJacobian();

// Fonction à tester
auto inverse_func = [&R](const Vector3& delta) {
    return (SO3d::exp(delta) * R).inverse().log();
};

// Calculer le jacobien numérique
Vector3 zero = Vector3::Zero();
auto J_numerical = TestUtils::centralDifferenceJacobian<3, 3>(
    inverse_func, zero
);

// Comparer
bool passed = TestUtils::compareMatrices<3, 3>(
    J_analytical, J_numerical, 1e-5, true  // verbose
);
```

### Utiliser pour l'Optimisation

```cpp
// Problème: minimiser ||g.translation() - target||²
SE3d g = SE3d::Identity();
Vector3 target(1.0, 0.0, 0.0);
Vector6 strain = /* paramètres à optimiser */;

// Fonction coût
auto cost = [](const Vector6& xi) {
    SE3d g = SE3d::exp(xi);
    return (g.translation() - target).squaredNorm();
};

// Gradient numérique
Vector6 gradient = TestUtils::centralDifferenceGradient<6>(cost, strain);

// Utiliser pour descente de gradient
strain -= learning_rate * gradient;
```

## Conventions Mathématiques

### Perturbation à Gauche

Tous les jacobiens utilisent la convention de **perturbation à gauche** :
```
g_δ = exp(δ) * g
```

Cela signifie que pour une variation `δ` dans l'espace tangent :
```
J = ∂f(g_δ)/∂δ |_{δ=0}
```

### Espace Tangent

Les jacobiens sont exprimés dans l'**espace tangent** (algèbre de Lie), pas dans l'espace ambiant.

Pour SO3 : tangent ∈ so(3) ≅ ℝ³  
Pour SE3 : tangent ∈ se(3) ≅ ℝ⁶

### Formules Clés

**Composition SO3** :
- ∂(R\*S)/∂R = S^T (adjoint de S⁻¹)
- ∂(R\*S)/∂S = I

**Inverse SO3** :
- ∂(R⁻¹)/∂R = -R^T

**Action SO3** :
- ∂(R\*p)/∂R = -[R\*p]× (négatif du skew-symmetric de R\*p)
- ∂(R\*p)/∂p = R

**Composition SE3** :
- ∂(g\*h)/∂g = Ad_{h⁻¹}
- ∂(g\*h)/∂h = I

**Action SE3** :
- ∂(g\*p)/∂g = [R | -[R\*p]×] (matrice 3×6)
- ∂(g\*p)/∂p = R

## Tests de Validation

Chaque jacobien analytique est validé contre des différences finies :

1. **Précision** : Erreur relative < 10⁻⁵ (ou 10⁻⁴ pour cas relaxés)
2. **Robustesse** : Testé sur plusieurs configurations
3. **Cohérence** : Propriétés mathématiques vérifiées (e.g., (R⁻¹)⁻¹ = R)

### Résultats Attendus

Tous les tests doivent **PASS** ✓ avec les tolérances par défaut.

Si un test échoue :
- Vérifier l'implémentation analytique
- Augmenter le pas de différences finies si près d'une singularité
- Utiliser `verbose=true` pour diagnostiquer

## Prochaines Étapes

- [ ] Ajouter jacobiens pour SO2 et SE2
- [ ] Implémenter dexp() et dexpInv() (différentielle de l'exponentielle)
- [ ] Support pour autodiff (optionnel)
- [ ] Benchmarks de performance
- [ ] Exemples d'optimisation de trajectoires

## Références

1. **Solà et al. (2018)** - "A micro Lie theory for state estimation in robotics"
2. **Eade (2017)** - "Lie Groups for Computer Vision"
3. Fichier `../../DIFFERENTIATION.md` - Guide complet d'utilisation
4. Fichier `../../IMPLEMENTATION_PLAN.md` - Plan d'implémentation

---

**Dernière mise à jour** : Décembre 2025  
**Statut** : Phase 2 - 60% Complete ✅
