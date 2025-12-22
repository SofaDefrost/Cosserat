# Guide de Différentiabilité - Liegroups

## Vue d'Ensemble

Ce guide explique comment utiliser les capacités de différentiabilité de la librairie `liegroups` pour l'optimisation, le contrôle optimal, et l'apprentissage différentiable.

## Table des Matières

1. [Introduction](#introduction)
2. [Tests de Différences Finies](#tests-de-différences-finies)
3. [Jacobiens Disponibles](#jacobiens-disponibles)
4. [Exemples d'Utilisation](#exemples-dutilisation)
5. [Intégration avec Autodiff](#intégration-avec-autodiff)
6. [Bonnes Pratiques](#bonnes-pratiques)

---

## Introduction

Les groupes de Lie sont des variétés différentiables, ce qui signifie que toutes leurs opérations (exp, log, composition, action) sont différentiables. Cette librairie fournit :

- **Jacobiens analytiques** pour les opérations principales
- **Utilitaires de test** pour valider les implémentations
- **Support optionnel** pour la différentiation automatique

### Pourquoi la Différentiabilité ?

La différentiabilité permet :
- **Optimisation basée sur gradient** de trajectoires
- **Contrôle optimal** avec méthodes du premier/second ordre
- **Apprentissage différentiable** (ML/RL sur groupes de Lie)
- **Calibration automatique** de paramètres de modèles

---

## Tests de Différences Finies

### Utilitaire de Test

La classe `DifferentiationTestUtils` fournit des méthodes pour tester les jacobiens :

```cpp
#include "Tests/differentiation/DifferentiationTestUtils.h"

using namespace sofa::component::cosserat::liegroups::testing;
using TestUtils = DifferentiationTestUtils<double>;

// Test d'un jacobien analytique
bool test_passed = TestUtils::testJacobian<3, 3>(
    my_function,           // Fonction à tester
    analytical_jacobian,   // Jacobien analytique
    test_point,            // Point de test
    1e-5,                  // Tolérance
    true,                  // Utiliser différences centrales
    true                   // Mode verbose
);
```

### Méthodes Disponibles

#### 1. Gradient (R^n → R)

```cpp
auto gradient = TestUtils::centralDifferenceGradient<3>(
    [](const Vector3& x) { return x.squaredNorm(); },
    test_point
);
```

#### 2. Jacobien (R^n → R^m)

```cpp
auto jacobian = TestUtils::centralDifferenceJacobian<3, 6>(
    [](const Vector6& xi) { return SE3d::exp(xi).translation(); },
    test_tangent
);
```

#### 3. Comparaison de Matrices

```cpp
bool is_close = TestUtils::compareMatrices<6, 6>(
    analytical_jacobian,
    numerical_jacobian,
    1e-5,     // Tolérance
    true      // Verbose
);
```

---

## Jacobiens Disponibles

### SO3 (Rotations 3D)

```cpp
SO3d R = SO3d::exp(omega);

// Adjoint (déjà implémenté)
Matrix3 Ad = R.adjoint();  // Ad_R : so(3) → so(3)

// Action sur un point
Vector3 Rp = R.act(p);

// TODO: Jacobien de composition
// auto [J_left, J_right] = R.composeJacobians(S);

// TODO: Jacobien de l'inverse
// Matrix3 J_inv = R.inverseJacobian();
```

### SE3 (Transformations Rigides 3D)

```cpp
SE3d g = SE3d::exp(xi);

// Adjoint (déjà implémenté)
Matrix6 Ad = g.adjoint();  // Ad_g : se(3) → se(3)

// Action sur un point
Vector3 gp = g.act(p);

// TODO: Jacobiens de composition et inverse
```

### Cartes Exp/Log

```cpp
// Différentielle de l'exponentielle (déclarée, à implémenter)
Matrix6 dexp_xi = SE3d::dexp(xi);

// Différentielle de l'inverse de exp
Matrix6 dexpInv_xi = SE3d::dexpInv(xi);

// Différentielle du logarithme
Matrix6 dlog_g = g.dlog();
```

---

## Exemples d'Utilisation

### Exemple 1 : Optimisation de Trajectoire

```cpp
#include "SE3.h"
#include "Tests/differentiation/DifferentiationTestUtils.h"

using SE3d = SE3<double>;
using Vector6 = SE3d::TangentVector;
using TestUtils = DifferentiationTestUtils<double>;

// Objectif : Minimiser l'énergie d'une trajectoire de Cosserat
double trajectory_energy(const std::vector<Vector6>& strains, double length) {
    SE3d g = SE3d::Identity();
    double energy = 0.0;
    
    for (const auto& strain : strains) {
        g = g * SE3d::expCosserat(strain, length);
        energy += strain.squaredNorm();  // Terme de régularisation
    }
    
    // Terme d'erreur de position finale
    Vector3 target(1.0, 0.0, 0.0);
    energy += (g.translation() - target).squaredNorm();
    
    return energy;
}

// Gradient numérique pour optimisation
VectorXd compute_gradient(const VectorXd& params) {
    // Utiliser les utilitaires pour calculer le gradient
    // ...
}
```

### Exemple 2 : Validation de Jacobien

```cpp
#include "SO3.h"

// Test du jacobien d'action de SO3
void test_SO3_action_jacobian() {
    using SO3d = SO3<double>;
    using Vector3 = SO3d::Vector;
    using TestUtils = DifferentiationTestUtils<double>;
    
    Vector3 omega(0.3, 0.2, -0.1);
    SO3d R = SO3d::exp(omega);
    Vector3 p(1.0, 2.0, 3.0);
    
    // Fonction : action de R(delta) * R sur p
    auto action_perturbed = [&R, &p](const Vector3& delta) -> Vector3 {
        return (SO3d::exp(delta) * R).act(p);
    };
    
    // Jacobien numérique
    Vector3 zero = Vector3::Zero();
    auto J_numerical = TestUtils::centralDifferenceJacobian<3, 3>(
        action_perturbed, zero
    );
    
    // Jacobien analytique (à implémenter)
    Vector3 Rp = R.act(p);
    Matrix3 J_analytical = -SO3d::buildAntisymmetric(Rp);
    
    // Comparaison
    bool passed = TestUtils::compareMatrices<3, 3>(
        J_analytical, J_numerical, 1e-5, true
    );
    
    if (passed) {
        std::cout << "✓ Jacobian test passed!" << std::endl;
    } else {
        std::cout << "✗ Jacobian test failed!" << std::endl;
    }
}
```

### Exemple 3 : Contrôle Optimal

```cpp
// Problème de contrôle optimal : atteindre une cible
struct OptimalControlProblem {
    SE3d target;
    int n_steps;
    double dt;
    
    // État initial
    SE3d initial_state = SE3d::Identity();
    
    // Coût : distance au carré + régularisation de contrôle
    double cost(const std::vector<Vector6>& controls) {
        SE3d state = initial_state;
        double total_cost = 0.0;
        
        for (const auto& u : controls) {
            // Dynamique : intégration avec expo
            state = state * SE3d::exp(u * dt);
            
            // Régularisation du contrôle
            total_cost += u.squaredNorm();
        }
        
        // Coût terminal : distance à la cible
        Vector6 error = (state.inverse() * target).log();
        total_cost += 100.0 * error.squaredNorm();
        
        return total_cost;
    }
    
    // Gradient (à calculer avec différences finies ou autodiff)
    VectorXd gradient(const VectorXd& controls_flat);
};
```

---

## Intégration avec Autodiff

### Option 1 : autodiff (Recommandé)

```cpp
// À venir : Support pour autodiff
#ifdef COSSERAT_WITH_AUTODIFF
#include <autodiff/forward/dual.hpp>

using dual = autodiff::dual;
using SE3dual = SE3<dual>;

// Fonction avec calcul automatique du gradient
auto compute_with_gradient(const VectorXd& params) {
    VectorXdual params_dual = params.cast<dual>();
    
    // Votre fonction coût
    dual cost = trajectory_energy(params_dual);
    
    // Gradient automatique
    return autodiff::gradient(cost, params_dual);
}
#endif
```

### Option 2 : CppAD

```cpp
// Support CppAD à venir
```

---

## Bonnes Pratiques

### 1. Tester Systématiquement les Jacobiens

Toujours valider les jacobiens analytiques avec des différences finies :

```cpp
// Dans vos tests
TEST(MyLieGroupTest, JacobianValidation) {
    // Implémenter la fonction
    auto f = [](const VectorN& x) { return my_function(x); };
    
    // Calculer jacobien analytique
    MatrixMN J_analytical = my_analytical_jacobian(test_point);
    
    // Valider
    EXPECT_TRUE(TestUtils::testJacobian<M, N>(
        f, J_analytical, test_point, 1e-5, true
    ));
}
```

### 2. Choix de la Méthode

- **Différences centrales** : Plus précises mais 2× plus lentes
- **Différences forward** : Plus rapides mais moins précises
- **Autodiff** : Exact et rapide (quand disponible)

```cpp
// Pour les tests : utiliser central
auto J_test = TestUtils::centralDifferenceJacobian<M, N>(f, x);

// Pour l'optimisation : utiliser autodiff ou jacobiens analytiques
```

### 3. Gestion de la Stabilité Numérique

```cpp
// Tester près des singularités
std::vector<Vector3> critical_points = {
    Vector3::Zero(),           // Identité
    Vector3(M_PI, 0, 0),      // Rotation de 180°
    Vector3(0.001, 0, 0),     // Près de l'identité
};

for (const auto& omega : critical_points) {
    // Tester avec tolérance adaptée
    double tolerance = (omega.norm() < 0.01) ? 1e-4 : 1e-5;
    // ...
}
```

### 4. Optimisation de Performance

```cpp
// Réutiliser les calculs intermédiaires
struct CachedJacobian {
    SE3d g;
    Matrix6 Ad;  // Pré-calculé
    
    CachedJacobian(const Vector6& xi) : g(SE3d::exp(xi)) {
        Ad = g.adjoint();
    }
    
    Vector6 apply(const Vector6& v) const {
        return Ad * v;  // Réutilisation
    }
};
```

---

## Application : Optimisation de Trajectoires

### Utilisation de `CosseratTrajectoryOptimizer`

La classe `CosseratTrajectoryOptimizer` utilise les jacobiens analytiques pour optimiser des configurations de poutres Cosserat.

#### Exemple Simple

```cpp
#include "optimization/CosseratTrajectoryOptimizer.h"

using namespace sofa::component::cosserat::liegroups::optimization;

// Configuration
const int n_sections = 10;
const double section_length = 0.1;  // 10cm par section

// Strains initiaux (poutre droite)
std::vector<Vector6> initial_strains(n_sections, Vector6::Zero());

// Cible : 80cm en X, 20cm en Z
SE3d target = SE3d::Identity();
target.translation() = Vector3(0.8, 0.0, 0.2);

// Paramètres d'optimisation
CosseratTrajectoryOptimizer<double>::Parameters params;
params.learning_rate = 0.05;
params.max_iterations = 500;
params.tolerance = 1e-6;
params.regularization = 0.001;
params.use_line_search = true;
params.verbose = true;

// Optimiser
CosseratTrajectoryOptimizer<double> optimizer;
auto result = optimizer.optimizeToTarget(
    initial_strains, target, section_length, params
);

// Résultats
std::cout << "Convergé: " << result.converged << std::endl;
std::cout << "Erreur finale: " 
          << (result.final_transform.translation() - target.translation()).norm() 
          << " m" << std::endl;
```

#### Fonctionnalités

1. **Backpropagation automatique** : Les gradients sont calculés en propageant en arrière à travers la chaîne cinématique

2. **Line search adaptatif** : Recherche linéaire d'Armijo pour déterminer le pas optimal

3. **Régularisation** : Pénalisation L2 sur les strains pour éviter les solutions extrêmes

4. **Fonction coût personnalisée** :
```cpp
auto custom_cost = [](const std::vector<Vector6>& strains, 
                      std::vector<Vector6>& gradient) -> double {
    // Implémentez votre fonction coût ici
    // gradient doit être rempli
    return cost_value;
};

auto result = optimizer.optimizeCustom(initial_strains, custom_cost, params);
```

#### Exemple Complet

Voir `examples/simple_trajectory_optimization.cpp` pour un exemple complet avec :
- Affichage progressif
- Analyse de convergence
- Validation des résultats

Compiler avec :
```bash
cmake -DCOSSERAT_BUILD_EXAMPLES=ON ..
make simple_trajectory_optimization
./bin/examples/simple_trajectory_optimization
```

---

## Références

1. **Lie Groups for Computer Vision** - Eade (2017)
2. **A micro Lie theory for state estimation** - Solà et al. (2018)
3. **Numerical Recipes** - Press et al. (2007) - Chapitre sur différentiation numérique
4. **autodiff Documentation** - https://autodiff.github.io
5. **Optimization on Manifolds** - Absil, Mahony, Sepulchre (2008)

---

## Contribuer

Pour ajouter de nouveaux jacobiens :

1. Implémenter la méthode analytique dans la classe du groupe
2. Ajouter un test dans `Tests/differentiation/`
3. Valider avec différences finies
4. Documenter dans ce fichier

### Template de Test

```cpp
TEST_F(DifferentiationTest, MyNewJacobian) {
    // Setup
    MyGroup g = /* ... */;
    
    // Fonction à tester
    auto f = [&g](const TangentVector& delta) {
        return /* ... */;
    };
    
    // Jacobien analytique
    Matrix J_analytical = g.myNewJacobian();
    
    // Validation
    EXPECT_TRUE(TestUtils::testJacobian<M, N>(
        f, J_analytical, test_point, 1e-5, true, true
    ));
}
```

---

**Dernière mise à jour** : Décembre 2025  
**Contact** : Voir GitHub Issues du projet
