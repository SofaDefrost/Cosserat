# Plan de différentiabilité pour la bibliothèque Lie Groups

## Analyse de différentiabilité actuelle

Après analyse approfondie du code source dans `src/liegroups/`, la bibliothèque **n'est pas entièrement différentiable** pour l'automatic differentiation (AD). Voici les problèmes identifiés :

### Problèmes majeurs pour la différentiabilité

1. **Branching conditionnel sur les données d'entrée** :
   - Les fonctions exponentielle (`exp`) et logarithmique (`log`) utilisent des instructions `if` basées sur des seuils comme `norm < epsilon`
   - Exemples :
     - `SO3::expImpl()`: `if (theta < Types<Scalar>::epsilon())`
     - `SE3::computeExp()`: `if (angle < Types<Scalar>::epsilon())`
     - `SE3::computeLog()`: `if (angle < Types<Scalar>::epsilon())`
   - Ces branches créent des discontinuités dans les gradients lorsque la norme d'entrée franchit le seuil

2. **Fonctions non différentiables** :
   - `std::abs()` utilisé dans les validations (`SO2`, `Types.h`)
   - `std::atan2()` dans `SO3::distance()` - discontinu aux points critiques (0,0) et le long des coupures de branche
   - Ces fonctions introduisent des points non-lisses où les gradients sont indéfinis ou infinis

3. **Mécanismes de secours dans les opérations** :
   - Méthodes comme `distance()` et `slerp()` incluent des blocs try-catch et des fallbacks conditionnels
   - Ces mécanismes peuvent mener à un comportement discontinu pendant l'optimisation

4. **Problèmes potentiels avec exceptions et limites** :
   - La formule BCH lance `NumericalInstabilityException` si les entrées sont trop grandes
   - Utilisation de `std::numeric_limits<Scalar>::max()` ou `quiet_NaN()` comme fallbacks

### Points forts pour la stabilité numérique

- Utilisation de seuils et d'approximations pour petits angles pour améliorer la stabilité
- Fonctions stables comme `sinc`, `cosc`, et `sinc3` dans `Types.h`
- Expansions en série de Taylor pour petits angles

### État actuel de compatibilité AD

La bibliothèque n'est **pas compatible AD** dans sa forme actuelle. Elle nécessiterait une refactorisation importante pour être utilisée dans des frameworks d'optimisation basés sur gradients.

## Plan étape par étape pour rendre la bibliothèque différentiable

### Phase 1 : Préparation et analyse (1-2 jours)

#### Étape 1.1 : Créer des tests de différentiabilité
- **Objectif** : Établir une base de tests pour vérifier la différentiabilité
- **Tâches** :
  - Intégrer un framework AD simple (CppAD, autodiff, ou CppNumericalSolvers)
  - Créer des tests unitaires vérifiant la correction des gradients
  - Tester les fonctions de base : `exp`, `log`, `compose`, `inverse`
- **Critères** : Tests passent avec scalaires AD et comparaisons avec différences finies

#### Étape 1.2 : Documenter les fonctions problématiques
- **Objectif** : Identifier et prioriser les points à refactoriser
- **Tâches** :
  - Lister toutes les fonctions avec branching ou fonctions non-différentiables
  - Classer par fréquence d'usage (exp/log prioritaires)
  - Créer une matrice de compatibilité AD pour chaque groupe
- **Délivrable** : Document `ad_compatibility_analysis.md`

### Phase 2 : Refactorisation des fonctions de base (3-5 jours)

#### Étape 2.1 : Remplacer les branches epsilon par des transitions smooth
- **Objectif** : Éliminer les discontinuités des seuils
- **Méthode** : Utiliser des fonctions de transition smooth
- **Exemples d'implémentation** :
  ```cpp
  // Au lieu de : if (norm < eps) { approx } else { exact }
  // Utiliser :
  Scalar smooth_factor = 1.0 / (1.0 + exp(-k * (norm - eps)));
  result = smooth_factor * approx + (1.0 - smooth_factor) * exact;
  ```
- **Fonctions à modifier** :
  - `SO3::exp()`, `SO3::log()`
  - `SE3::computeExp()`, `SE3::computeLog()`
  - `SE2::exp()`, `SE2::log()`
  - `SE23::computeExp()`, `SE23::computeLog()`

#### Étape 2.2 : Éliminer les fonctions non-différentiables
- **Objectif** : Remplacer abs et atan2 par approximations différentiables
- **Solutions** :
  - Remplacer `std::abs(x)` par `x * tanh(k * x)` (k grand pour approximation)
  - Remplacer `std::atan2(y, x)` par une implémentation smooth :
    ```cpp
    Scalar smooth_atan2(Scalar y, Scalar x) {
        Scalar r = sqrt(x*x + y*y);
        Scalar angle = acos(x / (r + 1e-8));  // Éviter division par zéro
        return (y >= 0) ? angle : -angle;
    }
    ```
- **Supprimer les exceptions** : Remplacer par des pénalités smooth ou clamps

#### Étape 2.3 : Refactoriser les conversions matrice-quaternion
- **Problème** : `SO3::matrixToQuaternion()` utilise du branching complexe
- **Solution** : Implémenter une version avec poids continus :
  ```cpp
  // Au lieu de if-else basé sur trace, utiliser des poids smooth
  Scalar w1 = smooth_max(0, trace + 1) / 4;
  Scalar w2 = smooth_max(0, 1 + R(0,0) - R(1,1) - R(2,2)) / 4;
  // etc.
  ```

### Phase 3 : Optimisation et tests (2-3 jours)

#### Étape 3.1 : Optimiser les performances
- **Challenge** : Les transitions smooth peuvent être plus coûteuses
- **Solutions** :
  - Utiliser des approximations polynomiales pour petits angles
  - Pré-calculer les transitions quand possible
  - Maintenir les seuils epsilon mais de manière différentiable

#### Étape 3.2 : Tests exhaustifs
- **Couverture** :
  - Gradients corrects pour tous les cas (petits angles, grands angles, cas limites)
  - Comparaison précision numérique vs différentiabilité
  - Tests de stabilité numérique
- **Performance** : Mesurer l'overhead AD (objectif : max 2-3x)

#### Étape 3.3 : Compatibilité frameworks AD
- **Frameworks cibles** :
  - CppAD (C++)
  - Stan Math Library (C++)
  - PyTorch C++ extensions
  - JAX (via bindings)
- **Validation** : S'assurer que le templating fonctionne avec différents types scalaires AD

### Phase 4 : Documentation et migration (1-2 jours)

#### Étape 4.1 : Mettre à jour la documentation
- **Contenu** :
  - Section dédiée à la compatibilité AD dans `README.md`
  - Guide d'usage avec frameworks AD
  - Limitations restantes documentées
  - Exemples d'intégration AD

#### Étape 4.2 : Migration graduelle
- **Stratégie** :
  - Créer des versions "AD-safe" des fonctions critiques
  - Maintenir la compatibilité backward avec macros/flags
  - Flag de compilation : `COSSAERAT_AD_SAFE` vs `COSSAERAT_PERFORMANCE`

### Alternatives plus simples

#### Option A : Wrappers AD-compatible (recommandée si refactorisation complète trop coûteuse)
- Créer des wrappers qui utilisent uniquement les parties différentiables
- Contourner les fonctions problématiques dans le code utilisateur
- Avantages : Développement rapide, compatibilité backward

#### Option B : Version spécialisée
- Fork des fonctions critiques pour mode AD
- Utiliser des directives de précompilation `#ifdef COSSAERAT_AD_MODE`
- Avantages : Performance optimale dans les deux modes

### Estimation des ressources

- **Développement** : 1-2 développeurs pendant 1-2 semaines
- **Tests** : Environnement de test AD configuré
- **Validation** : Tests de non-régression complets

### Critères de succès

1. **Correctitude** : Tous les tests AD passent (gradients corrects ± tolérance)
2. **Performance** : Overhead AD acceptable (max 2-3x vs version non-AD)
3. **Compatibilité** : Fonctionne avec au moins 2 frameworks AD majeurs
4. **Stabilité** : Maintenir la précision numérique existante
5. **Maintenabilité** : Code lisible et documenté

### Risques et mitigation

- **Performance** : Monitorer l'impact sur les simulations existantes
- **Précision** : Tests de régression numériques approfondis
- **Complexité** : Revue de code et documentation détaillée
- **Compatibilité** : Tests avec différentes versions de Eigen et compilateurs

### Métriques de suivi

- Nombre de fonctions refactorisées
- Couverture des tests AD (%)
- Performance relative (benchmarks)
- Taille du code ajouté/modifié

Cette refactorisation permettra d'utiliser la bibliothèque dans des contextes d'optimisation gradient-based et d'apprentissage automatique, tout en préservant ses qualités numériques actuelles pour les simulations traditionnelles.