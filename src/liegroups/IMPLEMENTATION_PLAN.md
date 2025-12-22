# Plan d'Implémentation - Librairie Différentiable

## Phase 1 : Préparation de l'Infrastructure ✅ TERMINÉE

### Objectifs
1. Créer une structure de tests dédiée à la différentiabilité
2. Ajouter le support optionnel pour autodiff
3. Préparer les utilitaires pour les tests de jacobiens

### Tâches

#### 1.1 Structure de Tests ✅
- [x] Créer `Tests/differentiation/` pour les tests de différentiabilité
- [x] Créer `DifferentiationTestUtils.h` avec utilitaires de test
- [x] Créer `test_finite_differences.cpp` pour validation numérique
- [ ] Créer `test_autodiff_integration.cpp` (si autodiff activé)

#### 1.2 Support Autodiff (Optionnel)
- [ ] Créer `AutodiffSupport.h` pour détecter et supporter autodiff
- [ ] Ajouter option CMake `COSSERAT_WITH_AUTODIFF`
- [ ] Créer des exemples d'utilisation avec autodiff

#### 1.3 Utilitaires de Différentiation
- [ ] Créer `Differentiation.h` avec :
  - Calcul de différences finies
  - Vérification de jacobiens
  - Comparaison de gradients
  - Tests de cohérence

#### 1.4 Documentation
- [ ] Ajouter `DIFFERENTIATION.md` expliquant l'usage
- [ ] Exemples d'optimisation de trajectoires
- [ ] Guide d'intégration avec autodiff

---

## Phase 2 : Implémentation des Jacobiens ✅ EN COURS

### Tâches Planifiées
- [x] Ajouter `composeJacobians()` pour SO3 et SE3
- [x] Ajouter `inverseJacobian()` pour SO3 et SE3
- [x] Implémenter `actionJacobians()` complet pour SO3 et SE3
- [x] Créer `test_analytical_jacobians.cpp` avec tests exhaustifs
- [ ] Ajouter jacobiens pour SO2 et SE2
- [ ] Valider tous les tests

---

## Phase 3 : Validation et Optimisation (À venir)

### Tâches Planifiées
- [ ] Benchmarks de performance
- [ ] Tests de précision numérique
- [ ] Documentation complète
- [ ] Exemples d'applications

---

## Statut Actuel

**Branche** : `feature/differentiable-liegroups`  
**Phase Actuelle** : Phase 2 - Implémentation Jacobiens  
**Progression** : 60%

### Accompli Récemment
- ✅ Infrastructure complète de tests
- ✅ Jacobiens analytiques pour SO3 (compose, inverse, action)
- ✅ Jacobiens analytiques pour SE3 (compose, inverse, action)
- ✅ Suite de tests compréhensive avec validation numérique

### Prochaines Actions Immédiates
1. Compiler et valider les tests
2. Ajouter jacobiens pour SO2 et SE2 (optionnel)
3. Documenter les nouvelles API
