# Lie Theory Improvements Based on Solà et al. (2021)

**Branch:** `feature/liegroups-sola-improvements`  
**Reference:** "A micro Lie theory for state estimation in robotics" (Solà et al., 2021, arXiv:1812.01537)

## 🎯 Objectif

Enrichir la bibliothèque liegroups avec les concepts modernes de la théorie de Lie appliquée à l'estimation d'état en robotique, rendant le code plus intuitif et conforme à la littérature récente.

## ✅ Implémentations Complétées (Phase 1)

### 1. Opérateurs ⊕ et ⊖ (Plus/Minus)

**Fichiers modifiés:** `LieGroupBase.h`

#### Avant
```cpp
SE3d X_new = X.compose(SE3d::exp(delta));
TangentVector error = X_measured.inverse().compose(X_estimated).log();
```

#### Après
```cpp
SE3d X_new = X + delta;              // ou X.plus(delta)
TangentVector error = X_measured - X_estimated;  // ou X_measured.minus(X_estimated)
```

**Avantages :**
- ✅ Notation intuitive cohérente avec la littérature
- ✅ Code plus lisible
- ✅ Réduction des erreurs de programmation
- ✅ Facilite l'écriture d'algorithmes d'estimation

### 2. Jacobiennes Droites et Gauches (Jr, Jl)

**Fichiers modifiés:** `LieGroupBase.h`, `SO3.inl`

#### API Ajoutée

```cpp
// Jacobiennes droites (variations locales)
AdjointMatrix Jr = SO3d::rightJacobian(omega);
AdjointMatrix Jr_inv = SO3d::rightJacobianInverse(omega);

// Jacobiennes gauches (variations globales)
AdjointMatrix Jl = SO3d::leftJacobian(omega);
AdjointMatrix Jl_inv = SO3d::leftJacobianInverse(omega);
```

#### Formules Implémentées (SO(3))

**Right Jacobian** (Équation 143 du papier) :
```
Jr(θ) = I - (1-cos θ)/θ² [θ]× + (θ-sin θ)/θ³ [θ]²×
```

**Right Jacobian Inverse** (Équation 144 du papier) :
```
Jr⁻¹(θ) = I + ½[θ]× + (1/θ² - (1+cos θ)/(2θ sin θ))[θ]²×
```

**Relations :**
- `Jl(θ) = Jr(-θ)`
- `Jl⁻¹(θ) = Jr⁻¹(-θ)`

**Gestion des petits angles :**
- Pour `‖θ‖ < ε` : approximations au premier ordre
- `Jr(θ) ≈ I - ½[θ]×`
- `Jr⁻¹(θ) ≈ I + ½[θ]×`

## 📊 Cas d'Usage Implémentés

### ESKF (Error-State Kalman Filter)

L'exemple `example_sola_operators.cpp` montre l'implémentation complète d'une étape de prédiction ESKF :

```cpp
// État : X ∈ SO(3)
// Commande : u ∈ ℝ³
// Covariance : P ∈ ℝ³ˣ³

// Prédiction de l'état
SO3d X_pred = X + u;  // Notation intuitive !

// Prédiction de la covariance
auto F = SO3d::exp(u).inverse().adjoint();  // Ad_Exp(u)⁻¹
auto G = SO3d::rightJacobian(u);             // Jr(u)
Matrix3d P_pred = F * P * F.transpose() + G * Q * G.transpose();
```

### Propagation d'Incertitude

```cpp
// État avec incertitude
SO3d R_mean = SO3d::exp(omega);
Matrix3d Sigma = /* covariance locale */;

// Appliquer un incrément
SO3d R_new = R_mean + delta;

// Propager la covariance
auto F = SO3d::exp(delta).inverse().adjoint();
auto G = SO3d::rightJacobian(delta);
Matrix3d Sigma_new = F * Sigma * F.transpose() + G * Sigma_delta * G.transpose();
```

## 🧪 Exemples et Tests

### Fichier d'exemple complet
**Emplacement :** `examples/liegroups/example_sola_operators.cpp`

**Contenu :**
1. **Example 1** : Démonstration des opérateurs ⊕/⊖
2. **Example 2** : Calcul et vérification de Jr
3. **Example 3** : Relation entre Jl et Jr
4. **Example 4** : Propagation d'incertitude complète
5. **Example 5** : Gestion des petits angles
6. **Example 6** : Étape ESKF complète

### Compilation et exécution
```bash
cd build
make example_sola_operators
./examples/liegroups/example_sola_operators
```

## 📈 Impact sur les Performances

- ⚡ **Aucune pénalité de performance** : les nouvelles méthodes sont des wrappers inline
- 🎯 **Meilleure stabilité numérique** : approximations petits angles intégrées
- 📖 **Lisibilité améliorée** : code plus court et plus clair

## 🔜 Prochaines Étapes (voir AMELIORATIONS_PROPOSEES.md)

### Phase 2 : Jacobiens des Opérations (2 semaines)
- [ ] Jacobien de l'inversion : `J_X⁻¹_X = -Ad_X`
- [ ] Jacobien de la composition : `J_XY_X = Ad_Y⁻¹`
- [ ] Jacobiens de ⊕ et ⊖
- [ ] Jacobiens de l'action de groupe

### Phase 3 : Support d'Incertitude (2 semaines)
- [ ] Classe `GaussianOnManifold<LieGroupType>`
- [ ] Transformation local ↔ global frame
- [ ] Propagation automatique de covariance
- [ ] Exemples ESKF complets

### Phase 4 : SE(3) et au-delà (2 semaines)
- [ ] Implémentation des jacobiens pour SE(3)
- [ ] Amélioration des groupes composites (Bundle)
- [ ] BCH ordre supérieur
- [ ] Intégration avec composants SOFA

## 🔗 Références

### Papier Principal
```bibtex
@article{sola2021micro,
  title={A micro Lie theory for state estimation in robotics},
  author={Sol{\`a}, Joan and Deray, Jeremie and Atchuthan, Dinesh},
  journal={arXiv preprint arXiv:1812.01537},
  year={2021}
}
```

### Bibliothèque Associée
- **manif** : https://github.com/artivis/manif
- Bibliothèque C++ template-only implémentant les mêmes concepts
- Notre implémentation s'en inspire mais est adaptée aux besoins Cosserat

## 📝 Notes de Développement

### Conventions de Nommage
- **Right/Left** : explicite la frame de référence (local vs global)
- **plus/minus** : noms de méthodes explicites
- **operator+/-** : surcharges pour syntaxe intuitive

### Compatibilité
- ✅ Rétro-compatible : anciennes méthodes toujours disponibles
- ✅ `dexp()` et `dexpInv()` maintenant deprecated mais fonctionnels
- ✅ Redirection vers `rightJacobian()` pour clarté

### Tests
TODO : Ajouter tests unitaires
```cpp
TEST(JacobianTests, RightJacobianIdentity) {
    Vector3d omega(0.5, 0.3, 0.2);
    auto Jr = SO3d::rightJacobian(omega);
    auto Jr_inv = SO3d::rightJacobianInverse(omega);
    EXPECT_TRUE((Jr * Jr_inv).isApprox(Matrix3d::Identity(), 1e-10));
}

TEST(OperatorTests, PlusMinusConsistency) {
    SO3d X = SO3d::exp(Vector3d(0.3, 0.1, 0.2));
    Vector3d tau(0.1, 0.05, 0.02);
    SO3d Y = X + tau;
    Vector3d tau_recovered = Y - X;
    EXPECT_TRUE(tau.isApprox(tau_recovered, 1e-10));
}
```

## 🤝 Contribution

Pour continuer le développement :

1. Lire `AMELIORATIONS_PROPOSEES.md` pour la feuille de route complète
2. Suivre l'ordre de priorité indiqué
3. Tester avec l'exemple fourni
4. Ajouter des tests unitaires pour toute nouvelle fonctionnalité

## 📊 État du Projet

| Feature | Status | Priorité | Fichiers |
|---------|--------|----------|----------|
| Opérateurs ⊕/⊖ | ✅ Complet | HAUTE | `LieGroupBase.h` |
| Jr/Jl pour SO(3) | ✅ Complet | HAUTE | `SO3.inl` |
| Jr/Jl pour SE(3) | ⏳ À faire | HAUTE | `SE3.inl` |
| Jacobiens ops élémentaires | ⏳ À faire | MOYENNE | `LieGroupBase.h` |
| GaussianOnManifold | ⏳ À faire | MOYENNE | `Uncertainty.h` |
| Groupes composites | ⏳ À faire | BASSE | `Bundle.h` |

**Temps estimé pour complétion totale :** 6-8 semaines

---

**Auteur :** Basé sur l'analyse du papier de Solà et al. (2021)  
**Date :** Janvier 2025  
**Branch :** `feature/liegroups-sola-improvements`
