# Rapport d'Analyse : Différentiabilité de la Librairie `src/liegroups/`

**Date**: 22 décembre 2025  
**Projet**: plugin.Cosserat  
**Auteur**: Assistant AI - Analyse de code

---

## 1. Résumé Exécutif

La librairie `src/liegroups/` implémente une bibliothèque moderne de groupes de Lie en C++20 pour la simulation de poutres de Cosserat. L'analyse révèle que **la librairie possède déjà une architecture partiellement différentiable**, mais nécessite des extensions significatives pour devenir pleinement différentiable au sens de la différentiation automatique (AD).

### Verdict Principal
✅ **OUI, il est possible de rendre cette librairie différentiable**, avec un effort modéré à important selon le niveau de différentiabilité souhaité.

---

## 2. Architecture Actuelle

### 2.1 Structure de la Librairie

La librairie est organisée autour d'un patron CRTP (Curiously Recurring Template Pattern) avec :

- **Classe de base** : `LieGroupBase<Derived, Scalar, Dim, AlgebraDim, ActionDim>`
- **Groupes implémentés** :
  - `SO2<Scalar>` : Rotations 2D (dim algèbre = 1)
  - `SO3<Scalar>` : Rotations 3D (dim algèbre = 3)
  - `SE2<Scalar>` : Transformations rigides 2D (dim algèbre = 3)
  - `SE3<Scalar>` : Transformations rigides 3D (dim algèbre = 6)
  - `SE23<Scalar>` : Groupe semi-direct SE(2)⋉ℝ
  - `SGal3<Scalar>` : Groupe de Galilée spécial 3D
  - `RealSpace<Scalar, N>` : Espaces euclidiens
  - `Bundle<Groups...>` : Produit de groupes de Lie

### 2.2 Dépendances
- **Eigen3** : Algèbre linéaire (types `Matrix`, `Vector`, `Quaternion`)
- **STL C++20** : `concepts`, `type_traits`, `ranges`
- **SOFA Framework** : Framework de simulation physique

---

## 3. Éléments Déjà Différentiables

### 3.1 Cartes Exponentielles et Logarithmiques ✅

Toutes les classes implémentent :
```cpp
static Derived computeExp(const TangentVector &xi);
TangentVector computeLog() const;
```

Ces fonctions sont **mathématiquement différentiables** et constituent la base de la géométrie différentielle des groupes de Lie.

**Exemple (SE3)** :
```cpp
static SE3 computeExp(const TangentVector &xi) {
    const Vector3 rho = xi.template tail<3>();
    const Vector3 phi = xi.template head<3>();
    const Scalar angle = phi.norm();
    const SO3Type R = SO3Type::exp(phi);
    // Formule de Rodrigues avec traitement des petits angles
    Matrix3 V = /* ... calcul stable ... */;
    return SE3(R, V * rho);
}
```

### 3.2 Représentation Adjointe ✅

Implémentation de l'adjoint `Ad_g: 𝔤 → 𝔤` :
```cpp
AdjointMatrix computeAdjoint() const noexcept;
```

Cette fonction est déjà différentiable car elle représente une application linéaire.

### 3.3 Différentielles de l'Exponentielle ✅

Des méthodes comme `dexp` et `dexpInv` sont déclarées dans `LieGroupBase.h` :
```cpp
static AdjointMatrix dexp(const TangentVector &v);
static AdjointMatrix dexpInv(const TangentVector &v);
AdjointMatrix dlog() const;
```

### 3.4 Interpolation et Distance ✅

Fonctions géométriques différentiables :
```cpp
Derived interpolate(const Derived &other, const Scalar &t) const;
Scalar distance(const Derived &other) const;
```

---

## 4. Lacunes pour la Différentiabilité Complète

### 4.1 Pas de Support pour la Différentiation Automatique ❌

**Problème principal** : La librairie utilise des scalaires standards (`float`, `double`) sans support pour les types duaux nécessaires à l'AD.

**Solutions possibles** :
1. **Templating sur le type scalaire** : Déjà partiellement fait, mais pas exploité
2. **Intégration avec des librairies AD** :
   - [autodiff](https://github.com/autodiff/autodiff) (C++)
   - [CppAD](https://github.com/coin-or/CppAD)
   - [Enzyme](https://enzyme.mit.edu/) (LLVM-based)

### 4.2 Jacobiens Analytiques Incomplets ❌

Certaines fonctions cruciales n'ont pas de Jacobiens explicites :

| Fonction | Jacobien Requis | État |
|----------|----------------|------|
| `exp(ξ)` | ∂exp/∂ξ | ⚠️ Partiel (dexp) |
| `log(g)` | ∂log/∂g | ⚠️ Partiel (dlog) |
| `compose(g, h)` | ∂(g∘h)/∂g, ∂(g∘h)/∂h | ❌ Manquant |
| `inverse(g)` | ∂g⁻¹/∂g | ❌ Manquant |
| `act(g, p)` | ∂(g·p)/∂g, ∂(g·p)/∂p | ⚠️ Partiel (actionJacobian) |

### 4.3 Gestion de la Stabilité Numérique ⚠️

Les approximations de Taylor pour petits angles sont bien gérées :
```cpp
if (theta < Types<Scalar>::epsilon()) {
    // Approximation de premier ordre
    return Vector(m_quat.x() * Scalar(2), ...);
}
```

Mais **pas de garantie de différentiabilité** au point de singularité (θ = 0).

### 4.4 Représentation Non Minimale ⚠️

**SO3** : Utilise des quaternions (4 paramètres) pour représenter 3 DoF
- Avantage : Pas de singularités gimbal lock
- Inconvénient : Contrainte de normalisation `‖q‖ = 1` non différentiable au sens classique

---

## 5. Recommandations pour Rendre la Librairie Différentiable

### 5.1 Approche Minimale (Effort : Modéré)

#### Objectif
Permettre le calcul de gradients pour l'optimisation (type PyTorch/TensorFlow).

#### Actions
1. **Templatiser complètement sur le type scalaire**
   ```cpp
   template<typename ADScalar = double>
   class SE3 : public LieGroupBase<SE3<ADScalar>, ADScalar, ...> {
       // Toutes les opérations deviennent AD-compatibles
   };
   ```

2. **Ajouter les Jacobiens analytiques manquants**
   - Composition : 
     ```cpp
     std::pair<AdjointMatrix, AdjointMatrix> 
     computeComposeJacobians(const Derived &other) const;
     ```
   - Inverse :
     ```cpp
     AdjointMatrix computeInverseJacobian() const;
     ```

3. **Intégrer une librairie AD légère**
   - Recommandation : **autodiff** (header-only, compatible Eigen)
   - Exemple d'intégration :
     ```cpp
     #include <autodiff/forward/dual.hpp>
     using dual = autodiff::dual;
     
     SE3<dual> g = SE3<dual>::exp(xi);
     auto [value, jacobian] = autodiff::gradient(
         [](const VectorXdual& x) { return g.log().norm(); },
         xi
     );
     ```

### 5.2 Approche Complète (Effort : Important)

#### Objectif
Support complet pour l'apprentissage différentiable (type JAX, PyTorch Geometric).

#### Actions Additionnelles
1. **Implémenter les règles de différentiation personnalisées**
   ```cpp
   template<typename Scalar>
   struct ExpBackward {
       static TangentVector apply(
           const TangentVector &grad_output,
           const TangentVector &xi
       ) {
           // Règle de chaîne pour exp
           return dexp(xi).transpose() * grad_output;
       }
   };
   ```

2. **Gérer les singularités**
   - Implémenter des versions "smooth" des fonctions :
     ```cpp
     // Au lieu de if/else, utiliser des fonctions lisses
     Scalar alpha = smoothstep(theta, 0, epsilon);
     return alpha * approx_small + (1-alpha) * exact;
     ```

3. **Support pour le calcul différentiel d'ordre supérieur**
   - Hessiens pour l'optimisation de second ordre
   - Utilisable pour le contrôle optimal

4. **Bindings Python avec différentiabilité**
   ```python
   import torch
   from cosserat_lie import SE3
   
   class SE3Function(torch.autograd.Function):
       @staticmethod
       def forward(ctx, xi):
           g = SE3.exp(xi)
           ctx.save_for_backward(xi)
           return g
       
       @staticmethod
       def backward(ctx, grad_output):
           xi, = ctx.saved_tensors
           return SE3.dexp(xi).T @ grad_output
   ```

### 5.3 Approche Hybride (Recommandée)

**Stratégie** : Différentiation automatique pour le prototypage + Jacobiens analytiques pour la performance.

1. **Développement** : Utiliser AD pour vérifier les implémentations
2. **Production** : Basculer sur les Jacobiens analytiques optimisés
3. **Tests** : Comparaison numérique systématique

```cpp
#ifdef USE_AUTODIFF
    return computeJacobianAD();
#else
    return computeJacobianAnalytic();
#endif
```

---

## 6. Cas d'Usage pour la Différentiabilité

### 6.1 Optimisation de Trajectoires
```cpp
// Minimiser l'énergie d'une courbe de Cosserat
auto energy = [&](const VectorXd& strains) {
    SE3<dual> g = SE3<dual>::identity();
    for (int i = 0; i < n_sections; ++i) {
        auto strain_i = strains.segment<6>(6*i);
        g = g * SE3<dual>::expCosserat(strain_i, ds);
    }
    return (g.translation() - target).squaredNorm();
};

auto [E, grad_E] = autodiff::gradient(energy, initial_strains);
```

### 6.2 Apprentissage de Modèles Inverses
- Calibration de paramètres de raideur
- Identification de forces externes
- Contrôle optimal basé sur gradient

### 6.3 Simulation Différentiable
- Backpropagation à travers les pas de temps
- Optimisation de design de robots souples
- Co-optimisation contrôle/design

---

## 7. Estimation de l'Effort

| Tâche | Difficulté | Temps Estimé |
|-------|-----------|--------------|
| Templatiser sur ADScalar | Facile | 1-2 jours |
| Ajouter jacobiens manquants (analytique) | Moyen | 1 semaine |
| Intégrer autodiff | Facile | 2-3 jours |
| Gérer les singularités | Difficile | 2 semaines |
| Tests et validation | Moyen | 1 semaine |
| Bindings Python | Facile | 3-4 jours |
| **TOTAL (Approche Minimale)** | - | **3-4 semaines** |
| **TOTAL (Approche Complète)** | - | **6-8 semaines** |

---

## 8. Exemples de Code à Ajouter

### 8.1 Jacobien de la Composition

```cpp
template<typename _Scalar>
std::pair<typename SE3<_Scalar>::AdjointMatrix, 
          typename SE3<_Scalar>::AdjointMatrix>
SE3<_Scalar>::composeJacobians(const SE3 &h) const {
    // ∂(g∘h)/∂g = Ad_h⁻¹  (Jacobien gauche)
    AdjointMatrix J_left = h.inverse().adjoint();
    
    // ∂(g∘h)/∂h = I  (Jacobien droit)
    AdjointMatrix J_right = AdjointMatrix::Identity();
    
    return {J_left, J_right};
}
```

### 8.2 Jacobien de l'Inverse

```cpp
template<typename _Scalar>
typename SE3<_Scalar>::AdjointMatrix
SE3<_Scalar>::inverseJacobian() const {
    // ∂g⁻¹/∂g = -Ad_{g⁻¹}
    return -inverse().adjoint();
}
```

### 8.3 Support pour autodiff

```cpp
// Dans Types.h
template<typename _Scalar>
struct is_autodiff_type : std::false_type {};

#ifdef WITH_AUTODIFF
template<typename T>
struct is_autodiff_type<autodiff::dual<T>> : std::true_type {};

template<typename T>
struct is_autodiff_type<autodiff::dual2nd<T>> : std::true_type {};
#endif

// Adapter les fonctions selon le type
template<typename Scalar>
Scalar safeSqrt(const Scalar &x) {
    if constexpr (is_autodiff_type<Scalar>::value) {
        return autodiff::sqrt(autodiff::max(Scalar(0), x));
    } else {
        return std::sqrt(std::max(Scalar(0), x));
    }
}
```

---

## 9. Risques et Limitations

### 9.1 Performance
- **AD Forward** : Coût = O(n) avec n = nombre de variables
- **AD Reverse** : Coût = O(m) avec m = nombre de sorties
- **Jacobiens Analytiques** : Coût = O(1) mais effort d'implémentation

**Recommandation** : Mode hybride avec flag de compilation.

### 9.2 Compatibilité
- **Eigen et AD** : Généralement compatible mais attention aux types
- **SOFA Framework** : Peut nécessiter des adaptateurs pour passer des types AD

### 9.3 Singularités
Les groupes de Lie ont des **points singuliers** :
- SO3 : θ = π (rotation de 180°)
- SE3 : Près de l'identité
- Requiert des **branch cuts** bien définis

---

## 10. Conclusion

### Points Forts Actuels ✅
1. Architecture CRTP flexible et moderne
2. Implémentation mathématiquement correcte
3. Gestion des petits angles pour la stabilité numérique
4. Code bien structuré et documenté
5. Fondations différentielles présentes (exp, log, adjoint)

### Améliorations Nécessaires 🔧
1. Templating complet sur le type scalaire
2. Jacobiens analytiques pour toutes les opérations clés
3. Intégration d'une librairie AD
4. Gestion robuste des singularités
5. Tests exhaustifs des dérivées

### Faisabilité Globale
**Score : 8/10** - La librairie est déjà bien conçue pour devenir différentiable. L'effort principal réside dans :
- L'ajout des Jacobiens manquants (effort moyen)
- L'intégration d'une librairie AD (effort faible)
- La validation extensive (effort moyen)

### Recommandation Finale
**Procéder avec l'Approche Hybride** :
1. Phase 1 (1 mois) : Intégrer autodiff + jacobiens critiques
2. Phase 2 (1 mois) : Optimisation + tests exhaustifs
3. Phase 3 (optionnel) : Bindings Python pour ML/RL

Cette approche permettra d'obtenir une librairie **production-ready** pour :
- L'optimisation de trajectoires
- Le contrôle optimal
- L'apprentissage différentiable
- La calibration automatique

---

## Références

1. **Lie Groups for Computer Vision** - Eade, E. (2017)
2. **A micro Lie theory for state estimation in robotics** - Sola, J. et al. (2018)
3. **autodiff Documentation** - https://autodiff.github.io
4. **Eigen Documentation** - https://eigen.tuxfamily.org
5. **Discrete Cosserat approach** - IEEE Paper référencé dans le README

---

**Contact pour questions** : Voir issues GitHub du projet plugin.Cosserat
