# Notes de Debug - CosseratTrajectoryOptimizer

## Statut Actuel

L'optimiseur compile et s'exécute, mais **ne converge pas** vers la cible.

### Symptômes Observés

1. **Position reste à (0,0,0)** tout au long de l'optimisation
2. **Coût constant** : 0.34 à chaque itération
3. **Strains quasi-nuls** : ~1e-6, pratiquement pas de changement
4. **Gradients presque nuls** : La backpropagation ne produit pas de gradients significatifs

### Diagnostic

Le test `test_gradients.cpp` montre que les gradients **existent bien** numériquement :
```
Position actuelle: 0.01    0    0
Cible: 0.3   0   0

Section 0:
  strain[3]: grad = -0.029  <- ρx (élongation) affecte la position en X !
Section 1:
  strain[3]: grad = -0.029
Section 2:
  strain[3]: grad = -0.029
```

**NOTE** : strain[3] correspond à ρx (elongation), pas une rotation.
Convention : strain = [φx, φy, φz, ρx, ρy, ρz]
Voir `STRAIN_CONVENTION.md` pour les détails.

**Conclusion** : Le problème est dans l'implémentation de `backpropagateThroughChain()`.

---

## Problème Identifié

### Dans `CosseratTrajectoryOptimizer.h` lignes 407-487

La fonction `backpropagateThroughChain()` utilise une approximation **trop simplifiée** :

```cpp
// Ligne 467 : Traitement de la translation
strain_gradients[i].template head<3>() = local_grad * section_length;

// Lignes 469-482 : Traitement de la rotation
// Utilise position_in_local = Vector3::Zero() ou g_section.translation()
// Ceci ne propage PAS correctement à travers toute la chaîne !
```

### Pourquoi ça ne marche pas

1. **Propagation incomplète** : Le gradient ne se propage pas correctement à travers la chaîne de transformations composées
2. **Simplification excessive** : Le code suppose que `position_in_local` est soit zero soit la translation locale, mais en réalité il faut propager la position de tous les segments suivants
3. **Pas d'utilisation des jacobiens SE3** : Les méthodes `composeJacobians()` et `actionJacobians()` implémentées en Phase 2 ne sont pas utilisées !

---

## Solutions Proposées

### Solution 1 : Utiliser les Différences Finies (RAPIDE ⚡)

La plus simple pour avoir un optimiseur fonctionnel immédiatement :

```cpp
Cost computeCostAndGradient(...) const {
    Cost cost;
    // ... calcul du coût ...
    
    // Gradient par différences finies
    const double h = 1e-7;
    for (int i = 0; i < n_sections; ++i) {
        for (int j = 0; j < 6; ++j) {
            auto strains_plus = strains;
            strains_plus[i](j) += h;
            double cost_plus = computeCost(strains_plus, target, ...);
            
            auto strains_minus = strains;
            strains_minus[i](j) -= h;
            double cost_minus = computeCost(strains_minus, target, ...);
            
            cost.gradient[i](j) = (cost_plus - cost_minus) / (2.0 * h);
        }
    }
    
    return cost;
}
```

**Avantages** :
- ✅ Simple à implémenter (30 lignes)
- ✅ Marche à coup sûr
- ✅ Pas besoin de debug complexe

**Inconvénients** :
- ❌ Plus lent (2 × n_sections × 6 évaluations du coût par itération)
- ❌ Moins précis numériquement

---

### Solution 2 : Backpropagation Correcte avec Jacobiens Analytiques (OPTIMAL 🎯)

Utiliser les jacobiens implémentés en Phase 2 :

```cpp
void backpropagateThroughChain(...) const {
    const int n_sections = static_cast<int>(strains.size());
    
    // Gradient accumulé (dans l'espace tangent)
    Vector6 grad_se3 = Vector6::Zero();
    grad_se3.template head<3>() = position_gradient;  // Partie translation
    
    // Backprop en ordre inverse
    for (int i = n_sections - 1; i >= 0; --i) {
        SE3Type g_i = transforms[i];
        SE3Type g_section = SE3Type::exp(strains[i] * section_length);
        
        // Utiliser les jacobiens de composition
        auto [J_left, J_right] = g_i.composeJacobians(g_section);
        
        // Propager le gradient à travers la composition
        // grad_i = J_right^T * grad_{i+1}
        Vector6 grad_local = J_right.transpose() * grad_se3;
        
        // Propager à travers l'exponentielle
        // TODO: Implémenter dexp^{-1} ou utiliser approximation
        strain_gradients[i] = grad_local * section_length;
        
        // Propager au niveau précédent
        grad_se3 = J_left.transpose() * grad_se3;
    }
}
```

**Avantages** :
- ✅ Rapide (1 évaluation par itération)
- ✅ Précis
- ✅ Utilise les jacobiens déjà implémentés

**Inconvénients** :
- ❌ Plus complexe à implémenter correctement
- ❌ Nécessite potentiellement `dexp^{-1}` (Jacobien inverse de l'exponentielle)

---

### Solution 3 : Optimisation sur le Manifold avec Rétraction (AVANCÉ 🚀)

Utiliser des méthodes d'optimisation géométrique qui respectent la structure du manifold :

```cpp
// Au lieu de : strains[i] -= learning_rate * gradient[i]
// Utiliser une rétraction sur SE(3)

for (int i = 0; i < n_sections; ++i) {
    Vector6 descent_direction = -learning_rate * gradient[i];
    SE3Type g_i = SE3Type::exp(strains[i] * section_length);
    SE3Type g_updated = g_i * SE3Type::exp(descent_direction);
    strains[i] = g_updated.log() / section_length;
}
```

**Avantages** :
- ✅ Géométriquement correct
- ✅ Meilleure convergence théorique

**Inconvénients** :
- ❌ Très complexe
- ❌ Nécessite bibliothèque d'optimisation sur manifolds (ex: Manopt)

---

## Recommandation Immédiate

**Pour avoir un optimiseur fonctionnel rapidement** : Implémenter **Solution 1** (différences finies).

Temps estimé : 15-30 minutes

### Code à Modifier

Dans `CosseratTrajectoryOptimizer.h`, remplacer la méthode `computeCostAndGradient()` par :

```cpp
Cost computeCostAndGradient(...) const {
    Cost cost;
    const int n_sections = static_cast<int>(strains.size());
    cost.gradient.resize(n_sections, Vector6::Zero());
    
    // Fonction coût auxiliaire
    auto eval_cost = [&](const std::vector<Vector6>& s) -> double {
        SE3Type g = SE3Type::Identity();
        for (const auto& strain : s) {
            g = g * SE3Type::exp(strain * section_length);
        }
        Vector3 error = g.translation() - target.translation();
        double c = 0.5 * error.squaredNorm();
        for (const auto& strain : s) {
            c += 0.5 * regularization * strain.squaredNorm();
        }
        return c;
    };
    
    // Coût actuel
    cost.value = eval_cost(strains);
    
    // Gradients par différences finies centrales
    const double h = 1e-7;
    for (int i = 0; i < n_sections; ++i) {
        for (int j = 0; j < 6; ++j) {
            auto strains_plus = strains;
            strains_plus[i](j) += h;
            
            auto strains_minus = strains;
            strains_minus[i](j) -= h;
            
            cost.gradient[i](j) = (eval_cost(strains_plus) - eval_cost(strains_minus)) / (2.0 * h);
        }
    }
    
    return cost;
}
```

Ensuite supprimer la méthode `backpropagateThroughChain()` qui n'est plus utilisée.

---

## Tests Après Correction

Après modification, l'exemple devrait :
- ✅ Converger vers la cible
- ✅ Coût décroissant à chaque itération
- ✅ Strains non-nuls et significatifs
- ✅ Position finale proche de (0.8, 0, 0.2)

---

## Pour Plus Tard : Solution 2

Une fois l'optimiseur fonctionnel avec les différences finies, on pourra implémenter la Solution 2 pour de meilleures performances, en utilisant les jacobiens analytiques correctement.

Cela nécessitera :
1. Implémenter `dexp()` et `dexpInv()` pour SE3
2. Propager correctement à travers la chaîne de compositions
3. Valider avec les différences finies

---

**Créé** : 23 Décembre 2025  
**Auteur** : Warp AI Agent
