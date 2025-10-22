# Améliorations proposées pour la bibliothèque liegroups

**Date :** 2025-01-22  
**Basé sur :** "A micro Lie theory for state estimation in robotics" (Solà et al., 2021)  
**Analyse du code :** ~/travail/plugin/plugin.Cosserat/src/liegroups/

## Vue d'ensemble

Votre implémentation actuelle est solide et bien structurée, utilisant le CRTP pattern moderne. Cependant, le papier de Solà propose plusieurs concepts qui pourraient enrichir considérablement votre bibliothèque, notamment pour les applications d'estimation d'état et de propagation d'incertitude.

---

## 1. ⭐ PRIORITÉ HAUTE : Opérateurs ⊕ et ⊖ (Plus/Minus)

### Problème actuel
Votre code implémente `compose()` et `log()` séparément, mais n'offre pas l'interface unifiée ⊕/⊖ qui simplifie grandement l'écriture d'algorithmes d'estimation.

### Solution proposée

Ajouter à `LieGroupBase.h` :

```cpp
/**
 * @brief Right-plus operator: X ⊕ τ = X ◦ Exp(τ)
 * @param tau Tangent vector in local frame (TₓM)
 * @return Composed element
 */
[[nodiscard]] Derived plus(const TangentVector &tau) const noexcept {
    return derived().compose(Derived::exp(tau));
}

/**
 * @brief Right-minus operator: Y ⊖ X = Log(X⁻¹ ◦ Y)
 * @param other Another group element
 * @return Tangent vector in local frame at X
 */
[[nodiscard]] TangentVector minus(const Derived &other) const noexcept {
    return derived().inverse().compose(other).log();
}

// Surcharge d'opérateurs pour usage intuitif
[[nodiscard]] Derived operator+(const TangentVector &tau) const noexcept {
    return plus(tau);
}

[[nodiscard]] TangentVector operator-(const Derived &other) const noexcept {
    return minus(other);
}
```

**Impact :**
- Simplifie l'écriture d'algorithmes ESKF
- Cohérence avec la notation du papier
- Facilite la propagation d'incertitudes

**Exemple d'utilisation :**
```cpp
// Avant
SE3d X_updated = X.compose(SE3d::exp(delta));

// Après (plus intuitif)
SE3d X_updated = X + delta;  // ou X.plus(delta)

// Pour les erreurs
TangentVector error = X_measured - X_estimated;  // ou X_measured.minus(X_estimated)
```

---

## 2. ⭐ PRIORITÉ HAUTE : Jacobiennes droites et gauches (Jr, Jl)

### Problème actuel
Votre code implémente `dexp()` mais ne fournit pas explicitement les jacobiennes **Jr(τ)** et **Jl(τ)** qui sont essentielles pour la propagation d'incertitude.

### Solution proposée

Ajouter à `LieGroupBase.h` :

```cpp
/**
 * @brief Right Jacobian of M: Jr(τ) = τD Exp(τ)/Dτ
 * Maps variations of τ to local tangent space at Exp(τ)
 * @param tau Tangent vector
 * @return Right Jacobian matrix
 */
[[nodiscard]] static AdjointMatrix rightJacobian(const TangentVector &tau) noexcept {
    return Derived::computeRightJacobian(tau);
}

/**
 * @brief Left Jacobian of M: Jl(τ) = ᴱD Exp(τ)/Dτ  
 * Maps variations of τ to global tangent space (Lie algebra)
 * @param tau Tangent vector
 * @return Left Jacobian matrix
 */
[[nodiscard]] static AdjointMatrix leftJacobian(const TangentVector &tau) noexcept {
    return Derived::computeLeftJacobian(tau);
}

/**
 * @brief Inverse of right Jacobian
 */
[[nodiscard]] static AdjointMatrix rightJacobianInverse(const TangentVector &tau) noexcept {
    return Derived::computeRightJacobianInverse(tau);
}

/**
 * @brief Inverse of left Jacobian
 */
[[nodiscard]] static AdjointMatrix leftJacobianInverse(const TangentVector &tau) noexcept {
    return Derived::computeLeftJacobianInverse(tau);
}
```

**Implémentation pour SO(3)** (ajouter à `SO3.h`) :

```cpp
static AdjointMatrix computeRightJacobian(const TangentVector &omega) noexcept {
    const Scalar theta = omega.norm();
    
    if (theta < Types<Scalar>::epsilon()) {
        return Matrix::Identity() - Scalar(0.5) * hat(omega);
    }
    
    const Matrix omega_hat = hat(omega);
    const Matrix omega_hat2 = omega_hat * omega_hat;
    const Scalar theta2 = theta * theta;
    const Scalar theta3 = theta2 * theta;
    
    // Jr(θ) = I - (1-cos θ)/θ² [θ]× + (θ-sin θ)/θ³ [θ]²×
    return Matrix::Identity() 
           - ((Scalar(1) - std::cos(theta)) / theta2) * omega_hat
           + ((theta - std::sin(theta)) / theta3) * omega_hat2;
}

static AdjointMatrix computeLeftJacobian(const TangentVector &omega) noexcept {
    // Jl(θ) = Jr(-θ)
    return computeRightJacobian(-omega);
}

static AdjointMatrix computeRightJacobianInverse(const TangentVector &omega) noexcept {
    const Scalar theta = omega.norm();
    
    if (theta < Types<Scalar>::epsilon()) {
        return Matrix::Identity() + Scalar(0.5) * hat(omega);
    }
    
    const Matrix omega_hat = hat(omega);
    const Matrix omega_hat2 = omega_hat * omega_hat;
    const Scalar theta2 = theta * theta;
    
    // Jr⁻¹(θ) = I + ½[θ]× + (1/θ² - (1+cos θ)/(2θ sin θ))[θ]²×
    const Scalar coeff = Scalar(1) / theta2 
                        - (Scalar(1) + std::cos(theta)) / (Scalar(2) * theta * std::sin(theta));
    
    return Matrix::Identity() 
           + Scalar(0.5) * omega_hat 
           + coeff * omega_hat2;
}
```

**Impact :**
- Propagation d'incertitude correcte dans les filtres de Kalman
- Approximations linéaires précises autour des éléments du groupe
- Support des algorithmes d'optimisation sur variétés

---

## 3. ⭐ PRIORITÉ MOYENNE : Jacobiens des opérations élémentaires

### Problème actuel
Manque les jacobiens pour l'inversion, la composition, et l'action de groupe.

### Solution proposée

Ajouter à `LieGroupBase.h` :

```cpp
/**
 * @brief Jacobian of inverse: J_X⁻¹_X = -Ad_X
 * @return Jacobian matrix
 */
[[nodiscard]] AdjointMatrix inverseJacobian() const noexcept {
    return -adjoint();
}

/**
 * @brief Jacobian of composition wrt first argument: J_XY_X = Ad_Y⁻¹
 * @param other Second element in composition
 * @return Jacobian matrix
 */
[[nodiscard]] AdjointMatrix composeJacobianFirst(const Derived &other) const noexcept {
    return other.inverse().adjoint();
}

/**
 * @brief Jacobian of composition wrt second argument: J_XY_Y = I
 * @return Identity matrix
 */
[[nodiscard]] static AdjointMatrix composeJacobianSecond() noexcept {
    return AdjointMatrix::Identity();
}

/**
 * @brief Jacobian of plus operator wrt base point: J_X⊕τ_X = Ad_Exp(τ)⁻¹
 * @param tau Tangent increment
 * @return Jacobian matrix
 */
[[nodiscard]] AdjointMatrix plusJacobianPoint(const TangentVector &tau) const noexcept {
    return Derived::exp(tau).inverse().adjoint();
}

/**
 * @brief Jacobian of plus operator wrt tangent: J_X⊕τ_τ = Jr(τ)
 * @param tau Tangent increment
 * @return Jacobian matrix
 */
[[nodiscard]] static AdjointMatrix plusJacobianTangent(const TangentVector &tau) noexcept {
    return rightJacobian(tau);
}

/**
 * @brief Jacobian of minus operator wrt first argument: J_Y⊖X_Y = Jr⁻¹(τ)
 * where τ = Y ⊖ X
 * @param other Second element
 * @return Jacobian matrix
 */
[[nodiscard]] AdjointMatrix minusJacobianFirst(const Derived &other) const noexcept {
    TangentVector tau = minus(other);
    return rightJacobianInverse(tau);
}

/**
 * @brief Jacobian of minus operator wrt second argument: J_Y⊖X_X = -Jl⁻¹(τ)
 * where τ = Y ⊖ X
 * @param other Second element
 * @return Jacobian matrix
 */
[[nodiscard]] AdjointMatrix minusJacobianSecond(const Derived &other) const noexcept {
    TangentVector tau = minus(other);
    return -leftJacobianInverse(tau);
}
```

**Implémentation pour SO(3)** (ajouter à `SO3.h`) :

```cpp
/**
 * @brief Jacobian of rotation action on point: J_Rv_R = -R[v]×
 * @param point Point to rotate
 * @return Jacobian matrix 3×3
 */
JacobianMatrix actionJacobianRotation(const Vector &point) const noexcept {
    return -matrix() * hat(point);
}

/**
 * @brief Jacobian of rotation action wrt point: J_Rv_v = R
 * @return Rotation matrix 3×3
 */
JacobianMatrix actionJacobianPoint() const noexcept {
    return matrix();
}
```

---

## 4. ⭐ PRIORITÉ MOYENNE : Support pour la gestion d'incertitude

### Solution proposée

Créer un nouveau fichier `Uncertainty.h` :

```cpp
#pragma once

#include "LieGroupBase.h"
#include <Eigen/Dense>

namespace sofa::component::cosserat::liegroups {

/**
 * @brief Gaussian distribution on a Lie group manifold
 * 
 * Represents X ~ N(X̄, Σ) where X̄ ∈ M is the mean and 
 * Σ ∈ ℝᵐˣᵐ is the covariance in the tangent space at X̄
 */
template<typename LieGroupType>
class GaussianOnManifold {
public:
    using Scalar = typename LieGroupType::Scalar;
    using TangentVector = typename LieGroupType::TangentVector;
    using Covariance = typename LieGroupType::AdjointMatrix;
    
    GaussianOnManifold(const LieGroupType &mean, const Covariance &cov)
        : m_mean(mean), m_covariance(cov) {}
    
    const LieGroupType& mean() const { return m_mean; }
    const Covariance& covariance() const { return m_covariance; }
    
    /**
     * @brief Transform covariance from local to global frame
     * Using: ᴱΣ = Ad_X ᵡΣ Ad_X^T
     */
    Covariance toGlobalFrame() const {
        auto Ad = m_mean.adjoint();
        return Ad * m_covariance * Ad.transpose();
    }
    
    /**
     * @brief Transform covariance from global to local frame
     * Using: ᵡΣ = Ad_X⁻¹ ᴱΣ (Ad_X⁻¹)^T
     */
    static Covariance toLocalFrame(const LieGroupType &X, const Covariance &global_cov) {
        auto Ad_inv = X.inverse().adjoint();
        return Ad_inv * global_cov * Ad_inv.transpose();
    }
    
    /**
     * @brief Propagate uncertainty through a function f: M → N
     * Using: Σ_Y ≈ J Σ_X J^T where J = Df(X)/DX
     */
    template<typename FunctionType, typename OutputLieGroupType>
    static GaussianOnManifold<OutputLieGroupType> 
    propagate(const GaussianOnManifold<LieGroupType> &input,
              FunctionType &&func,
              const typename OutputLieGroupType::AdjointMatrix &jacobian) {
        
        auto output_mean = func(input.mean());
        auto output_cov = jacobian * input.covariance() * jacobian.transpose();
        
        return GaussianOnManifold<OutputLieGroupType>(output_mean, output_cov);
    }

private:
    LieGroupType m_mean;
    Covariance m_covariance;
};

// Alias pour usage pratique
template<typename Scalar>
using GaussianSO3 = GaussianOnManifold<SO3<Scalar>>;

template<typename Scalar>
using GaussianSE3 = GaussianOnManifold<SE3<Scalar>>;

} // namespace
```

**Exemple d'utilisation :**
```cpp
// Créer une distribution gaussienne sur SO(3)
SO3d mean_rotation = SO3d::exp(Eigen::Vector3d(0.1, 0.2, 0.0));
Eigen::Matrix3d covariance = Eigen::Matrix3d::Identity() * 0.01;
GaussianSO3<double> gaussian(mean_rotation, covariance);

// Propager à travers une composition
SO3d delta = SO3d::exp(Eigen::Vector3d(0.05, 0.0, 0.0));
auto jacobian = mean_rotation.composeJacobianFirst(delta);
auto propagated = GaussianSO3<double>::propagate(
    gaussian,
    [&](const SO3d& R) { return R * delta; },
    jacobian
);
```

---

## 5. ⭐ PRIORITÉ BASSE : Groupes composites (Bundles)

### Problème actuel
Votre fichier `Bundle.h` existe mais pourrait être amélioré avec les opérateurs ⊞/⊟.

### Solution proposée

Améliorer `Bundle.h` avec :

```cpp
/**
 * @brief Composite group with diamond operators
 * Allows heterogeneous state vectors: X = ⟨X₁, X₂, ..., Xₘ⟩
 */
template<typename... Groups>
class Bundle {
public:
    using TangentVector = /* Concatenated tangent vectors */;
    
    /**
     * @brief Diamond-plus: X ⊞ τ with non-interacting blocks
     */
    Bundle diamondPlus(const TangentVector &tau) const {
        return applyToEach([](auto& g, auto& t) { 
            return g.plus(t); 
        }, tau);
    }
    
    /**
     * @brief Diamond-minus: Y ⊟ X with non-interacting blocks
     */
    TangentVector diamondMinus(const Bundle &other) const {
        return concatenate([](auto& g1, auto& g2) { 
            return g1.minus(g2); 
        }, other);
    }
    
    /**
     * @brief Block-wise Jacobian computation
     */
    template<typename Function>
    auto jacobian(Function&& f) const {
        // Compute Jacobian block by block
        // Returns matrix [∂f₁/∂X₁ ... ∂f₁/∂Xₘ]
        //                [   ⋮      ⋱     ⋮    ]
        //                [∂fₙ/∂X₁ ... ∂fₙ/∂Xₘ]
    }
};
```

**Cas d'usage :**
```cpp
// État hétérogène pour calibration
using State = Bundle<RealSpace<double,2>, SO3<double>, SE3<double>>;
// ⟨bias_calib, orientation, pose⟩

State X = /* ... */;
TangentVector perturbation = /* ... */;
State X_updated = X.diamondPlus(perturbation);
```

---

## 6. ⭐ PRIORITÉ BASSE : Approximations BCH améliorées

### Problème actuel
Votre `BCH()` existe mais pourrait être plus précis.

### Solution proposée

Améliorer dans `LieGroupBase.inl` :

```cpp
template<typename Derived, FloatingPoint Scalar, int Dim, int AlgDim, int ActDim>
typename LieGroupBase<Derived, Scalar, Dim, AlgDim, ActDim>::TangentVector
LieGroupBase<Derived, Scalar, Dim, AlgDim, ActDim>::BCH(
    const TangentVector &v, const TangentVector &w, int order) {
    
    TangentVector result = v + w;
    
    if (order < 2) return result;
    
    // 2nd order: + ½[v,w]
    auto lie_bracket = [](const TangentVector& a, const TangentVector& b) {
        return Derived::vee(Derived::hat(a) * Derived::hat(b) 
                          - Derived::hat(b) * Derived::hat(a));
    };
    
    result += Scalar(0.5) * lie_bracket(v, w);
    
    if (order < 3) return result;
    
    // 3rd order: + 1/12([v,[v,w]] + [w,[w,v]])
    result += (Scalar(1)/Scalar(12)) * 
              (lie_bracket(v, lie_bracket(v, w)) + 
               lie_bracket(w, lie_bracket(w, v)));
    
    if (order < 4) return result;
    
    // 4th order terms (more complex)
    // ... (voir papier Chirikjian pour formules complètes)
    
    return result;
}
```

---

## 7. Améliorations de documentation

### Ajouter dans `docs/` :

1. **`jacobians.md`** : Guide complet sur l'usage des jacobiens
2. **`uncertainty_propagation.md`** : Exemples ESKF complets
3. **`notation_guide.md`** : Correspondance avec notation du papier de Solà

Exemple de contenu pour `jacobians.md` :

```markdown
# Guide des Jacobiens

## Notation
- **Right Jacobian** : ᵡDf(X)/DX (variations locales)
- **Left Jacobian** : ᴱDf(X)/DX (variations globales)

## Règle de la chaîne
Pour Z = g(f(X)) :
```
DZ/DX = DZ/DY · DY/DX
```

## Exemples pratiques

### ESKF Prédiction
```cpp
// État : X ∈ SE(3)
// Commande : u ∈ ℝ⁶
// Dynamique : X_j = X_i ⊕ u

// Jacobiens
F = X_i.plusJacobianPoint(u);  // = Ad_Exp(u)⁻¹
G = SE3d::plusJacobianTangent(u);  // = Jr(u)

// Propagation covariance
P_j = F * P_i * F.transpose() + G * Q * G.transpose();
```
```

---

## 8. Tests unitaires à ajouter

Créer `tests/JacobianTests.cpp` :

```cpp
#include <gtest/gtest.h>
#include "liegroups/SO3.h"
#include "liegroups/SE3.h"

using namespace sofa::component::cosserat::liegroups;

TEST(JacobianTests, RightJacobianSO3_SmallAngle) {
    Eigen::Vector3d omega(0.001, 0.002, 0.001);
    auto Jr = SO3d::rightJacobian(omega);
    auto Jr_numerical = numericalJacobian([](auto w) {
        return SO3d::exp(w);
    }, omega);
    
    EXPECT_TRUE(Jr.isApprox(Jr_numerical, 1e-6));
}

TEST(JacobianTests, PlusMinusConsistency) {
    SO3d X = SO3d::exp(Eigen::Vector3d(0.5, 0.2, 0.1));
    Eigen::Vector3d tau(0.1, 0.05, 0.02);
    
    SO3d Y = X.plus(tau);
    Eigen::Vector3d tau_recovered = Y.minus(X);
    
    EXPECT_TRUE(tau.isApprox(tau_recovered, 1e-10));
}

TEST(UncertaintyTests, CovarianceTransform) {
    SO3d X = SO3d::exp(Eigen::Vector3d(0.3, 0.1, 0.2));
    Eigen::Matrix3d local_cov = Eigen::Matrix3d::Identity() * 0.01;
    
    auto global_cov = X.adjoint() * local_cov * X.adjoint().transpose();
    auto recovered = X.inverse().adjoint() * global_cov * 
                     X.inverse().adjoint().transpose();
    
    EXPECT_TRUE(local_cov.isApprox(recovered, 1e-10));
}
```

---

## 9. Intégration avec SOFA

### Créer des composants SOFA utilisant les nouveaux outils

`components/LieGroupStateComponent.h` :

```cpp
#include <sofa/core/State.h>
#include "liegroups/SE3.h"
#include "liegroups/Uncertainty.h"

template<typename LieGroupType>
class LieGroupStateComponent : public sofa::core::State {
public:
    using Gaussian = GaussianOnManifold<LieGroupType>;
    
    Data<std::vector<LieGroupType>> d_positions;
    Data<std::vector<typename LieGroupType::Covariance>> d_covariances;
    
    /**
     * @brief Apply velocity increment with uncertainty propagation
     */
    void applyVelocity(double dt, const std::vector<TangentVector>& velocities) {
        auto& pos = *d_positions.beginEdit();
        auto& cov = *d_covariances.beginEdit();
        
        for (size_t i = 0; i < pos.size(); ++i) {
            TangentVector delta = velocities[i] * dt;
            
            // Update position
            pos[i] = pos[i].plus(delta);
            
            // Update covariance
            auto F = pos[i].plusJacobianPoint(delta);
            cov[i] = F * cov[i] * F.transpose();
        }
        
        d_positions.endEdit();
        d_covariances.endEdit();
    }
};
```

---

## 10. Checklist de migration

### Phase 1 : Fondations (2-3 semaines)
- [ ] Ajouter opérateurs `plus()` / `minus()` à `LieGroupBase.h`
- [ ] Implémenter `rightJacobian()` / `leftJacobian()` pour SO(3)
- [ ] Implémenter `rightJacobian()` / `leftJacobian()` pour SE(3)
- [ ] Tests unitaires pour ⊕/⊖ et jacobiens
- [ ] Documentation des nouveaux opérateurs

### Phase 2 : Propagation d'incertitude (2 semaines)
- [ ] Créer `Uncertainty.h` avec `GaussianOnManifold`
- [ ] Ajouter tous les jacobiens d'opérations élémentaires
- [ ] Exemples ESKF dans `examples/`
- [ ] Tests de propagation de covariance

### Phase 3 : Avancé (2-3 semaines)
- [ ] Améliorer `Bundle.h` avec opérateurs ⊞/⊟
- [ ] BCH amélioré avec ordres supérieurs
- [ ] Intégration composants SOFA
- [ ] Documentation complète et tutoriels

---

## Conclusion

**Points forts actuels :**
- ✅ Architecture CRTP bien pensée
- ✅ Implémentations correctes de exp/log
- ✅ Support quaternions pour SO(3)

**Gains attendus des améliorations :**
1. **Notation unifiée** : Code plus lisible et conforme à la littérature
2. **Jacobiens complets** : Support natif pour ESKF, optimisation sur variétés
3. **Gestion d'incertitude** : Propagation rigoureuse dans les algorithmes d'estimation
4. **Interopérabilité** : Meilleure intégration avec frameworks d'optimisation (Ceres, g2o)

**Effort estimé :** 6-8 semaines pour implémentation complète + tests + documentation

**Ordre recommandé :**
1. Opérateurs ⊕/⊖ (impact immédiat sur lisibilité)
2. Jacobiens Jr/Jl (nécessaires pour estimation)
3. Support incertitude (valeur ajoutée pour applications robotiques)
4. Reste selon besoins spécifiques du projet Cosserat

N'hésite pas si tu veux que je détaille l'implémentation de certaines parties !
