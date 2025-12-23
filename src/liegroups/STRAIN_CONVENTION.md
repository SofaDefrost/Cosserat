# Convention des Strains Cosserat

## Vue d'Ensemble

Dans la librairie liegroups, les strains Cosserat sont représentés par des vecteurs 6D appartenant à l'algèbre de Lie se(3). Ce document clarifie la convention utilisée et les subtilités importantes.

---

## Convention du Vector6

Un strain Cosserat est représenté par un `Vector6` avec la convention suivante :

```
strain = [φx, φy, φz, ρx, ρy, ρz]ᵀ
         └────┬────┘  └────┬────┘
           rotation    translation
```

### Partie Rotation (indices 0-2) : φ = [φx, φy, φz]

| Indice | Symbole | Nom | Description Physique |
|--------|---------|-----|---------------------|
| 0 | φx | **Torsion** | Rotation autour de l'axe X (axe principal de la poutre) |
| 1 | φy | **Bending Y** | Rotation autour de l'axe Y (flexion dans le plan XZ) |
| 2 | φz | **Bending Z** | Rotation autour de l'axe Z (flexion dans le plan XY) |

### Partie Translation (indices 3-5) : ρ = [ρx, ρy, ρz]

| Indice | Symbole | Nom | Description Physique |
|--------|---------|-----|---------------------|
| 3 | ρx | **Elongation** | Étirement/compression le long de X (déviation de 1.0) |
| 4 | ρy | **Shearing Y** | Cisaillement dans la direction Y |
| 5 | ρz | **Shearing Z** | Cisaillement dans la direction Z |

---

## Subtilité Importante : Configuration de Repos

⚠️ **Point Clé** : Dans la théorie de Cosserat, la configuration de repos a une **élongation nominale de 1.0** le long de l'axe X.

### Implémentation dans `SE3::buildXiHat()` (lignes 686-688)

```cpp
xi_hat(0, 3) = 1.0 + rho.x();  // Translation en X = 1 + strain_x
xi_hat(1, 3) = rho.y();         // Translation en Y = strain_y  
xi_hat(2, 3) = rho.z();         // Translation en Z = strain_z
```

Cela signifie :
- **ρx = 0** → élongation nominale (poutre rectiligne naturelle)
- **ρx > 0** → étirement (poutre plus longue)
- **ρx < 0** → compression (poutre plus courte)

La translation effective dans `buildXiHat` est donc `[1+ρx, ρy, ρz]`, pas simplement `[ρx, ρy, ρz]`.

### Pourquoi cette Convention ?

Dans la modélisation de poutres Cosserat :
1. L'axe X est l'**axe curviligne** de la poutre
2. La configuration de référence est une poutre **rectiligne le long de X**
3. Le paramètre arc-length s varie de 0 à L
4. La dérivée spatiale du frame `d/ds` a une composante nominale de 1 en X
5. Les strains représentent des **déviations** par rapport à cette configuration

Donc : `v = v_rest + strain` où `v_rest = [1, 0, 0]`

---

## Correspondance avec se(3)

Dans l'algèbre de Lie se(3), un élément ξ est souvent écrit :

```
ξ = [v, ω]ᵀ  où v ∈ ℝ³ (vélocité linéaire), ω ∈ ℝ³ (vélocité angulaire)
```

**Attention** : Notre convention est **inversée** !

| Notre Convention | Convention se(3) Standard |
|------------------|---------------------------|
| `[φ, ρ]ᵀ` | `[v, ω]ᵀ` |
| head<3>() = φ (rotation) | head<3>() = v (translation) |
| tail<3>() = ρ (translation) | tail<3>() = ω (rotation) |

### Dans le Code

Voir `SE3.h` lignes 89-90 :
```cpp
const Vector3 rho = strain.template tail<3>(); // Linear strain (translation rate)
const Vector3 phi = strain.template head<3>(); // Angular strain (rotation rate)
```

Et lignes 103-104 pour `computeExp` :
```cpp
const Vector3 rho = xi.template tail<3>();
const Vector3 phi = xi.template head<3>();
```

---

## Exemples d'Utilisation

### Exemple 1 : Torsion Pure

```cpp
Eigen::Matrix<double, 6, 1> strain_torsion;
strain_torsion << 0.1,  // φx = torsion de 0.1 rad
                  0.0,  // φy = pas de bending Y
                  0.0,  // φz = pas de bending Z
                  0.0,  // ρx = pas d'élongation
                  0.0,  // ρy = pas de shearing Y
                  0.0;  // ρz = pas de shearing Z
```

### Exemple 2 : Flexion (Bending) dans le plan XZ

```cpp
Eigen::Matrix<double, 6, 1> strain_bending;
strain_bending << 0.0,   // φx = pas de torsion
                  0.2,   // φy = flexion autour Y (courbe dans XZ)
                  0.0,   // φz = pas de flexion autour Z
                  0.0,   // ρx = pas d'élongation
                  0.0,   // ρy = pas de shearing Y
                  0.0;   // ρz = pas de shearing Z
```

### Exemple 3 : Élongation + Cisaillement

```cpp
Eigen::Matrix<double, 6, 1> strain_elongation_shear;
strain_elongation_shear << 0.0,   // φx = pas de torsion
                           0.0,   // φy = pas de bending Y
                           0.0,   // φz = pas de bending Z
                           0.05,  // ρx = élongation 5% (1.0 → 1.05)
                           0.02,  // ρy = cisaillement Y
                           0.0;   // ρz = pas de cisaillement Z
```

### Exemple 4 : Configuration Complexe

```cpp
Eigen::Matrix<double, 6, 1> strain_complex;
strain_complex << 0.1,    // Torsion
                  0.15,   // Bending Y
                  -0.05,  // Bending Z (sens opposé)
                  -0.02,  // Compression (1.0 → 0.98)
                  0.01,   // Shearing Y
                  0.01;   // Shearing Z
```

---

## Utilisation dans l'Optimisation

### CosseratTrajectoryOptimizer

Quand on optimise des strains avec `CosseratTrajectoryOptimizer`, les gradients sont calculés par rapport à ce Vector6 :

```cpp
// Gradient du coût par rapport aux strains
std::vector<Vector6> gradients;  // gradient[i] est ∂cost/∂strain_i

// Interprétation :
// gradient[i][0] → influence de la torsion sur le coût
// gradient[i][1] → influence du bending Y sur le coût
// gradient[i][2] → influence du bending Z sur le coût
// gradient[i][3] → influence de l'élongation sur le coût
// gradient[i][4] → influence du shearing Y sur le coût
// gradient[i][5] → influence du shearing Z sur le coût
```

### Forward Kinematics

Pour calculer la pose finale d'une poutre avec N sections :

```cpp
SE3d g = SE3d::Identity();
for (int i = 0; i < n_sections; ++i) {
    SE3d g_section = SE3d::expCosserat(strains[i], section_length);
    g = g * g_section;  // Composition
}
Vector3 tip_position = g.translation();
```

---

## Lien avec les Propriétés Matérielles

Les strains sont liés aux forces/moments internes via la loi constitutive de Hooke :

### Dans `HookeSeratBaseForceField`

```cpp
// Relation contrainte-déformation
Moment_x = G * J * φx        // Torsion
Moment_y = E * Iy * φy       // Bending Y
Moment_z = E * Iz * φz       // Bending Z
Force_x  = E * A * ρx        // Traction/compression
Force_y  = k * G * A * ρy    // Shearing Y
Force_z  = k * G * A * ρz    // Shearing Z
```

Où :
- **E** : Module de Young
- **G** : Module de cisaillement
- **A** : Aire de section
- **Iy, Iz** : Moments d'inertie (bending)
- **J** : Constante de torsion
- **k** : Coefficient de cisaillement (Timoshenko)

---

## Validation et Tests

### Test de Cohérence

```cpp
// Vérifier que head/tail sont corrects
Vector6 strain;
strain << 1, 2, 3, 4, 5, 6;

Vector3 phi = strain.template head<3>();  // [1, 2, 3] rotation
Vector3 rho = strain.template tail<3>();  // [4, 5, 6] translation

assert(phi[0] == 1.0);  // φx (torsion)
assert(phi[1] == 2.0);  // φy (bending Y)
assert(phi[2] == 3.0);  // φz (bending Z)
assert(rho[0] == 4.0);  // ρx (elongation)
assert(rho[1] == 5.0);  // ρy (shearing Y)
assert(rho[2] == 6.0);  // ρz (shearing Z)
```

### Tests Unitaires

Voir :
- `test_gradients.cpp` : Validation des gradients numériques
- `test_trajectory_optimization.cpp` : Tests d'optimisation
- `test_HookeSerat_*.cpp` : Tests du modèle physique

---

## Références

### Dans le Code Source

1. **SE3.h** lignes 86-100 : `expCosserat()` - Exponentielle avec convention Cosserat
2. **SE3.h** lignes 644-692 : `buildXiHat()` - Construction de la matrice se(3) avec élongation nominale
3. **HookeSeratBaseForceField.inl** lignes 104-174 : Calcul des propriétés de section
4. **CosseratTrajectoryOptimizer.h** : Optimisation utilisant les strains

### Documentation Externe

1. **Cosserat Theory** : Antman, S. S. (2005). "Nonlinear Problems of Elasticity"
2. **Lie Groups for Robotics** : Sola et al. (2018). "A micro Lie theory"
3. **Beam Theory** : Timoshenko, S. P. (1921). "On the correction for shear"

---

## Résumé Visuel

```
         STRAIN VECTOR (6D)
    ┌──────────────────────────┐
    │  φx  │ Torsion           │
    │  φy  │ Bending Y         │ ← head<3>() = Rotation
    │  φz  │ Bending Z         │
    ├──────┼───────────────────┤
    │  ρx  │ Elongation        │
    │  ρy  │ Shearing Y        │ ← tail<3>() = Translation
    │  ρz  │ Shearing Z        │
    └──────┴───────────────────┘
    
    REPÈRE LOCAL POUTRE
         Z
         ↑
         │
         │
         └────→ X (axe poutre)
        ╱
       ╱
      Y
```

---

**Créé** : 23 Décembre 2025  
**Auteur** : Warp AI Agent  
**Version** : 1.0
