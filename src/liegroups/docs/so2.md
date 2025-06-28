# SO(2) Implementation

## Overview

`SO2<Scalar>` implements the Special Orthogonal group in 2D, which represents rotations in a plane. SO(2) is a 1-dimensional Lie group, where the dimension corresponds to the angle of rotation.

## Mathematical Properties

- **Dimension**: 1 (the angle θ)
- **Group operation**: composition of rotations
- **Identity element**: rotation by 0 radians
- **Inverse**: rotation by -θ
- **Lie algebra**: so(2), which is isomorphic

En algèbre de Lie, **SO(2)** et **SE(2)** sont deux groupes de Lie fondamentaux pour les transformations dans le plan 2D :

## **SO(2) - Special Orthogonal Group**

**SO(2)** représente le groupe des **rotations** dans le plan 2D.

- **Définition mathématique** : Matrices 2×2 orthogonales de déterminant 1
- **Forme matricielle** :
  ```
  R(θ) = [cos θ  -sin θ]
         [sin θ   cos θ]
  ```
- **Paramètre** : Un seul angle θ ∈ [0, 2π)
- **Dimension** : 1 (une seule dimension de liberté)
- **Algèbre de Lie so(2)** : Matrices antisymétriques 2×2
  ```
  ξ = [0   -θ]
      [θ    0]
  ```

**Applications** :
- Rotations d'objets 2D
- Orientation d'un robot mobile
- Rotations de caméras autour de l'axe optique

## **SE(2) - Special Euclidean Group**

**SE(2)** représente le groupe des **transformations rigides** dans le plan (rotation + translation).

- **Définition mathématique** : Transformations qui préservent les distances et les angles
- **Forme matricielle homogène** :
  ```
  T = [R(θ)  t]  = [cos θ  -sin θ  tx]
      [0     1]    [sin θ   cos θ  ty]
                   [0      0      1 ]
  ```
- **Paramètres** : θ (rotation) + (tx, ty) (translation)
- **Dimension** : 3 (trois degrés de liberté)
- **Algèbre de Lie se(2)** : 
  ```
  ξ = [ω×  ρ]  = [0   -θ  tx]
      [0   0]    [θ    0  ty]
                 [0    0   0]
  ```

**Applications** :
- Position et orientation d'un robot mobile
- Transformations d'images 2D
- Mouvements planaires en robotique

## **Relation entre SO(2) et SE(2)**

- **SO(2) ⊂ SE(2)** : SO(2) est un sous-groupe de SE(2)
- **SE(2) = SO(2) ⋉ ℝ²** : SE(2) est le produit semi-direct de SO(2) et des translations ℝ²

## **Dans votre code**

Dans le contexte de votre plugin Cosserat (qui traite des poutres déformables), ces groupes sont utilisés pour :
- **SO(2)** : Représenter les rotations de sections transversales
- **SE(2)** : Représenter les transformations complètes (position + orientation) des sections

C'est pourquoi ce code définit des classes `SO2<Scalar>` et `SE2<Scalar>` qui héritent de `LieGroupBase` pour implémenter ces structures mathématiques avec les opérations appropriées (exp, log, composition, etc.).