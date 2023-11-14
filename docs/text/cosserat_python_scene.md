# Principes de base

<aside> 💡 _Toutes les informations concernant le logiciel SOFA peuvent être retrouver sur le site : [SOFA – Documentation (sofa-framework.org)](https://www.sofa-framework.org/community/doc/), pour y comprendre le concept, les principes, ... ._
</aside>

## Théorie de Cosserat

Cosserat sur base sur la notion de _Piecewise Constant Strain_. Elle repose sur le fait d’estimer la **deformation** de manière très grossière, comme ci-dessous.

![PCS.png](images/PCS.png)

### Base rigide

![rigidbase.png](images/rigidbase.png)

Comme son nom l’indique ici il s’agit de la base de la poutre. Il s’agit d’un repère qui pouvant être assigné à une translation et/ou rotations désirée, et qui résultera d’une position finale. A noter, que le repère est toujours associé à un ressort afin d’être fixe dans le repère global de SOFA.

$$ Rotation =[x , y, z] $$

$$ Translation = [x,y,z] $$

$$ Position = [x , y, z, xq, yq, zq, wq] $$

### Exemple - Base rigide

![exemple rigid translation.png](images/exemple_rigid_translation.png)

$$ Rotation =[0 , 0, 0] $$

$$ Translation = [0,1,0] $$

$$ Position = [0 , 1, 0, 0, 0, 0, 1] $$

### Section

![noeud curv.png](images/noeud_curv.png)

Les sections sont situées dans le repère local de Cosserat. Une section représente la distance entre deux nœuds :

$$ L0 /L1 = Section 1 $$

$$ L1/L2 = Section 2 $$

Une cross-section est une partie de la poutre qu’on vient couper de manière transversale. De manière générale, on suppose que les cross sections sont rondes et ont le même rayon tout le long de la poutre. Mais il est possible que la forme de la poutre ne soit pas circulaire mais rectangulaire, et que le rayon change d’une distance à une autre. Les coordonnées des sections et des noeuds se retrouvent dans un tableau nommé “curv_abs_input” (:Curviligne abscisse inputs).

<aside> ❗ _**Attention à ne pas confondre : section et cross-section. Ceux ne sont pas les mêmes notions**_

</aside>

Il est possible d’appliquer une déformation à partir d’un certain nœud, les deux déformations possibles sont la torsion (selon x), et la flexion (selon y ou z).

$$ positionS = [x,y,z] $$

$$ Flexion = -+(1/LenghtSection) * pi $$

Pi : Angle de rotation

(1/LengthSection) : Vecteur unitaire (normalisation)

-/+ : Sens de rotation

### Exemple - Section (1)

Pour une poutre de longueur de 8cm avec un nombre de sections de 6 on a :

![exemple curv abs input.png](images/exemple_curv_abs_input.png)

**Tableau curv_abs_input :**

![Untitled](images/Untitled.png)

### Exemple - Section (2)

Pour une poutre de longueur de 8cm avec un nombre de sections de 6 on a :

![position 90degres.png](images/position_90degres.png)

**Tableau curv_abs_input :**

||X|Y|Z|
|---|---|---|---|
|0|0|0|0|
|1|0|0|0|
|2|0|0|0|
|3|0|0|0|
|4|0|0|0|
|5|0|0|-2.35|

### Frame

![frame1.png](images/frame1.png)

La frame est représenté dans le repère locale de Cosserat, puis elle est intégrée dans le repère global de SOFA. Une frame est représentée comme les bases rigides, on retrouve les coordonnées de translation et rotation :

$$ R :[x , y, z, 0, 0, 0, 1] $$

Les frames ne sont pas obligatoirement à équidistances les unes par rapport aux autres. Il est possible de concentrer un certain nombre de frames à un certain endroit, ce qui permet une rapidité de calcul sur l’endroit d’étude désirée.

<aside> 💡 Plus le nombre de frame augmentent plus la précision de la courbe sera haute. ******

</aside>

## Plugin Cosserat - SOFA

### RigidBases

- MechanicalObject
    
    - Coordonnées du système :
    
    $$ R :[x , y, z, 0, 0, 0, 1] $$
    
- Spring
    - Rigidité
    - Rigidité angulaire

### CosseratCoordinate

- MechanicalObject
    - Nombre de sections
    - Taille totale de la poutre
- BeamHookLaw
    - Forme de la cross-section
    - Taille de la section
    - Rayon

### CosseratFrame

- MechanicalObject
    - Nombre de frame
    - Taille totale de la poutre
- CosseratMapping
    - MechanicalObject de RigidBases
    - MechanicalObject de Cosserat Coordinate
    - Curv in
    - Curv out