---
tutorial: Introduction
tags:
  - tuto/cosserat
cssclass:
  - dashboard
---
_Welcome to this tutorial on SOFA-Cosserat Plugin._  
  
[Formation plugin : Cosserat - CodiMD](https://notes.inria.fr/gcfuFPDYSeeAlG4gzfDJwA#)

Yinoussa Adagolodjo

---
#### SOFA

- SOFA : Simulation Open Framework Architecture
	- Dedicated to research
- Use for prototyping and development of physics-based simulation
	- FEM / Rigid (articulated) Bodies / Contacts / Other models...
	- Projections = MAPPING !
---
#### SOFA
- Free & Open Source: [https://www.sofa-framework.org/](https://www.sofa-framework.org/)
	- Download, use, modify, cite and contribute to SOFA!
- International community and two major events per year 
- SOFA Week from $13-17^{th}$  of November 
	- Onsite with online access
	- Free – all you need is to register
---
#### Tutorial Roadmap

- What to expect ?
	- An overview of the theory of Discrete Cosserat Model
	- Numerics & Hands-on examples
	    - Reduce coordinates
	    - Reduce coordinate (Cosserat) to Global coordinate (SOFA) State 
	    - Boundary conditions and interaction forces
- Ask your questions and actively participate throughout this tutorial.
	- The first time I am doing this, so I need your feedback
---
#### Introduction to Soft Robotics (SR)

- Soft robotics is an emerging and innovative field of robotics 
- Focuses on the design and development of robots made from :  
	- *Flexible*, 
	- *Deformable*,
	- *Compliant materials*.
---
#### Introduction to Soft Robotics (SR)
- Numerous advantages
	- *Adaptability* : Their ability to deform and adapt to their environment
	- *Safety* : Ideal for interactions with humans, delicate objects, or unstructured surroundings
	- *Versatility* : Various applications

---
####  **Key Applications**

- *Healthcare* : The gentle and non-invasive nature of SR is necessary
- *Industrial automation* : The high [Compliance](../../../Soft_Robot/Compliance.md), reduces the risk of damage during interactions with products
- *Search and rescue* : they can navigate through tight spaces and uneven terrains
- *Space exploration*, and *extreme environments*

---
#### **Challenges in Soft Robotics**

- *Modeling* : due to deformable materials
- [*Control*](../../../Soft_Robot/kinematics_dynamics_control.md) : due to the non-linear, multi-body dynamics of deformable materials
- [*Multi-dimension*](../../../../../../Projects/ANR/brainStormingANR2023.md) : due to a wide range of shape, volume (3D), surface (2D) and cable (1D)
- [*Multi-physics*](../../../../../../Projects/ANR/RobotMulti-physics.md) : due a wide range of physical behaviors, including mechanical deformation, thermal effects, electrical responses,

---
#### **Challenges in Soft Robotics**

- Addressing these challenges is crucial for the widespread adoption of soft robots in various fields.
(To go further on the introduction of deformable robotics ⇾ [Introduction_General](../../../Soft_Robot/Introduction_General.md) )

---
##### Why combine different models ?
- Introducing a versatile modeling tool for multidimensional Soft-Robot:
	- Animal locomotion (Bio-inspired robotics)
	- Flexible arms actuated by cables or similar mechanisms (Robotics)
	- Medical devices like endoscopes, needles, and catheters (Medical robotics)
---
##### Why combine different models ?
- Existing design and modeling tools demand enhancement, particularly in the mechanical aspect.
![](../images/Pasted%20image%2020231108201224.png)
- Offer a wide range of design possibilities:
	- Selection of actuators, force transmission within the structure, utilization of materials with varying stiffness profiles, and more.
---
##### Cosserat Theory
Choose strain as generalized coordinates, defined in global (or local) frame!
![400](../images/Pasted%20image%2020231108233708.png)
[Lazarus et al. 2013]
- The configuration of the Cosserat rod is defined by its centerline r(s).
- The orientation of each mass point of the rod is represented by an orthonormal basis ($d_1(s), d_2(s), d_3(s)$)
- The three local modes of deformation of the elastic rod : - (b1) material curvature $κ_1$ related to the direction $d_1$ of the cross-section, (b2) material curvature $κ_2$ related to $d_2$, and $κ_3$ twist.
---
##### Discrete Cosserat Modeling: DCM
- Serial rigid (6DoFs) bodies with reduced coordinates
![700](https://lh7-us.googleusercontent.com/q9X7GJY2GxkPSOqb7w_jzbsf5sIiUglQTnJDySqUer-mVAPQPr-ENDjkSMxFlB0LuyXX0DKcuMR3rKqdMmWGJiBBoXu9zM7sbXbgCrZZhKiD0mkcY7Pwru_J7JvyVSCD6o4cYXsV7L65TSJprRY_=nw)
---
##### Discrete Cosserat Modeling: DCM
- Piece-wise Constant Strain (PCS, treats rigid, soft, or hybrid robots uniformly)
![500](https://lh7-us.googleusercontent.com/hzJA2pzS-naPfFkf-98bPkGsQ86ZGdGwqGW3-un56s3ZcVkdOB0_Jus4a9W_nqO7jU7Tt_FDzCrFIbfA9XFqaPBmmq-do-TIJkFn6NX-RimX2UlWBTis_7bKzAp7fEmIeiuOZ1FueZ5yxijFSJls=nw)
---
##### Discrete Cosserat Modeling: DCM
- Models the deformation of a soft manipulator arm using a finite number of sections
- Assumes constant strains along each section
- Accounts for shears and torsion
- Simulates the inextensible behavior of a rod or cable
- Reduces model size, resulting in faster calculation times through the use of reduced coordinates.
---
##### DCM (Kinematics)
- Configuration $g= \begin{pmatrix}  \mathcal{R} & u \\  0 & 1  \end{pmatrix} ∈SE(3)$
- Velocity $\begin{align}\eta(s,t) &= \begin{bmatrix}\mathcal{w} \\ \mathcal{v} \end{bmatrix}\end{align} \in R^6$
- Strain $\begin{align}\xi(s,t) &= \begin{bmatrix}\mathcal{k} \\ \mathcal{p} \end{bmatrix}\end{align} \in R^6$
- Kinematics : $g'=g\hat{\xi}(X)$ ; $\dot{g} = g\hat{\eta}$
- Differntial Kinematics : $\eta'=\dot{\xi}(X)-ad_{\xi(X)}\eta(X)$
![](../images/Pasted%20image%2020231109002926.png)
---
![0](../images/Pasted%20image%2020231108234643.png)

---
##### DCM (Dynamics)
![800](../images/Pasted%20image%2020231109003349.png)

---
##### DCM Dynamic
![](../images/Pasted%20image%2020231109003734.png)

---
##### Approximation via PCS, VS and PLS
![300](../images/Pasted%20image%2020231109003934.png)
- **PCS**: A local approximation scheme employing a local constant strain assumption.
- **VS**: A global approximation method based on the chosen basis functions.
- **PLS**: A local approximation scheme with a linear strain function assumption. 
---
##### PCS Cosserat
![](../images/Pasted%20image%2020231109004616.png)

---
##### Discrete Cosserat Modeling: DCM
###### Limitations:
- Challenges in simulating truss structures, intricate geometries, or volumetric deformations
- The PCS parametrization of a manipulator is not rooted in the arm's intrinsic variables
---
##### Finite Element Modeling (FEM)
The Finite Element (FE) approach is typically represented by the position and velocity (global coordinates) of a system of interconnected nodes. 
- This approach offers several advantages:
	- **Versatility in object geometries:** It can be applied to a wide range of object geometries, including beams, shells, truss structures, and deformable volumes.
	- **Material law customization:** Material properties can be tailored to meet specific requirements.
	- **Ease of defining boundary conditions:** Boundary conditions for numerical models can be defined with ease.
	- **Flexibility for creating truss structures:** Beams can be interconnected freely to create truss structures.
---
##### Finite Element Modeling (FEM)
However, this approach also has limitations:
- **Time-consuming:** Simulations using the FE approach may be computationally intensive and time-consuming.
- **Additional constraints:** Additional constraints are often needed to prevent the extension of rod-like structures when modeling certain systems.
---
##### Modeling Soft Robots: Combining DCM with FEM

We combine Discrete Cosserat Modeling (DCM) with Finite Element Modeling (FEM) to harness the strengths of each model. This hybrid approach is particularly useful in scenarios like:
- Modeling the stiffness of cables used to actuate a soft robot with a deformable volumetric structure.
- This leads to a more realistic representation of the entire robot's behavior.
For example, it allows us to effectively simulate scenarios where inextensible cables are employed to control the motion of a soft robot.

----
  *Combined Accuracy*: By combining FEM and DCM, you can leverage the strengths of both methods. FEM provides fine-grained material modeling, while DCM captures the shape and motion of deformable elements accurately. This leads to a more realistic representation of the entire robot's behavior.

- *Unified Simulation Framework*:A combined FEM-DCM framework creates a unified simulation environment that can model both the deformable body of the robot (using FEM) and the actuation components (using DCM). This simplifies simulation setup and control algorithms.

---

[Hands on](Setting%20up%20the%20Environment.md)

---

[[Complement]]

---
##### [FEM and DCM](../../_docs/DCM_FEM.md) 

- *FEM's Material Modeling*: FEM excels at modeling the deformations and stress distributions in complex materials, including soft and deformable ones. It considers the local behavior of materials, making it more accurate for understanding the mechanical properties of soft robots.
- *Cosserat theory's Beam-Like Modeling*: DCM, on the other hand, is suitable for modeling the overall shape and bending of structures, making it a natural choice for cables, rods, and flexible elements in soft robots.
--- 