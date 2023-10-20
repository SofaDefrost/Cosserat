# Cosserat
<div class="highlight-content">
<strong>Overview</strong> 

An open-source plugin, designed to be compatible with the Sofa framework, facilitates the simulation of 1D objects. 
Specifically, it caters to the modeling of both rigid and flexible 1D entities, like rods, wires or needles, using the Cosserat beam theory. 
In this context, we have outlined a range of potential applications for this plugin. If you wish to explore its functionality, you have the flexibility to 
construct scenes using Python or XML, or you can take it a step further by developing new C++ components. 
We also welcome contributions from the community to enhance and expand the capabilities of this plugin.

</div>

## Description related to Soft-body modeling

The Cosserat model has found applications in the realm of continuum robotics, particularly for simulating the deformation of robot bodies with geometries and mechanical properties akin to rods. 
This model aligns closely with the dynamic deformation patterns exhibited by soft manipulators, as it can effectively replicate nonlinear deformations encompassing bending, torsion, extension, and shearing.

One distinctive feature of Cosserat's theory, within the domain of continuous media mechanics, lies in its conceptualization: 
it views each material point of an object as a rigid body with six degrees of freedom (three translations and three rotations). 
In contrast, many other models in continuum media mechanics tend to treat material points as particles with only three translation degrees of freedom.

When modeling linear structures, this framework enables the creation of a structure closely resembling articulated solids, consisting of a series of rigid bodies whose relative positions are defined by their strain states. 
Consequently, this model serves as a versatile tool for modeling and controlling a variety of systems, including concentric tube robots, continuum robots driven by cables, or pneumatic soft robots with constant cross-sections.

## Theory
Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et.

## Numerics
Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. 

## Some use cases
Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. 

### Modeling and control
Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. 

#### Direct control
Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. 

#### Modeling cochlear implant using Discret Cosserat Model (DCM)


| View 1                                                                                       | View 2                                                                                       | View 3                                                                                       |
|----------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------|
| ![333](docs/images/multiSectionWithColorMap1.png) | ![333](docs/images/multiSectionWithColorMap2.png) | ![333](docs/images/multiSectionWithColorMap3.png) |


## Utilizing the Discrete Cosserat Model for Cable Modeling in Deformable Robot Control:


| Direct simulation of a soft gripper                                                       | The study of the model convergence                                                            |
|-------------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------------|
| ![400](docs/images/cosseratgripper_2.png) | ![400](docs/images/tenCossseratSections.png) |


---

<strong> Actuation </strong>

| <img src="/docs/images/actuationConstraint_2.png" width="300" title="DCM Beam actuation using a given cable"> | <img src="docs/images/circleActuationConstraint.png" width="300" title="DCM Beam actuation using a given cable"> | <img src="/docs/images/actuationConstraint_1.png" width="300" title="DCM Beam actuation using a cable"> |
|---------------------------------------------------------------------------------------------------------------|:----------------------------------------------------------------------------------------------------------------:|--------------------------------------------------------------------------------------------------------:|
| DCM Beam actuation using a cable ```d =```                                                                    |                                    DCM Beam actuation using a cable ```d =```                                    |                                                                  Beam actuation using a cable ```d =``` |


| <img src="/docs/images/example1.gif" width="300" title="PCS_Example1.py "> | <img src="./docs/images/example2.gif" widt="300" title="PCS_Example2.py"> | <img src="./docs/images/example2.gif" width="300" title="PCS_Example3.PCS"> |
|----------------------------------------------------------------------------|:-------------------------------------------------------------------------:|----------------------------------------------------------------------------:|
| DCM Beam actuation using a cable ```d =```                                 |                DCM Beam actuation using a cable ```d =```                 |                                      Beam actuation using a cable ```d =``` |

<strong> Tripod using bending lab sensors </strong>
Format: ![Alt Text]


### Inverse Control


## Publications
1. Pieces-wise constant Strain PCS: This feature is based on the paper
- __Discrete cosserat approach for soft robot dynamics: A new piece-wise constant strain model with torsion and shears__ [Link to the paper](https://ieeexplore.ieee.org/document/7759808)
- __Coupling numerical deformable models in global and reduced coordinates for the simulation of the direct and the inverse kinematics of Soft Robots__ [Link to the paper](https://ieeexplore.ieee.org/abstract/document/9362217).
The link to download the Pdf of the Paper: __https://hal.archives-ouvertes.fr/hal-03192168/document__

<div align="center">
  <a href="https://www.youtube.com/watch?v=qwzKAgw31pU"><img src="https://img.youtube.com/vi/qwzKAgw31pU/0.jpg" alt="link to youtube"></a>
</div>

2. Pieces-wise Non-constant Strain:
3. DCM with Plastic model

