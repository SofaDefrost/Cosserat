
**Title:** Tutorial on Combining Deformable Models for Soft Robot Simulation

**Duration:** This tutorial can vary in duration based on the depth and complexity of the material, but a typical plan could span 2-4 hours.

**Target Audience:** Researchers, engineers, and students interested in soft robotics and the simulation of deformable robots.

**Prerequisites:** Participants should have a basic understanding of robotics, mechanics, and computer simulations. Knowledge of finite element methods (FEM) and rigid body dynamics can be helpful.

**Tutorial Outline:**

**1. [[../docs/text/Introduction]]**
   - Briefly introduce the field of soft robotics.
   - Explain the motivation for combining FEM and DCM models.

**2. [[../docs/text/Background Concepts]]**
   - An overview of Discrete Cosserat Model (DCM) and their relevance in soft robotics.
   - Finite Element Methods (FEM) and their use in simulating deformable structures ?
   - The concept of compliance matrices and their role in soft robot modeling ?

**3. [[../docs/text/Setting up the Environment]]**
   - The software tools and libraries used (SOFA, **Cosserat** plugin).
   - Instructions for setting up the simulation environment.

**4. [[../docs/text/Direct Simulation]]**
   - The direct simulation of a soft robot (e.g., the gripper example from the paper).
   - How to create the DCM model for cables and integrate it with FEM for the soft body ?
   - **Discussion** : the role of constraints, friction, and contact forces.

---
**5. [[../docs/text/Inverse Simulation]]**
   - Demonstrate an inverse simulation for soft robots (e.g., using the soft finger or tentacle as examples).
   - Explain how to compute actuation values to achieve desired end-effector positions.
   - Discuss the handling of sliding constraints in the inverse problem.

---
**6. [[../docs/text/Model Convergence and Sensitivity Analysis]]**
   - Here we will share insights on how to evaluate the convergence of your model with varying parameters 
	   - e.g., number of sections defining a cable
   - Discuss sensitivity to material properties and other simulation parameters.

---

**7. [[../docs/text/Advanced Topics and Future Work]]**
   - Explore potential applications of your method in soft robot design.
   - Discuss the ongoing and future research in this field, including real-world experiments and stress validation on cables.
---

**8. Performance Optimization**
   - How to **optimize the computation speed**, potentially using **parallelization**, model **order reduction**, or **recursive algorithms** for DCM.
---
**9. Conclusion and Q&A**
   - Summarize key takeaways from the tutorial.
   - Open the floor for questions and discussions.
---
**10. Practical Hands-On Session (Optional)**
   - If feasible, you can provide participants with exercises or hands-on practice to apply the concepts learned.

---
**Materials:**
- Presentation slides (see GitHub Repo)
- Repo folder tutorials.
- Code examples and simulations on the Repo/example.
- Reference: 
	- Paper RAL-SoRo :
	- Paper (Féderico) :
	- Paper (Flavie) :  
- Q&A session for participant engagement.
- Feedback : direct and mails
---

