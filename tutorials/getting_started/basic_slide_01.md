## **Introduction to SOFA**

- Have SOFA installed on your machine
- Install Cosserat plugin
  - In Tree
  - Out Tree

---

## **Step 1: Installing SOFA**

Before you begin with the specific Cosserat plugin, you need to install SOFA. Follow these steps:

1. Go to the official SOFA website (https://www.sofa-framework.org/) to download the latest version.
2. Choose the appropriate version for your operating system (Windows, Linux, or macOS).
3. Follow the installation instructions for your OS. Typically, this involves extracting the downloaded archive and setting environment variables.

---

## **Step 2: Setting Up the Cosserat Plugin**

Now, we'll dive into the essential part – configuring the Cosserat plugin within SOFA.

1. **Create plugins folder:**
   - Create folder externalPlugins
   - **sofa**
     - ├── **src**
     - ├── **build**
     - ├── **externalPlugins**

---

2. **Obtaining the Plugin:**

- GitHub : https://github.com/SofaDefrost/Cosserat
- Download the plugin :
  - git clone git@github.com:SofaDefrost/Cosserat.git (if you are using ssh-key)
  - git clone https://github.com/SofaDefrost/Cosserat.git
  - or Download the **Zip**

---

**3. Add _CMakeList.txt_ file inside the _externalPlugin_ folder**

```Cmake
	cmake_minimum_required(VERSION 3.1)
	sofa_find_package(SofaFramework)

	sofa_add_subdirectory(plugin SofaPython3 SofaPython3 ON) # If you want to use python
	sofa_add_subdirectory(plugin STLIB STLIB ON) # If you want to use python & Cosserat prefabs
	sofa_add_subdirectory(plugin Cosserat Cosserat ON)
```

---

**4. Activating the Plugin:** To activate the Cosserat plugin, follow these steps:

- Open your terminal and go to SOFA's **build-directory**
  - run
  ```bash
      cmake-gui .
  ```
  - In the _Search_ bar, type **external**,
  - In $SOFA\_EXTERNAL\_DIRECTORIES$
    - Fill in the empty box with:
      - **path-to-cosserat-directory**
  - Find the Cosserat plugin and enable it

---

5. **First Cosserat Scene: `tuto_1.py`**
   - As said previously, this component is based on the PCS (Piece-wise Constant Strain) formulation.
     ![](../images/Pasted%20image%2020231102173536.png)

---

## **Goals**:

- how to create a basic scene with the cosserat plugin
  - It is important to note the difference between :
    - **section** and **frames**
    - **section** and **cross-section**
- The notion of force-field : here **BeamHookeLawForceField**
- The notion of mapping: here **DiscreteCosseratMapping**
- Functions: **apply, applyJ**, **applyJT** for forces and **ApplyJ^T** for constraints

---

#### Start with the base

![600](../../docs/images/exemple_rigid_translation.png)>

---

[[./tutorial_00_basic_beam.py]]

---
