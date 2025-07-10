# Dependency Tree of Liegroups Folder

```mermaid
graph BT
    %% Styling
    classDef core fill:#d4f1f9,stroke:#333,stroke-width:1px
    classDef implementation fill:#d5f5e3,stroke:#333,stroke-width:1px
    classDef wrapper fill:#fadbd8,stroke:#333,stroke-width:1px
    classDef docs fill:#fef9e7,stroke:#333,stroke-width:1px

    %% Core Components
    subgraph Core["Core Components"]
        LieGroupBase["LieGroupBase.h/inl"]
        Types["Types.h"]
        Utils["Utils.h"]
    end

    %% Group Implementations
    subgraph Implementations["Group Implementations"]
        RealSpace["RealSpace.h"]
        SO2["SO2.h"]
        SO3["SO3.h"]
        SE2["SE2.h"]
        SE3["SE3.h"]
        SE23["SE23.h"]
        SGal3["SGal3.h"]
    end

    %% High-Level Wrapper
    subgraph Wrapper["High-Level Wrapper"]
        Bundle["Bundle.h"]
    end

    %% Documentation
    subgraph Documentation["Documentation"]
        Docs["Readme.md<br/>USAGE.md<br/>docs/<br/>tasks.md"]
    end

    %% Dependencies
    RealSpace --> LieGroupBase
    SO2 --> LieGroupBase
    SO3 --> LieGroupBase
    SE2 --> LieGroupBase
    SE2 --> SO2
    SE3 --> LieGroupBase
    SE3 --> SO3
    SE23 --> LieGroupBase
    SE23 --> SE3
    SGal3 --> LieGroupBase

    %% Bundle dependencies
    Bundle --> LieGroupBase
    Bundle --> Types

    %% Apply styling
    class LieGroupBase,Types,Utils core
    class RealSpace,SO2,SO3,SE2,SE3,SE23,SGal3 implementation
    class Bundle wrapper
    class Docs docs
```

## Explanation of the Dependency Tree

The diagram above illustrates the dependencies between different files in the Liegroups folder:

1. **Core Components**:
   - `LieGroupBase.h/inl`: Base class for all Lie group implementations
   - `Types.h`: Contains common type definitions
   - `Utils.h`: Provides utility functions

2. **Group Implementations** (all depend on LieGroupBase):
   - `RealSpace.h`: Simplest implementation 
   - `SO2.h`: 2D rotation group
   - `SO3.h`: 3D rotation group
   - `SE2.h`: 2D rigid transformation (depends on SO2)
   - `SE3.h`: 3D rigid transformation (depends on SO3)
   - `SE23.h`: 3D rigid transformation with scale (depends on SE3)
   - `SGal3.h`: 3D similarity transformation

3. **High-Level Wrapper**:
   - `Bundle.h`: Combines multiple Lie groups using std::tuple and type_traits

4. **Documentation**:
   - Various documentation files explaining usage and implementation details

The arrows represent "depends on" relationships, showing which files include other files.

