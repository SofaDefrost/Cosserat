project(Cosserat.src VERSION 21.06.99)

set(HEADER_FILES
    initCosserat.h

    mapping/BaseCosserat.h
    mapping/BaseCosserat.inl
    mapping/DiscreteCosseratMapping.h
    mapping/DiscreteCosseratMapping.inl
    mapping/DiscretDynamicCosseratMapping.h
    mapping/DiscreteDynamicCosseratMapping.inl
    mapping/ProjectionEngine.h
    mapping/ProjectionEngine.inl
    mapping/DifferenceMultiMapping.h
    mapping/DifferenceMultiMapping.inl
    mapping/RigidDistanceMapping.h
    mapping/RigidDistanceMapping.inl
    mapping/LegendrePolynomialsMapping.h
    mapping/LegendrePolynomialsMapping.inl


    forcefield/BeamHookeLawForceField.h
    forcefield/BeamHookeLawForceField.inl
    forcefield/CosseratInternalActuation.h
    forcefield/CosseratInternalActuation.inl
    forcefield/MyUniformVelocityDampingForceField.h
    forcefield/MyUniformVelocityDampingForceField.inl

    constraint/CosseratSlidingConstraint.h
    constraint/CosseratSlidingConstraint.inl
    constraint/QPSlidingConstraint.h
    constraint/QPSlidingConstraint.inl
    constraint/CosseratNeedleSlidingConstraint.h
    constraint/CosseratNeedleSlidingConstraint.inl
    )

set(SOURCE_FILES
    initCosserat.cpp

    mapping/BaseCosserat.cpp
    mapping/ProjectionEngine.cpp
        mapping/DiscreteCosseratMapping.cpp
    mapping/DiscreteDynamicCosseratMapping.cpp
    mapping/DifferenceMultiMapping.cpp
    mapping/RigidDistanceMapping.cpp
    mapping/LegendrePolynomialsMapping.cpp

    forcefield/BeamHookeLawForceField.cpp
    forcefield/CosseratInternalActuation.cpp
    forcefield/MyUniformVelocityDampingForceField.cpp

    constraint/CosseratSlidingConstraint.cpp
    constraint/QPSlidingConstraint.cpp
    constraint/CosseratNeedleSlidingConstraint.cpp
    )

# add_subdirectory(Binding)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})

target_include_directories(${PROJECT_NAME} PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../>")

target_link_libraries(${PROJECT_NAME} SofaTest SofaGTestMain SofaCore
    SofaConstraint SofaBaseMechanics SofaUserInteraction)
