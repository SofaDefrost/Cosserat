cmake_minimum_required(VERSION 3.12)

project(Cosserat_SRC)

set(HEADER_FILES
    config.h.in
    fwd.h
    types.h
    initCosserat.cpp
    engine/GeometricStiffnessEngine.h
    engine/GeometricStiffnessEngine.inl
    forcefield/BeamHookeLawForceField.h
    forcefield/BeamHookeLawForceField.inl
    mapping/CosseratNonLinearMapping2D.h
    mapping/CosseratNonLinearMapping2D.inl
    mapping/DifferenceMultiMapping.h
    mapping/DifferenceMultiMapping.inl
    constraint/UnilateralPlaneConstraint.h
    constraint/UnilateralPlaneConstraint.inl
)

set(SOURCE_FILES
    initCosserat.cpp
    engine/GeometricStiffnessEngine.cpp
    forcefield/BeamHookeLawForceField.cpp
    mapping/CosseratNonLinearMapping2D.cpp
    mapping/DifferenceMultiMapping.cpp
    constraint/UnilateralPlaneConstraint.cpp
)

find_package(SofaFramework REQUIRED)
find_package(SofaBase REQUIRED)
find_package(SofaOpenglVisual REQUIRED)
find_package(Geometric_stiffness REQUIRED)
find_package(Eigen3 REQUIRED)

add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} PUBLIC
    SofaBase
    SofaOpenglVisual
    Geometric_stiffness
)
target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>")
target_include_directories(${PROJECT_NAME} PUBLIC ${EIGEN3_INCLUDE_DIR})

set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "-DSOFA_BUILD_COSSERAT")
