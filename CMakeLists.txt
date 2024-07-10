cmake_minimum_required(VERSION 3.12)
project(Cosserat VERSION 21.12.0)

include(cmake/environment.cmake)

find_package(Sofa.Config REQUIRED)
sofa_find_package(Sofa.Component.Constraint.Lagrangian.Model REQUIRED)
sofa_find_package(Sofa.Component.StateContainer REQUIRED)
sofa_find_package(Sofa.Component.Mapping.NonLinear REQUIRED)
sofa_find_package(Sofa.GL REQUIRED)
sofa_find_package(Sofa.Component.Topology.Container.Dynamic REQUIRED)

sofa_find_package(STLIB QUIET)
if(STLIB_FOUND)
    message("-- Found dependency : 'STLIB' plugin .")
else()
    message("-- The highly recommended 'STLIB' plugin is missing. You can compile Cosserat but some of the provided python examples will not work. ")
endif()

set(SRC_ROOT_DIR src/${PROJECT_NAME})

set(HEADER_FILES
    ${SRC_ROOT_DIR}/config.h.in
    ${SRC_ROOT_DIR}/fwd.h
    ${SRC_ROOT_DIR}/types.h
    ${SRC_ROOT_DIR}/mapping/BaseCosseratMapping.h
    ${SRC_ROOT_DIR}/mapping/BaseCosseratMapping.inl
    ${SRC_ROOT_DIR}/mapping/DiscreteCosseratMapping.h
    ${SRC_ROOT_DIR}/mapping/DiscreteCosseratMapping.inl
    ${SRC_ROOT_DIR}/mapping/DiscreteDynamicCosseratMapping.h
    ${SRC_ROOT_DIR}/mapping/DiscreteDynamicCosseratMapping.inl
    ${SRC_ROOT_DIR}/engine/ProjectionEngine.h
    ${SRC_ROOT_DIR}/engine/ProjectionEngine.inl
    ${SRC_ROOT_DIR}/mapping/DifferenceMultiMapping.h
    ${SRC_ROOT_DIR}/mapping/DifferenceMultiMapping.inl
    ${SRC_ROOT_DIR}/mapping/RigidDistanceMapping.h
    ${SRC_ROOT_DIR}/mapping/RigidDistanceMapping.inl
    ${SRC_ROOT_DIR}/engine/PointsManager.h
    ${SRC_ROOT_DIR}/engine/PointsManager.inl
    ${SRC_ROOT_DIR}/forcefield/BeamHookeLawForceField.h
    ${SRC_ROOT_DIR}/forcefield/BeamHookeLawForceField.inl
    ${SRC_ROOT_DIR}/forcefield/BeamHookeLawForceFieldRigid.h
    ${SRC_ROOT_DIR}/forcefield/BeamHookeLawForceFieldRigid.inl
    ${SRC_ROOT_DIR}/forcefield/CosseratInternalActuation.h
    ${SRC_ROOT_DIR}/forcefield/CosseratInternalActuation.inl
    ${SRC_ROOT_DIR}/constraint/CosseratSlidingConstraint.h
    ${SRC_ROOT_DIR}/constraint/CosseratSlidingConstraint.inl
    ${SRC_ROOT_DIR}/mapping/LegendrePolynomialsMapping.h
    ${SRC_ROOT_DIR}/mapping/LegendrePolynomialsMapping.inl
    ${SRC_ROOT_DIR}/constraint/CosseratNeedleSlidingConstraint.h
    ${SRC_ROOT_DIR}/constraint/CosseratNeedleSlidingConstraint.inl
    )
set(SOURCE_FILES
    ${SRC_ROOT_DIR}/initCosserat.cpp
    ${SRC_ROOT_DIR}/mapping/BaseCosseratMapping.cpp
    ${SRC_ROOT_DIR}/mapping/DiscreteCosseratMapping.cpp
    ${SRC_ROOT_DIR}/mapping/DiscreteDynamicCosseratMapping.cpp
    ${SRC_ROOT_DIR}/engine/ProjectionEngine.cpp
    ${SRC_ROOT_DIR}/mapping/DifferenceMultiMapping.cpp
    ${SRC_ROOT_DIR}/mapping/RigidDistanceMapping.cpp
    ${SRC_ROOT_DIR}/engine/PointsManager.cpp
    ${SRC_ROOT_DIR}/forcefield/BeamHookeLawForceField.cpp
    ${SRC_ROOT_DIR}/forcefield/BeamHookeLawForceFieldRigid.cpp
    ${SRC_ROOT_DIR}/forcefield/CosseratInternalActuation.cpp
    ${SRC_ROOT_DIR}/constraint/CosseratSlidingConstraint.cpp
    ${SRC_ROOT_DIR}/mapping/LegendrePolynomialsMapping.cpp
    ${SRC_ROOT_DIR}/constraint/CosseratNeedleSlidingConstraint.cpp
    )

sofa_find_package(SoftRobots QUIET)
if(SoftRobots_FOUND)
    list(APPEND HEADER_FILES
        ${SRC_ROOT_DIR}/constraint/CosseratActuatorConstraint.h
        ${SRC_ROOT_DIR}/constraint/CosseratActuatorConstraint.inl
        
        ${SRC_ROOT_DIR}/constraint/QPSlidingConstraint.h
        ${SRC_ROOT_DIR}/constraint/QPSlidingConstraint.inl
        )
    list(APPEND SOURCE_FILES
        ${SRC_ROOT_DIR}/constraint/CosseratActuatorConstraint.cpp
        
        ${SRC_ROOT_DIR}/constraint/QPSlidingConstraint.cpp
        )
endif()

file(GLOB_RECURSE RESOURCE_FILES  "*.md" "*.psl" "*.py" "*.pyscn" "*.py3scn" "*.scn" "*.ah")

if(WIN32)
    add_definitions(-D_WINSOCKAPI_)
endif()

add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES} ${DOC_FILES} )
target_link_libraries(${PROJECT_NAME}
    Sofa.Component.Constraint.Lagrangian.Model
    Sofa.Component.StateContainer
    Sofa.Component.Mapping.NonLinear
    Sofa.GL
    Sofa.Component.Topology.Container.Dynamic
)

if(Sofa.GL_FOUND)
    target_link_libraries(${PROJECT_NAME} Sofa.GL)
endif()

if(SoftRobots_FOUND)
    target_link_libraries(${PROJECT_NAME} SoftRobots)
    message("-- Found dependency: 'SoftRobots' plugin.")
endif()

sofa_install_pythonscripts(PLUGIN_NAME ${PROJECT_NAME} PYTHONSCRIPTS_SOURCE_DIR "python")

find_file(SofaPython3Tools NAMES "SofaPython3/lib/cmake/SofaPython3/SofaPython3Tools.cmake")
if(SofaPython3Tools)
    message("-- Found SofaPython3Tools.")
    include(${SofaPython3Tools})
else()
    # try again with the find_package mechanism
    find_package(SofaPython3 QUIET)
endif()
if(SofaPython3Tools OR SofaPython3_FOUND)
    if(NOT SP3_PYTHON_PACKAGES_DIRECTORY)
        set(SP3_PYTHON_PACKAGES_DIRECTORY "python3/site-packages")
    endif()
    message("-- Python3 packages will be installed in ${SP3_PYTHON_PACKAGES_DIRECTORY}.")
    # SP3_add_python_package(
    #         SOURCE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python3/cosserat
    #         TARGET_DIRECTORY cosserat
    # )
    add_subdirectory(${SRC_ROOT_DIR}/Binding)
endif()



## Install rules for the library and headers; CMake package configurations files
sofa_create_package_with_targets(
    PACKAGE_NAME ${PROJECT_NAME}
    PACKAGE_VERSION ${PROJECT_VERSION}
    TARGETS ${PROJECT_NAME} AUTO_SET_TARGET_PROPERTIES
    INCLUDE_SOURCE_DIR "src"
    INCLUDE_INSTALL_DIR ${PROJECT_NAME}
    EXAMPLE_INSTALL_DIR "examples"
    RELOCATABLE "plugins"
    )


# Tests
# If SOFA_BUILD_TESTS exists and is OFF, then these tests will be auto-disabled
cmake_dependent_option(COSSERAT_BUILD_TESTS "Compile the tests" ON "SOFA_BUILD_TESTS OR NOT DEFINED SOFA_BUILD_TESTS" OFF)
if(COSSERAT_BUILD_TESTS)
    add_subdirectory(tests)
endif()

# Config files and install rules for pythons scripts
sofa_install_pythonscripts(PLUGIN_NAME ${PROJECT_NAME} PYTHONSCRIPTS_SOURCE_DIR "examples/python3/")

include(cmake/packaging.cmake)
