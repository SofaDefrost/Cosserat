cmake_minimum_required(VERSION 3.12)
project(CosseratPlugin VERSION 21.12.0)

# set (COSSERAT_VERSION ${PROJECT_VERSION})
include(cmake/environment.cmake)

find_package(SofaFramework REQUIRED)
find_package(SofaUserInteraction REQUIRED)
find_package(SofaBaseMechanics REQUIRED)
find_package(SofaConstraint REQUIRED)
find_package(SofaRigid REQUIRED)
find_package(SofaMiscMapping REQUIRED)
find_package(SofaEngine REQUIRED)
find_package(SofaSparseSolver REQUIRED)
find_package(SofaBaseLinearSolver REQUIRED)

sofa_find_package(Sofa.GL QUIET)
sofa_find_package(SofaBaseTopology REQUIRED)

# option(COSSERAT_IGNORE_ERRORS "Enable this option to ignore the recommendations and generate anyway." ON)

# if(NOT COSSERAT_IGNORE_ERRORS)
#     if(NOT PLUGIN_SOFAPYTHON AND NOT PLUGIN_SOFAPYTHON3)
#         message(SEND_ERROR "
#             You should not use the plugin Cosserat without enabling PLUGIN_SOFAPYTHON or PLUGIN_SOFAPYTHON3.
#             Since all the scenes to test this plugin are using one of this two plugins. 
#             To fix this error you can either follow the recommendation or enable COSSERAT_IGNORE_ERRORS.")
#     endif()    
# endif()

# set(COSSERAT_PYTHON ${SofaPython_FOUND}) # config.h.in

sofa_find_package(STLIB QUIET)
if(STLIB_FOUND)
    message("-- Found dependency : 'STLIB' plugin .")
else()
    message("-- The highly recommended 'STLIB' plugin is missing. You can compile Cosserat but some of the provided python examples will not work. ")
endif()

set(HEADER_FILES
    src/initCosserat.h
    src/mapping/BaseCosserat.h
    src/mapping/BaseCosserat.inl
    src/mapping/DiscreteCosseratMapping.h
    src/mapping/DiscreteCosseratMapping.inl
    src/mapping/DiscreteDynamicCosseratMapping.h
    src/mapping/DiscreteDynamicCosseratMapping.inl
    src/mapping/ProjectionEngine.h
    src/mapping/ProjectionEngine.inl
    src/mapping/DifferenceMultiMapping.h
    src/mapping/DifferenceMultiMapping.inl
    src/mapping/RigidDistanceMapping.h
    src/mapping/RigidDistanceMapping.inl
    src/forcefield/BeamHookeLawForceField.h
    src/forcefield/BeamHookeLawForceField.inl
    src/forcefield/BeamPlasticLawForceField.h
    src/forcefield/BeamPlasticLawForceField.inl
    src/forcefield/CosseratInternalActuation.h
    src/forcefield/CosseratInternalActuation.inl
    src/constraint/CosseratSlidingConstraint.h
    src/constraint/CosseratSlidingConstraint.inl    
    src/mapping/LegendrePolynomialsMapping.h
    src/mapping/LegendrePolynomialsMapping.inl
    src/engine/PhysicalState.h
    src/engine/PhysicalState.inl

    )

set(SOURCE_FILES
    src/initCosserat.cpp
    src/mapping/BaseCosserat.cpp
    src/mapping/DiscreteCosseratMapping.cpp
    src/mapping/DiscreteDynamicCosseratMapping.cpp
    src/mapping/ProjectionEngine.cpp
    src/mapping/DifferenceMultiMapping.cpp
    src/mapping/RigidDistanceMapping.cpp
    src/forcefield/BeamHookeLawForceField.cpp
    src/forcefield/BeamPlasticLawForceField.cpp
    src/forcefield/CosseratInternalActuation.cpp
    src/constraint/CosseratSlidingConstraint.cpp
    src/mapping/LegendrePolynomialsMapping.cpp
    src/engine/PhysicalState.cpp
    )


set(README_FILE "Cosserat.txt" )

find_package(SoftRobots QUIET)
if(SoftRobots_FOUND)
    list(APPEND HEADER_FILES
        src/constraint/CosseratActuatorConstraint.h
        src/constraint/CosseratActuatorConstraint.inl
#        src/constraint/CosseratUnilateralInteractionConstraint.h
#        src/constraint/CosseratUnilateralInteractionConstraint.inl

        src/constraint/CosseratNeedleSlidingConstraint.h
        src/constraint/CosseratNeedleSlidingConstraint.inl
        src/constraint/QPSlidingConstraint.h
        src/constraint/QPSlidingConstraint.inl
        )
    list(APPEND SOURCE_FILES
        src/constraint/CosseratActuatorConstraint.cpp
        src/constraint/CosseratNeedleSlidingConstraint.cpp
#         src/constraint/CosseratUnilateralInteractionConstraint.cpp
        src/constraint/QPSlidingConstraint.cpp    
        )
endif()

file(GLOB_RECURSE RESOURCE_FILES  "*.md" "*.psl" "*.py" "*.pyscn" "*.py3scn" "*.scn" "*.ah")

IF(WIN32)
	add_definitions(-D_WINSOCKAPI_)
ENDIF(WIN32)

add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES} ${DOC_FILES} ${RESOURCE_FILES} )

target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>")
target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>")
target_include_directories(${PROJECT_NAME} PUBLIC "$<INSTALL_INTERFACE:include>")

set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "-DSOFA_BUILD_COSSERAT")
set_target_properties(${PROJECT_NAME} PROPERTIES DEBUG_POSTFIX "_d")
target_link_libraries(${PROJECT_NAME} Sofa.Helper SofaCore SofaConstraint SofaBaseMechanics SofaUserInteraction SofaRigid SofaMiscMapping SofaSparseSolver SofaEngine)

if(Sofa.GL_FOUND)
    target_link_libraries(${PROJECT_NAME} Sofa.GL)
endif()

if(SoftRobots_FOUND)
    target_link_libraries(${PROJECT_NAME} SoftRobots)
    message("-- Found dependency : 'SoftRobots' plugin .")
endif()

sofa_install_pythonscripts(PLUGIN_NAME ${PROJECT_NAME} PYTHONSCRIPTS_SOURCE_DIR "python")
find_file(SofaPython3Tools NAMES "SofaPython3/lib/cmake/SofaPython3/SofaPython3Tools.cmake")

if(SofaPython3Tools)
    message("-- Found SofaPython3Tools. Python3 packages will be installed.")
    include(${SofaPython3Tools})
    set(SP3_PYTHON_PACKAGES_DIRECTORY "python3/site-packages")
else()
    # try again with the find_package mechanism
    find_package(SofaPython3 QUIET)
endif()

if(SofaPython3Tools OR SofaPython3_FOUND)
    SP3_add_python_package(
            SOURCE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python3/cosserat
            TARGET_DIRECTORY cosserat
    )
    message("-- Found SofaPython3_FOUND. Python3 packages will be installed.")
endif()

## Install rules for the library and headers; CMake package configurations files
sofa_create_package_with_targets(
    PACKAGE_NAME ${PROJECT_NAME}
    PACKAGE_VERSION ${PROJECT_VERSION}
    TARGETS ${PROJECT_NAME} AUTO_SET_TARGET_PROPERTIES
    INCLUDE_SOURCE_DIR "src"
    INCLUDE_INSTALL_DIR ${PROJECT_NAME}
    RELOCATABLE "plugins"
    )

# find_package(GTest CONFIG QUIET)
# if (NOT GTest_FOUND)
#     # This find_package needs to be executed at MOST once, else it will fails on the second call
#     find_package(GMock QUIET)
# endif ()

option(COSSERAT_BUILD_TESTS "Build unit tests" OFF )
if(COSSERAT_BUILD_TESTS)
    add_subdirectory(Tests)
endif()

include(cmake/packaging.cmake)
