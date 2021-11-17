cmake_minimum_required(VERSION 3.1)
project(CosseratPlugin VERSION 21.06.99)

# Policies
cmake_policy(SET CMP0072 NEW)
cmake_policy(SET CMP0074 NEW)
cmake_policy(SET CMP0076 NEW)
cmake_policy(SET CMP0094 NEW)

set (COSSERAT_VERSION ${PROJECT_VERSION})

## Default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()


include(cmake/environment.cmake)
find_package(SofaFramework REQUIRED)
#find_package(SOFA COMPONENTS SofaFramework REQUIRED)

option(COSSERAT_IGNORE_ERRORS "Enable this option to ignore the recommendations and generate anyway." ON)

if(NOT COSSERAT_IGNORE_ERRORS)
    if(NOT PLUGIN_SOFAPYTHON AND NOT PLUGIN_SOFAPYTHON3)
        message(SEND_ERROR "
            You should not use the plugin Cosserat without enabling PLUGIN_SOFAPYTHON or PLUGIN_SOFAPYTHON3.
            Since all the scenes to test this plugin are using one of this two plugins. 
            To fix this error you can either follow the recommendation or enable COSSERAT_IGNORE_ERRORS.")
    endif()    
endif()

# set(SofaPython_VERSION 3.8.5)
# find_package(SofaPython3 QUIET)
set(COSSERAT_PYTHON ${SofaPython_FOUND}) # config.h.in
find_package(STLIB QUIET)
if(STLIB_FOUND)
    message("-- Found dependency : 'STLIB' plugin .")
else()
    message("-- The highly recommanded 'STLIB' plugin is missing. You can compile Cosserat but some of the provided python examples will not work. ")
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
    src/forcefield/CosseratInternalActuation.h
    src/forcefield/CosseratInternalActuation.inl
    src/forcefield/MyUniformVelocityDampingForceField.h
    src/forcefield/MyUniformVelocityDampingForceField.inl
    
    src/constraint/CosseratSlidingConstraint.h
    src/constraint/CosseratSlidingConstraint.inl
    src/constraint/QPSlidingConstraint.h
    src/constraint/QPSlidingConstraint.inl
    src/constraint/CosseratNeedleSlidingConstraint.h
    src/constraint/CosseratNeedleSlidingConstraint.inl
    src/constraint/CosseratUnilateralInteractionConstraint.h
    src/constraint/CosseratUnilateralInteractionConstraint.inl
    # src/mapping/LegendrePolynomialsMapping.h
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
    src/forcefield/CosseratInternalActuation.cpp
    src/forcefield/MyUniformVelocityDampingForceField.cpp

    src/constraint/CosseratSlidingConstraint.cpp
    src/constraint/QPSlidingConstraint.cpp
    src/constraint/CosseratNeedleSlidingConstraint.cpp
    src/constraint/CosseratUnilateralInteractionConstraint.cpp
    # src/mapping/LegendrePolynomialsMapping.cpp
    )

set(README_FILE "Cosserat.txt" )

find_package(SoftRobots QUIET)
if(SoftRobots_FOUND)
    list(APPEND HEADER_FILES
        src/constraint/CosseratActuatorConstraint.h
        src/constraint/CosseratActuatorConstraint.inl
        )
    list(APPEND SOURCE_FILES
        src/constraint/CosseratActuatorConstraint.cpp
        )
endif()

file(GLOB_RECURSE RESOURCE_FILES  "*.md" "*.psl" "*.py" "*.pyscn" "*.scn" "*.ah")

add_library(${PROJECT_NAME} SHARED ${HEADER_FILES} ${SOURCE_FILES} ${RESOURCE_FILES} )

target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>")
target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>")
target_include_directories(${PROJECT_NAME} PUBLIC "$<INSTALL_INTERFACE:include>")

set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "-DSOFA_BUILD_COSSERAT")
set_target_properties(${PROJECT_NAME} PROPERTIES DEBUG_POSTFIX "_d")
target_link_libraries(${PROJECT_NAME} SofaCore SofaConstraint SofaSimpleFem SofaBaseMechanics SofaRigid SofaBaseVisual SofaOpenglVisual)


# if(SofaPython_FOUND)
#     target_link_libraries(${PROJECT_NAME} SofaPython)
#     Config files and install rules for pythons scripts
#     sofa_install_pythonscripts(PLUGIN_NAME ${PROJECT_NAME} PYTHONSCRIPTS_SOURCE_DIR "python/cosserat")
# endif()

sofa_install_pythonscripts(PLUGIN_NAME ${PROJECT_NAME} PYTHONSCRIPTS_SOURCE_DIR "python")
if(SofaPython3_FOUND)
    get_filename_component(SP3_ROOT_DIR "${SofaPython3_DIR}/../../.." ABSOLUTE)
    message(STATUS
            "Cosserat find SofaPython3 support  \n"
            "SP3 version: ${SofaPython3_VERSION}\n"
            "SP3 location: ${SP3_ROOT_DIR}"
    )
    SP3_add_python_package(
            SOURCE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python3/cosserat
            TARGET_DIRECTORY cosserat
    )
endif()

if(SoftRobots_FOUND)
    target_link_libraries(${PROJECT_NAME} SoftRobots)
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


find_package(GTest CONFIG QUIET)
if (NOT GTest_FOUND)
    # This find_package needs to be executed at MOST once, else it will fails on the second call
    find_package(GMock QUIET)
endif ()


option(COSSERAT_BUILD_TESTS "Build unit tests" ON )
if(COSSERAT_BUILD_TESTS)
    add_subdirectory(Tests)
endif()

SET_PROPERTY(GLOBAL PROPERTY USE_FOLDERS ON)
SET_PROPERTY(TARGET ${PROJECT_NAME} PROPERTY FOLDER "plugins")

include(cmake/packaging.cmake)
