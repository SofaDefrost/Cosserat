# This CMakeLists.txt defines how the Python bindings for the Lie groups are built.

project(LieGroups.Python)

set(SOURCE_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/Module_LieGroupe.cpp
 # # ${CMAKE_CURRENT_SOURCE_DIR}/src/Binding_LieGroups.cpp
)

set(HEADER_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Binding_LieGroups.h
)

if (NOT TARGET SofaPython3::Plugin)
    find_package(SofaPython3 REQUIRED COMPONENTS Plugin Bindings.Sofa)
endif()

sofa_find_package(Sofa.GL QUIET)

SP3_add_python_module(
    TARGET       ${PROJECT_NAME}
    PACKAGE      LieGroups.Python
    MODULE       LieGroups
    DESTINATION  sofa
    SOURCES      ${SOURCE_FILES}
    HEADERS      ${HEADER_FILES}
    DEPENDS      Cosserat SofaPython3::Plugin SofaPython3::Bindings.Sofa
)

message("-- SofaPython3 binding for the Cosserat plugin will be created.")
