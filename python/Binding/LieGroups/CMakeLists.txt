# This CMakeLists.txt defines how the Python bindings for the Lie groups are built.

project(LieGroupsBindings)

set(SOURCE_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/Module_LieGroupe.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Binding_LieGroups.cpp
)

set(HEADER_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/Binding_LieGroups.h
)

if (NOT TARGET SofaPython3::Plugin)
    find_package(SofaPython3 REQUIRED COMPONENTS Plugin Bindings.Sofa)
endif()

sofa_find_package(Sofa.GL QUIET)

SP3_add_python_module(
    TARGET       ${PROJECT_NAME}
    PACKAGE      Binding
    MODULE       LieGroups
    DESTINATION  /
    SOURCES      ${SOURCE_FILES}
    HEADERS      ${HEADER_FILES}
    DEPENDS      SofaPython3::Plugin SofaPython3::Bindings.Sofa LieGroups
)

message("-- SofaPython3 binding for the Cosserat plugin will be created.")
