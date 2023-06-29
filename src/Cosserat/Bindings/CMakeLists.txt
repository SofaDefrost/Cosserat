project(CosseratBindings)

set(SOURCE_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/Module_Cosserat.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Binding_PointsManager.cpp
)

set(HEADER_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/Binding_PointsManager.h
)

if (NOT TARGET SofaPython3::Plugin)
    find_package(SofaPython3 REQUIRED COMPONENTS Plugin Bindings.Sofa)
endif()

sofa_find_package(Sofa.GL QUIET)

SP3_add_python_module(
    TARGET       ${PROJECT_NAME}
    PACKAGE      Binding
    MODULE       Cosserat
    DESTINATION  Sofa
    SOURCES      ${SOURCE_FILES}
    HEADERS      ${HEADER_FILES}
    DEPENDS      SofaPython3::Plugin SofaPython3::Bindings.Sofa Sofa.GL Cosserat

)
message("-- SofaPython3 binding for the Cosserat plugin will be created.")
