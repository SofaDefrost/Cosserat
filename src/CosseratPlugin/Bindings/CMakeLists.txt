project(Bindings.Cosserat)

message("------> ${CMAKE_CURRENT_SOURCE_DIR}")
set(HEADER_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/Binding_PointManager.h
)

set(SOURCE_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/Binding_PointManager.cpp
)


find_package(SofaPython3 REQUIRED COMPONENTS Plugin)

find_package(Sofa.Simulation.Core REQUIRED)
find_package(Sofa.Component.Collision.Response.Contact REQUIRED) # for Binding_ContactListener
find_package(Sofa.Component.Visual REQUIRED) # for Binding_BaseCamera

SP3_add_python_module(
    TARGET       ${PROJECT_NAME}
    MODULE       BindingsCosserat
    DESTINATION  .
    SOURCES      ${SOURCE_FILES}
    HEADERS      ${HEADER_FILES}
    DEPENDS      SofaPython3::Plugin SofaPython3::Bindings.Sofa.Core
)
