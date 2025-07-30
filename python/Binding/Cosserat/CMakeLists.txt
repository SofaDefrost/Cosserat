project(Cosserat.Python VERSION 21.12.0)

set(SOURCE_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Module_Cosserat.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Binding_PointsManager.cpp
)

set(HEADER_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Binding_PointsManager.h
)

if (NOT TARGET SofaPython3::Plugin)
    find_package(SofaPython3 REQUIRED COMPONENTS Plugin Bindings.Sofa)
endif()
message("t===> " ${SOURCE_FILES}, " <===", )
message("t===> " ${HEADER_FILES}, " <===", )
SP3_add_python_module(
    TARGET       ${PROJECT_NAME}
    PACKAGE      Cosserat.Python
    MODULE       Cosserat
    DESTINATION  Sofa
    SOURCES      ${SOURCE_FILES}
    HEADERS      ${HEADER_FILES}
    DEPENDS      Cosserat SofaPython3::Plugin SofaPython3::Bindings.Sofa.Core
)

#target_compile_features(${PROJECT_NAME} PRIVATE cxx_std_20)
#target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/src)
message("-- SofaPython3 binding for the Cosserat plugin will be created.")
