add_subdirectory(Cosserat)
add_subdirectory(LieGroups)
#
# project(CoBindings)
#
# set(HEADER_FILES
#     ${CMAKE_CURRENT_SOURCE_DIR}/Binding_QPInverseProblemSolver.h
#     ${CMAKE_CURRENT_SOURCE_DIR}/Binding_QPInverseProblemSolver_doc.h
# )
#
# set(SOURCE_FILES
#     ${CMAKE_CURRENT_SOURCE_DIR}/Binding_QPInverseProblemSolver.cpp
#     ${CMAKE_CURRENT_SOURCE_DIR}/Module_SoftRobotsInverse.cpp
# )
#
# if (NOT TARGET SofaPython3::Plugin)
#     find_package(SofaPython3 REQUIRED COMPONENTS Plugin Bindings.Sofa)
# endif()
#
# SP3_add_python_module(
#     TARGET       ${PROJECT_NAME}
#     PACKAGE      Bindings.Modules
#     MODULE       SoftRobotsInverse
#     DESTINATION  /
#     SOURCES      ${SOURCE_FILES}
#     HEADERS      ${HEADER_FILES}
#     DEPENDS      SofaPython3::Plugin SofaPython3::Bindings.Sofa Sofa.GL SoftRobots.Inverse
# )
# message("-- SofaPython3 bindings for SoftRobots.Inverse will be created.")
