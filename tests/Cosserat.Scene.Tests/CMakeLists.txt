cmake_minimum_required(VERSION 3.12)
project(Cosserat.Scene.Tests C CXX)

find_package(SofaPython3 REQUIRED)
find_package(Sofa.Testing REQUIRED)

enable_testing()

set(HEADER_FILES
        Cosserat_Python_Tests.h
        )

set(SOURCE_FILES
        Cosserat_Python_Tests.cpp
    )

set(PYTHON_FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/Components/BaseCosseratMapping.py
    ${CMAKE_CURRENT_SOURCE_DIR}/Components/DiscreteCosseratMapping.py
    )

add_executable(${PROJECT_NAME} ${SOURCE_FILES} ${HEADER_FILES} ${PYTHON_FILES})

target_link_libraries(${PROJECT_NAME} 
        SofaPython3Testing
        Sofa.Testing
        Cosserat
)

target_include_directories(${PROJECT_NAME}
        PUBLIC 
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../>"
)

add_test(
        NAME ${PROJECT_NAME}
        COMMAND $${PROJECT_NAME}
)


find_package(Sofa.Testing REQUIRED)

# This will set rpaths relative to all SOFA core AND relocatable dependencies
sofa_auto_set_target_rpath(
    TARGETS ${PROJECT_NAME}
    RELOCATABLE "plugins"
)


set(DIR_BINDING_LIST
    Components
    )
foreach(dir_binding ${DIR_BINDING_LIST})
    if (_isMultiConfig) # MSVC
        foreach(config_type ${CMAKE_CONFIGURATION_TYPES})
            SP3_configure_directory(${CMAKE_CURRENT_SOURCE_DIR}/${dir_binding} ${CMAKE_BINARY_DIR}/${RUNTIME_OUTPUT_DIRECTORY}/${config_type}/${PROJECT_NAME}.d/${dir_binding})
            install(DIRECTORY ${CMAKE_BINARY_DIR}/${RUNTIME_OUTPUT_DIRECTORY}/${config_type}/${PROJECT_NAME}.d/${dir_binding} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY}/${PROJECT_NAME}.d)
        endforeach()
    else()
        SP3_configure_directory(${CMAKE_CURRENT_SOURCE_DIR}/${dir_binding} ${CMAKE_BINARY_DIR}/${RUNTIME_OUTPUT_DIRECTORY}/${PROJECT_NAME}.d/${dir_binding})
        install(DIRECTORY ${CMAKE_BINARY_DIR}/${RUNTIME_OUTPUT_DIRECTORY}/${PROJECT_NAME}.d/${dir_binding} DESTINATION ${RUNTIME_OUTPUT_DIRECTORY}/${PROJECT_NAME}.d)
    endif()
endforeach()
