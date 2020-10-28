#######################
# CPack configuration #
#######################
set(CPACK_GENERATOR ZIP)
set(CPACK_PACKAGE_VERSION "${PROJECT_VERSION}")
set(CPACK_PACKAGE_NAME "${PROJECT_NAME} v${CPACK_PACKAGE_VERSION}")
set(CPACK_PACKAGE_VENDOR "Defrost team")
set(CPACK_PACKAGE_CONTACT "https://project.inria.fr/softrobot/contact/")
set(CPACK_PACKAGE_DESCRIPTION "Ease the modeling, the simulation and the control of 1D objecyts in SOFA.")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.txt")
set(CPACK_PACKAGE_HOMEPAGE_URL "https://project.inria.fr/cosserat")
set(CPACK_PACKAGE_FILE_NAME "${PROJECT_NAME}_v${CPACK_PACKAGE_VERSION}")
if(WIN32)
#    set(CPACK_PACKAGE_INSTALL_DIRECTORY "SOFA\\\\v${CPACK_PACKAGE_VERSION}\\\\plugins\\\\${PROJECT_NAME}")
    if("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "AMD64")
        set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_FILE_NAME}_Win64")
    else()
        set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_FILE_NAME}_Win32")
    endif()
elseif(UNIX)
    set(CPACK_PACKAGE_INSTALL_DIRECTORY "SOFA/v${CPACK_PACKAGE_VERSION}/plugins/${PROJECT_NAME}")
    if(APPLE)
        set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_FILE_NAME}_MacOS")
    else()
        set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_FILE_NAME}_Linux")
    endif()
endif()
# set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.txt")
set(CPACK_RESOURCE_FILE_WELCOME "${CMAKE_CURRENT_SOURCE_DIR}/README.txt")

include(CPack)
