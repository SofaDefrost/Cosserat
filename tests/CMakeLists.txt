cmake_minimum_required(VERSION 3.1)

project(DiscretCosseratMapping_test VERSION 1.0)


find_package(CosseratPlugin REQUIRED)
find_package(SofaTest REQUIRED)
find_package(SofaGTestMain REQUIRED)

set(SOURCE_FILES

    mapping/DiscretCosseratMappingTest.cpp
    
    )

add_executable(${PROJECT_NAME} ${SOURCE_FILES})

target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>")

target_link_libraries(${PROJECT_NAME} CosseratPlugin SofaTest SofaGTestMain SofaCore SceneCreator SofaConstraint SofaBaseMechanics SofaUserInteraction SofaSimpleFem SofaRigid)

add_definitions("-DBEAMADAPTER_TEST_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/\"")
 
add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})
