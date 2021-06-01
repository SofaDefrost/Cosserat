cmake_minimum_required(VERSION 3.1)

project(CosseratPlugin_test VERSION 1.0)


find_package(CosseratPlugin REQUIRED)
find_package(SofaGTestMain REQUIRED)
find_package(SofaTest REQUIRED)

set(HEADER_FILES
        Tests/constraint/ClassName.h
        )
set(SOURCE_FILES

        Tests/constraint/first-test.cpp
    )

add_executable(${PROJECT_NAME} ${SOURCE_FILES})

target_include_directories(${PROJECT_NAME} PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>")

target_link_libraries(${PROJECT_NAME} SofaTest SofaGTestMain CosseratPlugin SofaCore SofaConstraint SofaBaseMechanics SofaUserInteraction)

add_definitions("-DBEAMADAPTER_TEST_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}/\"")
 
add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})

