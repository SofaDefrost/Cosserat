cmake_minimum_required(VERSION 3.12)

project(Cosserat_test)

if (UNIT_TEST)

    enable_testing()

    set(HEADER_FILES
    )

    set(SOURCE_FILES
        unit/Example.cpp
        unit/SO2Test.cpp
        unit/SO3Test.cpp
        unit/SE3Test.cpp
        unit/SE23Test.cpp
        unit/SGal3Test.cpp
        unit/BundleTest.cpp
        unit/BeamHookeLawForceFieldTest.cpp
        unit/DiscretCosseratMappingTest.cpp
        unit/CosseratUnilateralInteractionConstraintTest.cpp
    )

    add_executable(${PROJECT_NAME} ${HEADER_FILES} ${SOURCE_FILES})
    target_link_libraries(${PROJECT_NAME} Cosserat_SRC)
    target_link_libraries(${PROJECT_NAME} gtest gtest_main)
    target_link_libraries(${PROJECT_NAME} SofaGTestMain)

    if(APPLE)
    target_link_libraries(${PROJECT_NAME} "-framework Cocoa")
    endif()

    add_test(NAME ${PROJECT_NAME} COMMAND ${PROJECT_NAME})

endif()

# Add integration tests if enabled
if (INTEGRATION_TEST)
    # Set up integration tests
    set(INTEGRATION_SOURCE_FILES
        integration/LieGroupIntegrationTest.cpp
    )

    add_executable(${PROJECT_NAME}_integration ${INTEGRATION_SOURCE_FILES})
    target_link_libraries(${PROJECT_NAME}_integration Cosserat_SRC)
    target_link_libraries(${PROJECT_NAME}_integration gtest gtest_main)
    target_link_libraries(${PROJECT_NAME}_integration SofaGTestMain)

    if(APPLE)
        target_link_libraries(${PROJECT_NAME}_integration "-framework Cocoa")
    endif()

    add_test(NAME ${PROJECT_NAME}_integration COMMAND ${PROJECT_NAME}_integration)
endif()

# Add benchmarks subdirectory if benchmarks are enabled
if(COSSERAT_BUILD_BENCHMARKS)
    # Create benchmarks CMakeLists.txt if it doesn't exist
    if(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/CMakeLists.txt)
        file(WRITE ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/CMakeLists.txt
"cmake_minimum_required(VERSION 3.12)

project(Cosserat_benchmarks)

set(SOURCE_FILES
    LieGroupBenchmark.cpp
)

add_executable(${PROJECT_NAME} ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} Cosserat_SRC)
target_link_libraries(${PROJECT_NAME} benchmark::benchmark_main)

if(APPLE)
    target_link_libraries(${PROJECT_NAME} \"-framework Cocoa\")
endif()
")
    endif()

    add_subdirectory(benchmarks)
endif()
