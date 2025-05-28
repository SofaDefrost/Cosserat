cmake_minimum_required(VERSION 3.12)

project(Cosserat_test)

if (UNIT_TEST)

    enable_testing()

    set(HEADER_FILES
    )

    set(SOURCE_FILES
        engine/GeometricStiffnessEngineTest.cpp
        forcefield/BeamHookeLawForceFieldTest.cpp
        mapping/AdaptiveBeamMappingTest.cpp
        mapping/CosseratNonLinearMapping2DTest.cpp
        liegroups/RealSpaceTest.cpp
        liegroups/SO2Test.cpp
        liegroups/SE2Test.cpp
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

# Add benchmarks subdirectory if benchmarks are enabled
if(COSSERAT_BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()
