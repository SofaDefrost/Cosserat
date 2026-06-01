cmake_minimum_required(VERSION 3.12)

project(Cosserat_benchmarks)

set(SOURCE_FILES
    LieGroupBenchmark.cpp
)

add_executable(
  Cosserat_tests unit/Strain2RigidCosseratMappingTest.cpp;
  unit/LiegroupsIntegrationTest.cpp;
  unit/DiscreteCosseratMappingApplyDJTTest.cpp
)
target_link_libraries(Cosserat_tests Cosserat_SRC)
target_link_libraries(Cosserat_tests benchmark::benchmark_main)

if(APPLE)
    target_link_libraries(Cosserat_tests "-framework Cocoa")
endif()
