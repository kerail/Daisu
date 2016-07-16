set(SO3_TEST_SOURCE_FILES test/test_so3.cpp)
set(SO3_TEST_PERF_SOURCE_FILES test/test_so3_performance.cpp)
set(SO3_TEST_PREC_SOURCE_FILES test/test_so3_precision.cpp)

find_package(Boost 1.55 COMPONENTS unit_test_framework REQUIRED)
find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif ()

add_executable(test_so3 ${SO3_TEST_SOURCE_FILES})
target_link_libraries(test_so3
        ${Boost_FILESYSTEM_LIBRARY}
        ${Boost_SYSTEM_LIBRARY}
        ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})


add_executable(test_so3_performance ${SO3_TEST_PERF_SOURCE_FILES})
add_executable(test_so3_precision ${SO3_TEST_PREC_SOURCE_FILES})