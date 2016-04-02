cmake_minimum_required(VERSION 2.8)

project(Felt)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -ggdb3")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -ggdb3")

set(PROJECT_ROOT ${PROJECT_SOURCE_DIR})
set(SRC_ROOT ${PROJECT_SOURCE_DIR}/src)
set(TESTS_ROOT ${SRC_ROOT}/tests)
set(EIGEN3_INCLUDE_DIR /usr/include/eigen3 CACHE PATH "Eigen include directory")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_ROOT}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_ROOT}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_ROOT}/bin)

include_directories(
    ${PROJECT_ROOT}/include ${EIGEN3_INCLUDE_DIR}
)

add_definitions(-std=c++14 -fopenmp -Wfatal-errors)


# C++ unit tests
file(GLOB TEST_CPPS "${TESTS_ROOT}/*.cpp")
add_executable(tests ${TEST_CPPS})
target_link_libraries(tests boost_unit_test_framework gomp)
