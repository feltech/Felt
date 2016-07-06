cmake_minimum_required(VERSION 2.8)

project(Felt)

cmake_policy(SET CMP0015 NEW)

set(PROJECT_ROOT ${PROJECT_SOURCE_DIR})
set(SRC_ROOT ${PROJECT_SOURCE_DIR}/src)
set(TESTS_ROOT ${SRC_ROOT}/tests)
set(EIGEN3_INCLUDE_DIR /usr/include/eigen3 CACHE PATH "Eigen include directory")
#set(GTEST_INCLUDE_DIR ../googletest/googletest/include CACHE PATH "Google Test include directory")
#set(GTEST_LIBRARY_DIR ../googletest/build/Debug CACHE PATH "Google Test library directory")

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_ROOT}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_ROOT}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_ROOT}/bin)

include_directories(
    ${PROJECT_ROOT}/include ${EIGEN3_INCLUDE_DIR}
)

add_definitions(
    -std=c++14 -fopenmp -Wfatal-errors -fdiagnostics-color=auto
)


# C++ unit tests
#file(GLOB TEST_CPPS "${TESTS_ROOT}/test_*.cpp")

#add_executable(tests ${TEST_CPPS})
#target_link_libraries(tests boost_unit_test_framework gomp)

#link_directories(${GTEST_LIBRARY_DIR})
#file(GLOB GTEST_CPPS "${TESTS_ROOT}/*Test.cpp")
#add_executable(gtests ${GTEST_CPPS} "${TESTS_ROOT}/Utils.cpp")
#include_directories(${GTEST_INCLUDE_DIR} ${PROJECT_ROOT}/include ${EIGEN3_INCLUDE_DIR})
#target_link_libraries(gtests gtest gomp pthread)


file(GLOB CTEST_CPPS "${TESTS_ROOT}/catchtest*.cpp")
add_executable(ctests ${CTEST_CPPS} ${TESTS_ROOT}/Utils.cpp)
include_directories(${PROJECT_ROOT}/include ${EIGEN3_INCLUDE_DIR})
target_link_libraries(ctests gomp pthread)
