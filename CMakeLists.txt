project(Felt)

cmake_minimum_required(VERSION 2.8)
include(CMakeToolsHelpers OPTIONAL)
find_program(CTEST_MEMORYCHECK_COMMAND NAMES valgrind)

cmake_policy(SET CMP0015 NEW)

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

add_definitions(
	-std=c++14 -fopenmp -Wfatal-errors -fdiagnostics-color=auto -DCATCH_CONFIG_CONSOLE_WIDTH=200
	-DFELT_EXCEPTIONS -Wall -Wextra -Wconversion -Wsign-conversion
)

# C++ unit tests
enable_testing()
file(GLOB_RECURSE CTEST_CPPS "${TESTS_ROOT}/test*.cpp")
add_executable(ctests "${TESTS_ROOT}/main.cpp" ${CTEST_CPPS})
include_directories(${PROJECT_ROOT}/include ${EIGEN3_INCLUDE_DIR})
target_link_libraries(ctests gomp pthread)
add_test(
	NAME felt
	WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
	COMMAND ctests
)
