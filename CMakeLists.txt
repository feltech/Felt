cmake_minimum_required(VERSION 2.8)

project(Felt)

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

add_definitions(-std=c++14 -fopenmp)


# C++ unit tests

add_executable(
	tests 
	${TESTS_ROOT}/test_all.cpp ${TESTS_ROOT}/test_Grid.cpp 
	${TESTS_ROOT}/test_MappedGrid.cpp ${TESTS_ROOT}/test_PartitionedGrid.cpp 
	${TESTS_ROOT}/test_Poly.cpp ${TESTS_ROOT}/test_Surface.cpp
	${TESTS_ROOT}/test_PolyGrid.cpp ${TESTS_ROOT}/Utils.cpp 
)
target_link_libraries(tests boost_unit_test_framework gomp)

# Main Felt library
target_link_libraries(Felt gomp)

