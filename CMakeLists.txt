cmake_minimum_required(VERSION 2.8)

project(Felt)

set(PROJECT_ROOT ${PROJECT_SOURCE_DIR})
set(SRC_ROOT ${PROJECT_SOURCE_DIR}/src)
set(TESTS_ROOT ${SRC_ROOT}/tests)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_ROOT}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_ROOT}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_ROOT}/bin)

include_directories(
    ${PROJECT_ROOT}/include /usr/include/eigen3
)

add_definitions(-std=c++14 -fopenmp)


# C++ unit tests

add_executable(
	tests 
	${SRC_ROOT}/Poly.cpp
	${TESTS_ROOT}/test_all.cpp ${TESTS_ROOT}/test_Grid.cpp 
	${TESTS_ROOT}/test_MappedGrid.cpp ${TESTS_ROOT}/test_PartitionedGrid.cpp 
	${TESTS_ROOT}/test_Poly.cpp ${TESTS_ROOT}/test_Surface.cpp
	${TESTS_ROOT}/test_PolyGrid.cpp ${TESTS_ROOT}/Utils.cpp 
)
target_link_libraries(tests boost_unit_test_framework gomp)

# Main Felt library
add_library(Felt SHARED ${SRC_ROOT}/Poly.cpp)
target_link_libraries(Felt gomp)

