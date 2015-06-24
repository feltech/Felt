cmake_minimum_required(VERSION 2.8)

project(felt)

set(PROJECT_ROOT ${PROJECT_SOURCE_DIR})
set(TESTS_ROOT ${PROJECT_SOURCE_DIR}/src/tests)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_ROOT}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_ROOT}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_ROOT}/bin)

include_directories(
    ${PROJECT_ROOT}/include
)

add_definitions(-std=c++11)


# C++ unit tests

add_executable(
	tests 
	${TESTS_ROOT}/test_all.cpp ${TESTS_ROOT}/test_Grid.cpp 
	${TESTS_ROOT}/test_MappedGrid.cpp ${TESTS_ROOT}/test_PartitionedGrid.cpp 
	${TESTS_ROOT}/test_Poly.cpp ${TESTS_ROOT}/test_Surface.cpp
)

target_link_libraries(tests boost_unit_test_framework gomp)
