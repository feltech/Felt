cmake_minimum_required(VERSION 3.16)

project(Felt)

find_package(Eigen3 REQUIRED)
find_package(OpenMP REQUIRED)

add_library(felt INTERFACE)
target_include_directories(
	felt
	INTERFACE
	${PROJECT_SOURCE_DIR}/include
)
target_compile_features(felt INTERFACE cxx_std_14)
target_link_libraries(felt INTERFACE Eigen3::Eigen OpenMP::OpenMP_CXX)

# C++ unit tests
if (FELT_ENABLE_TESTS)
	enable_testing()
	set(TESTS_ROOT ${PROJECT_SOURCE_DIR}/src/tests)
	file(GLOB_RECURSE CTEST_CPPS "${TESTS_ROOT}/test*.cpp")
	add_executable(ctests "${TESTS_ROOT}/main.cpp" ${CTEST_CPPS})
	target_link_libraries(ctests PRIVATE felt)
	add_test(
		NAME felt
		WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
		COMMAND ctests
	)
endif ()
