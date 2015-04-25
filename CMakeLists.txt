cmake_minimum_required(VERSION 2.8)

project(felt)

set(PROJECT_ROOT ${PROJECT_SOURCE_DIR})
set(TESTS_ROOT ${PROJECT_SOURCE_DIR}/src/tests)
set(GUI_ROOT ${PROJECT_SOURCE_DIR}/src/gui)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_ROOT}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_ROOT}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_ROOT}/bin)

set(
    PYTHON_2_7_INCLUDE "/usr/include/python2.7" 
    CACHE FILEPATH "Python 2.7 include directory"
)
set(
    CYTHON_INCLUDE "/usr/local/lib/python2.7/dist-packages/Cython/Includes" 
    CACHE FILEPATH "Cython include directory"
)

include_directories(
    ${PROJECT_ROOT}/include ${PROJECT_ROOT}/vendor/include ${PYTHON_2_7_INCLUDE}
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


# Python bindings

# RktSFMLInterface shared library
add_library(
    RktSFMLInterface SHARED
    ${GUI_ROOT}/RktFileInterface.cpp ${GUI_ROOT}/RktSFMLInterface.cpp 
    ${GUI_ROOT}/RktSFMLRenderer.cpp ${GUI_ROOT}/RktSFMLSystem.cpp
)

# PyRktSFMLInterface Python extension.
add_custom_target(
    PyRktSFMLInterface ALL
    PROJECT_ROOT=${PROJECT_ROOT} 
    PYTHONPATH=${CYTHON_INCLUDE}
    python setup.py build_ext -I${PYTHON_2_7_INCLUDE}
    WORKING_DIRECTORY ${PROJECT_ROOT}/scripts/felt
    COMMENT "\n\n<<<<< Building Python/C++ libRocket/SFML interface >>>>>"
)
add_dependencies(PyRktSFMLInterface RktSFMLInterface)

# pyfelt.c embedding pyfelt.py
add_custom_command(
    OUTPUT ${PROJECT_ROOT}/scripts/pyfelt.c
    COMMAND cython --embed -o ${PROJECT_ROOT}/scripts/pyfelt.c 
        ${PROJECT_ROOT}/scripts/pyfelt.py
    COMMENT "\n\n<<<<< Building embedded Python executable >>>>>"
    DEPENDS ${PROJECT_ROOT}/scripts/pyfelt.py
)
add_custom_target(cyfelt DEPENDS ${PROJECT_ROOT}/scripts/pyfelt.c)

# pyfelt executable.
add_executable(pyfelt ${PROJECT_ROOT}/scripts/pyfelt.c)
add_dependencies(pyfelt cyfelt)


target_link_libraries(pyfelt python2.7)