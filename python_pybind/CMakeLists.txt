# cmake_minimum_required(VERSION 3.1.0)
# project(opengv_python_pybind_standalone CXX)

# set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/modules/")
## find_package(opengv REQUIRED)
# include( C:/tools/cmake_install_libs/opengv/CMake/opengv-config.cmake )


set(LIBNAME pyopengv)

add_subdirectory( pybind11 )

find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

find_package(NumPy REQUIRED)
include_directories(${NUMPY_INCLUDE_DIRS})

add_library(${LIBNAME} MODULE pyopengv.cpp) # note the module for building libs
target_link_libraries(${LIBNAME} opengv pybind11::module)

if( MSVC )
set_target_properties(${LIBNAME} PROPERTIES
    SUFFIX ".pyd"
)
else()
set_target_properties(${LIBNAME} PROPERTIES
    SUFFIX ".so"
)
endif()

set(py_ver "${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}")
MESSAGE( STATUS "python version: " ${py_ver} )

set(PYTHON_INSTALL_DIR
    "${CMAKE_INSTALL_PREFIX}/lib/python${py_ver}")

install(TARGETS ${LIBNAME} DESTINATION "${PYTHON_INSTALL_DIR}")