
add_subdirectory(pybind11)

include_directories(${PYTHON_INCLUDE_DIRS})


pybind11_add_module(pyopengv pyopengv.cpp)
target_link_libraries(pyopengv PRIVATE opengv)


# Find whether to install python libs in site-packages or dist-packages
execute_process( COMMAND
    python -c "import distutils.sysconfig; print('dist-packages' if distutils.sysconfig.get_python_lib().endswith('dist-packages') else 'site-packages')"
    OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)

set(PYTHON_INSTALL_DIR
    "${CMAKE_INSTALL_PREFIX}/lib/python2.7/${PYTHON_SITE_PACKAGES}"
    CACHE PATH "Path where to install pyopengv")

install(TARGETS pyopengv DESTINATION "${PYTHON_INSTALL_DIR}")
