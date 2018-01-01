

find_package(Boost REQUIRED COMPONENTS python3)

if(${Boost_VERSION} LESS 106300)
    message("Not using boost/python/numpy.")
    set(USE_BOOST_PYTHON_NUMPY "off")
else()
    message("Using boost/python/numpy.")
    set(USE_BOOST_PYTHON_NUMPY "on")
    add_definitions(-DUSE_BOOST_PYTHON_NUMPY)

    find_package(Boost 1.63 REQUIRED COMPONENTS python3 numpy3)
endif()

include_directories(${Boost_INCLUDE_DIRS})

set(Python_ADDITIONAL_VERSIONS "3.6")
find_package(PythonLibs 3 REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

find_package(NumPy REQUIRED)
include_directories(${NUMPY_INCLUDE_DIRS})

add_library( pyopengv SHARED pyopengv.cpp )
target_link_libraries(pyopengv opengv)

set_target_properties(pyopengv PROPERTIES
    PREFIX ""
    SUFFIX ".so"
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED ON
)

if(APPLE)
    set_target_properties(pyopengv PROPERTIES
        LINK_FLAGS "-undefined dynamic_lookup"
    )
    target_link_libraries(pyopengv ${Boost_LIBRARIES})
else()
    target_link_libraries(pyopengv
        ${Boost_LIBRARIES}
        ${PYTHON_LIBRARIES}
    )
endif()

# Find whether to install python libs in site-packages or dist-packages
execute_process ( COMMAND
    python -c "import distutils.sysconfig; print('dist-packages' if distutils.sysconfig.get_python_lib().endswith('dist-packages') else 'site-packages')"
    OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)

set(PYTHON_INSTALL_DIR
    "${CMAKE_INSTALL_PREFIX}/lib/python2.7/${PYTHON_SITE_PACKAGES}"
    CACHE PATH "Path where to install pyopengv")

install(TARGETS pyopengv DESTINATION "${PYTHON_INSTALL_DIR}")
