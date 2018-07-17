
find_package(Boost REQUIRED)

option(BUILD_FOR_PYTHON3 "build pyopengv for Python 3" OFF)
set(BOOST_PYTHON3_COMPONENT "python3" CACHE STRING
    "name of the Boost Python 3 component (python3, python-py36 or python-py37)")

if(${BUILD_FOR_PYTHON3})
    set(BOOST_PYTHON_COMPONENT ${BOOST_PYTHON3_COMPONENT})
    set(BOOST_NUMPY_COMPONENT numpy3)
    set(Python_ADDITIONAL_VERSIONS 3.7 3.6)
    set(Python_VERSION_SUFFIX 3)
else()
    if(${Boost_VERSION} LESS 106700)
        set(BOOST_PYTHON_COMPONENT python)
        set(BOOST_NUMPY_COMPONENT numpy)
    else()
        set(BOOST_PYTHON_COMPONENT python27)
        set(BOOST_NUMPY_COMPONENT numpy27)
    endif()
    set(Python_ADDITIONAL_VERSIONS 2.7)
    set(Python_VERSION_SUFFIX "")
endif()

if(${Boost_VERSION} LESS 106300)
    message("Not using boost/python/numpy.")
    set(USE_BOOST_PYTHON_NUMPY "off")

    find_package(Boost REQUIRED COMPONENTS ${BOOST_PYTHON_COMPONENT})
else()
    message("Using boost/python/numpy.")
    set(USE_BOOST_PYTHON_NUMPY "on")
    add_definitions(-DUSE_BOOST_PYTHON_NUMPY)

    find_package(Boost 1.63 REQUIRED COMPONENTS
                 ${BOOST_PYTHON_COMPONENT} ${BOOST_NUMPY_COMPONENT})
endif()

include_directories(${Boost_INCLUDE_DIRS})

find_package(PythonLibs ${Python_VERSION_SUFFIX} REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

find_package(NumPy REQUIRED)
include_directories(${NUMPY_INCLUDE_DIRS})

add_library(pyopengv SHARED pyopengv.cpp)
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