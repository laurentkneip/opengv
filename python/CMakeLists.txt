
find_package(Boost COMPONENTS python REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(Python_ADDITIONAL_VERSIONS "2.7")
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})         

find_package(NumPy REQUIRED)
include_directories(${NUMPY_INCLUDE_DIRS})
               
        
add_library( pyopengv SHARED pyopengv.cpp )
target_link_libraries(pyopengv opengv)


set_target_properties(pyopengv PROPERTIES
    PREFIX ""
    SUFFIX ".so"
)
if(APPLE)
    set_target_properties(pyopengv PROPERTIES
        LINK_FLAGS "-undefined dynamic_lookup"
    )
    target_link_libraries(pyopengv "${Boost_LIBRARY_DIR}/libboost_python.a")
else()
    target_link_libraries(pyopengv
        ${Boost_LIBRARIES}
        ${PYTHON_LIBRARIES}
    )
endif()

# Find whether to install python libs in site-packages or dist-packages
execute_process ( COMMAND
    python -c "import distutils.sysconfig; print 'dist-packages' if distutils.sysconfig.get_python_lib().endswith('dist-packages') else 'site-packages'"
    OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)

set(PYTHON_INSTALL_DIR
    "${CMAKE_INSTALL_PREFIX}/lib/python2.7/${PYTHON_SITE_PACKAGES}"
    CACHE PATH "Path where to install pyopengv")

install(TARGETS pyopengv DESTINATION "${PYTHON_INSTALL_DIR}")
