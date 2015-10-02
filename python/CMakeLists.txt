
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

install(TARGETS pyopengv DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/python2.7/site-packages")
