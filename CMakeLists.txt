cmake_minimum_required(VERSION 3.12)
project(opengv VERSION 1.0 LANGUAGES CXX)

# Set C++ standard (required by Eigen 3.x)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set the build type.  Options are:
#
#  None (CMAKE_C_FLAGS or CMAKE_CXX_FLAGS used)
#  Debug (CMAKE_C_FLAGS_DEBUG or CMAKE_CXX_FLAGS_DEBUG)
#  Release (CMAKE_C_FLAGS_RELEASE or CMAKE_CXX_FLAGS_RELEASE)
#  RelWithDebInfo (CMAKE_C_FLAGS_RELWITHDEBINFO or CMAKE_CXX_FLAGS_RELWITHDEBINFO)
#  MinSizeRel (CMAKE_C_FLAGS_MINSIZEREL or CMAKE_CXX_FLAGS_MINSIZEREL)

set(CMAKE_BUILD_TYPE Release)

#set the default path for built executables to the "bin" directory
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)
#set the default path for built libraries to the "lib" directory
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)

OPTION(BUILD_TESTS "Build tests" ON)
OPTION(BUILD_PYTHON "Build Python extension" OFF)

IF(MSVC)
  set(BUILD_SHARED_LIBS OFF)
ELSE()
  OPTION(BUILD_SHARED_LIBS "Build shared libraries" OFF)
  OPTION(BUILD_POSITION_INDEPENDENT_CODE "Build position independent code (-fPIC)" ON)
ENDIF()

# Define M_PI for all Windows compilers (MSVC, Clang, MinGW)
IF(WIN32)
  add_definitions(-D_USE_MATH_DEFINES)
ENDIF()

IF(MSVC)
  add_compile_options(/wd4514 /wd4267 /bigobj)
  # Disable Eigen half-precision float to avoid MSVC ostream operator conflict
  add_definitions(-DEIGEN_HAS_HALF_FLOAT=0)
  # Increase stack size for large Groebner matrix allocations in fivept_kneip and sixpt algorithms
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /STACK:8388608")
ELSE()
  # Check if we're using Clang on Windows (more conservative flags needed)
  IF(CMAKE_CXX_COMPILER_ID MATCHES "Clang" AND WIN32)
    # Clang on Windows: use safer optimization flags
    add_definitions (
      -O2
      -Wall
      -Wextra
      -Wwrite-strings
      -Wno-unused-parameter
      -fno-strict-aliasing
    )
  ELSE()
    # GCC or Clang on Unix: use aggressive optimizations
    IF (CMAKE_SYSTEM_PROCESSOR MATCHES "(arm64)|(ARM64)|(aarch64)|(AARCH64)")
      add_definitions (-march=armv8-a)
    ELSEIF (CMAKE_SYSTEM_PROCESSOR MATCHES
            "(arm)|(ARM)|(armhf)|(ARMHF)|(armel)|(ARMEL)")
      add_definitions (-march=armv7-a)
    ELSE ()
      add_definitions (-march=native) #TODO use correct c++11 def once everybody has moved to gcc 4.7 # for now I even removed std=gnu++0x
    ENDIF()
    add_definitions (
      -O3
      -Wall
      -Wextra
      #-Werror
      -Wwrite-strings
      -Wno-unused-parameter
      -fno-strict-aliasing
    )
  ENDIF()
ENDIF()

IF(BUILD_POSITION_INDEPENDENT_CODE)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
ENDIF()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/modules/")

# Find Eigen - try multiple methods for better cross-platform support
# First try Eigen3 (modern CMake find_package)
find_package(Eigen3 QUIET)
if(Eigen3_FOUND)
  # Modern CMake target-based Eigen (e.g., vcpkg, modern installs)
  if(TARGET Eigen3::Eigen)
    # Get include directory from target properties
    get_target_property(_EIGEN3_INCLUDE_DIR Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
    if(_EIGEN3_INCLUDE_DIR)
      set(ADDITIONAL_INCLUDE_DIRS ${_EIGEN3_INCLUDE_DIR} ${_EIGEN3_INCLUDE_DIR}/unsupported)
      message(STATUS "Found Eigen3 via target Eigen3::Eigen: ${_EIGEN3_INCLUDE_DIR}")
    endif()
  elseif(EIGEN3_INCLUDE_DIR)
    # Legacy variable-based Eigen
    set(ADDITIONAL_INCLUDE_DIRS ${EIGEN3_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR}/unsupported)
    message(STATUS "Found Eigen3 via find_package(Eigen3): ${EIGEN3_INCLUDE_DIR}")
  endif()
else()
  # Fall back to custom FindEigen module
  find_package(Eigen REQUIRED)
  if(EIGEN_INCLUDE_DIRS)
    set(ADDITIONAL_INCLUDE_DIRS ${EIGEN_INCLUDE_DIRS} ${EIGEN_INCLUDE_DIRS}/unsupported)
    message(STATUS "Found Eigen via FindEigen module: ${EIGEN_INCLUDE_DIRS}")
  elseif(EIGEN_INCLUDE_DIR)
    set(ADDITIONAL_INCLUDE_DIRS ${EIGEN_INCLUDE_DIR} ${EIGEN_INCLUDE_DIR}/unsupported)
    message(STATUS "Found Eigen via FindEigen module: ${EIGEN_INCLUDE_DIR}")
  else()
    message(FATAL_ERROR "Eigen not found. Please install Eigen3 or set EIGEN_INCLUDE_DIR. "
                        "On Windows with conda: conda install eigen")
  endif()
endif()

# Verify Eigen unsupported modules exist
set(_EIGEN_UNSUPPORTED_FOUND FALSE)
foreach(_DIR IN LISTS ADDITIONAL_INCLUDE_DIRS)
  if(EXISTS "${_DIR}/Eigen/NonLinearOptimization")
    set(_EIGEN_UNSUPPORTED_FOUND TRUE)
    message(STATUS "Eigen unsupported modules found at: ${_DIR}")
    break()
  endif()
endforeach()

if(NOT _EIGEN_UNSUPPORTED_FOUND)
  message(WARNING "Eigen unsupported modules (NonLinearOptimization, NumericalDiff) not found. "
                  "Some features may not compile. Please ensure Eigen is fully installed.")
endif()

set( OPENGV_SOURCE_FILES
  src/absolute_pose/modules/main.cpp
  src/absolute_pose/modules/gp3p/code.cpp
  src/absolute_pose/modules/gp3p/init.cpp
  src/absolute_pose/modules/gp3p/reductors.cpp
  src/absolute_pose/modules/gp3p/spolynomials.cpp
  src/absolute_pose/modules/Epnp.cpp
  src/absolute_pose/modules/Sqpnp.cpp
  src/absolute_pose/modules/gpnp1/code.cpp
  src/absolute_pose/modules/gpnp1/init.cpp
  src/absolute_pose/modules/gpnp1/reductors.cpp
  src/absolute_pose/modules/gpnp1/spolynomials.cpp
  src/absolute_pose/modules/gpnp2/code.cpp
  src/absolute_pose/modules/gpnp2/init.cpp
  src/absolute_pose/modules/gpnp2/reductors.cpp
  src/absolute_pose/modules/gpnp2/spolynomials.cpp
  src/absolute_pose/modules/gpnp3/code.cpp
  src/absolute_pose/modules/gpnp3/init.cpp
  src/absolute_pose/modules/gpnp3/reductors.cpp
  src/absolute_pose/modules/gpnp3/spolynomials.cpp
  src/absolute_pose/modules/gpnp4/code.cpp
  src/absolute_pose/modules/gpnp4/init.cpp
  src/absolute_pose/modules/gpnp4/reductors.cpp
  src/absolute_pose/modules/gpnp4/spolynomials.cpp
  src/absolute_pose/modules/gpnp5/code.cpp
  src/absolute_pose/modules/gpnp5/init.cpp
  src/absolute_pose/modules/gpnp5/reductors.cpp
  src/absolute_pose/modules/gpnp5/spolynomials.cpp
  src/absolute_pose/modules/upnp2.cpp
  src/absolute_pose/modules/upnp4.cpp
  src/relative_pose/modules/main.cpp
  src/relative_pose/modules/fivept_nister/modules.cpp
  src/relative_pose/modules/fivept_stewenius/modules.cpp
  src/relative_pose/modules/fivept_kneip/code.cpp
  src/relative_pose/modules/fivept_kneip/init.cpp
  src/relative_pose/modules/fivept_kneip/reductors.cpp
  src/relative_pose/modules/fivept_kneip/spolynomials.cpp
  src/relative_pose/modules/sixpt/modules2.cpp
  src/relative_pose/modules/eigensolver/modules.cpp
  src/relative_pose/modules/ge/modules.cpp
  src/math/cayley.cpp
  src/math/quaternion.cpp
  src/math/arun.cpp
  src/math/Sturm.cpp
  src/math/roots.cpp
  src/math/gauss_jordan.cpp
  src/absolute_pose/methods.cpp
  src/absolute_pose/CentralAbsoluteAdapter.cpp
  src/absolute_pose/NoncentralAbsoluteAdapter.cpp
  src/absolute_pose/NoncentralAbsoluteMultiAdapter.cpp
  src/relative_pose/methods.cpp
  src/relative_pose/CentralRelativeAdapter.cpp
  src/relative_pose/CentralRelativeWeightingAdapter.cpp
  src/relative_pose/NoncentralRelativeAdapter.cpp
  src/relative_pose/CentralRelativeMultiAdapter.cpp
  src/relative_pose/NoncentralRelativeMultiAdapter.cpp
  src/triangulation/methods.cpp
  src/point_cloud/methods.cpp
  src/point_cloud/PointCloudAdapter.cpp
  src/sac_problems/absolute_pose/AbsolutePoseSacProblem.cpp
  src/sac_problems/absolute_pose/MultiNoncentralAbsolutePoseSacProblem.cpp
  src/sac_problems/relative_pose/CentralRelativePoseSacProblem.cpp
  src/sac_problems/relative_pose/NoncentralRelativePoseSacProblem.cpp
  src/sac_problems/relative_pose/RotationOnlySacProblem.cpp
  src/sac_problems/relative_pose/TranslationOnlySacProblem.cpp
  src/sac_problems/relative_pose/EigensolverSacProblem.cpp
  src/sac_problems/relative_pose/MultiCentralRelativePoseSacProblem.cpp
  src/sac_problems/relative_pose/MultiNoncentralRelativePoseSacProblem.cpp
  src/sac_problems/point_cloud/PointCloudSacProblem.cpp
  src/absolute_pose/MACentralAbsolute.cpp
  src/absolute_pose/MANoncentralAbsolute.cpp
  src/relative_pose/MACentralRelative.cpp
  src/relative_pose/MANoncentralRelative.cpp
  src/relative_pose/MANoncentralRelativeMulti.cpp
  src/point_cloud/MAPointCloud.cpp )

set( OPENGV_HEADER_FILES
  include/opengv/types.hpp
  include/opengv/OptimizationFunctor.hpp
  include/opengv/absolute_pose/methods.hpp
  include/opengv/relative_pose/methods.hpp
  include/opengv/triangulation/methods.hpp
  include/opengv/point_cloud/methods.hpp
  include/opengv/math/cayley.hpp
  include/opengv/math/quaternion.hpp
  include/opengv/math/arun.hpp
  include/opengv/math/Sturm.hpp
  include/opengv/math/roots.hpp
  include/opengv/math/gauss_jordan.hpp
  include/opengv/absolute_pose/AbsoluteAdapterBase.hpp
  include/opengv/absolute_pose/CentralAbsoluteAdapter.hpp
  include/opengv/absolute_pose/NoncentralAbsoluteAdapter.hpp
  include/opengv/absolute_pose/NoncentralAbsoluteMultiAdapter.hpp
  include/opengv/absolute_pose/AbsoluteMultiAdapterBase.hpp
  include/opengv/relative_pose/RelativeAdapterBase.hpp
  include/opengv/relative_pose/RelativeMultiAdapterBase.hpp
  include/opengv/relative_pose/CentralRelativeAdapter.hpp
  include/opengv/relative_pose/CentralRelativeWeightingAdapter.hpp
  include/opengv/relative_pose/NoncentralRelativeAdapter.hpp
  include/opengv/relative_pose/CentralRelativeMultiAdapter.hpp
  include/opengv/relative_pose/NoncentralRelativeMultiAdapter.hpp
  include/opengv/point_cloud/PointCloudAdapterBase.hpp
  include/opengv/point_cloud/PointCloudAdapter.hpp
  include/opengv/sac_problems/absolute_pose/AbsolutePoseSacProblem.hpp
  include/opengv/sac_problems/absolute_pose/MultiNoncentralAbsolutePoseSacProblem.hpp
  include/opengv/sac_problems/relative_pose/CentralRelativePoseSacProblem.hpp
  include/opengv/sac_problems/relative_pose/NoncentralRelativePoseSacProblem.hpp
  include/opengv/sac_problems/relative_pose/MultiCentralRelativePoseSacProblem.hpp
  include/opengv/sac_problems/relative_pose/MultiNoncentralRelativePoseSacProblem.hpp
  include/opengv/sac_problems/relative_pose/EigensolverSacProblem.hpp
  include/opengv/sac_problems/relative_pose/RotationOnlySacProblem.hpp
  include/opengv/sac_problems/relative_pose/TranslationOnlySacProblem.hpp
  include/opengv/sac_problems/point_cloud/PointCloudSacProblem.hpp
  include/opengv/absolute_pose/MACentralAbsolute.hpp
  include/opengv/absolute_pose/MANoncentralAbsolute.hpp
  include/opengv/relative_pose/MACentralRelative.hpp
  include/opengv/relative_pose/MANoncentralRelative.hpp
  include/opengv/relative_pose/MANoncentralRelativeMulti.hpp
  include/opengv/point_cloud/MAPointCloud.hpp )

add_library( opengv ${OPENGV_SOURCE_FILES} ${OPENGV_HEADER_FILES} )
add_library( random_generators test/random_generators.cpp test/random_generators.hpp test/experiment_helpers.cpp test/experiment_helpers.hpp test/time_measurement.cpp test/time_measurement.hpp )
set_target_properties( opengv random_generators PROPERTIES
                    SOVERSION ${PROJECT_VERSION}
                    VERSION ${PROJECT_VERSION}
                    CXX_STANDARD 14
                    CXX_STANDARD_REQUIRED ON
                    DEBUG_POSTFIX d )

target_include_directories( opengv PUBLIC 
                # only when building from the source tree
                $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
                # only when using the lib from the install path
                $<INSTALL_INTERFACE:include>
                ${ADDITIONAL_INCLUDE_DIRS} )

target_include_directories( random_generators PRIVATE ${ADDITIONAL_INCLUDE_DIRS} )

target_link_libraries( random_generators opengv )

IF (BUILD_TESTS)
  enable_testing()
  include_directories( ${ADDITIONAL_INCLUDE_DIRS} )

  add_executable( test_absolute_pose test/test_absolute_pose.cpp )
  target_link_libraries( test_absolute_pose opengv random_generators )
  add_test(NAME test_absolute_pose
    WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}
    COMMAND test_absolute_pose)

  add_executable( test_sqpnp_omnidirectional test/test_sqpnp_omnidirectional.cpp )
  target_link_libraries( test_sqpnp_omnidirectional opengv random_generators )
  add_test(NAME test_sqpnp_omnidirectional
    WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}
    COMMAND test_sqpnp_omnidirectional)

  add_executable( test_panorama_sqpnp test/test_panorama_sqpnp.cpp )
  target_link_libraries( test_panorama_sqpnp opengv random_generators )
  add_test(NAME test_panorama_sqpnp
    WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH}
    COMMAND test_panorama_sqpnp)

  # Benchmark executable (not a test)
  add_executable( benchmark_panorama_algorithms benchmark/benchmark_panorama_algorithms.cpp )
  target_include_directories( benchmark_panorama_algorithms PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/test )
  target_link_libraries( benchmark_panorama_algorithms opengv random_generators )

  add_executable( test_p3p_backward_vectors test/test_p3p_backward_vectors.cpp )
  target_link_libraries( test_p3p_backward_vectors opengv )
  add_test(NAME test_p3p_backward_vectors
    WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
    COMMAND test_p3p_backward_vectors)

  add_executable( test_absolute_pose_sac test/test_absolute_pose_sac.cpp )
  target_link_libraries( test_absolute_pose_sac opengv random_generators )
  add_test(NAME test_absolute_pose_sac
    WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
    COMMAND test_absolute_pose_sac)

  add_executable( test_noncentral_absolute_pose test/test_noncentral_absolute_pose.cpp )
  target_link_libraries( test_noncentral_absolute_pose opengv random_generators )
  add_test(NAME test_noncentral_absolute_pose
    WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
    COMMAND test_noncentral_absolute_pose)

  add_executable( test_noncentral_absolute_pose_sac test/test_noncentral_absolute_pose_sac.cpp )
  target_link_libraries( test_noncentral_absolute_pose_sac opengv random_generators )
  add_test(NAME test_noncentral_absolute_pose_sac
    WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
    COMMAND test_noncentral_absolute_pose_sac)

  add_executable( test_multi_noncentral_absolute_pose_sac test/test_multi_noncentral_absolute_pose_sac.cpp )
  target_link_libraries( test_multi_noncentral_absolute_pose_sac opengv random_generators )
  add_test(NAME test_multi_noncentral_absolute_pose_sac
    WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
    COMMAND test_multi_noncentral_absolute_pose_sac)

  add_executable( test_relative_pose test/test_relative_pose.cpp )
  target_link_libraries( test_relative_pose opengv random_generators )
  add_test(NAME test_relative_pose
    WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
    COMMAND test_relative_pose)

  add_executable( test_relative_pose_rotationOnly test/test_relative_pose_rotationOnly.cpp )
  target_link_libraries( test_relative_pose_rotationOnly opengv random_generators )
  add_test(NAME test_relative_pose_rotationOnly
    WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
    COMMAND test_relative_pose_rotationOnly)

  add_executable( test_relative_pose_rotationOnly_sac test/test_relative_pose_rotationOnly_sac.cpp )
  target_link_libraries( test_relative_pose_rotationOnly_sac opengv random_generators )
  add_test(NAME test_relative_pose_rotationOnly_sac
    WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
    COMMAND test_relative_pose_rotationOnly_sac)

  add_executable( test_relative_pose_sac test/test_relative_pose_sac.cpp )
  target_link_libraries( test_relative_pose_sac opengv random_generators )
  add_test(NAME test_relative_pose_sac
    WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
    COMMAND test_relative_pose_sac)

  add_executable( test_noncentral_relative_pose test/test_noncentral_relative_pose.cpp )
  target_link_libraries( test_noncentral_relative_pose opengv random_generators )
  add_test(NAME test_noncentral_relative_pose
    WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
    COMMAND test_noncentral_relative_pose)

  add_executable( test_noncentral_relative_pose_sac test/test_noncentral_relative_pose_sac.cpp )
  target_link_libraries( test_noncentral_relative_pose_sac opengv random_generators )
  add_test(NAME test_noncentral_relative_pose_sac
    WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
    COMMAND test_noncentral_relative_pose_sac)

  add_executable( test_multi_noncentral_relative_pose_sac test/test_multi_noncentral_relative_pose_sac.cpp )
  target_link_libraries( test_multi_noncentral_relative_pose_sac opengv random_generators )
  add_test(NAME test_multi_noncentral_relative_pose_sac
    WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
    COMMAND test_multi_noncentral_relative_pose_sac)

  add_executable( test_eigensolver_sac test/test_eigensolver_sac.cpp )
  target_link_libraries( test_eigensolver_sac opengv random_generators )
  add_test(NAME test_eigensolver_sac
    WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
    COMMAND test_eigensolver_sac)

  add_executable( test_triangulation test/test_triangulation.cpp )
  target_link_libraries( test_triangulation opengv random_generators )
  add_test(NAME test_triangulation
    WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
    COMMAND test_triangulation)

  add_executable( test_eigensolver test/test_eigensolver.cpp )
  target_link_libraries( test_eigensolver opengv random_generators )
  add_test(NAME test_eigensolver
    WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
    COMMAND test_eigensolver)

  add_executable( test_point_cloud test/test_point_cloud.cpp )
  target_link_libraries( test_point_cloud opengv random_generators )
  add_test(NAME  test_point_cloud
    WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
    COMMAND test_point_cloud)

  add_executable( test_point_cloud_sac test/test_point_cloud_sac.cpp )
  target_link_libraries( test_point_cloud_sac opengv random_generators )
  add_test(NAME test_point_cloud_sac
    WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
    COMMAND test_point_cloud_sac)

  add_executable( test_Sturm test/test_Sturm.cpp )
  target_link_libraries( test_Sturm opengv random_generators )
  add_test(NAME test_Sturm
    WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}
    COMMAND test_Sturm)

  set_target_properties(
    test_absolute_pose
    test_absolute_pose_sac
    test_noncentral_absolute_pose
    test_noncentral_absolute_pose_sac
    test_multi_noncentral_absolute_pose_sac
    test_relative_pose
    test_relative_pose_rotationOnly
    test_relative_pose_rotationOnly_sac
    test_relative_pose_sac
    test_noncentral_relative_pose
    test_noncentral_relative_pose_sac
    test_multi_noncentral_relative_pose_sac
    test_eigensolver_sac
    test_triangulation
    test_eigensolver
    test_point_cloud
    test_point_cloud_sac
    test_Sturm
    test_sqpnp_omnidirectional
    test_panorama_sqpnp
    benchmark_panorama_algorithms
    # test_panorama_algorithms  # Skipped - fails on MSVC build
    # test_sqpnp_edge_cases  # Skipped - uses normal_distribution with zero stddev
    # test_backward_vector_analysis  # Skipped - SQPnP not suitable for backward vectors
    # test_sqpnp_hybrid_benchmark  # Skipped - uses normal_distribution with zero stddev

    PROPERTIES
    CXX_STANDARD 14
    CXX_STANDARD_REQUIRED ON
    DEBUG_POSTFIX d )


ENDIF()

# Configuration
set(config_install_dir "lib/cmake/${PROJECT_NAME}-${PROJECT_VERSION}")
set(include_install_dir "include")
set(version_config "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake")
set(project_config "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake")
set(targets_export_name "${PROJECT_NAME}Targets")

# Include module with fuction 'write_basic_package_version_file'
include(CMakePackageConfigHelpers)

# Configure '<PROJECT-NAME>ConfigVersion.cmake'
# Note: PROJECT_VERSION is used as a VERSION
write_basic_package_version_file(
    "${version_config}" COMPATIBILITY SameMajorVersion )

# Configure '<PROJECT-NAME>Config.cmake'
# Use variables:
#   * targets_export_name
#   * PROJECT_NAME
configure_package_config_file(
    "modules/Config.cmake.in"
    "${project_config}"
    INSTALL_DESTINATION "${config_install_dir}")

# Export '<PROJECT-NAME>Targets.cmake' to build dir (to find package in build dir without install)
export(TARGETS opengv FILE "${CMAKE_CURRENT_BINARY_DIR}/${targets_export_name}.cmake")

# Targets:
install(TARGETS opengv
        EXPORT "${targets_export_name}"
        LIBRARY DESTINATION "lib"
        ARCHIVE DESTINATION "lib"
        RUNTIME DESTINATION "bin"
        INCLUDES DESTINATION "${include_install_dir}")

# Config
#   * <prefix>/lib/cmake/opengv/opengvConfig.cmake
#   * <prefix>/lib/cmake/opengv/opengvConfigVersion.cmake
install(FILES "${project_config}" "${version_config}"
        DESTINATION "${config_install_dir}")

# Config
#   * <prefix>/lib/cmake/opengv/opengvTargets.cmake
install(EXPORT "${targets_export_name}"
        NAMESPACE "opengv::"
        DESTINATION "${config_install_dir}")

# Headers
install(DIRECTORY include/ 
        DESTINATION ${include_install_dir} 
        FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")


if (BUILD_PYTHON)
  # Find Python components - required for Python bindings
  # Priority: 1) User-specified PYTHON_EXECUTABLE, 2) Conda Python, 3) System Python
  
  # Check if user explicitly set PYTHON_EXECUTABLE
  if(NOT PYTHON_EXECUTABLE)
    # Prefer conda Python if available (common on Windows with conda environments)
    if(DEFINED ENV{CONDA_PREFIX})
      set(_CONDA_PYTHON "$ENV{CONDA_PREFIX}/python.exe")
      if(EXISTS "${_CONDA_PYTHON}")
        set(PYTHON_EXECUTABLE "${_CONDA_PYTHON}" CACHE PATH "Python executable")
        message(STATUS "Using conda Python: ${PYTHON_EXECUTABLE}")
      endif()
    endif()
  endif()
  
  # Try modern find_package (CMake 3.12+)
  # If PYTHON_EXECUTABLE is set, use it as hint
  if(NOT PYTHON_EXECUTABLE)
    find_package(Python3 COMPONENTS Interpreter Development QUIET)
    if(Python3_FOUND)
      set(PYTHON_EXECUTABLE ${Python3_EXECUTABLE} CACHE PATH "Python executable")
      if(Python3_INCLUDE_DIRS)
        set(Python_INCLUDE_DIRS ${Python3_INCLUDE_DIRS})
      endif()
      if(Python3_LIBRARIES)
        set(Python_LIBRARIES ${Python3_LIBRARIES})
      endif()
    endif()
  else()
    # PYTHON_EXECUTABLE is set, find components for it
    find_package(Python3 COMPONENTS Interpreter Development QUIET)
    if(Python3_FOUND AND Python3_EXECUTABLE STREQUAL PYTHON_EXECUTABLE)
      if(Python3_INCLUDE_DIRS)
        set(Python_INCLUDE_DIRS ${Python3_INCLUDE_DIRS})
      endif()
      if(Python3_LIBRARIES)
        set(Python_LIBRARIES ${Python3_LIBRARIES})
      endif()
    else()
      # Fallback: manually find components for specified executable
      execute_process(
        COMMAND "${PYTHON_EXECUTABLE}" -c "import sysconfig; print(sysconfig.get_path('include'))"
        OUTPUT_VARIABLE _PYTHON_INCLUDE_DIR
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
      )
      if(_PYTHON_INCLUDE_DIR)
        set(Python_INCLUDE_DIRS ${_PYTHON_INCLUDE_DIR})
      endif()
    endif()
  endif()
  
  # Fallback: use legacy find_package if needed
  if(NOT PYTHON_EXECUTABLE)
    find_package(PythonInterp REQUIRED)
    find_package(PythonLibs QUIET)
    if(PYTHON_INCLUDE_DIRS)
      set(Python_INCLUDE_DIRS ${PYTHON_INCLUDE_DIRS})
    endif()
    if(PYTHON_LIBRARIES)
      set(Python_LIBRARIES ${PYTHON_LIBRARIES})
    endif()
  endif()
  
  # Ensure Python executable is set
  if(NOT PYTHON_EXECUTABLE)
    message(FATAL_ERROR "Python executable not found. Please install Python or set PYTHON_EXECUTABLE")
  endif()
  
  message(STATUS "Building Python bindings with Python: ${PYTHON_EXECUTABLE}")
  if(Python_INCLUDE_DIRS)
    message(STATUS "Python include dirs: ${Python_INCLUDE_DIRS}")
  endif()
  if(Python_LIBRARIES)
    message(STATUS "Python libraries: ${Python_LIBRARIES}")
  endif()
  
  add_subdirectory( python )
endif()
