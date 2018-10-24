cmake_minimum_required(VERSION 3.1.3)
project(opengv VERSION 1.0 LANGUAGES CXX)

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

IF(MSVC)
  add_compile_options(/wd4514 /wd4267 /bigobj)
  add_definitions(-D_USE_MATH_DEFINES)
ELSE()
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

IF(BUILD_POSITION_INDEPENDENT_CODE)
  add_definitions( -fPIC )
ENDIF()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/modules/")
find_package(Eigen REQUIRED)
set(ADDITIONAL_INCLUDE_DIRS ${EIGEN_INCLUDE_DIRS} ${EIGEN_INCLUDE_DIR}/unsupported)

set( OPENGV_SOURCE_FILES
  src/absolute_pose/modules/main.cpp
  src/absolute_pose/modules/gp3p/code.cpp
  src/absolute_pose/modules/gp3p/init.cpp
  src/absolute_pose/modules/gp3p/reductors.cpp
  src/absolute_pose/modules/gp3p/spolynomials.cpp
  src/absolute_pose/modules/Epnp.cpp
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
                    CXX_STANDARD 11
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

    PROPERTIES
    CXX_STANDARD 11
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
        NAMESPACE "${namespace}"
        DESTINATION "${config_install_dir}")

# Headers
install(DIRECTORY include/ 
        DESTINATION ${include_install_dir} 
        FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")


if (BUILD_PYTHON)
  add_subdirectory( python )
endif()
