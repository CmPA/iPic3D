cmake_minimum_required(VERSION 2.8.5)

#
# Set project name
#

project(iPic3D)

#
# Find packages (already managed by CMake)
#

find_package(Git QUIET REQUIRED)

execute_process(
    COMMAND "${GIT_EXECUTABLE}" describe --always HEAD
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
    RESULT_VARIABLE res
    OUTPUT_VARIABLE MYAPP_VERSION
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE)

set_property(GLOBAL APPEND
    PROPERTY CMAKE_CONFIGURE_DEPENDS
    "${CMAKE_SOURCE_DIR}/.git/index")

string(REGEX REPLACE "^([0-9]+)\\.([0-9]+)\\.([0-9]+).*$"
"\\1;\\2;\\3" _ver_parts "${MYAPP_VERSION}")
#[[list(GET _ver_parts 0 MYAPP_VERSION_MAJOR)
list(GET _ver_parts 1 MYAPP_VERSION_MINOR)
list(GET _ver_parts 2 MYAPP_VERSION_PATCH)]]

if("${MYAPP_VERSION}" MATCHES "^.*-(.*)-g.*$")
    string(REGEX REPLACE "^.*-(.*)-g.*$" "\\1" MYAPP_VERSION_MICRO
"${MYAPP_VERSION}")
else()
    set(MYAPP_VERSION_MICRO "0")
endif()

message("Code revision: ${MYAPP_VERSION}")

find_package(HDF5 COMPONENTS HL C REQUIRED)
find_package(MPI REQUIRED)

#
# Options
#

option(IPIC_TESTS "Set up the code tests" OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Chose the type of build, options
are: None, Debug, Release, RelWithDebInfo, MinSizeRel" FORCE)
endif()

#
# Include sub-directories
#

add_subdirectory(main)

include_directories(
    include
    ${MPI_INCLUDE_PATH}
    ${HDF5_INCLUDE_DIRS}
)

add_definitions("-DMPICH_IGNORE_CXX_SEEK")

#
# Executable compilation
#

add_executable(
    iPic3D-${MYAPP_VERSION}
    iPic3D.cpp
)

#
# Libraries linked
#

target_link_libraries(
    iPic3D-${MYAPP_VERSION}
    iPic3Dlib
    ${MPI_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${HDF5_HL_LIBRARIES}	
)

SET( legacy_link iPic3D)
SET( legacy_target iPic3D-${MYAPP_VERSION})
ADD_CUSTOM_TARGET(link_target ALL
                  COMMAND ${CMAKE_COMMAND} -E create_symlink ${legacy_target} ${legacy_link})

#
# Code testing
#

if(IPIC_TESTS)
  enable_testing()
  set(IPIC_TESTS_DIR "${CMAKE_BINARY_DIR}/tests" CACHE STRING "Location of the source files for iPic3D")

  add_test(NAME GEM-test
           COMMAND ${CMAKE_COMMAND}
           -DIPIC_TESTS_DIR=${IPIC_TESTS_DIR}
           -DIPIC_SOURCE_DIR=${CMAKE_SOURCE_DIR}
           -DIPICEXEC=$<TARGET_FILE:iPic3D>
           -DMPIEXEC=${MPIEXEC}
           -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
           -DMPIEXEC_POSTFLAGS=${MPIEXEC_POSTFLAGS}
           -DIPIC_TESTS_DIR=${IPIC_TESTS_DIR}
           -P ${CMAKE_SOURCE_DIR}/testfiles/CMakeRunTest-GEM.txt)

  add_test(NAME uname-test
           COMMAND ${CMAKE_COMMAND}
           -P ${CMAKE_SOURCE_DIR}/testfiles/CMakeRunTest-uname.txt)

endif(IPIC_TESTS)
