cmake_minimum_required(VERSION 2.8.5)

#
# Set project name
#

project(iPic3D)

#
# Find packages (already managed by CMake)
#

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

#
# Executable compilation
#

add_executable(
    iPic3D
    iPic3D.cpp
)

#
# Libraries linked
#

target_link_libraries(
    iPic3D
    iPic3Dlib
    ${MPI_LIBRARIES}
    ${HDF5_LIBRARIES}
)

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
