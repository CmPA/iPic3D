cmake_minimum_required(VERSION 2.8.8)

#
# Set project name
#

project(h5hut_test)

#
# Find packages (already managed by CMake)
#

find_package(HDF5 COMPONENTS HL C REQUIRED)
find_package(MPI REQUIRED)

#
# Set the location of iPic3D
#

set(IPIC_HOME $ENV{IPIC_HOME} CACHE STRING "Path to iPic3D source directory")
if("${IPIC_HOME}" STREQUAL "")
  message(FATAL_ERROR "ERROR: Please enter the path to the iPic3D source directory in the IPIC_HOME varialbe.")
endif()

#
# Set include/library paths to libraries (not managed by CMake)
#

set(IPIC_INCLUDE_DIRS    ${IPIC_HOME}/include)
set(H5HUT_INCLUDE_DIRS   $ENV{H5HUT_DIR}/include)
set(H5HUT_LIBRARY_DIRS   $ENV{H5HUT_DIR}/lib)

#
# Include directories
#

add_subdirectory(H5hut-io)
set(H5HUTIO_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/H5hut-io/include)

include_directories(
    include
    ${IPIC_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}
    ${H5HUT_INCLUDE_DIRS}
    ${H5HUTIO_INCLUDE_DIRS}
)

#
# External libraries location (not managed by CMake)
#

link_directories(
    ${H5HUT_LIBRARY_DIRS}
)

#
# Executable compilation
#

add_executable(
    h5hut-test
    src/h5hut-test.cpp
)

#
# Libraries linked
#

target_link_libraries(
    h5hut-test
    H5hut-io
    H5hut
    ${HDF5_LIBRARIES}
    ${MPI_LIBRARIES}
)
