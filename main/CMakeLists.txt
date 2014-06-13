cmake_minimum_required(VERSION 2.8.5)

#
# Project declaration
#

project(iPic3Dlib)

#
# Set exec path
#

set(EXECUTABLE_OUTPUT_PATH exec/${CMAKE_BUILD_TYPE})
set(LIBRARY_OUTPUT_PATH lib)
set(IPIC_HOME ${CMAKE_CURRENT_SOURCE_DIR}/..)

#
# Find third class libraries
#

find_package(HDF5 COMPONENTS HL C REQUIRED)
find_package(MPI REQUIRED)

#
# Header file list
#

file(
        GLOB
        inc_files
        ${IPIC_HOME}/include/*.h
)

#
# Source file list
#

file(
        GLOB
        src_files
        ${IPIC_HOME}/ConfigFile/src/*.cpp
        ${IPIC_HOME}/PSKOutput3D/*.cpp
        ${IPIC_HOME}/bc/*.cpp
        ${IPIC_HOME}/communication/*.cpp
        ${IPIC_HOME}/fields/*.cpp
        ${IPIC_HOME}/grids/*.cpp
        ${IPIC_HOME}/iPIC3D.cpp/*.cpp
        ${IPIC_HOME}/inputoutput/*.cpp
        ${IPIC_HOME}/mathlib/*.cpp
        ${IPIC_HOME}/mpidata/*.cpp
        ${IPIC_HOME}/particles/*.cpp
        ${IPIC_HOME}/performances/*.cpp
        ${IPIC_HOME}/processtopology/*.cpp
        ${IPIC_HOME}/solvers/*.cpp
        ${IPIC_HOME}/utility/*.cpp
        ${IPIC_HOME}/main/iPic3Dlib.cpp
)

#
# Compilation options
#

option(IPIC_H5HUT_OUTPUT "Output is done using Parallel HDF5 with the H5hut library")
if(IPIC_H5HUT_OUTPUT)
  set(H5HUT_HOME $ENV{H5HUT_DIR} CACHE STRING "H5hut home directory")
  set(H5HUT_INCLUDE_DIRS ${H5HUT_HOME}/include)

  set(H5HUTIO_HOME ${IPIC_HOME}/H5hut-io CACHE STRING "H5hut-io library home")
  set(H5HUTIO_INCLUDE_DIRS ${H5HUTIO_HOME}/include)
  set(H5HUTIO_LIBRARY_DIRS ${H5HUTIO_HOME}/lib)
  set(H5HUTIO_LIBRARIES H5hut-io)

  add_definitions(-DUSEH5HUT)
endif()

option(IPIC_BATSRUS "Compile library with coupling code for BATS-R-US" OFF)
if(IPIC_BATSRUS)
  add_definitions( -DBATSRUS )
endif()

option(IPIC_XEONPHI "Compile options for Xeon Phi" OFF)
if(IPIC_XEONPHI)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endif()

option(IPIC_BUILD_SHARED "Compile shared library" OFF)
if(IPIC_BUILD_SHARED)
  set(IPIC_BUILD_TYPE SHARED)
else()
  set(IPIC_BUILD_TYPE STATIC)
endif()

#
# include and lib directories
#

include_directories(
        ${IPIC_HOME}/include
        ${H5HUTIO_INCLUDE_DIRS}
        ${H5HUT_INCLUDE_DIRS}
        ${HDF5_INCLUDE_DIRS}
        ${MPI_INCLUDE_PATH}
)

link_directories(
        ${H5HUTIO_LIBRARY_DIRS}
)

#
# Create the library
#

add_library(
        iPic3Dlib
        ${IPIC_BUILD_TYPE}
        ${inc_files}
        ${src_files}
)

#
# Link external libraries
#

target_link_libraries(
         iPic3Dlib
         ${H5HUTIO_LIBRARIES}
         ${HDF5_LIBRARIES}
         ${MPI_LIBRARIES}
         ${CMAKE_DL_LIBS}
)
