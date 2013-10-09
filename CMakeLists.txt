cmake_minimum_required(VERSION 2.8.8)

#
# Project declaration
#

project(iPic3D)

#
# Set exec path
#

set(EXECUTABLE_OUTPUT_PATH exec/${CMAKE_BUILD_TYPE})
set(LIBRARY_OUTPUT_PATH lib)

#
# Find third class libraries
#

find_package(HDF5 COMPONENTS HL C REQUIRED)
find_package(MPI REQUIRED)

#
# include and lib directories
#

include_directories(
        include
        ${HDF5_INCLUDE_DIRS}
        ${MPI_INCLUDE_PATH}
)

link_directories(
        ${HDF5_LIBRARY_DIRS}
        ${MPI_LIBRARY_DIRS}
)

#
# Header file list
#

file(
        GLOB
        inc_files
        include/*.h
)

#
# Source file list
#

file(
        GLOB
        src_files
        ConfigFile/src/*.cpp
        PSKOutput3D/*.cpp
        bc/*.cpp
        communication/*.cpp
        fields/*.cpp
        grids/*.cpp
        iPIC3D.cpp/*.cpp
        inputoutput/*.cpp
        mathlib/*.cpp
        mpidata/*.cpp
        particles/*.cpp
        performances/*.cpp
        processtopology/*.cpp
        solvers/*.cpp
        utility/*.cpp
        main/iPic3Dlib.cpp
)

#
# Compilation options
#

option(IPIC_PARALLEL_HDF5 "Output is done using the parallel HDF5 library" OFF)
if(IPIC_PARALLEL_HDF5)
  add_definitions(-DUSE_PHDF5)
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
# Executable declaration
#

add_executable(
        iPic3D
        iPic3D.cpp
)

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
         ${HDF5_LIBRARIES}
         ${HDF5_HL_LIBRARIES}
         ${MPI_LIBRARIES}
)

target_link_libraries(
         iPic3D
         iPic3Dlib
)
