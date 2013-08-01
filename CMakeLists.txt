cmake_minimum_required(VERSION 2.8.8)

#
# Project declaration
#

project(iPic3D)

#
# Set exec path
#

set(EXECUTABLE_OUTPUT_PATH work/${CMAKE_BUILD_TYPE})
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
# Macro definitions
#

set(TEST_PH5 $ENV{USE_PHDF5})
if(DEFINED TEST_PH5)
  add_definitions( -DPHDF5 )
  message(" WARNING: Compiling the code using the Parallel HDF5 functions.")
else(DEFINED TEST_PH5)
  message(" INFO: Compiling using the sequantial HDF5 functions.")
endif(DEFINED TEST_PH5)

set(TEST_B $ENV{BATSRUS})
if(DEFINED TEST_B)
  add_definitions( -DBATSRUS )
  message(" WARNING: BATSRUS flag is active.")
else(DEFINED TEST_B)
  message(" INFO: BATSRUS is not active.")
endif(DEFINED TEST_B)

#
# Executable declaration
#

add_executable(
        iPic3D
        iPic3D.cpp
)

add_library(
        iPic3Dlib
        SHARED
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
