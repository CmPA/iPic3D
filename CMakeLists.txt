cmake_minimum_required(VERSION 2.8.10)

#
# Project declaration
#

project(iPic3D-clean)

#
# Set exec path
#

set(EXECUTABLE_OUTPUT_PATH bin/${CMAKE_BUILD_TYPE})
set(LIBRARY_OUTPUT_PATH lib)

#
# Find third class libraries
#

find_package(HDF5 COMPONENTS C HL REQUIRED)
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
	ConfigFile/src/*.h
        PSKOutput3D/*.h
        bc/*.h
        communication/*.h
        fields/*.h
        grids/*.h
        iPIC3D.h/*.h
        inputoutput/*.h
        mathlib/*.h
        mpidata/*.h
        particles/*.h
        performances/*.h
        processtopology/*.h
        solvers/*.h
        utility/*.h
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
        iPIC3D.cpp
)

#
# Executable declaration
#

add_executable(
        iPIC3D
        ${inc_files}
        ${src_files}
)

#
# Link external libraries
#

target_link_libraries(
         iPIC3D
         ${HDF5_LIBRARIES}
         ${HDF5_HL_LIBRARIES}
         ${MPI_LIBRARIES}
)
