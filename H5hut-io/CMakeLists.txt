cmake_minimum_required(VERSION 2.8.8)

#
# Project declaration
#

project(H5hut-io)

#
# Set exec path
#

set(LIBRARY_OUTPUT_PATH lib)

#
# Find third class libraries
#

find_package(HDF5 COMPONENTS C REQUIRED)
find_package(MPI REQUIRED)

set(H5HUT_HOME $ENV{H5HUT_DIR} CACHE STRING "H5hut home directory")
set(H5HUT_INCLUDE_DIRS ${H5HUT_HOME}/include)
set(H5HUT_LIBRARY_DIRS ${H5HUT_HOME}/lib)
set(H5HUT_LIBRARIES ${H5HUT_LIBRARY_DIRS}/libH5hut.a)

#
# include and lib directories
#

include_directories(
        include
        ${H5HUT_INCLUDE_DIRS}
        ${HDF5_INCLUDE_DIRS}
        ${MPI_INCLUDE_PATH}
)

link_directories(
        ${H5HUT_LIBRARY_DIRS}
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
        src/*.cpp
)

#
# Compilation options
#

add_library(
        H5hut-io
        STATIC
        ${src_files}
)

#
# Link external libraries
#

target_link_libraries(
         H5hut-io
         ${H5HUT_LIBRARIES}
)

