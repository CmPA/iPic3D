cmake_minimum_required(VERSION 2.8.8)

#
# Project declaration
#

project(h5Xdfm)

#
# Set exec path
#

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BUILD_TYPE})

#
# Set executable
#

add_executable(
        h5Xdmf
        h5Xdmf.cpp
)
