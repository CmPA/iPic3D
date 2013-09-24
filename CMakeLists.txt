cmake_minimum_required(VERSION 2.8.8) 
# compiler set in ../cmake/cmake_template.cmake.XeonPhi
message ("for Xeon Phi:")
message ("cmake .. -DCMAKE_TOOLCHAIN_FILE=../cmake/cmake_template.cmake.XeonPhi")
#message ("for Xeon:")
#message ("cmake .. -DCMAKE_TOOLCHAIN_FILE=../cmake/cmake_template.cmake.Xeon")
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
# Set compiler flags per system
#
if (${CMAKE_SYSTEM_PROCESSOR} STREQUAL "k1om") ## Xeon Phi
   set(CMAKE_CXX_FLAGS "-O3 -openmp -fno-exceptions -vec-report -mmic")
   #set(CMAKE_CXX_FLAGS "-openmp -g -mmic") # set flags for Xeon Phi, totalview
elseif (${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64") ## Xeon
   set(CMAKE_CXX_COMPILER "icpc")  
   set(CMAKE_CXX_FLAGS "-O3 -openmp -fno-exceptions -vec-report")
else()
   set(CMAKE_CXX_FLAGS "-O3")
endif()

#
# Find third class libraries
#
if (${CMAKE_SYSTEM_PROCESSOR} STREQUAL "k1om") ## Xeon Phi
   set(CMAKE_PREFIX_PATH /opt/hdf5/1.8.10-patch1-mic)			
   set(VARIOUS_LIB /opt/intel/lib/mic)
   set (EXTRA_LIBS ${VARIOUS_LIB}/libimf.so ${VARIOUS_LIB}/libsvml.so ${VARIOUS_LIB}/libirng.so ${VARIOUS_LIB}/libintlc.so ${MPELIB})
else()
   find_package(MPI REQUIRED)
   set(EXTRA_LIBS "")   
endif()

find_package(HDF5 COMPONENTS HL C REQUIRED)
message("HDF5_INCLUDE_DIRS is ${HDF5_INCLUDE_DIRS}")


#
# include and lib directories
#

# include_directories: files there are accessible to the project
include_directories(
        include
        ${HDF5_INCLUDE_DIRS}
        ${MPI_INCLUDE_PATH}
)

link_directories(
        ${HDF5_LIBRARY_DIRS}
        ${MPI_LIBRARY_DIRS}
	${EXTRA_LIBS}
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

#build iPic as a library also
#libiPic3Dlib.so in folder lib
add_library(
        iPic3Dlib  #name of the library
        SHARED	   #type of the library
        ${inc_files}    # stuff to build the library
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
	 ${EXTRA_LIBS}
)

target_link_libraries(
         iPic3D
         iPic3Dlib
)


## to save the executable in the folder where the CMakeLists.txt file is, i.e. CMAKE_CURRENT_SOURCE_DIR
set_target_properties(iPic3D PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

## debug releases have a _d appended to the executable
set_target_properties(iPic3D PROPERTIES DEBUG_POSTFIX "_d")


message("Which system am I compiling for:")
message("CMAKE_SYSTEM_PROCESSOR is ${CMAKE_SYSTEM_PROCESSOR}")

message("Compiler & compiler flags:")
message("CMAKE_CXX_COMPILER is ${CMAKE_CXX_COMPILER}")
message("CMAKE_CXX_FLAGS is ${CMAKE_CXX_FLAGS}")

