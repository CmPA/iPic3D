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

set(EXECUTABLE_OUTPUT_PATH exec/${CMAKE_BUILD_TYPE})
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
        main/*.cpp
)

#
# Compilation options
#

set(IPIC_TESTS_DIR "${CMAKE_BINARY_DIR}/tests" CACHE STRING "Location of the source files for iPic3D")

option(IPIC_PARALLEL_HDF5 "Output is done using the parallel HDF5 library" OFF)
if(IPIC_PARALLEL_HDF5)
  add_definitions(-DPHDF5)
endif()

option(IPIC_BATSRUS "Compile library with coupling code for BATS-R-US" OFF)
if(IPIC_BATSRUS)
  add_definitions( -DBATSRUS )
endif()

option(IPIC_XEONPHI "Compile options for Xeon Phi" OFF)
if(IPIC_XEONPHI)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endif()

option(IPIC_TESTS "Set up the code tests" OFF)

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

#build iPic as a library also
#libiPic3Dlib.so in folder lib
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
	 ${EXTRA_LIBS}
)

target_link_libraries(
         iPic3D
         iPic3Dlib
)

## ## to save the executable in the folder where the CMakeLists.txt file is, i.e. CMAKE_CURRENT_SOURCE_DIR
## #set_target_properties(iPic3D PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
## 
## ## debug releases have a _d appended to the executable
## set_target_properties(iPic3D PROPERTIES DEBUG_POSTFIX "_d")
## 
## 
## message("Which system am I compiling for:")
## message("CMAKE_SYSTEM_PROCESSOR is ${CMAKE_SYSTEM_PROCESSOR}")
## 
## message("Compiler & compiler flags:")
## message("CMAKE_CXX_COMPILER is ${CMAKE_CXX_COMPILER}")
## message("CMAKE_CXX_FLAGS is ${CMAKE_CXX_FLAGS}")

#
# Code testing
#

if(IPIC_TESTS)
  enable_testing()

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
