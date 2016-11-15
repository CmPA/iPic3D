cmake_minimum_required(VERSION 2.8.5)

#
# Set project name
#

project(iPic3D CXX Fortran)

set (FFLAGS -std="f2003" -fbounds-check -Wall) #not sure this is enough for f03; may be needed to un-comment the line down as well
#set (CMAKE_Fortran_FLAGS -std="f2003") #to have it compile for f2003
set (CMAKE_Fortran_FLAGS -std="f2003" -fbounds-check -Wall)

if (CMAKE_Fortran_COMPILER MATCHES "gfortran.*") #for same reason, ${Fortran_COMPILER_NAME} is empty
  # gfortran
  #set (CMAKE_Fortran_FLAGS "-Wtabs -ffree-form -ffree-line-length-none") #-Wtabs: to disable annoying tab warnings # -ffree-form to eliminate line length limit
  set (CMAKE_Fortran_FLAGS " -ffree-form -ffree-line-length-none") #-Wtabs: to disable annoying tab warnings # -ffree-form to eliminate line length limit 
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -Wall")
endif (CMAKE_Fortran_COMPILER MATCHES "gfortran.*")


message ("CMAKE _Fortran _FLAGS_RELEASE: " ${CMAKE_Fortran_FLAGS_RELEASE})
message ("CMAKE_Fortran_FLAGS: " ${CMAKE_Fortran_FLAGS})
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

# to see comments
#add_definitions ("-Wall")


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
