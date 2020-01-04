###################################################################################
#
# CBrick
#
# Copyright (c) 2017-2019 Research Institute for Information Technology(RIIT),
#                    Kyushu University.  All rights reserved.
#
####################################################################################
#
#-*- CMAKE -*-#
# Process this file with cmake to produce a configure script.

# Following parameters must be given at cmake command line
# INSTALL_DIR, CXX, CC, F90

# Build options, see Readme.md for details.
#
# -D INSTALL_DIR=install_directory
#
# -D enable_OPENMP=(OFF | ON)
#
# -D with_example= {OFF | ON}
#
# -D with_Diagonal={OFF | ON}
#

cmake_minimum_required(VERSION 2.6)


#######
# Build options
#######

option(enable_OPENMP "Enable OpenMP" "OFF")
option(with_Diagonal "Enable Diagonal communication" "OFF")
option(with_example "Build Example" "OFF")

# Default
set(with_MPI "ON")



#######
# Project setting
#######

if (with_example)
  project(CBrick C CXX Fortran)
else()
  project(CBrick C CXX)
endif()

if(INSTALL_DIR STREQUAL "")
  set(INSTALL_DIR "/usr/local/CBrick")
endif()

# project()でCMAKE_INSTALL_PREFIX=/usr/localとなるので、再定義
set (CMAKE_INSTALL_PREFIX "${INSTALL_DIR}")

#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)


#######
# Policy setting
#######

cmake_policy(SET CMP0012 NEW)


#######
set(PROJECT_VERSION "1.4.3")
set(LIB_REVISION "20200104_1500")
#######



#######
# Cmake module
#######

set(CMAKE_MODULE_PATH  ${PROJECT_SOURCE_DIR}/cmake)
include(${PROJECT_SOURCE_DIR}/cmake/CompileOptionSelector.cmake)
include(${PROJECT_SOURCE_DIR}/cmake/CBrick_package.cmake)



#######
# Options
#######

AddOptimizeOption()


# Fortran free form
if (with_example)
  FreeForm()
endif()


# OpennMP
if (enable_OPENMP STREQUAL "ON")
  checkOpenMP()
endif()

# Diagonal communication
if (with_Diagonal)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_DIAGONAL_COMM")
endif()


#######
# Display options
#######

message(" ")
message("${PROJECT_NAME} Version : ${PROJECT_VERSION} : ${LIB_REVISION}")
message( STATUS "Destination PATH      : "      ${CMAKE_INSTALL_PREFIX})
message( STATUS "CMAKE_CXX_COMPILER    : "      ${CMAKE_CXX_COMPILER})
message( STATUS "CMAKE_CXX_FLAGS       : "      ${CMAKE_CXX_FLAGS})
if (with_example)
  message( STATUS "CMAKE_Fortran_COMPILER: "      ${CMAKE_Fortran_COMPILER})
  message( STATUS "CMAKE_Fortran_FLAGS   : "      ${CMAKE_Fortran_FLAGS})
endif()

message( STATUS "OpenMP support        : "      ${enable_OPENMP})
message( STATUS "Example               : "      ${with_example})
message( STATUS "Diagonal comm         : "      ${with_Diagonal})
message(" ")


if(USE_F_TCS STREQUAL "YES")
else()
  if(with_MPI)
    include(FindMPI)

    if(MPI_CXX_FOUND)
      message("MPI-CXX-CMP  : ${MPI_CXX_COMPILER}")
      message("MPI-CXX-FLAGS: ${MPI_CXX_COMPILE_FLAGS}")
      message("MPI-CXX-INC  : ${MPI_CXX_INCLUDE_PATH}")
      message("MPI-LNK-FLAGS: ${MPI_CXX_LINK_FLAGS}")
      message("MPI-LIBS     : ${MPI_CXX_LIBRARIES}")
    endif(MPI_CXX_FOUND)
  endif()
endif()



#######
# find Libraries to depend on
#######




#######
# Sub Directories
#######

add_subdirectory(src)
add_subdirectory(doc)

if (with_example STREQUAL "ON")
  enable_testing()
  add_subdirectory(example)
endif()


#######
# configure files
#######

configure_file( src/CB_Version.h.in ${PROJECT_BINARY_DIR}/src/CB_Version.h @ONLY)



#######
# Install files
#######

install(FILES  ChangeLog.md
               Readme.md
               License.txt
               ${PROJECT_BINARY_DIR}/CMakeCache.txt
        DESTINATION share)



#######
# Installer
#######

CBrick_package()
