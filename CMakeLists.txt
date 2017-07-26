###################################################################################
#
# CBrick
#
# Copyright (c) 2017 Research Institute for Information Technology(RIIT),
#                    Kyushu University.  All rights reserved.
#
####################################################################################
#
#-*- CMAKE -*-#
# Process this file with cmake to produce a configure script.

# Following parameters must be given at cmake command line
# INSTALL_DIR, CXX, CC

# Build options, see Readme.md for details.
#
# -D INSTALL_DIR=install_directory
#
# -D real_type=(float|double)
#
# -D enable_OPENMP=(no|yes)
#
# -D with_example= {no|yes}
#
# -D with_DEBUG= {no|yes}

cmake_minimum_required(VERSION 2.6)


#######
# Build options
#######

option(enable_OPENMP "Enable OpenMP" "OFF")
option(real_type "Type of floating point" "OFF")
option(enable_DEBUG "Debug mode" "OFF")


#######
# Project setting
#######

project(CBrick C CXX)

if(INSTALL_DIR STREQUAL "")
  set(INSTALL_DIR "/usr/local/CBrick")
endif()

# project()でCMAKE_INSTALL_PREFIX=/usr/localとなるので、再定義
set (CMAKE_INSTALL_PREFIX "${INSTALL_DIR}")

if (enable_DEBUG)
  set(CMAKE_BUILD_TYPE Debug)
else()
  set(CMAKE_BUILD_TYPE Release)
endif()

#######
set(PROJECT_VERSION "0.3.0")
set(LIB_REVISION "20170726_1200")
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

# OpennMP
if (enable_OPENMP STREQUAL "yes")
  checkOpenMP()
endif()

# Real type
precision()


#######
# Display options
#######

message(" ")
message("${PROJECT_NAME} Version : ${PROJECT_VERSION} : ${LIB_REVISION}")
message( STATUS "Destination PATH      : "      ${CMAKE_INSTALL_PREFIX})
message( STATUS "CMAKE_CXX_COMPILER    : "      ${CMAKE_CXX_COMPILER})
message( STATUS "CMAKE_CXX_FLAGS       : "      ${CMAKE_CXX_FLAGS})

message( STATUS "OpenMP support        : "      ${enable_OPENMP})
message( STATUS "Type of floating point: "      ${real_type})
message( STATUS "Debug mode            : "      ${enable_DEBUG})
message( STATUS "Example               : "      ${with_example})
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

# Debug
IF(enable_DEBUG)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_DEBUG")
ENDIF()



#######
# Sub Directories
#######

add_subdirectory(src)
add_subdirectory(doc)

if (with_example STREQUAL "yes")
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
