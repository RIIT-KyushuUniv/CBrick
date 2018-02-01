###################################################################################
#
# CBrick
#
# Copyright (c) 2017-2018 Research Institute for Information Technology(RIIT),
#                    Kyushu University.  All rights reserved.
#
####################################################################################


include_directories(
      ${PROJECT_SOURCE_DIR}/src
      ${PROJECT_BINARY_DIR}/src
      ${PROJECT_SOURCE_DIR}/example/diff3d
)

link_directories(
      ${PROJECT_BINARY_DIR}/src
)

## float

add_executable(diff3d_f main.cpp fcore.f90)
set_target_properties(diff3d_f PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(diff3d_f -lCBrick_f)

# In the case of Fujitsu--linkfortran is mandatory
if(USE_F_TCS STREQUAL "YES")
  target_link_libraries(diff3d_f ${CMAKE_CXX_IMPLICIT_LINK_LIBRARIES} "--linkfortran")
else()
  target_link_libraries(diff3d_f ${CMAKE_CXX_IMPLICIT_LINK_LIBRARIES})
endif()


## double

#add_executable(diff3d_d main.cpp fcore.f90)

#if(CMAKE_Fortran_COMPILER MATCHES ".*frtpx$")
#  set_target_properties(diff3d_d PROPERTIES COMPILE_FLAGS "-D_REAL_IS_DOUBLE_ -CcdRR8")
#elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
#  set_target_properties(diff3d_d PROPERTIES COMPILE_FLAGS "-D_REAL_IS_DOUBLE_ -fdefault-real-8")
#elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
#  set_target_properties(diff3d_d PROPERTIES COMPILE_FLAGS "-D_REAL_IS_DOUBLE_ -r8")
#elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
#  set_target_properties(diff3d_d PROPERTIES COMPILE_FLAGS "-D_REAL_IS_DOUBLE_")
#  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8")
#endif()

#set_target_properties(diff3d_d PROPERTIES LINKER_LANGUAGE CXX)
#target_link_libraries(diff3d_d -lCBrick_d)

# In the case of Fujitsu--linkfortran is mandatory
#if(USE_F_TCS STREQUAL "YES")
#  target_link_libraries(diff3d_d ${CMAKE_CXX_IMPLICIT_LINK_LIBRARIES} "--linkfortran")
#else()
#  target_link_libraries(diff3d_d ${CMAKE_CXX_IMPLICIT_LINK_LIBRARIES})
#endif()
