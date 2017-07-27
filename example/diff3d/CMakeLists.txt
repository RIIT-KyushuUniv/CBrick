###################################################################################
#
# CBrick
#
# Copyright (c) 2017 Research Institute for Information Technology(RIIT),
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


add_executable(diff3d main.cpp fcore.f90)

set_target_properties(diff3d PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries(diff3d -lCBrick)


# In the case of Fujitsu--linkfortran is mandatory
if(USE_F_TCS STREQUAL "YES")
  target_link_libraries(diff3d ${CMAKE_CXX_IMPLICIT_LINK_LIBRARIES} "--linkfortran")
else()
  target_link_libraries(diff3d ${CMAKE_CXX_IMPLICIT_LINK_LIBRARIES})
endif()
