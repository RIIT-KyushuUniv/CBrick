###################################################################################
#
# CBrick
#
# Copyright (c) 2017-2020 Research Institute for Information Technology(RIIT),
#                    Kyushu University.  All rights reserved.
#
####################################################################################


include_directories(
      ${PROJECT_SOURCE_DIR}/src
      ${PROJECT_BINARY_DIR}/src
)

link_directories(
      ${PROJECT_BINARY_DIR}/src
)


add_executable(commtest commtest.cpp)
target_link_libraries(commtest -lCBrick)
set (test_parameters -np 4 "./commtest" "64" "64" "64" "2" "cell" "F")
add_test(NAME test1 COMMAND "mpirun" ${test_parameters})
