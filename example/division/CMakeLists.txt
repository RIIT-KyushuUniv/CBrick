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


## Cell: G_size, numProc

add_executable(div1 div1.cpp)
target_link_libraries(div1 -lCBrick)
set (test_parameters -np 16 "./div1" "30" "24" "24" "0")
add_test(NAME test1 COMMAND "mpirun" ${test_parameters})


## Cell: G_size, G_div, numProc

add_executable(div2 div2.cpp)
target_link_libraries(div2 -lCBrick)
set (test_parameters -np 16 "./div2" "30" "24" "24" "2" "4" "2" "0")
add_test(NAME test2 COMMAND "mpirun" ${test_parameters})


## Node: G_size, numProc

add_executable(div3 div3.cpp)
target_link_libraries(div3 -lCBrick)
set (test_parameters -np 16 "./div3" "30" "24" "24" "0")
add_test(NAME test3 COMMAND "mpirun" ${test_parameters})


## Node: G_size, G_div, numProc

add_executable(div4 div4.cpp)
target_link_libraries(div4 -lCBrick)
set (test_parameters -np 16 "./div4" "30" "24" "24" "2" "4" "2" "0")
add_test(NAME test4 COMMAND "mpirun" ${test_parameters})
