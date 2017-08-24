###################################################################################
#
# CBrick
#
# Copyright (c) 2017 Research Institute for Information Technology(RIIT),
#                    Kyushu University.  All rights reserved.
#
####################################################################################

include_directories(
       ${PROJECT_BINARY_DIR}/src  # CB_Version.h
       ${PROJECT_SOURCE_DIR}/src
)

link_directories()

set(cb_files CB_SubDomain.cpp
             CB_CommS.cpp
             CB_CommV.cpp 
             CB_PackingScalar.cpp
             CB_PackingVector.cpp)


## float

add_library(CBrick_f STATIC ${cb_files})
target_link_libraries(CBrick_f)
install(TARGETS CBrick_f DESTINATION lib)


## double

add_library(CBrick_d STATIC ${cb_files})
set_target_properties(CBrick_d PROPERTIES COMPILE_FLAGS "-D_REAL_IS_DOUBLE_")
target_link_libraries(CBrick_d)
install(TARGETS CBrick_d DESTINATION lib)

install(FILES
        ${PROJECT_BINARY_DIR}/src/CB_Version.h
        ${PROJECT_SOURCE_DIR}/src/CB_SubDomain.h
        ${PROJECT_SOURCE_DIR}/src/CB_Define.h
        DESTINATION include)
###
