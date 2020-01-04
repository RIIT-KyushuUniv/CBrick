###################################################################################
#
# CBrick
#
# Copyright (c) 2017-2020 Research Institute for Information Technology(RIIT),
#                    Kyushu University.  All rights reserved.
#
####################################################################################

include_directories(
       ${PROJECT_BINARY_DIR}/src  # CB_Version.h
       ${PROJECT_SOURCE_DIR}/src
)

link_directories()

set(cb_files CB_SubDomain.cpp
             CB_Comm.cpp
   )


add_library(CBrick STATIC ${cb_files})
target_link_libraries(CBrick)
install(TARGETS CBrick DESTINATION lib)

install(FILES
        ${PROJECT_BINARY_DIR}/src/CB_Version.h
        ${PROJECT_SOURCE_DIR}/src/CB_SubDomain.h
        ${PROJECT_SOURCE_DIR}/src/CB_Define.h
        ${PROJECT_SOURCE_DIR}/src/CB_Comm.h
        ${PROJECT_SOURCE_DIR}/src/CB_Comm_inline.h
        ${PROJECT_SOURCE_DIR}/src/CB_PackingScalarCell.h
        ${PROJECT_SOURCE_DIR}/src/CB_PackingVectorCell.h
        ${PROJECT_SOURCE_DIR}/src/CB_PackingScalarNode.h
        ${PROJECT_SOURCE_DIR}/src/CB_PackingVectorNode.h
        ${PROJECT_SOURCE_DIR}/src/CB_Pack.h
        DESTINATION include)
###
