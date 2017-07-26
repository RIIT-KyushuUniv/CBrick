###################################################################################
#
# RIAM-COMPACT
#
# Copyright (C) 2015-2017 Research Institute for Applied Mechanics(RIAM)
#                       / Research Institute for Information Technology(RIIT), Kyushu University.
# All rights reserved.
#
# Copyright (C) 2015-2016 Advanced Institute for Computational Science, RIKEN.
# All rights reserved.
#
###################################################################################


include_directories(
       ${PROJECT_BINARY_DIR}/src/cpp_src  # rcVersion.h
       ${PROJECT_SOURCE_DIR}/src/cpp_src
       ${TP_INC}
       ${CPM_INC}
       ${CDM_INC}
)

link_directories(
      ${TP_LIB}
      ${CPM_LIB}
      ${CDM_LIB}
)

IF(with_PM)
  include_directories(${PM_INC})
  link_directories(${PM_LIB})
endif()


set(rc_files
       riamc_allocateArray.cpp
       riamc_Initialize.cpp
       riamc.cpp
       riamc_util.cpp
       #IterationControl.cpp
       #riamc_LS.cpp
)

add_library(RC STATIC ${rc_files})


if(with_MPI)

  IF(with_CDM)
    if(with_PM)
      SET(rc_libs "-lCDMmpi -lCPMmpi -lPMmpi -lTPmpi")
    else()
      SET(rc_libs "-lCDMmpi -lCPMmpi -lTPmpi")
    endif()
  ELSE()
    if(with_PM)
      SET(rc_libs "-lCPMmpi -lPMmpi -lTPmpi")
    else()
      SET(rc_libs "-lCPMmpi -lTPmpi")
    endif()
  ENDIF()

else()

  IF(with_CDM)
    if(with_PM)
      SET(rc_libs "-lCDM -lCPM -lPM -lTP")
    else()
      SET(rc_libs "-lCDM -lCPM -lTP")
    endif()
  ELSE()
    if(with_PM)
      SET(rc_libs "-lCPM -lPM -lTP")
    else()
      SET(rc_libs "-lCPM -lTP")
    endif()
  ENDIF()

endif()

target_link_libraries(RC "${rc_libs}")

###
