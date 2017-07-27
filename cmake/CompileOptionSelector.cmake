###################################################################################
#
# CBrick
#
# Copyright (c) 2017 Research Institute for Information Technology(RIIT),
#                    Kyushu University.  All rights reserved.
#
####################################################################################

##
## Compile option selector
##


macro (AddOptimizeOption)
  if (USE_F_TCS STREQUAL "YES")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Kfast,ocl,preex,simd=2,array_private,parallel,optmsg=2 -V -Nsrc -x0 -Xg")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Kfast,ocl,preex,simd=2,array_private,parallel,optmsg=2 -V -Nsrc -x0 -Xg")
    # -Xg   : gcc compatible flag
    # -fPIC : PIC flag
    # -Nfjcex : to link PMlib
    if (with_example)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Cpp -Kfast,ocl,preex,simd=2,array_private,parallel,optmsg=2 -V -Qt")
    endif()

#  elseif (TARGET_ARCH STREQUAL "K")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Kfast,ocl,preex,simd=2,array_private,parallel,optmsg=2 -V -Nsrc -x0 -Xg")
#    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Kfast,ocl,preex,simd=2,array_private,parallel,optmsg=2 -V -Nsrc -x0 -Xg")

  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -Wall")
    if (with_example)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -cpp -Wall --free-line-length-none")
    endif()

  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -qopt-report=3 -DMPICH_IGNORE_CXX_SEEK")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -qopt-report=3")
    if (with_example)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -Warn unused -fpp -qopt-report=5")
    endif()

  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "PGI")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fastsse")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fastsse")
    if (with_example)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3")
    endif()

  else()
    message("using default option")
  endif()
endmacro()


macro (AddSSE)
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    if(CMAKE_SYSTEM_NAME MATCHES "Darwin")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse3")
    else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
    endif()
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -xHost")

  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "PGI")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fastsse")
  endif()
endmacro()


macro (FreeForm)
  if(CMAKE_Fortran_COMPILER MATCHES ".*frtpx$")
    #set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")

  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-form")

  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -free")

  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Mfree")

  endif()
endmacro()


macro(C99)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")
endmacro()


macro(CPP11)
  include(CheckCXXCompilerFlag)
  CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
  CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
  if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
  endif()
endmacro()


macro(checkOpenMP)
  if(enable_OPENMP)
    if(CMAKE_CXX_COMPILER MATCHES ".*FCCpx$")
      set(OpenMP_C_FLAGS "-Kopenmp")
      set(OpenMP_CXX_FLAGS "-Kopenmp")
      if (with_example)
        set(OpenMP_Fortran_FLAGS "-Kopenmp")
      endif()
    elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
      set(OpenMP_C_FLAGS "-qopenmp")
      set(OpenMP_CXX_FLAGS "-qopenmp")
      if (with_example)
        set(OpenMP_Fortran_FLAGS "-qopenmp")
      endif()
    else()
      find_package(OpenMP REQUIRED)
    endif()

    # OpenMP_*_FLAGSにはfind_package(OpenMP REQUIRED)でオプションフラグが設定される
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    if (with_example)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
    endif()
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif()
endmacro()


macro(precision)
  if(real_type STREQUAL "OFF")
  # nothing, default is float
  set(real_type "float")

  elseif(real_type STREQUAL "float")
  # nothing

  elseif(real_type STREQUAL "double")
    ADD_DEFINITIONS(-D_REAL_IS_DOUBLE_)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_REAL_IS_DOUBLE_")

    if (with_example)
      if(CMAKE_Fortran_COMPILER MATCHES ".*frtpx$")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -CcdRR8")

      elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8")

      elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8")

      elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
      endif()

    endif()

  else() # neither 'float' nor 'double'
    message("@@@@@@@@@@@")
    message("FATAL ERROR : Invalid floating type : ${real_type}")
    message("@@@@@@@@@@@")
  ENDIF()
endmacro()
