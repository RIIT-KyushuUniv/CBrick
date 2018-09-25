###################################################################################
#
# CBrick
#
# Copyright (c) 2017-2018 Research Institute for Information Technology(RIIT),
#                    Kyushu University.  All rights reserved.
#
####################################################################################

##
## Compile option selector
##


macro (AddOptimizeOption)
  if (TARGET_ARCH STREQUAL "INTEL_F_TCS")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Kfast,parallel,optmsg=2 -V -Xg")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Kfast,parallel,optmsg=2 -V -Xg")
    if (with_example)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Cpp -Kfast,optmsg=2 -V -Qt")
    endif()

  elseif (USE_F_TCS STREQUAL "YES")
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
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fastsse -Mipa=fast,inline -O4 -Minfo=intensity,mp,vect")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fastsse -Mipa=fast,inline -O4 -Minfo=intensity,mp,vect")
    if (with_example)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fastsse -Mipa=fast,inline -O4 -Minfo=intensity,mp,vect")
    endif()

  else()
    message("using default option")
  endif()
endmacro()


macro (FreeForm)
  if(CMAKE_Fortran_COMPILER MATCHES ".*frtpx$")
    #set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")

  elseif(TARGET_ARCH STREQUAL "INTEL_F_TCS")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Free")

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
    if(USE_F_TCS STREQUAL "YES")
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
