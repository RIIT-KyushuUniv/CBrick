set(CMAKE_SYSTEM_NAME Linux)

include(CMakeForceCompiler)

CMAKE_FORCE_C_COMPILER(mpifccpx GNU)
CMAKE_FORCE_CXX_COMPILER(mpiFCCpx GNU)

set(CMAKE_FIND_ROOT_PATH /opt/FJSVfxlang/1.2.1)   # RIIT fx10, hayaka
set(CMAKE_INCLUDE_PATH /opt/FJSVfxlang/1.2.1/include)
set(CMAKE_LIBRARY_PATH /opt/FJSVfxlang/1.2.1/lib64)

set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

set(TARGET_ARCH "FX100")
set(USE_F_TCS "YES")
