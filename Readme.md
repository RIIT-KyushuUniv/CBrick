
# CBrick

A brick for building Cartesian applications

## OUTLINE

A class library to help making applications on Cartesian grid, including simple domain management, communication wrapper, and utilities.


## Copyright
- Copyright (c) 2017 Research Institute for Information Technology(RIIT), Kyushu University. All rights reserved.


## Software Requirement

- Cmake
- MPI library


## How to build

### Build

~~~
$ export CB_HOME=/hogehoge
$ mkdir build
$ cd build
$ cmake [options] ..
$ make
$ sudo make install
~~~


#### Note

This Cmake generates CBrick_f for float, and CBrick_d for double precision.


### Options

`-D INSTALL_DIR=` *Install_directory*

>  Specify the directory that this library will be installed. Built library is
   installed at `install_directory/lib` and the header files are placed at
   `install_directory/include`.

`-D enable_OPENMP=` {yes | no}

>  This option makes OpenMP directives effect. Default is yes.


`-D with_example=` {no|yes}



## Configure Examples

`$ export CB_HOME=hogehoge`

### INTEL/GNU compiler

~~~
$ cmake -DINSTALL_DIR=${CB_HOME}/CBrick \
        -Denable_OPENMP=yes \
        -Dwith_example=yes ..
~~~

#### Note
In case of some Intel compiler environment, please specify environemnt variables
`export CC=icc CXX=icpc F90=ifort FC=ifort` before compiling.


### FUJITSU compiler / FX10, FX100, K computer on login nodes (Cross compilation) and Fujitsu TCS environment for intel PC

~~~
$ cmake -DINSTALL_DIR=${CB_HOME}/CBrick \
        -DCMAKE_TOOLCHAIN_FILE=../cmake/Toolchain_fx10.cmake \
        -Denable_OPENMP=yes \
        -Dwith_example=yes ..

$ cmake -DINSTALL_DIR=${CB_HOME}/CBrick \
        -DCMAKE_TOOLCHAIN_FILE=../cmake/Toolchain_fx100.cmake \
        -Denable_OPENMP=yes \
        -Dwith_example=yes ..

$ cmake -DINSTALL_DIR=${CB_HOME}/CBrick \
        -DCMAKE_TOOLCHAIN_FILE=../cmake/Toolchain_K.cmake \
        -Denable_OPENMP=yes \
        -Dwith_example=yes ..

$ cmake -DINSTALL_DIR=${CB_HOME}/CBrick \
        -DCMAKE_TOOLCHAIN_FILE=../cmake/Toolchain_intel_F_TCS.cmake \
        -Denable_OPENMP=yes \
        -Dwith_example=yes ..

~~~

##### Note
- On Fujitsu machines(fx10, K, fx100), confirm appropriate directrory path for compiler environment.
- Before building, execute following command to clean for sure. `$ make distclean`


## Contributors

- Kenji Ono


### comment for RIIT CX400

~~~
$ module load intel2017 intelMPI
$ export CC=mpiicc CXX=mpiicpc F90=mpiifort
~~~

see [here](https://www.cc.kyushu-u.ac.jp/scp/system/library/intel/intel_7_3.html) and [here](https://www.cc.kyushu-u.ac.jp/scp/system/library/intel/intel.html).
