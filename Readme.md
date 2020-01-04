
# CBrick

A brick for building Cartesian applications

## OUTLINE

A class library to help making applications on Cartesian grid, including simple domain management, communication wrapper, and utilities.


## Copyright
- Copyright (c) 2017-2020 Research Institute for Information Technology(RIIT), Kyushu University. All rights reserved.


## Software Requirement

- Cmake
- MPI library


## How to build

### Build

~~~
$ export HOME=hogehoge
$ mkdir build
$ cd build
$ export CC=mpiicc CXX=mpiicpc F90=mpiifort // if needed like this line.
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

`-D enable_OPENMP=` {OFF | ON}

>  This option makes OpenMP directives effect. Default is OFF.


`-D with_example=` {OFF | ON}

> Specify the build of test modules.


`-D with_Diagonal=` {OFF | ON}

> Specify diagonal communication.


## Configure Examples

`$ export HOME=hogehoge`

### INTEL/GNU/PGI compiler

~~~
$ cmake -DINSTALL_DIR=${HOME}/CBrick -Denable_OPENMP=yes -Dwith_example=yes -Dwith_Diagonal=yes ..
~~~

#### Note
In case of some Intel compiler environment, please specify environemnt variables
`export CC=icc CXX=icpc F90=ifort FC=ifort` before compiling.
Also, for PGI compiler, `export CC=mpicc CXX=mpic++ F90=mpif90 FC=mpif90`


### Mac

#### gnu

~~~
$ module load gcc/9.2.0 openmpi/4.0.2-gcc9
$ cmake -DINSTALL_DIR=${HOME}/CBrick -Denable_OPENMP=ON ..
~~~


### FUJITSU compiler / FX10, FX100, K computer on login nodes (Cross compilation) and Fujitsu TCS environment for intel PC

~~~
$ cmake -DINSTALL_DIR=${HOME}/CBrick \
        -DCMAKE_TOOLCHAIN_FILE=../cmake/Toolchain_fx10.cmake \
        -Denable_OPENMP=yes \
        -Dwith_example=yes ..

$ cmake -DINSTALL_DIR=${HOME}/CBrick \
        -DCMAKE_TOOLCHAIN_FILE=../cmake/Toolchain_fx100.cmake \
        -Denable_OPENMP=yes \
        -Dwith_example=yes ..

$ cmake -DINSTALL_DIR=${HOME}/CBrick \
        -DCMAKE_TOOLCHAIN_FILE=../cmake/Toolchain_K.cmake \
        -Denable_OPENMP=yes \
        -Dwith_example=yes ..
~~~


### ITO susbsystem A/B

#### Intel compiler

~~~
$ module load intel/2018.3 openmpi/3.1.3-nocuda-intel18.3

$ cmake -DINSTALL_DIR=${HOME}/CBrick -Denable_OPENMP=ON ..
~~~

#### Fujitsu TCS

~~~
$ cmake -DINSTALL_DIR=${HOME}/CBrick \
        -DCMAKE_TOOLCHAIN_FILE=../cmake/Toolchain_ITO_TCS.cmake \
        -Denable_OPENMP=ON ..
~~~

##### Note
- On Fujitsu machines(fx10, K, fx100), confirm appropriate directrory path for compiler environment.
- Before building, execute following command to clean for sure. `$ make distclean`


## Contributors

- Kenji Ono

