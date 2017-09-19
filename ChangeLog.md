# CBrick

## TODO
- padding for inner most loop

## REVISION HISTORY

---
- 2017-09-20 Version 0.9.6
  - bug fix : registerCandidates4JK_Cell() / _Node() >>  tbl[odr].dsz[0]=G_size[0]


---
- 2017-09-04 Version 0.9.5
  - bug fix : registerCandidates() >> registerCandidates_Cell() / _Node()
  - registerCandidates4JK_Cell() / _Node()


---
- 2017-09-04 Version 0.9.4
  - virtual instruction modifier for destructor
  - JK division

---
- 2017-08-31 Version 0.9.3
  - use CMAKE_BUILD_TYPE for Debug > > define -DNDEBUG

---
- 2017-08-30 Version 0.9.2
  - change implementation of createRankTable() for OpenMP
  - add omp single

---
- 2017-08-30 Version 0.9.1
  - change instance method of cntl_tbl class
  - add buf_flag for delete pointer

---
- 2017-08-24 Version 0.9.0
  - Vector pack/unpack
  - remove blocking communication


---
- 2017-08-08 Version 0.8.0
  - build both float and double libraries


---
- 2017-08-03 Version 0.7.4
  - add stub for MPI_Comm and MPI_Request


---
- 2017-08-03 Version 0.7.3
  - add #ifndef DISABLE_MPI for serial build of applications


---
- 2017-08-03 Version 0.7.2
  - found bug in blocking communication method, open issue.


---
- 2017-08-03 Version 0.7.1
  - expire use_NB


---
- 2017-08-03 Version 0.7.0
  - modify G2L_index()


---
- 2017-08-03 Version 0.6.9
  - modify getHeadIndex()


---
- 2017-08-02 Version 0.6.8
  - write size, head info into div_process.txt


---
- 2017-08-02 Version 0.6.7
  - bug fix of moemory release in findOptimalDivision() and findParameter()


---
- 2017-07-28 Version 0.6.6
  - change scope of Comm_V_blocking() to public


---
- 2017-07-28 Version 0.6.5
  - setSubDomain()


---
- 2017-07-28 Version 0.6.4
  - useNonblocking()


---
- 2017-07-28 Version 0.6.3
  - CB_Version.h


---
- 2017-07-28 Version 0.6.2
  - bug fix : if ( !f ) delete[] f; => if (f) ...


---
- 2017-07-28 Version 0.6.1
  - fix install doc file


---
- 2017-07-27 Version 0.6.0
  - install header files


---
- 2017-07-27 Version 0.5.0
  - convert to Fortran index >> f_index


---
- 2017-07-27 Version 0.4.0
  - example/diff3d
  - add Fortran env in case of turning on example
  - introduce div_mode, and merge findParameter() into findOptimalDivision()


---
- 2017-07-26 Version 0.3.0
  - pack/unpack
  - CB_CommV.cpp
  - sortLenX(), sortCube()
  - findParameter()
  - example/div3.cpp, div4.cpp


---
- 2017-07-23 Version 0.2.1
  - CB_Comm.cpp


---
- 2017-07-22 Version 0.2.0
  - createRankTable()

---
- 2017-07-20 Version 0.1.0
  - Cmake compilation
  - Node and Cell


---
- 2017-06-29 Version 0.0.1
  - initial material
  - serial compile by makefile
    - `CB_SubDomain.cpp` and `example/division/div.cpp`
