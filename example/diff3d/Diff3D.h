//
//  diff3d.h
//  diff3d
//
//  Created by keno on 2016/06/26.
//  Copyright © 2016年 keno. All rights reserved.
//


#ifndef Diff3D_h
#define Diff3D_h

#include <math.h>
#include <stdio.h>
#include <CB_Define.h>

// Physical parameters
typedef struct {
  REAL_TYPE dh;    ///< Mesh width
  REAL_TYPE dt;    ///< Time increment
  REAL_TYPE alpha; ///< coefficient
//  REAL_TYPE org[3];///< origin of subdomain
} Phys_Param;


// Cotrol parameters
typedef struct {
//  int sz[3];    ///< size of subdomain
//  int div[3];   ///< Number of division for each direction
  int laststep; ///< Number of steps to calculate
  int fileout;  ///< Interval to output to a file
  int blocking; ///< 0-blocking, 1-nonblocking
} Cntl_Param;


// Linear Solver, do not use zero. Zero indicates invalid.
enum Linear_Solver {
  Jacobi=1,
  SOR,
  SOR2SMA
};

class Diff3D : public DomainInfo {

public:


  Diff3D() {

  }

  ~Diff3D() {}

public:

  // @fn alloc_int
  // @brief allocatin of integer array
  // @param [in]     sz Size of array
  // @retval pointer
  int* alloc_int(const size_t sz)
  {
    int* p=NULL;
    if ( !(p=(int*)malloc(sizeof(int)*sz)) ) {
      printf("fail to allocate memory\n");
    }
    return p;
  }


  // @fn alloc_real
  // @brief allocatin of real array
  // @param [in]     sz Size of array
  // @retval pointer
  // @note Type of real is defined by REAL_TYPE
  REAL_TYPE* alloc_real(const size_t sz)
  {
    REAL_TYPE* p=NULL;
    if ( !(p=(REAL_TYPE*)malloc(sizeof(REAL_TYPE)*sz)) ) {
      printf("fail to allocate memory\n");
    }
    return p;
  }
};

#endif /* Diff3D_h */
