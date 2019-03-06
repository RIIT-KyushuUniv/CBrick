//
//  diff3d.h
//  diff3d
//
//  Created by keno on 2016/06/26.
//  Copyright © 2016年 keno. All rights reserved.
//


#ifndef Diff3D_h
#define Diff3D_h

// precision
#ifdef _REAL_IS_DOUBLE_
#define REAL_TYPE double
#else
/** 実数型の指定
 * - デフォルトでは、REAL_TYPE=float
 * - コンパイル時オプション-D_REAL_IS_DOUBLE_を付与することで
 *   REAL_TYPE=doubleになる
 */
#define REAL_TYPE float
#endif


#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <CB_Define.h>
#include "Ffunc.h"

// Physical parameters
typedef struct {
  REAL_TYPE dh;    ///< Mesh width
  REAL_TYPE dt;    ///< Time increment
  REAL_TYPE alpha; ///< coefficient
  REAL_TYPE org[3];///< origin of subdomain
} Phys_Param;


// Cotrol parameters
typedef struct {
  int laststep; ///< Number of steps to calculate
  int fileout;  ///< Interval to output to a file
  int div_mode; ///< 0-IJK, 1-JK
} Cntl_Param;


// Linear Solver, do not use zero. Zero indicates invalid.
enum Linear_Solver {
  Jacobi=1,
  SOR,
  SOR2SMA
};

#endif /* Diff3D_h */
