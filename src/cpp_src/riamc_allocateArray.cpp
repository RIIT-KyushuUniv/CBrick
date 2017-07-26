/*
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
*/

#include "riamc.h"


// #################################################################
/**
 * @brief V3D配列のアロケート
 * @param [in] sz 配列サイズ
 * @ret pointer
 */
REAL_TYPE* RIAMC::Alloc_Real_V3D(const int* sz)
{
  if ( !sz ) return NULL;

  size_t dims[3], nx;

  dims[0] = (size_t)(sz[0] + 2*GUIDE);
  dims[1] = (size_t)(sz[1] + 2*GUIDE);
  dims[2] = (size_t)(sz[2] + 2*GUIDE);

  nx = dims[0] * dims[1] * dims[2] * 3;

  REAL_TYPE* var = new REAL_TYPE[nx];

#pragma omp parallel for
  for (int i=0; i<nx; i++) var[i]=0.0;

  return var;
}



// #################################################################
/**
 * @brief S3D配列のアロケート
 * @param [in] sz 配列サイズ
 * @ret pointer
 */
REAL_TYPE* RIAMC::Alloc_Real_S3D(const int* sz)
{
  if ( !sz ) return NULL;

  size_t dims[3], nx;

  dims[0] = (size_t)(sz[0] + 2*GUIDE);
  dims[1] = (size_t)(sz[1] + 2*GUIDE);
  dims[2] = (size_t)(sz[2] + 2*GUIDE);

  nx = dims[0] * dims[1] * dims[2];

  REAL_TYPE* var = new REAL_TYPE[nx];

#pragma omp parallel for
  for (int i=0; i<nx; i++) var[i]=0.0;

  return var;
}



// #################################################################
/**
 * @brief 配列のアロケート
 * @param [in,out] total ソルバーに使用するメモリ量
 */
bool RIAMC::allocateArray(double &total)
{
  double array_size = (size[0]+2*GUIDE) * (size[1]+2*GUIDE) * (size[2]+2*GUIDE);

  total += ( array_size * 19 ) * (double)sizeof(REAL_TYPE);

  // Poisson の右辺項
  if( (RHS = Alloc_Real_S3D(size)) == NULL ) return false;

  // 圧力
  if( (P   = Alloc_Real_S3D(size)) == NULL ) return false;

  // SGS
  if( (SGS = Alloc_Real_S3D(size)) == NULL ) return false;


  // 座標値
  if( (Z   = Alloc_Real_S3D(size)) == NULL ) return false;


  // V3D配列
  if( (V   = Alloc_Real_V3D(size)) == NULL ) return false;
  if( (VV  = Alloc_Real_V3D(size)) == NULL ) return false;
  if( (VD  = Alloc_Real_V3D(size)) == NULL ) return false;
  if( (CV  = Alloc_Real_V3D(size)) == NULL ) return false;

  // ワーク配列 V3DEx配列
  if( (uvw = Alloc_Real_V3D(size)) == NULL ) return false;


  return true;
}
