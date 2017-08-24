/*
###################################################################################
#
# CBrick
#
# Copyright (c) 2017 Research Institute for Information Technology(RIIT),
#                    Kyushu University.  All rights reserved.
#
####################################################################################
*/

/*
 * @file   CB_PackingVector.cpp
 * @brief  SubDomain class
 */

#include "CB_SubDomain.h"


/* バッファへのインデクス変換 (I方向)
 *  @param [in] _I  i方向インデクス
 *  @param [in] _J  j方向インデクス
 *  @param [in] _K  k方向インデクス
 *  @param [in] _L  ベクトル成分インデクス {0,1,2}
 *  @param [in] _IS i方向の開始点インデクス
 *  @param [in] _NJ j方向インデクスサイズ
 *  @param [in] _NK k方向インデクスサイズ
 *  @param [in] _VC 実際に送受信する仮想セル数
 *  @return 1次元インデクス
 */
#define _IDX_VX(_I,_J,_K,_L,_IS,_NJ,_NK,_VC) \
( _L * _VC * (_NJ+2*_VC) * (_NK+2*_VC) \
+ (_K+_VC) * _VC * (_NJ+2*_VC) \
+ (_J+_VC) * _VC \
+ (_I-_IS) \
)

/* バッファへのインデクス変換 (J方向)
 *  @param [in] _I  i方向インデクス
 *  @param [in] _J  j方向インデクス
 *  @param [in] _K  k方向インデクス
 *  @param [in] _L  ベクトル成分インデクス {0,1,2}
 *  @param [in] _NI i方向インデクスサイズ
 *  @param [in] _JS j方向の開始点インデクス
 *  @param [in] _NK k方向インデクスサイズ
 *  @param [in] _VC 実際に送受信する仮想セル数
 *  @return 1次元インデクス
 */
#define _IDX_VY(_I,_J,_K,_L,_NI,_JS,_NK,_VC) \
( _L * (_NI+2*_VC) * _VC * (_NK+2*_VC) \
+ (_K+_VC) * (_NI+2*_VC) * _VC \
+ (_J-_JS) * (_NI+2*_VC) \
+ (_I+_VC) \
)

/* バッファへのインデクス変換 (K方向)
 *  @param [in] _I  i方向インデクス
 *  @param [in] _J  j方向インデクス
 *  @param [in] _K  k方向インデクス
 *  @param [in] _L  ベクトル成分インデクス {0,1,2}
 *  @param [in] _NI i方向インデクスサイズ
 *  @param [in] _NJ j方向インデクスサイズ
 *  @param [in] _KS ｋ方向の開始点インデクス
 *  @param [in] _VC 実際に送受信する仮想セル数
 *  @return 1次元インデクス
 */
#define _IDX_VZ(_I,_J,_K,_L,_NI,_NJ,_KS,_VC) \
( _L * (_NI+2*_VC) * (_NJ+2*_VC) * _VC \
+ (_K-_KS) * (_NI+2*_VC) * (_NJ+2*_VC) \
+ (_J+_VC) * (_NI+2*_VC) \
+ (_I+_VC) \
)


/*
 * @brief pack send data for X direction
 * @param [in]  array   source array
 * @param [in]  vc_comm number of guide cell layer to be sent
 * @param [out] sendm   send buffer of X- direction
 * @param [out] sendp   send buffer of X+ direction
 * @param [in]  nIDm    Rank number of X- direction
 * @param [in]  nIDp    Rank number of X+ direction
 */
void SubDomain::pack_VX(const REAL_TYPE *array,
                        const int vc_comm,
                        REAL_TYPE *sendm,
                        REAL_TYPE *sendp,
                        const int nIDm,
                        const int nIDp)
{
  int imax = size[0];
  int jmax = size[1];
  int kmax = size[2];
  int vc = halo_width;

  if( nIDm >= 0 )
  {
#pragma omp parallel for collapse(3)
    for (int l=0; l<3; l++) {
      for( int k=0-vc_comm; k<kmax+vc_comm; k++ ){
        for( int j=0-vc_comm; j<jmax+vc_comm; j++ ){
          for( int i=0; i<vc_comm; i++ ){
            sendm[_IDX_VX(i,j,k,l,0,jmax,kmax,vc_comm)] = array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)];
          }
        }
      }
    }
  }

  if( nIDp >= 0 )
  {
#pragma omp parallel for collapse(3)
    for (int l=0; l<3; l++) {
      for( int k=0-vc_comm; k<kmax+vc_comm; k++ ){
        for( int j=0-vc_comm; j<jmax+vc_comm; j++ ){
          for( int i=imax-vc_comm; i<imax; i++ ){
            sendp[_IDX_VX(i,j,k,l,imax-vc_comm,jmax,kmax,vc_comm)] = array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)];
          }
        }
      }
    }
  }
}


/*
 * @brief unpack send data for X direction
 * @param [in,out]  array   dest array
 * @param [in]  vc_comm number of guide cell layer to be sent
 * @param [in]  recvm   recv buffer of X- direction
 * @param [in]  recvp   recv buffer of X+ direction
 * @param [in]  nIDm    Rank number of X- direction
 * @param [in]  nIDp    Rank number of X+ direction
 */
void SubDomain::unpack_VX(REAL_TYPE *array,
                          const int vc_comm,
                          const REAL_TYPE *recvm,
                          const REAL_TYPE *recvp,
                          const int nIDm,
                          const int nIDp)
{
  int imax = size[0];
  int jmax = size[1];
  int kmax = size[2];
  int vc = halo_width;

  if( nIDm >= 0 )
  {
#pragma omp parallel for collapse(3)
    for (int l=0; l<3; l++) {
      for( int k=0-vc_comm; k<kmax+vc_comm; k++ ){
        for( int j=0-vc_comm; j<jmax+vc_comm; j++ ){
          for( int i=0-vc_comm; i<0; i++ ){
            array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)] = recvm[_IDX_VX(i,j,k,l,0-vc_comm,jmax,kmax,vc_comm)];
          }
        }
      }
    }
  }

  if( nIDp >= 0 )
  {
#pragma omp parallel for collapse(3)
    for (int l=0; l<3; l++) {
      for( int k=0-vc_comm; k<kmax+vc_comm; k++ ){
        for( int j=0-vc_comm; j<jmax+vc_comm; j++ ){
          for( int i=imax; i<imax+vc_comm; i++ ){
            array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)] = recvp[_IDX_VX(i,j,k,l,imax,jmax,kmax,vc_comm)];
          }
        }
      }
    }
  }
}

/*
 * @brief pack send data for Y direction
 * @param [in]  array   source array
 * @param [in]  vc_comm number of guide cell layer to be sent
 * @param [out] sendm   send buffer of Y- direction
 * @param [out] sendp   send buffer of Y+ direction
 * @param [in]  nIDm    Rank number of Y- direction
 * @param [in]  nIDp    Rank number of Y+ direction
 */
void SubDomain::pack_VY(const REAL_TYPE *array,
                        const int vc_comm,
                        REAL_TYPE *sendm,
                        REAL_TYPE *sendp,
                        const int nIDm,
                        const int nIDp)
{
  int imax = size[0];
  int jmax = size[1];
  int kmax = size[2];
  int vc = halo_width;

  if( nIDm >= 0 )
  {
#pragma omp parallel for collapse(3)
    for (int l=0; l<3; l++) {
      for( int k=0-vc_comm; k<kmax+vc_comm; k++ ){
        for( int j=0; j<vc_comm; j++ ){
          for( int i=0-vc_comm; i<imax+vc_comm; i++ ){
            sendm[_IDX_VY(i,j,k,l,imax,0,kmax,vc_comm)] = array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)];
          }
        }
      }
    }
  }

  if( nIDp >= 0 )
  {
#pragma omp parallel for collapse(3)
    for (int l=0; l<3; l++) {
      for( int k=0-vc_comm; k<kmax+vc_comm; k++ ){
        for( int j=jmax-vc_comm; j<jmax; j++ ){
          for( int i=0-vc_comm; i<imax+vc_comm; i++ ){
            sendp[_IDX_VY(i,j,k,l,imax,jmax-vc_comm,kmax,vc_comm)] = array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)];
          }
        }
      }
    }
  }
}


/*
 * @brief unpack send data for Y direction
 * @param [in,out]  array   dest array
 * @param [in]  vc_comm number of guide cell layer to be sent
 * @param [in]  recvm   recv buffer of Y- direction
 * @param [in]  recvp   recv buffer of Y+ direction
 * @param [in]  nIDm    Rank number of Y- direction
 * @param [in]  nIDp    Rank number of Y+ direction
 */
void SubDomain::unpack_VY(REAL_TYPE *array,
                          const int vc_comm,
                          const REAL_TYPE *recvm,
                          const REAL_TYPE *recvp,
                          const int nIDm,
                          const int nIDp)
{
  int imax = size[0];
  int jmax = size[1];
  int kmax = size[2];
  int vc = halo_width;

  if( nIDm >= 0 )
  {
#pragma omp parallel for collapse(3)
    for (int l=0; l<3; l++) {
      for( int k=0-vc_comm; k<kmax+vc_comm; k++ ){
        for( int j=0-vc_comm; j<0; j++ ){
          for( int i=0-vc_comm; i<imax+vc_comm; i++ ){
            array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)] = recvm[_IDX_VY(i,j,k,l,imax,0-vc_comm,kmax,vc_comm)];
          }
        }
      }
    }
  }

  if( nIDp >= 0 )
  {
#pragma omp parallel for collapse(3)
    for (int l=0; l<3; l++) {
      for( int k=0-vc_comm; k<kmax+vc_comm; k++ ){
        for( int j=jmax; j<jmax+vc_comm; j++ ){
          for( int i=0-vc_comm; i<imax+vc_comm; i++ ){
            array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)] = recvp[_IDX_VY(i,j,k,l,imax,jmax,kmax,vc_comm)];
          }
        }
      }
    }
  }
}


/*
 * @brief pack send data for Z direction
 * @param [in]  array   source array
 * @param [in]  vc_comm number of guide cell layer actually to be sent
 * @param [out] sendm   send buffer of Z- direction
 * @param [out] sendp   send buffer of Z+ direction
 * @param [in]  nIDm    Rank number of Z- direction
 * @param [in]  nIDp    Rank number of Z+ direction
 */
void SubDomain::pack_VZ(const REAL_TYPE *array,
                        const int vc_comm,
                        REAL_TYPE *sendm,
                        REAL_TYPE *sendp,
                        const int nIDm,
                        const int nIDp)
{
  int imax = size[0];
  int jmax = size[1];
  int kmax = size[2];
  int vc = halo_width;

  if( nIDm >= 0 )
  {
#pragma omp parallel for collapse(3)
    for (int l=0; l<3; l++) {
      for( int k=0; k<vc_comm; k++ ){
        for( int j=0-vc_comm; j<jmax+vc_comm; j++ ){
          for( int i=0-vc_comm; i<imax+vc_comm; i++ ){
            sendm[_IDX_VZ(i,j,k,l,imax,jmax,0,vc_comm)] = array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)];
          }
        }
      }
    }
  }

  if( nIDp >= 0 )
  {
#pragma omp parallel for collapse(3)
    for (int l=0; l<3; l++) {
      for( int k=kmax-vc_comm; k<kmax; k++ ){
        for( int j=0-vc_comm; j<jmax+vc_comm; j++ ){
          for( int i=0-vc_comm; i<imax+vc_comm; i++ ){
            sendp[_IDX_VZ(i,j,k,l,imax,jmax,kmax-vc_comm,vc_comm)] = array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)];
          }
        }
      }
    }
  }
}


/*
 * @brief unpack send data for Z direction
 * @param [in,out]  array   dest array
 * @param [in]  vc_comm number of guide cell layer to be sent
 * @param [in]  recvm   recv buffer of Z- direction
 * @param [in]  recvp   recv buffer of Z+ direction
 * @param [in]  nIDm    Rank number of Z- direction
 * @param [in]  nIDp    Rank number of Z+ direction
 */
void SubDomain::unpack_VZ(REAL_TYPE *array,
                          const int vc_comm,
                          const REAL_TYPE *recvm,
                          const REAL_TYPE *recvp,
                          const int nIDm,
                          const int nIDp)
{
  int imax = size[0];
  int jmax = size[1];
  int kmax = size[2];
  int vc = halo_width;

  if( nIDm >= 0 )
  {
#pragma omp parallel for collapse(3)
    for (int l=0; l<3; l++) {
      for( int k=0-vc_comm; k<0; k++ ){
        for( int j=0-vc_comm; j<jmax+vc_comm; j++ ){
          for( int i=0-vc_comm; i<imax+vc_comm; i++ ){
            array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)] = recvm[_IDX_VZ(i,j,k,l,imax,jmax,0-vc_comm,vc_comm)];
          }
        }
      }
    }
  }

  if( nIDp >= 0 )
  {
#pragma omp parallel for collapse(3)
    for (int l=0; l<3; l++) {
      for( int k=kmax; k<kmax+vc_comm; k++ ){
        for( int j=0-vc_comm; j<jmax+vc_comm; j++ ){
          for( int i=0-vc_comm; i<imax+vc_comm; i++ ){
            array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)] = recvp[_IDX_VZ(i,j,k,l,imax,jmax,kmax,vc_comm)];
          }
        }
      }
    }
  }
}
