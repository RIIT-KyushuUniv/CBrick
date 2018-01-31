/*
###################################################################################
#
# CBrick
#
# Copyright (c) 2017-2018 Research Institute for Information Technology(RIIT),
#                    Kyushu University.  All rights reserved.
#
####################################################################################
*/

/*
 * @file   CB_PackingVectorNode.cpp
 * @brief  SubDomain class
 */

#include "CB_SubDomain.h"
#include "CB_Pack.h"


/*
 * @brief pack send data for I direction
 * @param [in]  array   source array
 * @param [in]  vc_comm number of guide cell layer to be sent
 * @param [out] sendm   send buffer of I- direction
 * @param [out] sendp   send buffer of I+ direction
 * @param [in]  nIDm    Rank number of I- direction
 * @param [in]  nIDp    Rank number of I+ direction
 */
void SubDomain::pack_VIn(const REAL_TYPE *array,
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
          for( int i=1; i<=vc_comm; i++ ){
            sendm[_IDX_VI(i,j,k,l,1,jmax,kmax,vc_comm)] = array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)];
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
            sendp[_IDX_VI(i,j,k,l,imax-vc_comm,jmax,kmax,vc_comm)] = array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)];
          }
        }
      }
    }
  }
}


/*
 * @brief unpack send data for I direction
 * @param [in,out]  array   dest array
 * @param [in]  vc_comm number of guide cell layer to be sent
 * @param [in]  recvm   recv buffer of I- direction
 * @param [in]  recvp   recv buffer of I+ direction
 * @param [in]  nIDm    Rank number of I- direction
 * @param [in]  nIDp    Rank number of I+ direction
 */
void SubDomain::unpack_VIn(REAL_TYPE *array,
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
          for( int i=1-vc_comm; i<1; i++ ){
            array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)] = recvm[_IDX_VI(i,j,k,l,1-vc_comm,jmax,kmax,vc_comm)];
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
            array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)] = recvp[_IDX_VI(i,j,k,l,imax,jmax,kmax,vc_comm)];
          }
        }
      }
    }
  }
}

/*
 * @brief pack send data for J direction
 * @param [in]  array   source array
 * @param [in]  vc_comm number of guide cell layer to be sent
 * @param [out] sendm   send buffer of J- direction
 * @param [out] sendp   send buffer of J+ direction
 * @param [in]  nIDm    Rank number of J- direction
 * @param [in]  nIDp    Rank number of J+ direction
 */
void SubDomain::pack_VJn(const REAL_TYPE *array,
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
        for( int j=1; j<=vc_comm; j++ ){
          for( int i=0-vc_comm; i<imax+vc_comm; i++ ){
            sendm[_IDX_VJ(i,j,k,l,imax,1,kmax,vc_comm)] = array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)];
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
            sendp[_IDX_VJ(i,j,k,l,imax,jmax-vc_comm,kmax,vc_comm)] = array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)];
          }
        }
      }
    }
  }
}


/*
 * @brief unpack send data for J direction
 * @param [in,out]  array   dest array
 * @param [in]  vc_comm number of guide cell layer to be sent
 * @param [in]  recvm   recv buffer of J- direction
 * @param [in]  recvp   recv buffer of J+ direction
 * @param [in]  nIDm    Rank number of J- direction
 * @param [in]  nIDp    Rank number of J+ direction
 */
void SubDomain::unpack_VJn(REAL_TYPE *array,
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
        for( int j=1-vc_comm; j<1; j++ ){
          for( int i=0-vc_comm; i<imax+vc_comm; i++ ){
            array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)] = recvm[_IDX_VJ(i,j,k,l,imax,1-vc_comm,kmax,vc_comm)];
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
            array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)] = recvp[_IDX_VJ(i,j,k,l,imax,jmax,kmax,vc_comm)];
          }
        }
      }
    }
  }
}


/*
 * @brief pack send data for K direction
 * @param [in]  array   source array
 * @param [in]  vc_comm number of guide cell layer actually to be sent
 * @param [out] sendm   send buffer of K- direction
 * @param [out] sendp   send buffer of K+ direction
 * @param [in]  nIDm    Rank number of K- direction
 * @param [in]  nIDp    Rank number of K+ direction
 */
void SubDomain::pack_VKn(const REAL_TYPE *array,
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
      for( int k=1; k<=vc_comm; k++ ){
        for( int j=0-vc_comm; j<jmax+vc_comm; j++ ){
          for( int i=0-vc_comm; i<imax+vc_comm; i++ ){
            sendm[_IDX_VK(i,j,k,l,imax,jmax,1,vc_comm)] = array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)];
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
            sendp[_IDX_VK(i,j,k,l,imax,jmax,kmax-vc_comm,vc_comm)] = array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)];
          }
        }
      }
    }
  }
}


/*
 * @brief unpack send data for K direction
 * @param [in,out]  array   dest array
 * @param [in]  vc_comm number of guide cell layer to be sent
 * @param [in]  recvm   recv buffer of K- direction
 * @param [in]  recvp   recv buffer of K+ direction
 * @param [in]  nIDm    Rank number of K- direction
 * @param [in]  nIDp    Rank number of K+ direction
 */
void SubDomain::unpack_VKn(REAL_TYPE *array,
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
      for( int k=1-vc_comm; k<1; k++ ){
        for( int j=0-vc_comm; j<jmax+vc_comm; j++ ){
          for( int i=0-vc_comm; i<imax+vc_comm; i++ ){
            array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)] = recvm[_IDX_VK(i,j,k,l,imax,jmax,1-vc_comm,vc_comm)];
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
            array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)] = recvp[_IDX_VK(i,j,k,l,imax,jmax,kmax,vc_comm)];
          }
        }
      }
    }
  }
}
