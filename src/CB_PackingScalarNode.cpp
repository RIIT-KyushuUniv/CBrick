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
 * @file   CB_PackingScalarNode.cpp
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
void SubDomain::pack_SXnode(const REAL_TYPE *array,
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
#pragma omp parallel for collapse(2)
    for( int k=0; k<kmax; k++ ){
      for( int j=0; j<jmax; j++ ){
        for( int i=1; i<=vc_comm; i++ ){
          sendm[_IDX_SI(i,j,k,1,jmax,vc_comm)] = array[_IDX_S3D(i,j,k,imax,jmax,vc)];
        }
      }
    }
  }

  if( nIDp >= 0 )
  {
#pragma omp parallel for collapse(2)
    for( int k=0; k<kmax; k++ ){
      for( int j=0; j<jmax; j++ ){
        for( int i=imax-vc_comm; i<imax; i++ ){
          sendp[_IDX_SI(i,j,k,imax-vc_comm,jmax,vc_comm)] = array[_IDX_S3D(i,j,k,imax,jmax,vc)];
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
void SubDomain::unpack_SXnode(REAL_TYPE *array,
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
#pragma omp parallel for collapse(2)
    for( int k=0; k<kmax; k++ ){
      for( int j=0; j<jmax; j++ ){
        for( int i=1-vc_comm; i<1; i++ ){
          array[_IDX_S3D(i,j,k,imax,jmax,vc)] = recvm[_IDX_SI(i,j,k,1-vc_comm,jmax,vc_comm)];
        }
      }
    }
  }

  if( nIDp >= 0 )
  {
#pragma omp parallel for collapse(2)
    for( int k=0; k<kmax; k++ ){
      for( int j=0; j<jmax; j++ ){
        for( int i=imax; i<imax+vc_comm; i++ ){
          array[_IDX_S3D(i,j,k,imax,jmax,vc)] = recvp[_IDX_SI(i,j,k,imax,jmax,vc_comm)];
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
void SubDomain::pack_SYnode(const REAL_TYPE *array,
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
#pragma omp parallel for collapse(2)
    for( int k=0; k<kmax; k++ ){
      for( int j=1; j<=vc_comm; j++ ){
        for( int i=0; i<imax; i++ ){
          sendm[_IDX_SJ(i,j,k,imax,1,vc_comm)] = array[_IDX_S3D(i,j,k,imax,jmax,vc)];
        }
      }
    }
  }

  if( nIDp >= 0 )
  {
#pragma omp parallel for collapse(2)
    for( int k=0; k<kmax; k++ ){
      for( int j=jmax-vc_comm; j<jmax; j++ ){
        for( int i=0; i<imax; i++ ){
          sendp[_IDX_SJ(i,j,k,imax,jmax-vc_comm,vc_comm)] = array[_IDX_S3D(i,j,k,imax,jmax,vc)];
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
void SubDomain::unpack_SYnode(REAL_TYPE *array,
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
#pragma omp parallel for collapse(2)
    for( int k=0; k<kmax; k++ ){
      for( int j=1-vc_comm; j<1; j++ ){
        for( int i=0; i<imax; i++ ){
          array[_IDX_S3D(i,j,k,imax,jmax,vc)] = recvm[_IDX_SJ(i,j,k,imax,1-vc_comm,vc_comm)];
        }
      }
    }
  }

  if( nIDp >= 0 )
  {
#pragma omp parallel for collapse(2)
    for( int k=0; k<kmax; k++ ){
      for( int j=jmax; j<jmax+vc_comm; j++ ){
        for( int i=0; i<imax; i++ ){
          array[_IDX_S3D(i,j,k,imax,jmax,vc)] = recvp[_IDX_SJ(i,j,k,imax,jmax,vc_comm)];
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
void SubDomain::pack_SZnode(const REAL_TYPE *array,
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
#pragma omp parallel for collapse(2)
    for( int k=1; k<=vc_comm; k++ ){
      for( int j=0; j<jmax; j++ ){
        for( int i=0; i<imax; i++ ){
          sendm[_IDX_SK(i,j,k,imax,jmax,1,vc_comm)] = array[_IDX_S3D(i,j,k,imax,jmax,vc)];
        }
      }
    }
  }

  if( nIDp >= 0 )
  {
#pragma omp parallel for collapse(2)
    for( int k=kmax-vc_comm; k<kmax; k++ ){
      for( int j=0; j<jmax; j++ ){
        for( int i=0; i<imax; i++ ){
          sendp[_IDX_SK(i,j,k,imax,jmax,kmax-vc_comm,vc_comm)] = array[_IDX_S3D(i,j,k,imax,jmax,vc)];
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
void SubDomain::unpack_SZnode(REAL_TYPE *array,
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
#pragma omp parallel for collapse(2)
    for( int k=1-vc_comm; k<1; k++ ){
      for( int j=0; j<jmax; j++ ){
        for( int i=0; i<imax; i++ ){
          array[_IDX_S3D(i,j,k,imax,jmax,vc)] = recvm[_IDX_SK(i,j,k,imax,jmax,1-vc_comm,vc_comm)];
        }
      }
    }
  }

  if( nIDp >= 0 )
  {
#pragma omp parallel for collapse(2)
    for( int k=kmax; k<kmax+vc_comm; k++ ){
      for( int j=0; j<jmax; j++ ){
        for( int i=0; i<imax; i++ ){
          array[_IDX_S3D(i,j,k,imax,jmax,vc)] = recvp[_IDX_SK(i,j,k,imax,jmax,kmax,vc_comm)];
        }
      }
    }
  }
}

/*
 * @brief pack send data for diagonal edge
 * @param [in]  array    source array
 * @param [in]  vc_comm  number of guide cell layer to be sent
 * @param [out] sendbuf  send buffer
 * @param [out] recvbuf  recv buffer
 * @param [out] req      Array of MPI request
 * @retval true-success, false-fail
 */
bool SubDomain::pack_SEnode(REAL_TYPE *array,
                         const int vc_comm,
                         REAL_TYPE *sendbuf,
                         REAL_TYPE *recvbuf,
                         MPI_Request *req)
{
#ifdef _DIAGONAL_COMM

  int imax = size[0];
  int jmax = size[1];
  int kmax = size[2];
  int vc = halo_width;
  int tag = 0;
  MPI_Datatype dtype = MPI_FLOAT;
  if( sizeof(REAL_TYPE) == _SIZE_DOUBLE_ )
  {
    dtype = MPI_DOUBLE;
  }

  size_t ptr = 0;

  //// X edge ////
  for( int dir=int(E_mYmZ);dir<=int(E_pYpZ);dir++ )
  {
    if( comm_tbl[dir] >= 0 )
    {
      REAL_TYPE *sendptr = &sendbuf[ptr];
      REAL_TYPE *recvptr = &recvbuf[ptr];
      size_t sz = (imax-1) * vc_comm * vc_comm;

      // recv
      if ( MPI_SUCCESS != MPI_Irecv(recvptr,
                                    sz,
                                    dtype,
                                    comm_tbl[dir],
                                    tag,
                                    MPI_COMM_WORLD,
                                    &req[dir*2]) ) return false;

      // pack
      switch(dir)
      {
      case int(E_mYmZ):
#pragma omp parallel for collapse(3)
        for( int k=1; k<=vc_comm; k++ ){
          for( int j=1; j<=vc_comm; j++ ){
            for( int i=1; i<imax; i++ ){
              sendptr[_IDX_S3D(i-1,j-1,k-1,imax-1,vc_comm,0)] = array[_IDX_S3D(i,j,k,imax,jmax,vc)];
            }
          }
        }
        break;

      case int(E_pYmZ):
#pragma omp parallel for collapse(3)
        for( int k=1; k<=vc_comm; k++ ){
          for( int j=jmax-vc_comm; j<jmax; j++ ){
            for( int i=1; i<imax; i++ ){
              sendptr[_IDX_S3D(i-1,j-(jmax-vc_comm),k-1,imax-1,vc_comm,0)] = array[_IDX_S3D(i,j,k,imax,jmax,vc)];
            }
          }
        }
        break;

      case int(E_mYpZ):
#pragma omp parallel for collapse(3)
        for( int k=kmax-vc_comm; k<kmax; k++ ){
          for( int j=1; j<=vc_comm; j++ ){
            for( int i=1; i<imax; i++ ){
              sendptr[_IDX_S3D(i-1,j-1,k-(kmax-vc_comm),imax-1,vc_comm,0)] = array[_IDX_S3D(i,j,k,imax,jmax,vc)];
            }
          }
        }
        break;
      case int(E_pYpZ):

#pragma omp parallel for collapse(3)
        for( int k=kmax-vc_comm; k<kmax; k++ ){
          for( int j=jmax-vc_comm; j<jmax; j++ ){
            for( int i=1; i<imax; i++ ){
              sendptr[_IDX_S3D(i-1,j-(jmax-vc_comm),k-(kmax-vc_comm),imax-1,vc_comm,0)] = array[_IDX_S3D(i,j,k,imax,jmax,vc)];
            }
          }
        }
        break;
      }

      // send
      if ( MPI_SUCCESS != MPI_Isend(sendptr,
                                    sz,
                                    dtype,
                                    comm_tbl[dir],
                                    tag,
                                    MPI_COMM_WORLD,
                                    &req[dir*2+1]) ) return false;

      // pointer
      ptr += sz;
    }
  }

  //// Y edge ////
  for( int dir=int(E_mXmZ);dir<=int(E_pXpZ);dir++ )
  {
    if( comm_tbl[dir] >= 0 )
    {
      REAL_TYPE *sendptr = &sendbuf[ptr];
      REAL_TYPE *recvptr = &recvbuf[ptr];
      size_t sz = vc_comm * (jmax-1) * vc_comm;

      // recv
      if ( MPI_SUCCESS != MPI_Irecv(recvptr,
                                    sz,
                                    dtype,
                                    comm_tbl[dir],
                                    tag,
                                    MPI_COMM_WORLD,
                                    &req[dir*2]) ) return false;

      // pack
      switch(dir)
      {
      case int(E_mXmZ):
#pragma omp parallel for collapse(3)
        for( int k=1; k<=vc_comm; k++ ){
          for( int j=1; j<jmax; j++ ){
            for( int i=1; i<=vc_comm; i++ ){
              sendptr[_IDX_S3D(i-1,j-1,k-1,vc_comm,jmax-1,0)] = array[_IDX_S3D(i,j,k,imax,jmax,vc)];
            }
          }
        }
        break;

      case int(E_pXmZ):
#pragma omp parallel for collapse(3)
        for( int k=1; k<=vc_comm; k++ ){
          for( int j=1; j<jmax; j++ ){
            for( int i=imax-vc_comm; i<imax; i++ ){
              sendptr[_IDX_S3D(i-(imax-vc_comm),j-1,k-1,vc_comm,jmax-1,0)] = array[_IDX_S3D(i,j,k,imax,jmax,vc)];
            }
          }
        }
        break;

      case int(E_mXpZ):
#pragma omp parallel for collapse(3)
        for( int k=kmax-vc_comm; k<kmax; k++ ){
          for( int j=1; j<jmax; j++ ){
            for( int i=1; i<=vc_comm; i++ ){
              sendptr[_IDX_S3D(i-1,j-1,k-(kmax-vc_comm),vc_comm,jmax-1,0)] = array[_IDX_S3D(i,j,k,imax,jmax,vc)];
            }
          }
        }
        break;
      case int(E_pXpZ):

#pragma omp parallel for collapse(3)
        for( int k=kmax-vc_comm; k<kmax; k++ ){
          for( int j=1; j<jmax; j++ ){
            for( int i=imax-vc_comm; i<imax; i++ ){
              sendptr[_IDX_S3D(i-(imax-vc_comm),j-1,k-(kmax-vc_comm),vc_comm,jmax-1,0)] = array[_IDX_S3D(i,j,k,imax,jmax,vc)];
            }
          }
        }
        break;
      }

      // send
      if ( MPI_SUCCESS != MPI_Isend(sendptr,
                                    sz,
                                    dtype,
                                    comm_tbl[dir],
                                    tag,
                                    MPI_COMM_WORLD,
                                    &req[dir*2+1]) ) return false;

      // pointer
      ptr += sz;
    }
  }

  //// Z edge ////
  for( int dir=int(E_mXmY);dir<=int(E_pXpY);dir++ )
  {
    if( comm_tbl[dir] >= 0 )
    {
      REAL_TYPE *sendptr = &sendbuf[ptr];
      REAL_TYPE *recvptr = &recvbuf[ptr];
      size_t sz = vc_comm * vc_comm * (kmax-1);

      // recv
      if ( MPI_SUCCESS != MPI_Irecv(recvptr,
                                    sz,
                                    dtype,
                                    comm_tbl[dir],
                                    tag,
                                    MPI_COMM_WORLD,
                                    &req[dir*2]) ) return false;

      // pack
      switch(dir)
      {
      case int(E_mXmY):
#pragma omp parallel for collapse(3)
        for( int k=1; k<kmax; k++ ){
          for( int j=1; j<=vc_comm; j++ ){
            for( int i=1; i<=vc_comm; i++ ){
              sendptr[_IDX_S3D(i-1,j-1,k-1,vc_comm,vc_comm,0)] = array[_IDX_S3D(i,j,k,imax,jmax,vc)];
            }
          }
        }
        break;

      case int(E_pXmY):
#pragma omp parallel for collapse(3)
        for( int k=1; k<kmax; k++ ){
          for( int j=1; j<=vc_comm; j++ ){
            for( int i=imax-vc_comm; i<imax; i++ ){
              sendptr[_IDX_S3D(i-(imax-vc_comm),j-1,k-1,vc_comm,vc_comm,0)] = array[_IDX_S3D(i,j,k,imax,jmax,vc)];
            }
          }
        }
        break;

      case int(E_mXpY):
#pragma omp parallel for collapse(3)
        for( int k=1; k<kmax; k++ ){
          for( int j=jmax-vc_comm; j<jmax; j++ ){
            for( int i=1; i<=vc_comm; i++ ){
              sendptr[_IDX_S3D(i-1,j-(jmax-vc_comm),k-1,vc_comm,vc_comm,0)] = array[_IDX_S3D(i,j,k,imax,jmax,vc)];
            }
          }
        }
        break;
      case int(E_pXpY):

#pragma omp parallel for collapse(3)
        for( int k=1; k<kmax; k++ ){
          for( int j=jmax-vc_comm; j<jmax; j++ ){
            for( int i=imax-vc_comm; i<imax; i++ ){
              sendptr[_IDX_S3D(i-(imax-vc_comm),j-(jmax-vc_comm),k-1,vc_comm,vc_comm,0)] = array[_IDX_S3D(i,j,k,imax,jmax,vc)];
            }
          }
        }
        break;
      }

      // send
      if ( MPI_SUCCESS != MPI_Isend(sendptr,
                                    sz,
                                    dtype,
                                    comm_tbl[dir],
                                    tag,
                                    MPI_COMM_WORLD,
                                    &req[dir*2+1]) ) return false;

      // pointer
      ptr += sz;
    }
  }

#endif
  return true;
}

/*
 * @brief unpack send data for diagonal edge
 * @param [out] array    dest array
 * @param [in]  vc_comm  number of guide cell layer to be sent
 * @param [in]  recvbuf  recv buffer
 */
void SubDomain::unpack_SEnode(REAL_TYPE *array,
                           const int vc_comm,
                           const REAL_TYPE *recvbuf)
{
#ifdef _DIAGONAL_COMM

  int imax = size[0];
  int jmax = size[1];
  int kmax = size[2];
  int vc = halo_width;

  size_t ptr = 0;

  //// X edge ////
  for( int dir=int(E_mYmZ);dir<=int(E_pYpZ);dir++ )
  {
    if( comm_tbl[dir] >= 0 )
    {
      const REAL_TYPE *recvptr = &recvbuf[ptr];
      size_t sz = (imax-1) * vc_comm * vc_comm;

      // unpack
      switch(dir)
      {
      case int(E_mYmZ):
#pragma omp parallel for collapse(3)
        for( int k=1-vc_comm; k<=0; k++ ){
          for( int j=1-vc_comm; j<=0; j++ ){
            for( int i=1; i<imax; i++ ){
              array[_IDX_S3D(i,j,k,imax,jmax,vc)] = recvptr[_IDX_S3D(i-1,j-(1-vc_comm),k-(1-vc_comm),imax-1,vc_comm,0)];
            }
          }
        }
        break;

      case int(E_pYmZ):
#pragma omp parallel for collapse(3)
        for( int k=1-vc_comm; k<=0; k++ ){
          for( int j=jmax; j<jmax+vc_comm; j++ ){
            for( int i=1; i<imax; i++ ){
              array[_IDX_S3D(i,j,k,imax,jmax,vc)] = recvptr[_IDX_S3D(i-1,j-(jmax),k-(1-vc_comm),imax-1,vc_comm,0)];
            }
          }
        }
        break;

      case int(E_mYpZ):
#pragma omp parallel for collapse(3)
        for( int k=kmax; k<kmax+vc_comm; k++ ){
          for( int j=1-vc_comm; j<=0; j++ ){
            for( int i=1; i<imax; i++ ){
              array[_IDX_S3D(i,j,k,imax,jmax,vc)] = recvptr[_IDX_S3D(i-1,j-(1-vc_comm),k-(kmax),imax-1,vc_comm,0)];
            }
          }
        }
        break;

      case int(E_pYpZ):
#pragma omp parallel for collapse(3)
        for( int k=kmax; k<kmax+vc_comm; k++ ){
          for( int j=jmax; j<jmax+vc_comm; j++ ){
            for( int i=1; i<imax; i++ ){
              array[_IDX_S3D(i,j,k,imax,jmax,vc)] = recvptr[_IDX_S3D(i-1,j-(jmax),k-(kmax),imax-1,vc_comm,0)];
            }
          }
        }
        break;
      }

      ptr += sz;
    }
  }

  //// Y edge ////
  for( int dir=int(E_mXmZ);dir<=int(E_pXpZ);dir++ )
  {
    if( comm_tbl[dir] >= 0 )
    {
      const REAL_TYPE *recvptr = &recvbuf[ptr];
      size_t sz = vc_comm * (jmax-1) * vc_comm;

      // unpack
      switch(dir)
      {
      case int(E_mXmZ):
#pragma omp parallel for collapse(3)
        for( int k=1-vc_comm; k<=0; k++ ){
          for( int j=1; j<jmax; j++ ){
            for( int i=1-vc_comm; i<=0; i++ ){
              array[_IDX_S3D(i,j,k,imax,jmax,vc)] = recvptr[_IDX_S3D(i-(1-vc_comm),j-1,k-(1-vc_comm),vc_comm,jmax-1,0)];
            }
          }
        }
        break;

      case int(E_pXmZ):
#pragma omp parallel for collapse(3)
        for( int k=1-vc_comm; k<=0; k++ ){
          for( int j=1; j<jmax; j++ ){
            for( int i=imax; i<imax+vc_comm; i++ ){
              array[_IDX_S3D(i,j,k,imax,jmax,vc)] = recvptr[_IDX_S3D(i-(imax),j-1,k-(1-vc_comm),vc_comm,jmax-1,0)];
            }
          }
        }
        break;

      case int(E_mXpZ):
#pragma omp parallel for collapse(3)
        for( int k=kmax; k<kmax+vc_comm; k++ ){
          for( int j=1; j<jmax; j++ ){
            for( int i=1-vc_comm; i<=0; i++ ){
              array[_IDX_S3D(i,j,k,imax,jmax,vc)] = recvptr[_IDX_S3D(i-(1-vc_comm),j-1,k-(kmax),vc_comm,jmax-1,0)];
            }
          }
        }
        break;

      case int(E_pXpZ):
#pragma omp parallel for collapse(3)
        for( int k=kmax; k<kmax+vc_comm; k++ ){
          for( int j=1; j<jmax; j++ ){
            for( int i=imax; i<imax+vc_comm; i++ ){
              array[_IDX_S3D(i,j,k,imax,jmax,vc)] = recvptr[_IDX_S3D(i-(imax),j-1,k-(kmax),vc_comm,jmax-1,0)];
            }
          }
        }
        break;
      }

      ptr += sz;
    }
  }

  //// Z edge ////
  for( int dir=int(E_mXmY);dir<=int(E_pXpY);dir++ )
  {
    if( comm_tbl[dir] >= 0 )
    {
      const REAL_TYPE *recvptr = &recvbuf[ptr];
      size_t sz = vc_comm * vc_comm * (kmax-1);

      // unpack
      switch(dir)
      {
      case int(E_mXmY):
#pragma omp parallel for collapse(3)
        for( int k=1; k<kmax; k++ ){
          for( int j=1-vc_comm; j<=0; j++ ){
            for( int i=1-vc_comm; i<=0; i++ ){
              array[_IDX_S3D(i,j,k,imax,jmax,vc)] = recvptr[_IDX_S3D(i-(1-vc_comm),j-(1-vc_comm),k-1,vc_comm,vc_comm,0)];
            }
          }
        }
        break;

      case int(E_pXmY):
#pragma omp parallel for collapse(3)
        for( int k=1; k<kmax; k++ ){
          for( int j=1-vc_comm; j<=0; j++ ){
            for( int i=imax; i<imax+vc_comm; i++ ){
              array[_IDX_S3D(i,j,k,imax,jmax,vc)] = recvptr[_IDX_S3D(i-(imax),j-(1-vc_comm),k-1,vc_comm,vc_comm,0)];
            }
          }
        }
        break;

      case int(E_mXpY):
#pragma omp parallel for collapse(3)
        for( int k=1; k<kmax; k++ ){
          for( int j=jmax; j<jmax+vc_comm; j++ ){
            for( int i=1-vc_comm; i<=0; i++ ){
              array[_IDX_S3D(i,j,k,imax,jmax,vc)] = recvptr[_IDX_S3D(i-(1-vc_comm),j-(jmax),k-1,vc_comm,vc_comm,0)];
            }
          }
        }
        break;

      case int(E_pXpY):
#pragma omp parallel for collapse(3)
        for( int k=1; k<kmax; k++ ){
          for( int j=jmax; j<jmax+vc_comm; j++ ){
            for( int i=imax; i<imax+vc_comm; i++ ){
              array[_IDX_S3D(i,j,k,imax,jmax,vc)] = recvptr[_IDX_S3D(i-(imax),j-(jmax),k-1,vc_comm,vc_comm,0)];
            }
          }
        }
        break;
      }

      ptr += sz;
    }
  }

#endif
}

/*
 * @brief pack send data for diagonal corner
 * @param [in]  array    source array
 * @param [in]  vc_comm  number of guide cell layer to be sent
 * @param [out] sendbuf  send buffer
 * @param [out] recvbuf  recv buffer
 * @param [out] req      Array of MPI request
 * @retval true-success, false-fail
 */
bool SubDomain::pack_SCnode(REAL_TYPE *array,
                         const int vc_comm,
                         REAL_TYPE *sendbuf,
                         REAL_TYPE *recvbuf,
                         MPI_Request *req)
{
#ifdef _DIAGONAL_COMM

  int imax = size[0];
  int jmax = size[1];
  int kmax = size[2];
  int vc = halo_width;
  int tag = 0;
  MPI_Datatype dtype = MPI_FLOAT;
  if( sizeof(REAL_TYPE) == _SIZE_DOUBLE_ )
  {
    dtype = MPI_DOUBLE;
  }

  size_t ptr = 0;

  //// 8 corner ////
  for( int dir=int(C_mXmYmZ);dir<=int(C_pXpYpZ);dir++ )
  {
    if( comm_tbl[dir] >= 0 )
    {
      REAL_TYPE *sendptr = &sendbuf[ptr];
      REAL_TYPE *recvptr = &recvbuf[ptr];
      size_t sz = vc_comm * vc_comm * vc_comm;

      // recv
      if ( MPI_SUCCESS != MPI_Irecv(recvptr,
                                    sz,
                                    dtype,
                                    comm_tbl[dir],
                                    tag,
                                    MPI_COMM_WORLD,
                                    &req[dir*2]) ) return false;

      // pack
      switch(dir)
      {
      case int(C_mXmYmZ):
#pragma omp parallel for collapse(3)
        for( int k=1; k<=vc_comm; k++ ){
          for( int j=1; j<=vc_comm; j++ ){
            for( int i=1; i<=vc_comm; i++ ){
              sendptr[_IDX_S3D(i-1,j-1,k-1,vc_comm,vc_comm,0)] = array[_IDX_S3D(i,j,k,imax,jmax,vc)];
            }
          }
        }
        break;

      case int(C_pXmYmZ):
#pragma omp parallel for collapse(3)
        for( int k=1; k<=vc_comm; k++ ){
          for( int j=1; j<=vc_comm; j++ ){
            for( int i=imax-vc_comm; i<imax; i++ ){
              sendptr[_IDX_S3D(i-(imax-vc_comm),j-1,k-1,vc_comm,vc_comm,0)] = array[_IDX_S3D(i,j,k,imax,jmax,vc)];
            }
          }
        }
        break;

      case int(C_mXpYmZ):
#pragma omp parallel for collapse(3)
        for( int k=1; k<=vc_comm; k++ ){
          for( int j=jmax-vc_comm; j<jmax; j++ ){
            for( int i=1; i<=vc_comm; i++ ){
              sendptr[_IDX_S3D(i-1,j-(jmax-vc_comm),k-1,vc_comm,vc_comm,0)] = array[_IDX_S3D(i,j,k,imax,jmax,vc)];
            }
          }
        }
        break;

      case int(C_pXpYmZ):
#pragma omp parallel for collapse(3)
        for( int k=1; k<=vc_comm; k++ ){
          for( int j=jmax-vc_comm; j<jmax; j++ ){
            for( int i=imax-vc_comm; i<imax; i++ ){
              sendptr[_IDX_S3D(i-(imax-vc_comm),j-(jmax-vc_comm),k-1,vc_comm,vc_comm,0)] = array[_IDX_S3D(i,j,k,imax,jmax,vc)];
            }
          }
        }
        break;

      case int(C_mXmYpZ):
#pragma omp parallel for collapse(3)
        for( int k=kmax-vc_comm; k<kmax; k++ ){
          for( int j=1; j<=vc_comm; j++ ){
            for( int i=1; i<=vc_comm; i++ ){
              sendptr[_IDX_S3D(i-1,j-1,k-(kmax-vc_comm),vc_comm,vc_comm,0)] = array[_IDX_S3D(i,j,k,imax,jmax,vc)];
            }
          }
        }
        break;

      case int(C_pXmYpZ):
#pragma omp parallel for collapse(3)
        for( int k=kmax-vc_comm; k<kmax; k++ ){
          for( int j=1; j<=vc_comm; j++ ){
            for( int i=imax-vc_comm; i<imax; i++ ){
              sendptr[_IDX_S3D(i-(imax-vc_comm),j-1,k-(kmax-vc_comm),vc_comm,vc_comm,0)] = array[_IDX_S3D(i,j,k,imax,jmax,vc)];
            }
          }
        }
        break;

      case int(C_mXpYpZ):
#pragma omp parallel for collapse(3)
        for( int k=kmax-vc_comm; k<kmax; k++ ){
          for( int j=jmax-vc_comm; j<jmax; j++ ){
            for( int i=1; i<=vc_comm; i++ ){
              sendptr[_IDX_S3D(i-1,j-(jmax-vc_comm),k-(kmax-vc_comm),vc_comm,vc_comm,0)] = array[_IDX_S3D(i,j,k,imax,jmax,vc)];
            }
          }
        }
        break;

      case int(C_pXpYpZ):
#pragma omp parallel for collapse(3)
        for( int k=kmax-vc_comm; k<kmax; k++ ){
          for( int j=jmax-vc_comm; j<jmax; j++ ){
            for( int i=imax-vc_comm; i<imax; i++ ){
              sendptr[_IDX_S3D(i-(imax-vc_comm),j-(jmax-vc_comm),k-(kmax-vc_comm),vc_comm,vc_comm,0)] = array[_IDX_S3D(i,j,k,imax,jmax,vc)];
            }
          }
        }
        break;
      }

      // send
      if ( MPI_SUCCESS != MPI_Isend(sendptr,
                                    sz,
                                    dtype,
                                    comm_tbl[dir],
                                    tag,
                                    MPI_COMM_WORLD,
                                    &req[dir*2+1]) ) return false;

      // pointer
      ptr += sz;
    }
  }

#endif
  return true;
}

/*
 * @brief unpack send data for diagonal corner
 * @param [out] array    dest array
 * @param [in]  vc_comm  number of guide cell layer to be sent
 * @param [in]  recvbuf  recv buffer
 */
void SubDomain::unpack_SCnode(REAL_TYPE *array,
                           const int vc_comm,
                           const REAL_TYPE *recvbuf)
{
#ifdef _DIAGONAL_COMM

  int imax = size[0];
  int jmax = size[1];
  int kmax = size[2];
  int vc = halo_width;

  size_t ptr = 0;

  //// 8 corner ////
  for( int dir=int(C_mXmYmZ);dir<=int(C_pXpYpZ);dir++ )
  {
    if( comm_tbl[dir] >= 0 )
    {
      const REAL_TYPE *recvptr = &recvbuf[ptr];
      size_t sz = vc_comm * vc_comm * vc_comm;

      // unpack
      switch(dir)
      {
      case int(C_mXmYmZ):
#pragma omp parallel for collapse(3)
        for( int k=1-vc_comm; k<=0; k++ ){
          for( int j=1-vc_comm; j<=0; j++ ){
            for( int i=1-vc_comm; i<=0; i++ ){
              array[_IDX_S3D(i,j,k,imax,jmax,vc)] = recvptr[_IDX_S3D(i-(1-vc_comm),j-(1-vc_comm),k-(1-vc_comm),vc_comm,vc_comm,0)];
            }
          }
        }
        break;

      case int(C_pXmYmZ):
#pragma omp parallel for collapse(3)
        for( int k=1-vc_comm; k<=0; k++ ){
          for( int j=1-vc_comm; j<=0; j++ ){
            for( int i=imax; i<imax+vc_comm; i++ ){
              array[_IDX_S3D(i,j,k,imax,jmax,vc)] = recvptr[_IDX_S3D(i-(imax),j-(1-vc_comm),k-(1-vc_comm),vc_comm,vc_comm,0)];
            }
          }
        }
        break;

      case int(C_mXpYmZ):
#pragma omp parallel for collapse(3)
        for( int k=1-vc_comm; k<=0; k++ ){
          for( int j=jmax; j<jmax+vc_comm; j++ ){
            for( int i=1-vc_comm; i<=0; i++ ){
              array[_IDX_S3D(i,j,k,imax,jmax,vc)] = recvptr[_IDX_S3D(i-(1-vc_comm),j-(jmax),k-(1-vc_comm),vc_comm,vc_comm,0)];
            }
          }
        }
        break;

      case int(C_pXpYmZ):
#pragma omp parallel for collapse(3)
        for( int k=1-vc_comm; k<=0; k++ ){
          for( int j=jmax; j<jmax+vc_comm; j++ ){
            for( int i=imax; i<imax+vc_comm; i++ ){
              array[_IDX_S3D(i,j,k,imax,jmax,vc)] = recvptr[_IDX_S3D(i-(imax),j-(jmax),k-(1-vc_comm),vc_comm,vc_comm,0)];
            }
          }
        }
        break;

      case int(C_mXmYpZ):
#pragma omp parallel for collapse(3)
        for( int k=kmax; k<kmax+vc_comm; k++ ){
          for( int j=1-vc_comm; j<=0; j++ ){
            for( int i=1-vc_comm; i<=0; i++ ){
              array[_IDX_S3D(i,j,k,imax,jmax,vc)] = recvptr[_IDX_S3D(i-(1-vc_comm),j-(1-vc_comm),k-(kmax),vc_comm,vc_comm,0)];
            }
          }
        }
        break;

      case int(C_pXmYpZ):
#pragma omp parallel for collapse(3)
        for( int k=kmax; k<kmax+vc_comm; k++ ){
          for( int j=1-vc_comm; j<=0; j++ ){
            for( int i=imax; i<imax+vc_comm; i++ ){
              array[_IDX_S3D(i,j,k,imax,jmax,vc)] = recvptr[_IDX_S3D(i-(imax),j-(1-vc_comm),k-(kmax),vc_comm,vc_comm,0)];
            }
          }
        }
        break;

      case int(C_mXpYpZ):
#pragma omp parallel for collapse(3)
        for( int k=kmax; k<kmax+vc_comm; k++ ){
          for( int j=jmax; j<jmax+vc_comm; j++ ){
            for( int i=1-vc_comm; i<=0; i++ ){
              array[_IDX_S3D(i,j,k,imax,jmax,vc)] = recvptr[_IDX_S3D(i-(1-vc_comm),j-(jmax),k-(kmax),vc_comm,vc_comm,0)];
            }
          }
        }
        break;

      case int(C_pXpYpZ):
#pragma omp parallel for collapse(3)
        for( int k=kmax; k<kmax+vc_comm; k++ ){
          for( int j=jmax; j<jmax+vc_comm; j++ ){
            for( int i=imax; i<imax+vc_comm; i++ ){
              array[_IDX_S3D(i,j,k,imax,jmax,vc)] = recvptr[_IDX_S3D(i-(imax),j-(jmax),k-(kmax),vc_comm,vc_comm,0)];
            }
          }
        }
        break;
      }

      ptr += sz;
    }
  }

#endif
}
