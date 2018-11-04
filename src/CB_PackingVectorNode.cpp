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
 * @brief  BrickComm class
 */

#include "CB_Comm.h"


/*
 * @brief pack send data for I direction
 * @param [in]  array   source array
 * @param [in]  vc_comm number of guide cell layer to be sent
 * @param [out] sendm   send buffer of I- direction
 * @param [out] sendp   send buffer of I+ direction
 * @param [in]  nIDm    Rank number of I- direction
 * @param [in]  nIDp    Rank number of I+ direction
 */
void BrickComm::pack_VXnode(const REAL_TYPE *array,
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
      for( int k=0; k<kmax; k++ ){
        for( int j=0; j<jmax; j++ ){
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
      for( int k=0; k<kmax; k++ ){
        for( int j=0; j<jmax; j++ ){
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
void BrickComm::unpack_VXnode(REAL_TYPE *array,
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
      for( int k=0; k<kmax; k++ ){
        for( int j=0; j<jmax; j++ ){
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
      for( int k=0; k<kmax; k++ ){
        for( int j=0; j<jmax; j++ ){
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
void BrickComm::pack_VYnode(const REAL_TYPE *array,
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
      for( int k=0; k<kmax; k++ ){
        for( int j=1; j<=vc_comm; j++ ){
          for( int i=0; i<imax; i++ ){
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
      for( int k=0; k<kmax; k++ ){
        for( int j=jmax-vc_comm; j<jmax; j++ ){
          for( int i=0; i<imax; i++ ){
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
void BrickComm::unpack_VYnode(REAL_TYPE *array,
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
      for( int k=0; k<kmax; k++ ){
        for( int j=1-vc_comm; j<1; j++ ){
          for( int i=0; i<imax; i++ ){
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
      for( int k=0; k<kmax; k++ ){
        for( int j=jmax; j<jmax+vc_comm; j++ ){
          for( int i=0; i<imax; i++ ){
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
void BrickComm::pack_VZnode(const REAL_TYPE *array,
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
        for( int j=0; j<jmax; j++ ){
          for( int i=0; i<imax; i++ ){
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
        for( int j=0; j<jmax; j++ ){
          for( int i=0; i<imax; i++ ){
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
void BrickComm::unpack_VZnode(REAL_TYPE *array,
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
        for( int j=0; j<jmax; j++ ){
          for( int i=0; i<imax; i++ ){
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
        for( int j=0; j<jmax; j++ ){
          for( int i=0; i<imax; i++ ){
            array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)] = recvp[_IDX_VK(i,j,k,l,imax,jmax,kmax,vc_comm)];
          }
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
bool BrickComm::pack_VEnode(REAL_TYPE *array,
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
      size_t sz = (imax-1) * vc_comm * vc_comm * 3;

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
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=1; k<=vc_comm; k++ ){
            for( int j=1; j<=vc_comm; j++ ){
              for( int i=1; i<imax; i++ ){
                sendptr[_IDX_V3D(i-1,j-1,k-1,l,imax-1,vc_comm,vc_comm,0)] = array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)];
              }
            }
          }
        }
        break;

      case int(E_pYmZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=1; k<=vc_comm; k++ ){
            for( int j=jmax-vc_comm; j<jmax; j++ ){
              for( int i=1; i<imax; i++ ){
                sendptr[_IDX_V3D(i-1,j-(jmax-vc_comm),k-1,l,imax-1,vc_comm,vc_comm,0)] = array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)];
              }
            }
          }
        }
        break;

      case int(E_mYpZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=kmax-vc_comm; k<kmax; k++ ){
            for( int j=1; j<=vc_comm; j++ ){
              for( int i=1; i<imax; i++ ){
                sendptr[_IDX_V3D(i-1,j-1,k-(kmax-vc_comm),l,imax-1,vc_comm,vc_comm,0)] = array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)];
              }
            }
          }
        }
        break;
      case int(E_pYpZ):

#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=kmax-vc_comm; k<kmax; k++ ){
            for( int j=jmax-vc_comm; j<jmax; j++ ){
              for( int i=1; i<imax; i++ ){
                sendptr[_IDX_V3D(i-1,j-(jmax-vc_comm),k-(kmax-vc_comm),l,imax-1,vc_comm,vc_comm,0)] = array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)];
              }
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
      size_t sz = vc_comm * (jmax-1) * vc_comm * 3;

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
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=1; k<=vc_comm; k++ ){
            for( int j=1; j<jmax; j++ ){
              for( int i=1; i<=vc_comm; i++ ){
                sendptr[_IDX_V3D(i-1,j-1,k-1,l,vc_comm,jmax-1,vc_comm,0)] = array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)];
              }
            }
          }
        }
        break;

      case int(E_pXmZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=1; k<=vc_comm; k++ ){
            for( int j=1; j<jmax; j++ ){
              for( int i=imax-vc_comm; i<imax; i++ ){
                sendptr[_IDX_V3D(i-(imax-vc_comm),j-1,k-1,l,vc_comm,jmax-1,vc_comm,0)] = array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)];
              }
            }
          }
        }
        break;

      case int(E_mXpZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=kmax-vc_comm; k<kmax; k++ ){
            for( int j=1; j<jmax; j++ ){
              for( int i=1; i<=vc_comm; i++ ){
                sendptr[_IDX_V3D(i-1,j-1,k-(kmax-vc_comm),l,vc_comm,jmax-1,vc_comm,0)] = array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)];
              }
            }
          }
        }
        break;
      case int(E_pXpZ):

#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=kmax-vc_comm; k<kmax; k++ ){
            for( int j=1; j<jmax; j++ ){
              for( int i=imax-vc_comm; i<imax; i++ ){
                sendptr[_IDX_V3D(i-(imax-vc_comm),j-1,k-(kmax-vc_comm),l,vc_comm,jmax-1,vc_comm,0)] = array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)];
              }
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
      size_t sz = vc_comm * vc_comm * (kmax-1) * 3;

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
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=1; k<kmax; k++ ){
            for( int j=1; j<=vc_comm; j++ ){
              for( int i=1; i<=vc_comm; i++ ){
                sendptr[_IDX_V3D(i-1,j-1,k-1,l,vc_comm,vc_comm,kmax-1,0)] = array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)];
              }
            }
          }
        }
        break;

      case int(E_pXmY):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=1; k<kmax; k++ ){
            for( int j=1; j<=vc_comm; j++ ){
              for( int i=imax-vc_comm; i<imax; i++ ){
                sendptr[_IDX_V3D(i-(imax-vc_comm),j-1,k-1,l,vc_comm,vc_comm,kmax-1,0)] = array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)];
              }
            }
          }
        }
        break;

      case int(E_mXpY):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=1; k<kmax; k++ ){
            for( int j=jmax-vc_comm; j<jmax; j++ ){
              for( int i=1; i<=vc_comm; i++ ){
                sendptr[_IDX_V3D(i-1,j-(jmax-vc_comm),k-1,l,vc_comm,vc_comm,kmax-1,0)] = array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)];
              }
            }
          }
        }
        break;
      case int(E_pXpY):

#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=1; k<kmax; k++ ){
            for( int j=jmax-vc_comm; j<jmax; j++ ){
              for( int i=imax-vc_comm; i<imax; i++ ){
                sendptr[_IDX_V3D(i-(imax-vc_comm),j-(jmax-vc_comm),k-1,l,vc_comm,vc_comm,kmax-1,0)] = array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)];
              }
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
void BrickComm::unpack_VEnode(REAL_TYPE *array,
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
      size_t sz = (imax-1) * vc_comm * vc_comm * 3;

      // unpack
      switch(dir)
      {
      case int(E_mYmZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=1-vc_comm; k<=0; k++ ){
            for( int j=1-vc_comm; j<=0; j++ ){
              for( int i=1; i<imax; i++ ){
                array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)] = recvptr[_IDX_V3D(i-1,j-(1-vc_comm),k-(1-vc_comm),l,imax-1,vc_comm,vc_comm,0)];
              }
            }
          }
        }
        break;

      case int(E_pYmZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=1-vc_comm; k<=0; k++ ){
            for( int j=jmax; j<jmax+vc_comm; j++ ){
              for( int i=1; i<imax; i++ ){
                array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)] = recvptr[_IDX_V3D(i-1,j-(jmax),k-(1-vc_comm),l,imax-1,vc_comm,vc_comm,0)];
              }
            }
          }
        }
        break;

      case int(E_mYpZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=kmax; k<kmax+vc_comm; k++ ){
            for( int j=1-vc_comm; j<=0; j++ ){
              for( int i=1; i<imax; i++ ){
                array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)] = recvptr[_IDX_V3D(i-1,j-(1-vc_comm),k-(kmax),l,imax-1,vc_comm,vc_comm,0)];
              }
            }
          }
        }
        break;

      case int(E_pYpZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=kmax; k<kmax+vc_comm; k++ ){
            for( int j=jmax; j<jmax+vc_comm; j++ ){
              for( int i=1; i<imax; i++ ){
                array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)] = recvptr[_IDX_V3D(i-1,j-(jmax),k-(kmax),l,imax-1,vc_comm,vc_comm,0)];
              }
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
      size_t sz = vc_comm * (jmax-1) * vc_comm * 3;

      // unpack
      switch(dir)
      {
      case int(E_mXmZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=1-vc_comm; k<=0; k++ ){
            for( int j=1; j<jmax; j++ ){
              for( int i=1-vc_comm; i<=0; i++ ){
                array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)] = recvptr[_IDX_V3D(i-(1-vc_comm),j-1,k-(1-vc_comm),l,vc_comm,jmax-1,vc_comm,0)];
              }
            }
          }
        }
        break;

      case int(E_pXmZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=1-vc_comm; k<=0; k++ ){
            for( int j=1; j<jmax; j++ ){
              for( int i=imax; i<imax+vc_comm; i++ ){
                array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)] = recvptr[_IDX_V3D(i-(imax),j-1,k-(1-vc_comm),l,vc_comm,jmax-1,vc_comm,0)];
              }
            }
          }
        }
        break;

      case int(E_mXpZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=kmax; k<kmax+vc_comm; k++ ){
            for( int j=1; j<jmax; j++ ){
              for( int i=1-vc_comm; i<=0; i++ ){
                array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)] = recvptr[_IDX_V3D(i-(1-vc_comm),j-1,k-(kmax),l,vc_comm,jmax-1,vc_comm,0)];
              }
            }
          }
        }
        break;

      case int(E_pXpZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=kmax; k<kmax+vc_comm; k++ ){
            for( int j=1; j<jmax; j++ ){
              for( int i=imax; i<imax+vc_comm; i++ ){
                array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)] = recvptr[_IDX_V3D(i-(imax),j-1,k-(kmax),l,vc_comm,jmax-1,vc_comm,0)];
              }
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
      size_t sz = vc_comm * vc_comm * (kmax-1) * 3;

      // unpack
      switch(dir)
      {
      case int(E_mXmY):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=1; k<kmax; k++ ){
            for( int j=1-vc_comm; j<=0; j++ ){
              for( int i=1-vc_comm; i<=0; i++ ){
                array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)] = recvptr[_IDX_V3D(i-(1-vc_comm),j-(1-vc_comm),k-1,l,vc_comm,vc_comm,kmax-1,0)];
              }
            }
          }
        }
        break;

      case int(E_pXmY):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=1; k<kmax; k++ ){
            for( int j=1-vc_comm; j<=0; j++ ){
              for( int i=imax; i<imax+vc_comm; i++ ){
                array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)] = recvptr[_IDX_V3D(i-(imax),j-(1-vc_comm),k-1,l,vc_comm,vc_comm,kmax-1,0)];
              }
            }
          }
        }
        break;

      case int(E_mXpY):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=1; k<kmax; k++ ){
            for( int j=jmax; j<jmax+vc_comm; j++ ){
              for( int i=1-vc_comm; i<=0; i++ ){
                array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)] = recvptr[_IDX_V3D(i-(1-vc_comm),j-(jmax),k-1,l,vc_comm,vc_comm,kmax-1,0)];
              }
            }
          }
        }
        break;

      case int(E_pXpY):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=1; k<kmax; k++ ){
            for( int j=jmax; j<jmax+vc_comm; j++ ){
              for( int i=imax; i<imax+vc_comm; i++ ){
                array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)] = recvptr[_IDX_V3D(i-(imax),j-(jmax),k-1,l,vc_comm,vc_comm,kmax-1,0)];
              }
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
bool BrickComm::pack_VCnode(REAL_TYPE *array,
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
      size_t sz = vc_comm * vc_comm * vc_comm * 3;

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
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=1; k<=vc_comm; k++ ){
            for( int j=1; j<=vc_comm; j++ ){
              for( int i=1; i<=vc_comm; i++ ){
                sendptr[_IDX_V3D(i-1,j-1,k-1,l,vc_comm,vc_comm,vc_comm,0)] = array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)];
              }
            }
          }
        }
        break;

      case int(C_pXmYmZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=1; k<=vc_comm; k++ ){
            for( int j=1; j<=vc_comm; j++ ){
              for( int i=imax-vc_comm; i<imax; i++ ){
                sendptr[_IDX_V3D(i-(imax-vc_comm),j-1,k-1,l,vc_comm,vc_comm,vc_comm,0)] = array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)];
              }
            }
          }
        }
        break;

      case int(C_mXpYmZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=1; k<=vc_comm; k++ ){
            for( int j=jmax-vc_comm; j<jmax; j++ ){
              for( int i=1; i<=vc_comm; i++ ){
                sendptr[_IDX_V3D(i-1,j-(jmax-vc_comm),k-1,l,vc_comm,vc_comm,vc_comm,0)] = array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)];
              }
            }
          }
        }
        break;

      case int(C_pXpYmZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=1; k<=vc_comm; k++ ){
            for( int j=jmax-vc_comm; j<jmax; j++ ){
              for( int i=imax-vc_comm; i<imax; i++ ){
                sendptr[_IDX_V3D(i-(imax-vc_comm),j-(jmax-vc_comm),k-1,l,vc_comm,vc_comm,vc_comm,0)] = array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)];
              }
            }
          }
        }
        break;

      case int(C_mXmYpZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=kmax-vc_comm; k<kmax; k++ ){
            for( int j=1; j<=vc_comm; j++ ){
              for( int i=1; i<=vc_comm; i++ ){
                sendptr[_IDX_V3D(i-1,j-1,k-(kmax-vc_comm),l,vc_comm,vc_comm,vc_comm,0)] = array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)];
              }
            }
          }
        }
        break;

      case int(C_pXmYpZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=kmax-vc_comm; k<kmax; k++ ){
            for( int j=1; j<=vc_comm; j++ ){
              for( int i=imax-vc_comm; i<imax; i++ ){
                sendptr[_IDX_V3D(i-(imax-vc_comm),j-1,k-(kmax-vc_comm),l,vc_comm,vc_comm,vc_comm,0)] = array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)];
              }
            }
          }
        }
        break;

      case int(C_mXpYpZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=kmax-vc_comm; k<kmax; k++ ){
            for( int j=jmax-vc_comm; j<jmax; j++ ){
              for( int i=1; i<=vc_comm; i++ ){
                sendptr[_IDX_V3D(i-1,j-(jmax-vc_comm),k-(kmax-vc_comm),l,vc_comm,vc_comm,vc_comm,0)] = array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)];
              }
            }
          }
        }
        break;

      case int(C_pXpYpZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=kmax-vc_comm; k<kmax; k++ ){
            for( int j=jmax-vc_comm; j<jmax; j++ ){
              for( int i=imax-vc_comm; i<imax; i++ ){
                sendptr[_IDX_V3D(i-(imax-vc_comm),j-(jmax-vc_comm),k-(kmax-vc_comm),l,vc_comm,vc_comm,vc_comm,0)] = array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)];
              }
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
void BrickComm::unpack_VCnode(REAL_TYPE *array,
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
      size_t sz = vc_comm * vc_comm * vc_comm * 3;

      // unpack
      switch(dir)
      {
      case int(C_mXmYmZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=1-vc_comm; k<=0; k++ ){
            for( int j=1-vc_comm; j<=0; j++ ){
              for( int i=1-vc_comm; i<=0; i++ ){
                array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)] = recvptr[_IDX_V3D(i-(1-vc_comm),j-(1-vc_comm),k-(1-vc_comm),l,vc_comm,vc_comm,vc_comm,0)];
              }
            }
          }
        }
        break;

      case int(C_pXmYmZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=1-vc_comm; k<=0; k++ ){
            for( int j=1-vc_comm; j<=0; j++ ){
              for( int i=imax; i<imax+vc_comm; i++ ){
                array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)] = recvptr[_IDX_V3D(i-(imax),j-(1-vc_comm),k-(1-vc_comm),l,vc_comm,vc_comm,vc_comm,0)];
              }
            }
          }
        }
        break;

      case int(C_mXpYmZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=1-vc_comm; k<=0; k++ ){
            for( int j=jmax; j<jmax+vc_comm; j++ ){
              for( int i=1-vc_comm; i<=0; i++ ){
                array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)] = recvptr[_IDX_V3D(i-(1-vc_comm),j-(jmax),k-(1-vc_comm),l,vc_comm,vc_comm,vc_comm,0)];
              }
            }
          }
        }
        break;

      case int(C_pXpYmZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=1-vc_comm; k<=0; k++ ){
            for( int j=jmax; j<jmax+vc_comm; j++ ){
              for( int i=imax; i<imax+vc_comm; i++ ){
                array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)] = recvptr[_IDX_V3D(i-(imax),j-(jmax),k-(1-vc_comm),l,vc_comm,vc_comm,vc_comm,0)];
              }
            }
          }
        }
        break;

      case int(C_mXmYpZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=kmax; k<kmax+vc_comm; k++ ){
            for( int j=1-vc_comm; j<=0; j++ ){
              for( int i=1-vc_comm; i<=0; i++ ){
                array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)] = recvptr[_IDX_V3D(i-(1-vc_comm),j-(1-vc_comm),k-(kmax),l,vc_comm,vc_comm,vc_comm,0)];
              }
            }
          }
        }
        break;

      case int(C_pXmYpZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=kmax; k<kmax+vc_comm; k++ ){
            for( int j=1-vc_comm; j<=0; j++ ){
              for( int i=imax; i<imax+vc_comm; i++ ){
                array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)] = recvptr[_IDX_V3D(i-(imax),j-(1-vc_comm),k-(kmax),l,vc_comm,vc_comm,vc_comm,0)];
              }
            }
          }
        }
        break;

      case int(C_mXpYpZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=kmax; k<kmax+vc_comm; k++ ){
            for( int j=jmax; j<jmax+vc_comm; j++ ){
              for( int i=1-vc_comm; i<=0; i++ ){
                array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)] = recvptr[_IDX_V3D(i-(1-vc_comm),j-(jmax),k-(kmax),l,vc_comm,vc_comm,vc_comm,0)];
              }
            }
          }
        }
        break;

      case int(C_pXpYpZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=kmax; k<kmax+vc_comm; k++ ){
            for( int j=jmax; j<jmax+vc_comm; j++ ){
              for( int i=imax; i<imax+vc_comm; i++ ){
                array[_IDX_V3D(i,j,k,l,imax,jmax,kmax,vc)] = recvptr[_IDX_V3D(i-(imax),j-(jmax),k-(kmax),l,vc_comm,vc_comm,vc_comm,0)];
              }
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
