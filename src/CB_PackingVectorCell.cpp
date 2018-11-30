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
 * @file   CB_PackingVectorCell.cpp
 * @brief  BrickComm class
 */

#include "CB_Comm.h"


/*
 * @brief pack send data for I direction
 * @param [in]  array   source array
 * @param [in]  gc      number of guide cell layer to be sent
 * @param [out] sendm   send buffer of I- direction
 * @param [out] sendp   send buffer of I+ direction
 * @param [in]  nIDm    Rank number of I- direction
 * @param [in]  nIDp    Rank number of I+ direction
 */
void BrickComm::pack_VXcell(const REAL_TYPE *array,
                        const int gc,
                        REAL_TYPE *sendm,
                        REAL_TYPE *sendp,
                        const int nIDm,
                        const int nIDp)
{
  int NI = size[0];
  int NJ = size[1];
  int NK = size[2];
  int VC = halo_width;

  if( nIDm >= 0 )
  {
#pragma omp parallel for collapse(3)
    for (int l=0; l<3; l++) {
      for( int k=0; k<NK; k++ ){
        for( int j=0; j<NJ; j++ ){
          for( int i=0; i<gc; i++ ){
            sendm[_IDX_VI(i,j,k,l,NJ,NK,gc)] = array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)];
          }
        }
      }
    }
  }

  if( nIDp >= 0 )
  {
#pragma omp parallel for collapse(3)
    for (int l=0; l<3; l++) {
      for( int k=0; k<NK; k++ ){
        for( int j=0; j<NJ; j++ ){
          for( int i=0; i<gc; i++ ){
            sendp[_IDX_VI(i,j,k,l,NJ,NK,gc)] = array[_IDX_V3D(NI-gc+i,j,k,l,NI,NJ,NK,VC)];
          }
        }
      }
    }
  }
}


/*
 * @brief unpack send data for I direction
 * @param [in,out]  array   dest array
 * @param [in]  gc      number of guide cell layer to be sent
 * @param [in]  recvm   recv buffer of I- direction
 * @param [in]  recvp   recv buffer of I+ direction
 * @param [in]  nIDm    Rank number of I- direction
 * @param [in]  nIDp    Rank number of I+ direction
 */
void BrickComm::unpack_VXcell(REAL_TYPE *array,
                          const int gc,
                          const REAL_TYPE *recvm,
                          const REAL_TYPE *recvp,
                          const int nIDm,
                          const int nIDp)
{
  int NI = size[0];
  int NJ = size[1];
  int NK = size[2];
  int VC = halo_width;

  if( nIDm >= 0 )
  {
#pragma omp parallel for collapse(3)
    for (int l=0; l<3; l++) {
      for( int k=0; k<NK; k++ ){
        for( int j=0; j<NJ; j++ ){
          for( int i=0; i<gc; i++ ){
            array[_IDX_V3D(i-gc,j,k,l,NI,NJ,NK,VC)] = recvm[_IDX_VI(i,j,k,l,NJ,NK,gc)];
          }
        }
      }
    }
  }

  if( nIDp >= 0 )
  {
#pragma omp parallel for collapse(3)
    for (int l=0; l<3; l++) {
      for( int k=0; k<NK; k++ ){
        for( int j=0; j<NJ; j++ ){
          for( int i=0; i<gc; i++ ){
            array[_IDX_V3D(NI+i,j,k,l,NI,NJ,NK,VC)] = recvp[_IDX_VI(i,j,k,l,NJ,NK,gc)];
          }
        }
      }
    }
  }
}

/*
 * @brief pack send data for J direction
 * @param [in]  array   source array
 * @param [in]  gc number of guide cell layer to be sent
 * @param [out] sendm   send buffer of J- direction
 * @param [out] sendp   send buffer of J+ direction
 * @param [in]  nIDm    Rank number of J- direction
 * @param [in]  nIDp    Rank number of J+ direction
 */
void BrickComm::pack_VYcell(const REAL_TYPE *array,
                        const int gc,
                        REAL_TYPE *sendm,
                        REAL_TYPE *sendp,
                        const int nIDm,
                        const int nIDp)
{
  int NI = size[0];
  int NJ = size[1];
  int NK = size[2];
  int VC = halo_width;

  if( nIDm >= 0 )
  {
#pragma omp parallel for collapse(3)
    for (int l=0; l<3; l++) {
      for( int k=0; k<NK; k++ ){
        for( int j=0; j<gc; j++ ){
          for( int i=0; i<NI; i++ ){
            sendm[_IDX_VJ(i,j,k,l,NI,0,NK,gc)] = array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)];
          }
        }
      }
    }
  }

  if( nIDp >= 0 )
  {
#pragma omp parallel for collapse(3)
    for (int l=0; l<3; l++) {
      for( int k=0; k<NK; k++ ){
        for( int j=NJ-gc; j<NJ; j++ ){
          for( int i=0; i<NI; i++ ){
            sendp[_IDX_VJ(i,j,k,l,NI,NJ-gc,NK,gc)] = array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)];
          }
        }
      }
    }
  }
}


/*
 * @brief unpack send data for J direction
 * @param [in,out]  array   dest array
 * @param [in]  gc number of guide cell layer to be sent
 * @param [in]  recvm   recv buffer of J- direction
 * @param [in]  recvp   recv buffer of J+ direction
 * @param [in]  nIDm    Rank number of J- direction
 * @param [in]  nIDp    Rank number of J+ direction
 */
void BrickComm::unpack_VYcell(REAL_TYPE *array,
                          const int gc,
                          const REAL_TYPE *recvm,
                          const REAL_TYPE *recvp,
                          const int nIDm,
                          const int nIDp)
{
  int NI = size[0];
  int NJ = size[1];
  int NK = size[2];
  int VC = halo_width;

  if( nIDm >= 0 )
  {
#pragma omp parallel for collapse(3)
    for (int l=0; l<3; l++) {
      for( int k=0; k<NK; k++ ){
        for( int j=0-gc; j<0; j++ ){
          for( int i=0; i<NI; i++ ){
            array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)] = recvm[_IDX_VJ(i,j,k,l,NI,0-gc,NK,gc)];
          }
        }
      }
    }
  }

  if( nIDp >= 0 )
  {
#pragma omp parallel for collapse(3)
    for (int l=0; l<3; l++) {
      for( int k=0; k<NK; k++ ){
        for( int j=NJ; j<NJ+gc; j++ ){
          for( int i=0; i<NI; i++ ){
            array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)] = recvp[_IDX_VJ(i,j,k,l,NI,NJ,NK,gc)];
          }
        }
      }
    }
  }
}


/*
 * @brief pack send data for K direction
 * @param [in]  array   source array
 * @param [in]  gc number of guide cell layer actually to be sent
 * @param [out] sendm   send buffer of K- direction
 * @param [out] sendp   send buffer of K+ direction
 * @param [in]  nIDm    Rank number of K- direction
 * @param [in]  nIDp    Rank number of K+ direction
 */
void BrickComm::pack_VZcell(const REAL_TYPE *array,
                        const int gc,
                        REAL_TYPE *sendm,
                        REAL_TYPE *sendp,
                        const int nIDm,
                        const int nIDp)
{
  int NI = size[0];
  int NJ = size[1];
  int NK = size[2];
  int VC = halo_width;

  if( nIDm >= 0 )
  {
#pragma omp parallel for collapse(3)
    for (int l=0; l<3; l++) {
      for( int k=0; k<gc; k++ ){
        for( int j=0; j<NJ; j++ ){
          for( int i=0; i<NI; i++ ){
            sendm[_IDX_VK(i,j,k,l,NI,NJ,0,gc)] = array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)];
          }
        }
      }
    }
  }

  if( nIDp >= 0 )
  {
#pragma omp parallel for collapse(3)
    for (int l=0; l<3; l++) {
      for( int k=NK-gc; k<NK; k++ ){
        for( int j=0; j<NJ; j++ ){
          for( int i=0; i<NI; i++ ){
            sendp[_IDX_VK(i,j,k,l,NI,NJ,NK-gc,gc)] = array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)];
          }
        }
      }
    }
  }
}


/*
 * @brief unpack send data for K direction
 * @param [in,out]  array   dest array
 * @param [in]  gc number of guide cell layer to be sent
 * @param [in]  recvm   recv buffer of K- direction
 * @param [in]  recvp   recv buffer of K+ direction
 * @param [in]  nIDm    Rank number of K- direction
 * @param [in]  nIDp    Rank number of K+ direction
 */
void BrickComm::unpack_VZcell(REAL_TYPE *array,
                          const int gc,
                          const REAL_TYPE *recvm,
                          const REAL_TYPE *recvp,
                          const int nIDm,
                          const int nIDp)
{
  int NI = size[0];
  int NJ = size[1];
  int NK = size[2];
  int VC = halo_width;

  if( nIDm >= 0 )
  {
#pragma omp parallel for collapse(3)
    for (int l=0; l<3; l++) {
      for( int k=0-gc; k<0; k++ ){
        for( int j=0; j<NJ; j++ ){
          for( int i=0; i<NI; i++ ){
            array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)] = recvm[_IDX_VK(i,j,k,l,NI,NJ,0-gc,gc)];
          }
        }
      }
    }
  }

  if( nIDp >= 0 )
  {
#pragma omp parallel for collapse(3)
    for (int l=0; l<3; l++) {
      for( int k=NK; k<NK+gc; k++ ){
        for( int j=0; j<NJ; j++ ){
          for( int i=0; i<NI; i++ ){
            array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)] = recvp[_IDX_VK(i,j,k,l,NI,NJ,NK,gc)];
          }
        }
      }
    }
  }
}

/*
 * @brief pack send data for diagonal edge
 * @param [in]  array    source array
 * @param [in]  gc  number of guide cell layer to be sent
 * @param [out] sendbuf  send buffer
 * @param [out] recvbuf  recv buffer
 * @param [out] req      Array of MPI request
 * @retval true-success, false-fail
 */
bool BrickComm::pack_VEcell(REAL_TYPE *array,
                        const int gc,
                        REAL_TYPE *sendbuf,
                        REAL_TYPE *recvbuf,
                        MPI_Request *req)
{
#ifdef _DIAGONAL_COMM

  int NI = size[0];
  int NJ = size[1];
  int NK = size[2];
  int VC = halo_width;
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
      size_t sz = NI * gc * gc * 3;

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
          for( int k=0; k<gc; k++ ){
            for( int j=0; j<gc; j++ ){
              for( int i=0; i<NI; i++ ){
                sendptr[_IDX_V3D(i,j,k,l,NI,gc,gc,0)] = array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)];
              }
            }
          }
        }
        break;

      case int(E_pYmZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=0; k<gc; k++ ){
            for( int j=NJ-gc; j<NJ; j++ ){
              for( int i=0; i<NI; i++ ){
                sendptr[_IDX_V3D(i,j-(NJ-gc),k,l,NI,gc,gc,0)] = array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)];
              }
            }
          }
        }
        break;

      case int(E_mYpZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=NK-gc; k<NK; k++ ){
            for( int j=0; j<gc; j++ ){
              for( int i=0; i<NI; i++ ){
                sendptr[_IDX_V3D(i,j,k-(NK-gc),l,NI,gc,gc,0)] = array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)];
              }
            }
          }
        }
        break;
      case int(E_pYpZ):

#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=NK-gc; k<NK; k++ ){
            for( int j=NJ-gc; j<NJ; j++ ){
              for( int i=0; i<NI; i++ ){
                sendptr[_IDX_V3D(i,j-(NJ-gc),k-(NK-gc),l,NI,gc,gc,0)] = array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)];
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
      size_t sz = gc * NJ * gc * 3;

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
          for( int k=0; k<gc; k++ ){
            for( int j=0; j<NJ; j++ ){
              for( int i=0; i<gc; i++ ){
                sendptr[_IDX_V3D(i,j,k,l,gc,NJ,gc,0)] = array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)];
              }
            }
          }
        }
        break;

      case int(E_pXmZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=0; k<gc; k++ ){
            for( int j=0; j<NJ; j++ ){
              for( int i=NI-gc; i<NI; i++ ){
                sendptr[_IDX_V3D(i-(NI-gc),j,k,l,gc,NJ,gc,0)] = array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)];
              }
            }
          }
        }
        break;

      case int(E_mXpZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=NK-gc; k<NK; k++ ){
            for( int j=0; j<NJ; j++ ){
              for( int i=0; i<gc; i++ ){
                sendptr[_IDX_V3D(i,j,k-(NK-gc),l,gc,NJ,gc,0)] = array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)];
              }
            }
          }
        }
        break;
      case int(E_pXpZ):

#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=NK-gc; k<NK; k++ ){
            for( int j=0; j<NJ; j++ ){
              for( int i=NI-gc; i<NI; i++ ){
                sendptr[_IDX_V3D(i-(NI-gc),j,k-(NK-gc),l,gc,NJ,gc,0)] = array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)];
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
      size_t sz = gc * gc * NK * 3;

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
          for( int k=0; k<NK; k++ ){
            for( int j=0; j<gc; j++ ){
              for( int i=0; i<gc; i++ ){
                sendptr[_IDX_V3D(i,j,k,l,gc,gc,NK,0)] = array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)];
              }
            }
          }
        }
        break;

      case int(E_pXmY):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=0; k<NK; k++ ){
            for( int j=0; j<gc; j++ ){
              for( int i=NI-gc; i<NI; i++ ){
                sendptr[_IDX_V3D(i-(NI-gc),j,k,l,gc,gc,NK,0)] = array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)];
              }
            }
          }
        }
        break;

      case int(E_mXpY):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=0; k<NK; k++ ){
            for( int j=NJ-gc; j<NJ; j++ ){
              for( int i=0; i<gc; i++ ){
                sendptr[_IDX_V3D(i,j-(NJ-gc),k,l,gc,gc,NK,0)] = array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)];
              }
            }
          }
        }
        break;
      case int(E_pXpY):

#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=0; k<NK; k++ ){
            for( int j=NJ-gc; j<NJ; j++ ){
              for( int i=NI-gc; i<NI; i++ ){
                sendptr[_IDX_V3D(i-(NI-gc),j-(NJ-gc),k,l,gc,gc,NK,0)] = array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)];
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
 * @param [in]  gc  number of guide cell layer to be sent
 * @param [in]  recvbuf  recv buffer
 */
void BrickComm::unpack_VEcell(REAL_TYPE *array,
                          const int gc,
                          const REAL_TYPE *recvbuf)
{
#ifdef _DIAGONAL_COMM

  int NI = size[0];
  int NJ = size[1];
  int NK = size[2];
  int VC = halo_width;

  size_t ptr = 0;

  //// X edge ////
  for( int dir=int(E_mYmZ);dir<=int(E_pYpZ);dir++ )
  {
    if( comm_tbl[dir] >= 0 )
    {
      const REAL_TYPE *recvptr = &recvbuf[ptr];
      size_t sz = NI * gc * gc * 3;

      // unpack
      switch(dir)
      {
      case int(E_mYmZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=0-gc; k<0; k++ ){
            for( int j=0-gc; j<0; j++ ){
              for( int i=0; i<NI; i++ ){
                array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)] = recvptr[_IDX_V3D(i,j-(0-gc),k-(0-gc),l,NI,gc,gc,0)];
              }
            }
          }
        }
        break;

      case int(E_pYmZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=0-gc; k<0; k++ ){
            for( int j=NJ; j<NJ+gc; j++ ){
              for( int i=0; i<NI; i++ ){
                array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)] = recvptr[_IDX_V3D(i,j-(NJ),k-(0-gc),l,NI,gc,gc,0)];
              }
            }
          }
        }
        break;

      case int(E_mYpZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=NK; k<NK+gc; k++ ){
            for( int j=0-gc; j<0; j++ ){
              for( int i=0; i<NI; i++ ){
                array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)] = recvptr[_IDX_V3D(i,j-(0-gc),k-(NK),l,NI,gc,gc,0)];
              }
            }
          }
        }
        break;

      case int(E_pYpZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=NK; k<NK+gc; k++ ){
            for( int j=NJ; j<NJ+gc; j++ ){
              for( int i=0; i<NI; i++ ){
                array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)] = recvptr[_IDX_V3D(i,j-(NJ),k-(NK),l,NI,gc,gc,0)];
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
      size_t sz = gc * NJ * gc * 3;

      // unpack
      switch(dir)
      {
      case int(E_mXmZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=0-gc; k<0; k++ ){
            for( int j=0; j<NJ; j++ ){
              for( int i=0-gc; i<0; i++ ){
                array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)] = recvptr[_IDX_V3D(i-(0-gc),j,k-(0-gc),l,gc,NJ,gc,0)];
              }
            }
          }
        }
        break;

      case int(E_pXmZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=0-gc; k<0; k++ ){
            for( int j=0; j<NJ; j++ ){
              for( int i=NI; i<NI+gc; i++ ){
                array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)] = recvptr[_IDX_V3D(i-(NI),j,k-(0-gc),l,gc,NJ,gc,0)];
              }
            }
          }
        }
        break;

      case int(E_mXpZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=NK; k<NK+gc; k++ ){
            for( int j=0; j<NJ; j++ ){
              for( int i=0-gc; i<0; i++ ){
                array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)] = recvptr[_IDX_V3D(i-(0-gc),j,k-(NK),l,gc,NJ,gc,0)];
              }
            }
          }
        }
        break;

      case int(E_pXpZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=NK; k<NK+gc; k++ ){
            for( int j=0; j<NJ; j++ ){
              for( int i=NI; i<NI+gc; i++ ){
                array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)] = recvptr[_IDX_V3D(i-(NI),j,k-(NK),l,gc,NJ,gc,0)];
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
      size_t sz = gc * gc * NK * 3;

      // unpack
      switch(dir)
      {
      case int(E_mXmY):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=0; k<NK; k++ ){
            for( int j=0-gc; j<0; j++ ){
              for( int i=0-gc; i<0; i++ ){
                array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)] = recvptr[_IDX_V3D(i-(0-gc),j-(0-gc),k,l,gc,gc,NK,0)];
              }
            }
          }
        }
        break;

      case int(E_pXmY):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=0; k<NK; k++ ){
            for( int j=0-gc; j<0; j++ ){
              for( int i=NI; i<NI+gc; i++ ){
                array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)] = recvptr[_IDX_V3D(i-(NI),j-(0-gc),k,l,gc,gc,NK,0)];
              }
            }
          }
        }
        break;

      case int(E_mXpY):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=0; k<NK; k++ ){
            for( int j=NJ; j<NJ+gc; j++ ){
              for( int i=0-gc; i<0; i++ ){
                array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)] = recvptr[_IDX_V3D(i-(0-gc),j-(NJ),k,l,gc,gc,NK,0)];
              }
            }
          }
        }
        break;

      case int(E_pXpY):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=0; k<NK; k++ ){
            for( int j=NJ; j<NJ+gc; j++ ){
              for( int i=NI; i<NI+gc; i++ ){
                array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)] = recvptr[_IDX_V3D(i-(NI),j-(NJ),k,l,gc,gc,NK,0)];
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
 * @param [in]  gc  number of guide cell layer to be sent
 * @param [out] sendbuf  send buffer
 * @param [out] recvbuf  recv buffer
 * @param [out] req      Array of MPI request
 * @retval true-success, false-fail
 */
bool BrickComm::pack_VCcell(REAL_TYPE *array,
                        const int gc,
                        REAL_TYPE *sendbuf,
                        REAL_TYPE *recvbuf,
                        MPI_Request *req)
{
#ifdef _DIAGONAL_COMM

  int NI = size[0];
  int NJ = size[1];
  int NK = size[2];
  int VC = halo_width;
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
      size_t sz = gc * gc * gc * 3;

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
          for( int k=0; k<gc; k++ ){
            for( int j=0; j<gc; j++ ){
              for( int i=0; i<gc; i++ ){
                sendptr[_IDX_V3D(i,j,k,l,gc,gc,gc,0)] = array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)];
              }
            }
          }
        }
        break;

      case int(C_pXmYmZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=0; k<gc; k++ ){
            for( int j=0; j<gc; j++ ){
              for( int i=NI-gc; i<NI; i++ ){
                sendptr[_IDX_V3D(i-(NI-gc),j,k,l,gc,gc,gc,0)] = array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)];
              }
            }
          }
        }
        break;

      case int(C_mXpYmZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=0; k<gc; k++ ){
            for( int j=NJ-gc; j<NJ; j++ ){
              for( int i=0; i<gc; i++ ){
                sendptr[_IDX_V3D(i,j-(NJ-gc),k,l,gc,gc,gc,0)] = array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)];
              }
            }
          }
        }
        break;

      case int(C_pXpYmZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=0; k<gc; k++ ){
            for( int j=NJ-gc; j<NJ; j++ ){
              for( int i=NI-gc; i<NI; i++ ){
                sendptr[_IDX_V3D(i-(NI-gc),j-(NJ-gc),k,l,gc,gc,gc,0)] = array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)];
              }
            }
          }
        }
        break;

      case int(C_mXmYpZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=NK-gc; k<NK; k++ ){
            for( int j=0; j<gc; j++ ){
              for( int i=0; i<gc; i++ ){
                sendptr[_IDX_V3D(i,j,k-(NK-gc),l,gc,gc,gc,0)] = array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)];
              }
            }
          }
        }
        break;

      case int(C_pXmYpZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=NK-gc; k<NK; k++ ){
            for( int j=0; j<gc; j++ ){
              for( int i=NI-gc; i<NI; i++ ){
                sendptr[_IDX_V3D(i-(NI-gc),j,k-(NK-gc),l,gc,gc,gc,0)] = array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)];
              }
            }
          }
        }
        break;

      case int(C_mXpYpZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=NK-gc; k<NK; k++ ){
            for( int j=NJ-gc; j<NJ; j++ ){
              for( int i=0; i<gc; i++ ){
                sendptr[_IDX_V3D(i,j-(NJ-gc),k-(NK-gc),l,gc,gc,gc,0)] = array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)];
              }
            }
          }
        }
        break;

      case int(C_pXpYpZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=NK-gc; k<NK; k++ ){
            for( int j=NJ-gc; j<NJ; j++ ){
              for( int i=NI-gc; i<NI; i++ ){
                sendptr[_IDX_V3D(i-(NI-gc),j-(NJ-gc),k-(NK-gc),l,gc,gc,gc,0)] = array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)];
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
 * @param [in]  gc  number of guide cell layer to be sent
 * @param [in]  recvbuf  recv buffer
 */
void BrickComm::unpack_VCcell(REAL_TYPE *array,
                          const int gc,
                          const REAL_TYPE *recvbuf)
{
#ifdef _DIAGONAL_COMM

  int NI = size[0];
  int NJ = size[1];
  int NK = size[2];
  int VC = halo_width;

  size_t ptr = 0;

  //// 8 corner ////
  for( int dir=int(C_mXmYmZ);dir<=int(C_pXpYpZ);dir++ )
  {
    if( comm_tbl[dir] >= 0 )
    {
      const REAL_TYPE *recvptr = &recvbuf[ptr];
      size_t sz = gc * gc * gc * 3;

      // unpack
      switch(dir)
      {
      case int(C_mXmYmZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=0-gc; k<0; k++ ){
            for( int j=0-gc; j<0; j++ ){
              for( int i=0-gc; i<0; i++ ){
                array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)] = recvptr[_IDX_V3D(i-(0-gc),j-(0-gc),k-(0-gc),l,gc,gc,gc,0)];
              }
            }
          }
        }
        break;

      case int(C_pXmYmZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=0-gc; k<0; k++ ){
            for( int j=0-gc; j<0; j++ ){
              for( int i=NI; i<NI+gc; i++ ){
                array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)] = recvptr[_IDX_V3D(i-(NI),j-(0-gc),k-(0-gc),l,gc,gc,gc,0)];
              }
            }
          }
        }
        break;

      case int(C_mXpYmZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=0-gc; k<0; k++ ){
            for( int j=NJ; j<NJ+gc; j++ ){
              for( int i=0-gc; i<0; i++ ){
                array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)] = recvptr[_IDX_V3D(i-(0-gc),j-(NJ),k-(0-gc),l,gc,gc,gc,0)];
              }
            }
          }
        }
        break;

      case int(C_pXpYmZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=0-gc; k<0; k++ ){
            for( int j=NJ; j<NJ+gc; j++ ){
              for( int i=NI; i<NI+gc; i++ ){
                array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)] = recvptr[_IDX_V3D(i-(NI),j-(NJ),k-(0-gc),l,gc,gc,gc,0)];
              }
            }
          }
        }
        break;

      case int(C_mXmYpZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=NK; k<NK+gc; k++ ){
            for( int j=0-gc; j<0; j++ ){
              for( int i=0-gc; i<0; i++ ){
                array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)] = recvptr[_IDX_V3D(i-(0-gc),j-(0-gc),k-(NK),l,gc,gc,gc,0)];
              }
            }
          }
        }
        break;

      case int(C_pXmYpZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=NK; k<NK+gc; k++ ){
            for( int j=0-gc; j<0; j++ ){
              for( int i=NI; i<NI+gc; i++ ){
                array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)] = recvptr[_IDX_V3D(i-(NI),j-(0-gc),k-(NK),l,gc,gc,gc,0)];
              }
            }
          }
        }
        break;

      case int(C_mXpYpZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=NK; k<NK+gc; k++ ){
            for( int j=NJ; j<NJ+gc; j++ ){
              for( int i=0-gc; i<0; i++ ){
                array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)] = recvptr[_IDX_V3D(i-(0-gc),j-(NJ),k-(NK),l,gc,gc,gc,0)];
              }
            }
          }
        }
        break;

      case int(C_pXpYpZ):
#pragma omp parallel for collapse(4)
        for( int l=0; l<3; l++ ){
          for( int k=NK; k<NK+gc; k++ ){
            for( int j=NJ; j<NJ+gc; j++ ){
              for( int i=NI; i<NI+gc; i++ ){
                array[_IDX_V3D(i,j,k,l,NI,NJ,NK,VC)] = recvptr[_IDX_V3D(i-(NI),j-(NJ),k-(NK),l,gc,gc,gc,0)];
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
