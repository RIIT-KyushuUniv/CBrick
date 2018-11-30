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
void BrickComm::pack_SXnode(const REAL_TYPE *array,
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

/*
                 <--gc-->
rankA  [NI-3] [NI-2] [NI-1]  [NI]  [NI+1]
     -----+------+------|------+------+-------> i
rankB   [-2]   [-1]    [0]    [1]    [2]
                               <--gc-->
*/


  // 自領域のデータをマイナス側のランクに送る
  if( nIDm >= 0 )
  {
#pragma omp parallel for collapse(2)
    for( int k=0; k<NK; k++ ){
      for( int j=0; j<NJ; j++ ){
        for( int i=0; i<gc; i++ ){
          sendm[_IDX_SI(i,j,k,NJ,gc)] = array[_IDX_S3D(i+1,j,k,NI,NJ,VC)];
        }
      }
    }
  }

  // 自領域のデータをプラス側のランクに送る
  if( nIDp >= 0 )
  {
#pragma omp parallel for collapse(2)
    for( int k=0; k<NK; k++ ){
      for( int j=0; j<NJ; j++ ){
        for( int i=0; i<gc; i++ ){
          sendp[_IDX_SI(i,j,k,NJ,gc)] = array[_IDX_S3D(NI-2+i,j,k,NI,NJ,VC)];
        }
      }
    }
  }
}


/*
 * @brief unpack send data for I direction
 * @param [in,out]  array   dest array
 * @param [in]  gc number of guide cell layer to be sent
 * @param [in]  recvm   recv buffer of I- direction
 * @param [in]  recvp   recv buffer of I+ direction
 * @param [in]  nIDm    Rank number of I- direction
 * @param [in]  nIDp    Rank number of I+ direction
 */
void BrickComm::unpack_SXnode(REAL_TYPE *array,
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

/*
                 <--gc-->
rankA  [NI-3] [NI-2] [NI-1]  [NI]  [NI+1]
     -----+------+------|------+------+-------> i
rankB   [-2]   [-1]    [0]    [1]    [2]
                               <--gc-->
*/

  // マイナス側からのデータを自領域のガイドセルにコピー
  if( nIDm >= 0 )
  {
#pragma omp parallel for collapse(2)
    for( int k=0; k<NK; k++ ){
      for( int j=0; j<NJ; j++ ){
        for( int i=0; i<gc; i++ ){
          array[_IDX_S3D(i-1,j,k,NI,NJ,VC)] = recvm[_IDX_SI(i,j,k,NJ,gc)];
        }
      }
    }
  }

  // プラス側からのデータを自領域のガイドセルにコピー
  if( nIDp >= 0 )
  {
#pragma omp parallel for collapse(2)
    for( int k=0; k<NK; k++ ){
      for( int j=0; j<NJ; j++ ){
        for( int i=0; i<gc; i++ ){
          array[_IDX_S3D(NI+i,j,k,NI,NJ,VC)] = recvp[_IDX_SI(i,j,k,NJ,gc)];
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
void BrickComm::pack_SYnode(const REAL_TYPE *array,
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
#pragma omp parallel for collapse(2)
    for( int k=0; k<NK; k++ ){
      for( int j=1; j<=gc; j++ ){
        for( int i=0; i<NI; i++ ){
          sendm[_IDX_SJ(i,j,k,NI,1,gc)] = array[_IDX_S3D(i,j,k,NI,NJ,VC)];
        }
      }
    }
  }

  if( nIDp >= 0 )
  {
#pragma omp parallel for collapse(2)
    for( int k=0; k<NK; k++ ){
      for( int j=NJ-gc; j<NJ; j++ ){
        for( int i=0; i<NI; i++ ){
          sendp[_IDX_SJ(i,j,k,NI,NJ-gc,gc)] = array[_IDX_S3D(i,j,k,NI,NJ,VC)];
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
void BrickComm::unpack_SYnode(REAL_TYPE *array,
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
#pragma omp parallel for collapse(2)
    for( int k=0; k<NK; k++ ){
      for( int j=1-gc; j<1; j++ ){
        for( int i=0; i<NI; i++ ){
          array[_IDX_S3D(i,j,k,NI,NJ,VC)] = recvm[_IDX_SJ(i,j,k,NI,1-gc,gc)];
        }
      }
    }
  }

  if( nIDp >= 0 )
  {
#pragma omp parallel for collapse(2)
    for( int k=0; k<NK; k++ ){
      for( int j=NJ; j<NJ+gc; j++ ){
        for( int i=0; i<NI; i++ ){
          array[_IDX_S3D(i,j,k,NI,NJ,VC)] = recvp[_IDX_SJ(i,j,k,NI,NJ,gc)];
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
void BrickComm::pack_SZnode(const REAL_TYPE *array,
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
#pragma omp parallel for collapse(2)
    for( int k=1; k<=gc; k++ ){
      for( int j=0; j<NJ; j++ ){
        for( int i=0; i<NI; i++ ){
          sendm[_IDX_SK(i,j,k,NI,NJ,1,gc)] = array[_IDX_S3D(i,j,k,NI,NJ,VC)];
        }
      }
    }
  }

  if( nIDp >= 0 )
  {
#pragma omp parallel for collapse(2)
    for( int k=NK-gc; k<NK; k++ ){
      for( int j=0; j<NJ; j++ ){
        for( int i=0; i<NI; i++ ){
          sendp[_IDX_SK(i,j,k,NI,NJ,NK-gc,gc)] = array[_IDX_S3D(i,j,k,NI,NJ,VC)];
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
void BrickComm::unpack_SZnode(REAL_TYPE *array,
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
#pragma omp parallel for collapse(2)
    for( int k=1-gc; k<1; k++ ){
      for( int j=0; j<NJ; j++ ){
        for( int i=0; i<NI; i++ ){
          array[_IDX_S3D(i,j,k,NI,NJ,VC)] = recvm[_IDX_SK(i,j,k,NI,NJ,1-gc,gc)];
        }
      }
    }
  }

  if( nIDp >= 0 )
  {
#pragma omp parallel for collapse(2)
    for( int k=NK; k<NK+gc; k++ ){
      for( int j=0; j<NJ; j++ ){
        for( int i=0; i<NI; i++ ){
          array[_IDX_S3D(i,j,k,NI,NJ,VC)] = recvp[_IDX_SK(i,j,k,NI,NJ,NK,gc)];
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
bool BrickComm::pack_SEnode(REAL_TYPE *array,
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
      size_t sz = (NI-1) * gc * gc;

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
        for( int k=1; k<=gc; k++ ){
          for( int j=1; j<=gc; j++ ){
            for( int i=1; i<NI; i++ ){
              sendptr[_IDX_S3D(i-1,j-1,k-1,NI-1,gc,0)] = array[_IDX_S3D(i,j,k,NI,NJ,VC)];
            }
          }
        }
        break;

      case int(E_pYmZ):
#pragma omp parallel for collapse(3)
        for( int k=1; k<=gc; k++ ){
          for( int j=NJ-gc; j<NJ; j++ ){
            for( int i=1; i<NI; i++ ){
              sendptr[_IDX_S3D(i-1,j-(NJ-gc),k-1,NI-1,gc,0)] = array[_IDX_S3D(i,j,k,NI,NJ,VC)];
            }
          }
        }
        break;

      case int(E_mYpZ):
#pragma omp parallel for collapse(3)
        for( int k=NK-gc; k<NK; k++ ){
          for( int j=1; j<=gc; j++ ){
            for( int i=1; i<NI; i++ ){
              sendptr[_IDX_S3D(i-1,j-1,k-(NK-gc),NI-1,gc,0)] = array[_IDX_S3D(i,j,k,NI,NJ,VC)];
            }
          }
        }
        break;
      case int(E_pYpZ):

#pragma omp parallel for collapse(3)
        for( int k=NK-gc; k<NK; k++ ){
          for( int j=NJ-gc; j<NJ; j++ ){
            for( int i=1; i<NI; i++ ){
              sendptr[_IDX_S3D(i-1,j-(NJ-gc),k-(NK-gc),NI-1,gc,0)] = array[_IDX_S3D(i,j,k,NI,NJ,VC)];
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
      size_t sz = gc * (NJ-1) * gc;

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
        for( int k=1; k<=gc; k++ ){
          for( int j=1; j<NJ; j++ ){
            for( int i=1; i<=gc; i++ ){
              sendptr[_IDX_S3D(i-1,j-1,k-1,gc,NJ-1,0)] = array[_IDX_S3D(i,j,k,NI,NJ,VC)];
            }
          }
        }
        break;

      case int(E_pXmZ):
#pragma omp parallel for collapse(3)
        for( int k=1; k<=gc; k++ ){
          for( int j=1; j<NJ; j++ ){
            for( int i=NI-gc; i<NI; i++ ){
              sendptr[_IDX_S3D(i-(NI-gc),j-1,k-1,gc,NJ-1,0)] = array[_IDX_S3D(i,j,k,NI,NJ,VC)];
            }
          }
        }
        break;

      case int(E_mXpZ):
#pragma omp parallel for collapse(3)
        for( int k=NK-gc; k<NK; k++ ){
          for( int j=1; j<NJ; j++ ){
            for( int i=1; i<=gc; i++ ){
              sendptr[_IDX_S3D(i-1,j-1,k-(NK-gc),gc,NJ-1,0)] = array[_IDX_S3D(i,j,k,NI,NJ,VC)];
            }
          }
        }
        break;
      case int(E_pXpZ):

#pragma omp parallel for collapse(3)
        for( int k=NK-gc; k<NK; k++ ){
          for( int j=1; j<NJ; j++ ){
            for( int i=NI-gc; i<NI; i++ ){
              sendptr[_IDX_S3D(i-(NI-gc),j-1,k-(NK-gc),gc,NJ-1,0)] = array[_IDX_S3D(i,j,k,NI,NJ,VC)];
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
      size_t sz = gc * gc * (NK-1);

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
        for( int k=1; k<NK; k++ ){
          for( int j=1; j<=gc; j++ ){
            for( int i=1; i<=gc; i++ ){
              sendptr[_IDX_S3D(i-1,j-1,k-1,gc,gc,0)] = array[_IDX_S3D(i,j,k,NI,NJ,VC)];
            }
          }
        }
        break;

      case int(E_pXmY):
#pragma omp parallel for collapse(3)
        for( int k=1; k<NK; k++ ){
          for( int j=1; j<=gc; j++ ){
            for( int i=NI-gc; i<NI; i++ ){
              sendptr[_IDX_S3D(i-(NI-gc),j-1,k-1,gc,gc,0)] = array[_IDX_S3D(i,j,k,NI,NJ,VC)];
            }
          }
        }
        break;

      case int(E_mXpY):
#pragma omp parallel for collapse(3)
        for( int k=1; k<NK; k++ ){
          for( int j=NJ-gc; j<NJ; j++ ){
            for( int i=1; i<=gc; i++ ){
              sendptr[_IDX_S3D(i-1,j-(NJ-gc),k-1,gc,gc,0)] = array[_IDX_S3D(i,j,k,NI,NJ,VC)];
            }
          }
        }
        break;
      case int(E_pXpY):

#pragma omp parallel for collapse(3)
        for( int k=1; k<NK; k++ ){
          for( int j=NJ-gc; j<NJ; j++ ){
            for( int i=NI-gc; i<NI; i++ ){
              sendptr[_IDX_S3D(i-(NI-gc),j-(NJ-gc),k-1,gc,gc,0)] = array[_IDX_S3D(i,j,k,NI,NJ,VC)];
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
void BrickComm::unpack_SEnode(REAL_TYPE *array,
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
      size_t sz = (NI-1) * gc * gc;

      // unpack
      switch(dir)
      {
      case int(E_mYmZ):
#pragma omp parallel for collapse(3)
        for( int k=1-gc; k<=0; k++ ){
          for( int j=1-gc; j<=0; j++ ){
            for( int i=1; i<NI; i++ ){
              array[_IDX_S3D(i,j,k,NI,NJ,VC)] = recvptr[_IDX_S3D(i-1,j-(1-gc),k-(1-gc),NI-1,gc,0)];
            }
          }
        }
        break;

      case int(E_pYmZ):
#pragma omp parallel for collapse(3)
        for( int k=1-gc; k<=0; k++ ){
          for( int j=NJ; j<NJ+gc; j++ ){
            for( int i=1; i<NI; i++ ){
              array[_IDX_S3D(i,j,k,NI,NJ,VC)] = recvptr[_IDX_S3D(i-1,j-(NJ),k-(1-gc),NI-1,gc,0)];
            }
          }
        }
        break;

      case int(E_mYpZ):
#pragma omp parallel for collapse(3)
        for( int k=NK; k<NK+gc; k++ ){
          for( int j=1-gc; j<=0; j++ ){
            for( int i=1; i<NI; i++ ){
              array[_IDX_S3D(i,j,k,NI,NJ,VC)] = recvptr[_IDX_S3D(i-1,j-(1-gc),k-(NK),NI-1,gc,0)];
            }
          }
        }
        break;

      case int(E_pYpZ):
#pragma omp parallel for collapse(3)
        for( int k=NK; k<NK+gc; k++ ){
          for( int j=NJ; j<NJ+gc; j++ ){
            for( int i=1; i<NI; i++ ){
              array[_IDX_S3D(i,j,k,NI,NJ,VC)] = recvptr[_IDX_S3D(i-1,j-(NJ),k-(NK),NI-1,gc,0)];
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
      size_t sz = gc * (NJ-1) * gc;

      // unpack
      switch(dir)
      {
      case int(E_mXmZ):
#pragma omp parallel for collapse(3)
        for( int k=1-gc; k<=0; k++ ){
          for( int j=1; j<NJ; j++ ){
            for( int i=1-gc; i<=0; i++ ){
              array[_IDX_S3D(i,j,k,NI,NJ,VC)] = recvptr[_IDX_S3D(i-(1-gc),j-1,k-(1-gc),gc,NJ-1,0)];
            }
          }
        }
        break;

      case int(E_pXmZ):
#pragma omp parallel for collapse(3)
        for( int k=1-gc; k<=0; k++ ){
          for( int j=1; j<NJ; j++ ){
            for( int i=NI; i<NI+gc; i++ ){
              array[_IDX_S3D(i,j,k,NI,NJ,VC)] = recvptr[_IDX_S3D(i-(NI),j-1,k-(1-gc),gc,NJ-1,0)];
            }
          }
        }
        break;

      case int(E_mXpZ):
#pragma omp parallel for collapse(3)
        for( int k=NK; k<NK+gc; k++ ){
          for( int j=1; j<NJ; j++ ){
            for( int i=1-gc; i<=0; i++ ){
              array[_IDX_S3D(i,j,k,NI,NJ,VC)] = recvptr[_IDX_S3D(i-(1-gc),j-1,k-(NK),gc,NJ-1,0)];
            }
          }
        }
        break;

      case int(E_pXpZ):
#pragma omp parallel for collapse(3)
        for( int k=NK; k<NK+gc; k++ ){
          for( int j=1; j<NJ; j++ ){
            for( int i=NI; i<NI+gc; i++ ){
              array[_IDX_S3D(i,j,k,NI,NJ,VC)] = recvptr[_IDX_S3D(i-(NI),j-1,k-(NK),gc,NJ-1,0)];
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
      size_t sz = gc * gc * (NK-1);

      // unpack
      switch(dir)
      {
      case int(E_mXmY):
#pragma omp parallel for collapse(3)
        for( int k=1; k<NK; k++ ){
          for( int j=1-gc; j<=0; j++ ){
            for( int i=1-gc; i<=0; i++ ){
              array[_IDX_S3D(i,j,k,NI,NJ,VC)] = recvptr[_IDX_S3D(i-(1-gc),j-(1-gc),k-1,gc,gc,0)];
            }
          }
        }
        break;

      case int(E_pXmY):
#pragma omp parallel for collapse(3)
        for( int k=1; k<NK; k++ ){
          for( int j=1-gc; j<=0; j++ ){
            for( int i=NI; i<NI+gc; i++ ){
              array[_IDX_S3D(i,j,k,NI,NJ,VC)] = recvptr[_IDX_S3D(i-(NI),j-(1-gc),k-1,gc,gc,0)];
            }
          }
        }
        break;

      case int(E_mXpY):
#pragma omp parallel for collapse(3)
        for( int k=1; k<NK; k++ ){
          for( int j=NJ; j<NJ+gc; j++ ){
            for( int i=1-gc; i<=0; i++ ){
              array[_IDX_S3D(i,j,k,NI,NJ,VC)] = recvptr[_IDX_S3D(i-(1-gc),j-(NJ),k-1,gc,gc,0)];
            }
          }
        }
        break;

      case int(E_pXpY):
#pragma omp parallel for collapse(3)
        for( int k=1; k<NK; k++ ){
          for( int j=NJ; j<NJ+gc; j++ ){
            for( int i=NI; i<NI+gc; i++ ){
              array[_IDX_S3D(i,j,k,NI,NJ,VC)] = recvptr[_IDX_S3D(i-(NI),j-(NJ),k-1,gc,gc,0)];
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
bool BrickComm::pack_SCnode(REAL_TYPE *array,
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
      size_t sz = gc * gc * gc;

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
        for( int k=1; k<=gc; k++ ){
          for( int j=1; j<=gc; j++ ){
            for( int i=1; i<=gc; i++ ){
              sendptr[_IDX_S3D(i-1,j-1,k-1,gc,gc,0)] = array[_IDX_S3D(i,j,k,NI,NJ,VC)];
            }
          }
        }
        break;

      case int(C_pXmYmZ):
#pragma omp parallel for collapse(3)
        for( int k=1; k<=gc; k++ ){
          for( int j=1; j<=gc; j++ ){
            for( int i=NI-gc; i<NI; i++ ){
              sendptr[_IDX_S3D(i-(NI-gc),j-1,k-1,gc,gc,0)] = array[_IDX_S3D(i,j,k,NI,NJ,VC)];
            }
          }
        }
        break;

      case int(C_mXpYmZ):
#pragma omp parallel for collapse(3)
        for( int k=1; k<=gc; k++ ){
          for( int j=NJ-gc; j<NJ; j++ ){
            for( int i=1; i<=gc; i++ ){
              sendptr[_IDX_S3D(i-1,j-(NJ-gc),k-1,gc,gc,0)] = array[_IDX_S3D(i,j,k,NI,NJ,VC)];
            }
          }
        }
        break;

      case int(C_pXpYmZ):
#pragma omp parallel for collapse(3)
        for( int k=1; k<=gc; k++ ){
          for( int j=NJ-gc; j<NJ; j++ ){
            for( int i=NI-gc; i<NI; i++ ){
              sendptr[_IDX_S3D(i-(NI-gc),j-(NJ-gc),k-1,gc,gc,0)] = array[_IDX_S3D(i,j,k,NI,NJ,VC)];
            }
          }
        }
        break;

      case int(C_mXmYpZ):
#pragma omp parallel for collapse(3)
        for( int k=NK-gc; k<NK; k++ ){
          for( int j=1; j<=gc; j++ ){
            for( int i=1; i<=gc; i++ ){
              sendptr[_IDX_S3D(i-1,j-1,k-(NK-gc),gc,gc,0)] = array[_IDX_S3D(i,j,k,NI,NJ,VC)];
            }
          }
        }
        break;

      case int(C_pXmYpZ):
#pragma omp parallel for collapse(3)
        for( int k=NK-gc; k<NK; k++ ){
          for( int j=1; j<=gc; j++ ){
            for( int i=NI-gc; i<NI; i++ ){
              sendptr[_IDX_S3D(i-(NI-gc),j-1,k-(NK-gc),gc,gc,0)] = array[_IDX_S3D(i,j,k,NI,NJ,VC)];
            }
          }
        }
        break;

      case int(C_mXpYpZ):
#pragma omp parallel for collapse(3)
        for( int k=NK-gc; k<NK; k++ ){
          for( int j=NJ-gc; j<NJ; j++ ){
            for( int i=1; i<=gc; i++ ){
              sendptr[_IDX_S3D(i-1,j-(NJ-gc),k-(NK-gc),gc,gc,0)] = array[_IDX_S3D(i,j,k,NI,NJ,VC)];
            }
          }
        }
        break;

      case int(C_pXpYpZ):
#pragma omp parallel for collapse(3)
        for( int k=NK-gc; k<NK; k++ ){
          for( int j=NJ-gc; j<NJ; j++ ){
            for( int i=NI-gc; i<NI; i++ ){
              sendptr[_IDX_S3D(i-(NI-gc),j-(NJ-gc),k-(NK-gc),gc,gc,0)] = array[_IDX_S3D(i,j,k,NI,NJ,VC)];
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
void BrickComm::unpack_SCnode(REAL_TYPE *array,
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
      size_t sz = gc * gc * gc;

      // unpack
      switch(dir)
      {
      case int(C_mXmYmZ):
#pragma omp parallel for collapse(3)
        for( int k=1-gc; k<=0; k++ ){
          for( int j=1-gc; j<=0; j++ ){
            for( int i=1-gc; i<=0; i++ ){
              array[_IDX_S3D(i,j,k,NI,NJ,VC)] = recvptr[_IDX_S3D(i-(1-gc),j-(1-gc),k-(1-gc),gc,gc,0)];
            }
          }
        }
        break;

      case int(C_pXmYmZ):
#pragma omp parallel for collapse(3)
        for( int k=1-gc; k<=0; k++ ){
          for( int j=1-gc; j<=0; j++ ){
            for( int i=NI; i<NI+gc; i++ ){
              array[_IDX_S3D(i,j,k,NI,NJ,VC)] = recvptr[_IDX_S3D(i-(NI),j-(1-gc),k-(1-gc),gc,gc,0)];
            }
          }
        }
        break;

      case int(C_mXpYmZ):
#pragma omp parallel for collapse(3)
        for( int k=1-gc; k<=0; k++ ){
          for( int j=NJ; j<NJ+gc; j++ ){
            for( int i=1-gc; i<=0; i++ ){
              array[_IDX_S3D(i,j,k,NI,NJ,VC)] = recvptr[_IDX_S3D(i-(1-gc),j-(NJ),k-(1-gc),gc,gc,0)];
            }
          }
        }
        break;

      case int(C_pXpYmZ):
#pragma omp parallel for collapse(3)
        for( int k=1-gc; k<=0; k++ ){
          for( int j=NJ; j<NJ+gc; j++ ){
            for( int i=NI; i<NI+gc; i++ ){
              array[_IDX_S3D(i,j,k,NI,NJ,VC)] = recvptr[_IDX_S3D(i-(NI),j-(NJ),k-(1-gc),gc,gc,0)];
            }
          }
        }
        break;

      case int(C_mXmYpZ):
#pragma omp parallel for collapse(3)
        for( int k=NK; k<NK+gc; k++ ){
          for( int j=1-gc; j<=0; j++ ){
            for( int i=1-gc; i<=0; i++ ){
              array[_IDX_S3D(i,j,k,NI,NJ,VC)] = recvptr[_IDX_S3D(i-(1-gc),j-(1-gc),k-(NK),gc,gc,0)];
            }
          }
        }
        break;

      case int(C_pXmYpZ):
#pragma omp parallel for collapse(3)
        for( int k=NK; k<NK+gc; k++ ){
          for( int j=1-gc; j<=0; j++ ){
            for( int i=NI; i<NI+gc; i++ ){
              array[_IDX_S3D(i,j,k,NI,NJ,VC)] = recvptr[_IDX_S3D(i-(NI),j-(1-gc),k-(NK),gc,gc,0)];
            }
          }
        }
        break;

      case int(C_mXpYpZ):
#pragma omp parallel for collapse(3)
        for( int k=NK; k<NK+gc; k++ ){
          for( int j=NJ; j<NJ+gc; j++ ){
            for( int i=1-gc; i<=0; i++ ){
              array[_IDX_S3D(i,j,k,NI,NJ,VC)] = recvptr[_IDX_S3D(i-(1-gc),j-(NJ),k-(NK),gc,gc,0)];
            }
          }
        }
        break;

      case int(C_pXpYpZ):
#pragma omp parallel for collapse(3)
        for( int k=NK; k<NK+gc; k++ ){
          for( int j=NJ; j<NJ+gc; j++ ){
            for( int i=NI; i<NI+gc; i++ ){
              array[_IDX_S3D(i,j,k,NI,NJ,VC)] = recvptr[_IDX_S3D(i-(NI),j-(NJ),k-(NK),gc,gc,0)];
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
