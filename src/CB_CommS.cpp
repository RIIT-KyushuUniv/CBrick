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
 * @file   CB_CommS.cpp
 * @brief  SubDomain class
 */

#include "CB_SubDomain.h"

/*
 * @fn initComm()
 * @brief 通信バッファの確保
 */
bool SubDomain::initComm()
{
  int gc = halo_width;

  if (size[0]==0 || size[1]==0 || size[2]==0 || gc==0) {
    return false;
  }

  // バッファ領域としては、最大値で確保しておく
  int f_sz[3];
  f_sz[0] = (size[1]+2*gc) * (size[2]+2*gc) * gc;
  f_sz[1] = (size[0]+2*gc) * (size[2]+2*gc) * gc;
  f_sz[2] = (size[0]+2*gc) * (size[1]+2*gc) * gc;

  if ( !(f_xms = new REAL_TYPE [f_sz[0]]) ) {
    return false;
  }
  if ( !(f_xmr = new REAL_TYPE [f_sz[0]]) ) {
    return false;
  }
  if ( !(f_xps = new REAL_TYPE [f_sz[0]]) ) {
    return false;
  }
  if ( !(f_xpr = new REAL_TYPE [f_sz[0]]) ) {
    return false;
  }

  if ( !(f_yms = new REAL_TYPE [f_sz[1]]) ) {
    return false;
  }
  if ( !(f_ymr = new REAL_TYPE [f_sz[1]]) ) {
    return false;
  }
  if ( !(f_yps = new REAL_TYPE [f_sz[1]]) ) {
    return false;
  }
  if ( !(f_ypr = new REAL_TYPE [f_sz[1]]) ) {
    return false;
  }

  if ( !(f_zms = new REAL_TYPE [f_sz[2]]) ) {
    return false;
  }
  if ( !(f_zmr = new REAL_TYPE [f_sz[2]]) ) {
    return false;
  }
  if ( !(f_zps = new REAL_TYPE [f_sz[2]]) ) {
    return false;
  }
  if ( !(f_zpr = new REAL_TYPE [f_sz[2]]) ) {
    return false;
  }

  return true;
}

/*
 * @fn Comm_S_blocking
 * @brief スカラー変数のブロッキング通信
 * @param [in,out]  src     スカラー変数
 * @param [in]      gc_comm 実際に通信する通信面数
 * @retval true-success, false-fail
 */
bool SubDomain::Comm_S_blocking(REAL_TYPE* src, const int gc_comm)
{
  int imax = size[0];
  int jmax = size[1];
  int kmax = size[2];
  int gc = halo_width;

  // 実際に送受信するメッセージサイズ
  int msz[3];
  msz[0] = (size[1]+2*gc_comm) * (size[2]+2*gc_comm) * gc_comm;
  msz[1] = (size[0]+2*gc_comm) * (size[2]+2*gc_comm) * gc_comm;
  msz[2] = (size[0]+2*gc_comm) * (size[1]+2*gc_comm) * gc_comm;

  // X direction
  int nIDm = comm_tbl[X_minus];
  int nIDp = comm_tbl[X_plus];

  packX(src, gc_comm, f_xms, f_xps, nIDm, nIDp);
  if ( !send_and_recv(f_xms, f_xmr, f_xps, f_xpr, msz[0], nIDm, nIDp) ) return false;
  //if ( !sendrecv(f_xms, f_xmr, f_xps, f_xpr, msz[0], nIDm, nIDp) ) return false;
  unpackX(src, gc_comm, f_xmr, f_xpr, nIDm, nIDp);


  // Y direction
  nIDm = comm_tbl[Y_minus];
  nIDp = comm_tbl[Y_plus];

  packY(src, gc_comm, f_yms, f_yps, nIDm, nIDp);
  if ( !send_and_recv(f_yms, f_ymr, f_yps, f_ypr, msz[1], nIDm, nIDp) ) return false;
  //if ( !sendrecv(f_yms, f_ymr, f_yps, f_ypr, msz[1], nIDm, nIDp) ) return false;
  unpackY(src, gc_comm, f_ymr, f_ypr, nIDm, nIDp);


  // Z direction
  nIDm = comm_tbl[Z_minus];
  nIDp = comm_tbl[Z_plus];

  packZ(src, gc_comm, f_zms, f_zps, nIDm, nIDp);
  if ( !send_and_recv(f_zms, f_zmr, f_zps, f_zpr, msz[2], nIDm, nIDp) ) return false;
  //if ( !sendrecv(f_zms, f_zmr, f_zps, f_zpr, msz[2], nIDm, nIDp) ) return false;
  unpackZ(src, gc_comm, f_zmr, f_zpr, nIDm, nIDp);

  return true;
}


/*
 * @fn send_and_recv
 * @brief 隣接間通信
 * @param [in]  ms   Send buffer to Minus direction
 * @param [out] mr   Recieve buffer from Mminus direction
 * @param [in]  ps   Send buffer to Plus direction
 * @param [out] pr   Recieve buffer from Plus direction
 * @param [in]  msz  send/recieve size
 * @param [in]  nIDm Neighbor ID for Minus direction
 * @param [in]  nIDp Neighbor ID for Plus direction
 */
bool SubDomain::send_and_recv(REAL_TYPE* ms,
                              REAL_TYPE* mr,
                              REAL_TYPE* ps,
                              REAL_TYPE* pr,
                              int msz,
                              int nIDm,
                              int nIDp)
{
  // Plus side of subdomain
  int tag_p=0;
  MPI_Status *stat_p;

  if ( sizeof(REAL_TYPE) == _SIZE_DOUBLE_ ) {
    if ( nIDp >= 0 )
      if ( MPI_SUCCESS != MPI_Send(ps,
                                   msz,
                                   MPI_DOUBLE,
                                   nIDp,
                                   tag_p,
                                   MPI_COMM_WORLD) ) return false;
    if ( nIDp >= 0 )
      if ( MPI_SUCCESS != MPI_Recv(pr,
                                   msz,
                                   MPI_DOUBLE,
                                   nIDp,
                                   tag_p,
                                   MPI_COMM_WORLD,
                                   stat_p) ) return false;
  }
  else {
    if ( nIDp >= 0 )
      if ( MPI_SUCCESS != MPI_Send(ps,
                                   msz,
                                   MPI_FLOAT,
                                   nIDp,
                                   tag_p,
                                   MPI_COMM_WORLD) ) return false;
    if ( nIDp >= 0 )
      if ( MPI_SUCCESS != MPI_Recv(pr,
                                   msz,
                                   MPI_FLOAT,
                                   nIDp,
                                   tag_p,
                                   MPI_COMM_WORLD,
                                   stat_p) ) return false;
  }


  // Minus side of subdomain
  int tag_m=0;
  MPI_Status *stat_m;

  if ( sizeof(REAL_TYPE) == _SIZE_DOUBLE_ ) {
    if ( nIDm >= 0 )
      if ( MPI_SUCCESS != MPI_Recv(mr,
                                   msz,
                                   MPI_DOUBLE,
                                   nIDm,
                                   tag_m,
                                   MPI_COMM_WORLD,
                                   stat_m) ) return false;
    if ( nIDm >= 0 )
      if ( MPI_SUCCESS != MPI_Send(ms,
                                   msz,
                                   MPI_DOUBLE,
                                   nIDm,
                                   tag_m,
                                   MPI_COMM_WORLD) ) return false;
  }
  else {
    if ( nIDm >= 0 )
      if ( MPI_SUCCESS != MPI_Recv(mr,
                                   msz,
                                   MPI_FLOAT,
                                   nIDm,
                                   tag_m,
                                   MPI_COMM_WORLD,
                                   stat_m) ) return false;
    if ( nIDm >= 0 )
      if ( MPI_SUCCESS != MPI_Send(ms,
                                   msz,
                                   MPI_FLOAT,
                                   nIDm,
                                   tag_m,
                                   MPI_COMM_WORLD) ) return false;
  }

  return true;
}


/*
 * @fn sendrecv
 * @brief 隣接間通信
 * @param [in]  ms   Send buffer to Minus direction
 * @param [out] mr   Recieve buffer from Mminus direction
 * @param [in]  ps   Send buffer to Plus direction
 * @param [out] pr   Recieve buffer from Plus direction
 * @param [in]  msz  send/recieve size
 * @param [in]  nIDm Neighbor ID for Minus direction
 * @param [in]  nIDp Neighbor ID for Plus direction
 */
bool SubDomain::sendrecv(REAL_TYPE* ms,
                         REAL_TYPE* mr,
                         REAL_TYPE* ps,
                         REAL_TYPE* pr,
                         int msz,
                         int nIDm,
                         int nIDp)
{
  // Plus side of subdomain
  int tag_p=0;
  MPI_Status *stat_p;

  if ( sizeof(REAL_TYPE) == _SIZE_DOUBLE_ ) {
    if ( nIDp >= 0 ) {
      if ( MPI_SUCCESS != MPI_Sendrecv(ps, msz, MPI_DOUBLE, nIDp, tag_p,
                                       pr, msz, MPI_DOUBLE, nIDp, tag_p,
                                       MPI_COMM_WORLD,
                                       stat_p) ) return false;
    }
  }
  else {
    if ( nIDp >= 0 ) {
      if ( MPI_SUCCESS != MPI_Sendrecv(ps, msz, MPI_FLOAT, nIDp, tag_p,
                                       pr, msz, MPI_FLOAT, nIDp, tag_p,
                                       MPI_COMM_WORLD,
                                       stat_p) ) return false;
    }
  }


  // Minus side of subdomain
  int tag_m=0;
  MPI_Status *stat_m;

  if ( sizeof(REAL_TYPE) == _SIZE_DOUBLE_ ) {
    if ( nIDm >= 0 ) {
      if ( MPI_SUCCESS != MPI_Sendrecv(ms, msz, MPI_DOUBLE, nIDm, tag_m,
                                       mr, msz, MPI_DOUBLE, nIDm, tag_m,
                                       MPI_COMM_WORLD,
                                       stat_m) ) return false;
    }
  }
  else {
    if ( nIDm >= 0 ) {
      if ( MPI_SUCCESS != MPI_Sendrecv(ms, msz, MPI_FLOAT, nIDm, tag_m,
                                       mr, msz, MPI_FLOAT, nIDm, tag_m,
                                       MPI_COMM_WORLD,
                                       stat_m) ) return false;
      }
  }

  return true;
}


/*
 * @brief pack send data for X direction
 * @param [in]  array   source array
 * @param [in]  vc_comm number of guide cell layer to be sent
 * @param [out] sendm   send buffer of X- direction
 * @param [out] sendp   send buffer of X+ direction
 * @param [in]  nIDm    Rank number of X- direction
 * @param [in]  nIDp    Rank number of X+ direction
 */
void SubDomain::packX(const REAL_TYPE *array,
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
    for( int k=0-vc_comm; k<kmax+vc_comm; k++ ){
      for( int j=0-vc_comm; j<jmax+vc_comm; j++ ){
        for( int i=0; i<vc_comm; i++ ){
          sendm[_IDXFX(i,j,k,0,jmax,vc_comm)] = array[_IDX_S3D(i,j,k,imax,jmax,vc)];
        }
      }
    }
  }

  if( nIDp >= 0 )
  {
#pragma omp parallel for collapse(2)
    for( int k=0-vc_comm; k<kmax+vc_comm; k++ ){
      for( int j=0-vc_comm; j<jmax+vc_comm; j++ ){
        for( int i=imax-vc_comm; i<imax; i++ ){
          sendp[_IDXFX(i,j,k,imax-vc_comm,jmax,vc_comm)] = array[_IDX_S3D(i,j,k,imax,jmax,vc)];
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
void SubDomain::unpackX(REAL_TYPE *array,
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
    for( int k=0-vc_comm; k<kmax+vc_comm; k++ ){
      for( int j=0-vc_comm; j<jmax+vc_comm; j++ ){
        for( int i=0-vc_comm; i<0; i++ ){
          array[_IDX_S3D(i,j,k,imax,jmax,vc)] = recvm[_IDXFX(i,j,k,0-vc_comm,jmax,vc_comm)];
        }
      }
    }
  }

  if( nIDp >= 0 )
  {
#pragma omp parallel for collapse(2)
    for( int k=0-vc_comm; k<kmax+vc_comm; k++ ){
      for( int j=0-vc_comm; j<jmax+vc_comm; j++ ){
        for( int i=imax; i<imax+vc_comm; i++ ){
          array[_IDX_S3D(i,j,k,imax,jmax,vc)] = recvp[_IDXFX(i,j,k,imax,jmax,vc_comm)];
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
void SubDomain::packY(const REAL_TYPE *array,
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
    for( int k=0-vc_comm; k<kmax+vc_comm; k++ ){
      for( int j=0; j<vc_comm; j++ ){
        for( int i=0-vc_comm; i<imax+vc_comm; i++ ){
          sendm[_IDXFY(i,j,k,imax,0,vc_comm)] = array[_IDX_S3D(i,j,k,imax,jmax,vc)];
        }
      }
    }
  }

  if( nIDp >= 0 )
  {
#pragma omp parallel for collapse(2)
    for( int k=0-vc_comm; k<kmax+vc_comm; k++ ){
      for( int j=jmax-vc_comm; j<jmax; j++ ){
        for( int i=0-vc_comm; i<imax+vc_comm; i++ ){
          sendp[_IDXFY(i,j,k,imax,jmax-vc_comm,vc_comm)] = array[_IDX_S3D(i,j,k,imax,jmax,vc)];
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
void SubDomain::unpackY(REAL_TYPE *array,
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
    for( int k=0-vc_comm; k<kmax+vc_comm; k++ ){
      for( int j=0-vc_comm; j<0; j++ ){
        for( int i=0-vc_comm; i<imax+vc_comm; i++ ){
          array[_IDX_S3D(i,j,k,imax,jmax,vc)] = recvm[_IDXFY(i,j,k,imax,0-vc_comm,vc_comm)];
        }
      }
    }
  }

  if( nIDp >= 0 )
  {
#pragma omp parallel for collapse(2)
    for( int k=0-vc_comm; k<kmax+vc_comm; k++ ){
      for( int j=jmax; j<jmax+vc_comm; j++ ){
        for( int i=0-vc_comm; i<imax+vc_comm; i++ ){
          array[_IDX_S3D(i,j,k,imax,jmax,vc)] = recvp[_IDXFY(i,j,k,imax,jmax,vc_comm)];
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
void SubDomain::packZ(const REAL_TYPE *array,
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
    for( int k=0; k<vc_comm; k++ ){
      for( int j=0-vc_comm; j<jmax+vc_comm; j++ ){
        for( int i=0-vc_comm; i<imax+vc_comm; i++ ){
          sendm[_IDXFZ(i,j,k,imax,jmax,0,vc_comm)] = array[_IDX_S3D(i,j,k,imax,jmax,vc)];
        }
      }
    }
  }

  if( nIDp >= 0 )
  {
#pragma omp parallel for collapse(2)
    for( int k=kmax-vc_comm; k<kmax; k++ ){
      for( int j=0-vc_comm; j<jmax+vc_comm; j++ ){
        for( int i=0-vc_comm; i<imax+vc_comm; i++ ){
          sendp[_IDXFZ(i,j,k,imax,jmax,kmax-vc_comm,vc_comm)] = array[_IDX_S3D(i,j,k,imax,jmax,vc)];
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
void SubDomain::unpackZ(REAL_TYPE *array,
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
    for( int k=0-vc_comm; k<0; k++ ){
      for( int j=0-vc_comm; j<jmax+vc_comm; j++ ){
        for( int i=0-vc_comm; i<imax+vc_comm; i++ ){
          array[_IDX_S3D(i,j,k,imax,jmax,vc)] = recvm[_IDXFZ(i,j,k,imax,jmax,0-vc_comm,vc_comm)];
        }
      }
    }
  }

  if( nIDp >= 0 )
  {
#pragma omp parallel for collapse(2)
    for( int k=kmax; k<kmax+vc_comm; k++ ){
      for( int j=0-vc_comm; j<jmax+vc_comm; j++ ){
        for( int i=0-vc_comm; i<imax+vc_comm; i++ ){
          array[_IDX_S3D(i,j,k,imax,jmax,vc)] = recvp[_IDXFZ(i,j,k,imax,jmax,kmax,vc_comm)];
        }
      }
    }
  }
}


/*
 * @brief スカラー変数のノンブロッキング通信
 * @param [in,out]  src     スカラー変数
 * @param [in]      gc_comm 実際に通信する通信面数
 * @retval true-success, false-fail
 */
bool SubDomain::Comm_S_nonblocking(REAL_TYPE* src,
                                   const int gc_comm,
                                   MPI_Request *req)
{
  int imax = size[0];
  int jmax = size[1];
  int kmax = size[2];
  int gc = halo_width;

  // Communication identifier
  for (int i=0; i<12; i++) req[i] = MPI_REQUEST_NULL;

  // 実際に送受信するメッセージサイズ
  int msz[3];
  msz[0] = (size[1]+2*gc_comm) * (size[2]+2*gc_comm) * gc_comm;
  msz[1] = (size[0]+2*gc_comm) * (size[2]+2*gc_comm) * gc_comm;
  msz[2] = (size[0]+2*gc_comm) * (size[1]+2*gc_comm) * gc_comm;

  // X direction
  int nIDm = comm_tbl[X_minus];
  int nIDp = comm_tbl[X_plus];

  packX(src, gc_comm, f_xms, f_xps, nIDm, nIDp);
  if ( !IsendIrecv(f_xms, f_xmr, f_xps, f_xpr, msz[0], nIDm, nIDp, &req[0]) ) return false;

  // Y direction
  nIDm = comm_tbl[Y_minus];
  nIDp = comm_tbl[Y_plus];

  packY(src, gc_comm, f_yms, f_yps, nIDm, nIDp);
  if ( !IsendIrecv(f_yms, f_ymr, f_yps, f_ypr, msz[1], nIDm, nIDp, &req[4]) ) return false;

  // Z direction
  nIDm = comm_tbl[Z_minus];
  nIDp = comm_tbl[Z_plus];

  packZ(src, gc_comm, f_zms, f_zps, nIDm, nIDp);
  if ( !IsendIrecv(f_zms, f_zmr, f_zps, f_zpr, msz[2], nIDm, nIDp, &req[8]) ) return false;

  return true;
}


/*
 * @fn IsendIrecv
 * @brief 隣接間通信
 * @param [in]  ms   Send buffer to Minus direction
 * @param [out] mr   Recieve buffer from Mminus direction
 * @param [in]  ps   Send buffer to Plus direction
 * @param [out] pr   Recieve buffer from Plus direction
 * @param [in]  msz  send/recieve size
 * @param [in]  nIDm Neighbor ID for Minus direction
 * @param [in]  nIDp Neighbor ID for Plus direction
 * @param [out] req  Array of MPI request
 * @retval true-success, false-fail
 */
bool SubDomain::IsendIrecv(REAL_TYPE* ms,
                           REAL_TYPE* mr,
                           REAL_TYPE* ps,
                           REAL_TYPE* pr,
                           int msz,
                           int nIDm,
                           int nIDp,
                           MPI_Request *req)
{
  // Identifier
  MPI_Request r0 = MPI_REQUEST_NULL;
  MPI_Request r1 = MPI_REQUEST_NULL;
  MPI_Request r2 = MPI_REQUEST_NULL;
  MPI_Request r3 = MPI_REQUEST_NULL;

  int tag_p=0;
  int tag_m=0;

  // Recieve Minus
  if ( sizeof(REAL_TYPE) == _SIZE_DOUBLE_ ) {
    if ( nIDm >= 0 )
      if ( MPI_SUCCESS != MPI_Irecv(mr,
                                    msz,
                                    MPI_DOUBLE,
                                    nIDm,
                                    tag_m,
                                    MPI_COMM_WORLD,
                                    &r1) ) return false;
  }
  else {
    if ( nIDm >= 0 )
      if ( MPI_SUCCESS != MPI_Irecv(mr,
                                    msz,
                                    MPI_FLOAT,
                                    nIDm,
                                    tag_m,
                                    MPI_COMM_WORLD,
                                    &r1) ) return false;
  }

  // Recieve Plus
  if ( sizeof(REAL_TYPE) == _SIZE_DOUBLE_ ) {
    if ( nIDp >= 0 )
      if ( MPI_SUCCESS != MPI_Irecv(pr,
                                    msz,
                                    MPI_DOUBLE,
                                    nIDp,
                                    tag_p,
                                    MPI_COMM_WORLD,
                                    &r3) ) return false;
  }
  else {
    if ( nIDp >= 0 )
      if ( MPI_SUCCESS != MPI_Irecv(pr,
                                    msz,
                                    MPI_FLOAT,
                                    nIDp,
                                    tag_p,
                                    MPI_COMM_WORLD,
                                    &r3) ) return false;
  }

  // Send Plus
  if ( sizeof(REAL_TYPE) == _SIZE_DOUBLE_ ) {
    if ( nIDp >= 0 )
      if ( MPI_SUCCESS != MPI_Isend(ps,
                                    msz,
                                    MPI_DOUBLE,
                                    nIDp,
                                    tag_p,
                                    MPI_COMM_WORLD,
                                    &r2) ) return false;
  }
  else {
    if ( nIDp >= 0 )
      if ( MPI_SUCCESS != MPI_Isend(ps,
                                    msz,
                                    MPI_FLOAT,
                                    nIDp,
                                    tag_p,
                                    MPI_COMM_WORLD,
                                    &r2) ) return false;
  }

  // Send Minus
  if ( sizeof(REAL_TYPE) == _SIZE_DOUBLE_ ) {
    if ( nIDm >= 0 )
      if ( MPI_SUCCESS != MPI_Isend(ms,
                                    msz,
                                    MPI_DOUBLE,
                                    nIDm,
                                    tag_m,
                                    MPI_COMM_WORLD,
                                    &r0) ) return false;
  }
  else {
    if ( nIDm >= 0 )
      if ( MPI_SUCCESS != MPI_Isend(ms,
                                    msz,
                                    MPI_FLOAT,
                                    nIDm,
                                    tag_m,
                                    MPI_COMM_WORLD,
                                    &r0) ) return false;
  }

  req[0] = r0;
  req[1] = r1;
  req[2] = r2;
  req[3] = r3;

  return true;
}


/*
 * @brief スカラー変数のノンブロッキング通信
 * @param [in,out]  dest    スカラー変数
 * @param [in]      gc_comm 実際に通信する通信面数
 * @param [out]     req     Array of MPI request
 * @retval true-success, false-fail
 */
bool SubDomain::Comm_S_wait_nonblocking(REAL_TYPE* dest,
                                        const int gc_comm,
                                        MPI_Request *req)
{
  MPI_Status stat[4];

  //// X face ////
  int nIDm = comm_tbl[X_minus];
  int nIDp = comm_tbl[X_plus];
  if ( MPI_SUCCESS != MPI_Waitall( 4, &req[0], stat ) ) return false;
  unpackX(dest, gc_comm, f_xmr, f_xpr, nIDm, nIDp);


  //// Y face ////
  nIDm = comm_tbl[Y_minus];
  nIDp = comm_tbl[Y_plus];
  if ( MPI_SUCCESS != MPI_Waitall( 4, &req[4], stat ) ) return false;
  unpackY(dest, gc_comm, f_ymr, f_ypr, nIDm, nIDp);


  //// Z face ////
  nIDm = comm_tbl[Z_minus];
  nIDp = comm_tbl[Z_plus];
  if ( MPI_SUCCESS != MPI_Waitall( 4, &req[8], stat ) ) return false;
  unpackZ(dest, gc_comm, f_zmr, f_zpr, nIDm, nIDp);

  return true;
}
