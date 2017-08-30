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
bool SubDomain::initComm(const int num_compo)
{
  int gc = halo_width;

  if (size[0]==0 || size[1]==0 || size[2]==0 || gc==0 || num_compo==0) {
    return false;
  }

  // バッファ領域としては、最大値で確保しておく
  int f_sz[3];
  f_sz[0] = (size[1]+2*gc) * (size[2]+2*gc) * gc * num_compo;
  f_sz[1] = (size[0]+2*gc) * (size[2]+2*gc) * gc * num_compo;
  f_sz[2] = (size[0]+2*gc) * (size[1]+2*gc) * gc * num_compo;

  if ( !(f_xms = new REAL_TYPE [f_sz[0]]) ) return false;
  if ( !(f_xmr = new REAL_TYPE [f_sz[0]]) ) return false;
  if ( !(f_xps = new REAL_TYPE [f_sz[0]]) ) return false;
  if ( !(f_xpr = new REAL_TYPE [f_sz[0]]) ) return false;

  if ( !(f_yms = new REAL_TYPE [f_sz[1]]) ) return false;
  if ( !(f_ymr = new REAL_TYPE [f_sz[1]]) ) return false;
  if ( !(f_yps = new REAL_TYPE [f_sz[1]]) ) return false;
  if ( !(f_ypr = new REAL_TYPE [f_sz[1]]) ) return false;

  if ( !(f_zms = new REAL_TYPE [f_sz[2]]) ) return false;
  if ( !(f_zmr = new REAL_TYPE [f_sz[2]]) ) return false;
  if ( !(f_zps = new REAL_TYPE [f_sz[2]]) ) return false;
  if ( !(f_zpr = new REAL_TYPE [f_sz[2]]) ) return false;

  buf_flag = 1; // バッファ確保ずみ

  return true;
}


/*
 * @brief スカラー変数のノンブロッキング通信
 * @param [in,out]  src     スカラー変数
 * @param [in]      gc_comm 実際に通信する通信面数
 * @param [in,out]  req     MPI_Request
 * @retval true-success, false-fail
 */
bool SubDomain::Comm_S_nonblocking(REAL_TYPE* src,
                                   const int gc_comm,
                                   MPI_Request *req)
{
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

  pack_SX(src, gc_comm, f_xms, f_xps, nIDm, nIDp);
  if ( !IsendIrecv(f_xms, f_xmr, f_xps, f_xpr, msz[0], nIDm, nIDp, &req[0]) ) return false;

  // Y direction
  nIDm = comm_tbl[Y_minus];
  nIDp = comm_tbl[Y_plus];

  pack_SY(src, gc_comm, f_yms, f_yps, nIDm, nIDp);
  if ( !IsendIrecv(f_yms, f_ymr, f_yps, f_ypr, msz[1], nIDm, nIDp, &req[4]) ) return false;

  // Z direction
  nIDm = comm_tbl[Z_minus];
  nIDp = comm_tbl[Z_plus];

  pack_SZ(src, gc_comm, f_zms, f_zps, nIDm, nIDp);
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
  unpack_SX(dest, gc_comm, f_xmr, f_xpr, nIDm, nIDp);


  //// Y face ////
  nIDm = comm_tbl[Y_minus];
  nIDp = comm_tbl[Y_plus];
  if ( MPI_SUCCESS != MPI_Waitall( 4, &req[4], stat ) ) return false;
  unpack_SY(dest, gc_comm, f_ymr, f_ypr, nIDm, nIDp);


  //// Z face ////
  nIDm = comm_tbl[Z_minus];
  nIDp = comm_tbl[Z_plus];
  if ( MPI_SUCCESS != MPI_Waitall( 4, &req[8], stat ) ) return false;
  unpack_SZ(dest, gc_comm, f_zmr, f_zpr, nIDm, nIDp);

  return true;
}
