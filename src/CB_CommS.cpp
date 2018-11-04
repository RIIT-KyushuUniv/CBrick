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
 * @file   CB_CommS.cpp
 * @brief  BrickComm class
 */

#include "CB_Comm.h"


/*
 * @brief スカラー変数のノンブロッキング通信
 * @param [in,out]  src     スカラー変数
 * @param [in]      gc_comm 実際に通信する通信面数
 * @param [in,out]  req     MPI_Request
 * @retval true-success, false-fail
 */
bool BrickComm::Comm_S_nonblocking(REAL_TYPE* src,
                                   const int gc_comm,
                                   MPI_Request *req)
{
  // Communication identifier
  for (int i=0; i<NOFACE*2; i++) req[i] = MPI_REQUEST_NULL;

  // 実際に送受信するメッセージサイズ
  int msz[3];
  msz[0] = size[1] * size[2] * gc_comm;
  msz[1] = size[0] * size[2] * gc_comm;
  msz[2] = size[0] * size[1] * gc_comm;

  // X direction
  int nIDm = comm_tbl[I_minus];
  int nIDp = comm_tbl[I_plus];

  if (grid_type == "node")
  {
    pack_SXnode(src, gc_comm, f_ims, f_ips, nIDm, nIDp);
  }
  else
  {
    pack_SXcell(src, gc_comm, f_ims, f_ips, nIDm, nIDp);
  }

  if ( !IsendIrecv(f_ims, f_imr, f_ips, f_ipr, msz[0], nIDm, nIDp, &req[0]) ) return false;

  // Y direction
  nIDm = comm_tbl[J_minus];
  nIDp = comm_tbl[J_plus];

  if (grid_type == "node")
  {
    pack_SYnode(src, gc_comm, f_jms, f_jps, nIDm, nIDp);
  }
  else
  {
    pack_SYcell(src, gc_comm, f_jms, f_jps, nIDm, nIDp);
  }

  if ( !IsendIrecv(f_jms, f_jmr, f_jps, f_jpr, msz[1], nIDm, nIDp, &req[4]) ) return false;

  // Z direction
  nIDm = comm_tbl[K_minus];
  nIDp = comm_tbl[K_plus];

  if (grid_type == "node")
  {
    pack_SZnode(src, gc_comm, f_kms, f_kps, nIDm, nIDp);
  }
  else
  {
    pack_SZcell(src, gc_comm, f_kms, f_kps, nIDm, nIDp);
  }

  if ( !IsendIrecv(f_kms, f_kmr, f_kps, f_kpr, msz[2], nIDm, nIDp, &req[8]) ) return false;

#ifdef _DIAGONAL_COMM
  // edge
  if (grid_type == "node")
  {
    if( !pack_SEnode(src, gc_comm, f_es, f_er, req) ) return false;
  }
  else
  {
    if( !pack_SEcell(src, gc_comm, f_es, f_er, req) ) return false;
  }

  // corner
  if (grid_type == "node")
  {
    if( !pack_SCnode(src, gc_comm, f_cs, f_cr, req) ) return false;
  }
  else
  {
    if( !pack_SCcell(src, gc_comm, f_cs, f_cr, req) ) return false;
  }
#endif

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
bool BrickComm::IsendIrecv(REAL_TYPE* ms,
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
bool BrickComm::Comm_S_wait_nonblocking(REAL_TYPE* dest,
                                        const int gc_comm,
                                        MPI_Request *req)
{
#ifndef _DIAGONAL_COMM
  MPI_Status stat[4];
#else
  MPI_Status stat[26];
#endif

  //// X face ////
  int nIDm = comm_tbl[I_minus];
  int nIDp = comm_tbl[I_plus];
  if ( MPI_SUCCESS != MPI_Waitall( 4, &req[0], stat ) ) return false;
  if (grid_type == "node")
  {
    unpack_SXnode(dest, gc_comm, f_imr, f_ipr, nIDm, nIDp);
  }
  else
  {
    unpack_SXcell(dest, gc_comm, f_imr, f_ipr, nIDm, nIDp);
  }


  //// Y face ////
  nIDm = comm_tbl[J_minus];
  nIDp = comm_tbl[J_plus];
  if ( MPI_SUCCESS != MPI_Waitall( 4, &req[4], stat ) ) return false;
  if (grid_type == "node")
  {
    unpack_SYnode(dest, gc_comm, f_jmr, f_jpr, nIDm, nIDp);
  }
  else
  {
    unpack_SYcell(dest, gc_comm, f_jmr, f_jpr, nIDm, nIDp);
  }


  //// Z face ////
  nIDm = comm_tbl[K_minus];
  nIDp = comm_tbl[K_plus];
  if ( MPI_SUCCESS != MPI_Waitall( 4, &req[8], stat ) ) return false;
  if (grid_type == "node")
  {
    unpack_SZnode(dest, gc_comm, f_kmr, f_kpr, nIDm, nIDp);
  }
  else
  {
    unpack_SZcell(dest, gc_comm, f_kmr, f_kpr, nIDm, nIDp);
  }

#ifdef _DIAGONAL_COMM
  //// edge ////
  if ( MPI_SUCCESS != MPI_Waitall( 24, &req[12], stat ) ) return false;
  if (grid_type == "node")
  {
    unpack_SEnode(dest, gc_comm, f_er);
  }
  else
  {
    unpack_SEcell(dest, gc_comm, f_er);
  }

  //// corner ////
  if ( MPI_SUCCESS != MPI_Waitall( 16, &req[36], stat ) ) return false;
  if (grid_type == "node")
  {
    unpack_SCnode(dest, gc_comm, f_cr);
  }
  else
  {
    unpack_SCcell(dest, gc_comm, f_cr);
  }
#endif

  return true;
}
