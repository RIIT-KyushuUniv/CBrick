/*
###################################################################################
#
# CBrick
#
# Copyright (c) 2017-2019 Research Institute for Information Technology(RIIT),
#                    Kyushu University.  All rights reserved.
#
####################################################################################
*/

/*
 * @file   CB_Comm.cpp
 * @brief  BrickComm class
 */

#include "CB_Comm.h"


// #############################################################
// IsendIrecv Interface の実体
bool BrickComm::IsendIrecv(MPI_Datatype dtype,
                           void* ms,
                           void* mr,
                           void* ps,
                           void* pr,
                           int msz,
                           int nIDm,
                           int nIDp,
                           MPI_Request* req)
{
  // Identifier
  MPI_Request r0 = MPI_REQUEST_NULL;
  MPI_Request r1 = MPI_REQUEST_NULL;
  MPI_Request r2 = MPI_REQUEST_NULL;
  MPI_Request r3 = MPI_REQUEST_NULL;
  
  int tag_p=0;
  int tag_m=0;
  
  // Recieve Minus
  if ( nIDm >= 0 )
  {
    if ( MPI_SUCCESS != MPI_Irecv(mr,
                                  msz,
                                  dtype,
                                  nIDm,
                                  tag_m,
                                  MPI_COMM_WORLD,
                                  &r1) ) return false;
  }
  
  // Recieve Plus
  if ( nIDp >= 0 )
  {
    if ( MPI_SUCCESS != MPI_Irecv(pr,
                                  msz,
                                  dtype,
                                  nIDp,
                                  tag_p,
                                  MPI_COMM_WORLD,
                                  &r3) ) return false;
  }
  
  // Send Plus
  if ( nIDp >= 0 )
  {
    if ( MPI_SUCCESS != MPI_Isend(ps,
                                  msz,
                                  dtype,
                                  nIDp,
                                  tag_p,
                                  MPI_COMM_WORLD,
                                  &r2) ) return false;
  }
  
  // Send Minus
  if ( nIDm >= 0 )
  {
    if ( MPI_SUCCESS != MPI_Isend(ms,
                                  msz,
                                  dtype,
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


// #############################################################
// IrecvData Interface の実体
bool BrickComm::IrecvData(MPI_Datatype dtype,
                          void* ptr,
                          int sz,
                          int nID,
                          MPI_Request* req)
{  
  int tag = 0;
  
  if ( MPI_SUCCESS != MPI_Irecv(ptr,
                                sz,
                                dtype,
                                nID,
                                tag,
                                MPI_COMM_WORLD,
                                req) ) return false;
  return true;
}


// #############################################################
// IsendData Interface の実体
bool BrickComm::IsendData(MPI_Datatype dtype,
                          void* ptr,
                          int sz,
                          int nID,
                          MPI_Request* req)
{
  int tag = 0;
  
  if ( MPI_SUCCESS != MPI_Irecv(ptr,
                                sz,
                                dtype,
                                nID,
                                tag,
                                MPI_COMM_WORLD,
                                req) ) return false;
  return true;
}


// #############################################################
template
bool BrickComm::Comm_S_node(float* src, const int gc_comm, MPI_Request *req);

template
bool BrickComm::Comm_S_node(double* src, const int gc_comm, MPI_Request *req);

template
bool BrickComm::Comm_S_node(int* src, const int gc_comm, MPI_Request *req);

template
bool BrickComm::Comm_S_node(unsigned* src, const int gc_comm, MPI_Request *req);

template
bool BrickComm::Comm_S_node(long long* src, const int gc_comm, MPI_Request *req);


/* #########################################################
 * @brief スカラー変数 node
 * @param [in,out]  src     スカラー変数
 * @param [in]      gc_comm 実際に通信する通信面数
 * @param [in,out]  req     MPI_Request
 * @retval true-success, false-fail
 */
template <class T>
bool BrickComm::Comm_S_node(T* src,
                            const int gc_comm,
                            MPI_Request *req)
{
  T* b_ims = (T*)f_ims;  // I- direction send
  T* b_imr = (T*)f_imr;  // I- direction recv
  T* b_ips = (T*)f_ips;  // I+ direction send
  T* b_ipr = (T*)f_ipr;  // I+ direction recv
  T* b_jms = (T*)f_jms;  // J- direction send
  T* b_jmr = (T*)f_jmr;  // J- direction recv
  T* b_jps = (T*)f_jps;  // J+ direction send
  T* b_jpr = (T*)f_jpr;  // J+ direction recv
  T* b_kms = (T*)f_kms;  // K- direction send
  T* b_kmr = (T*)f_kmr;  // K- direction recv
  T* b_kps = (T*)f_kps;  // K+ direction send
  T* b_kpr = (T*)f_kpr;  // K+ direction recv
#ifdef _DIAGONAL_COMM
  T* b_es = (T*)f_es;   // edge send
  T* b_er = (T*)f_er;   // edge recv
  T* b_cs = (T*)f_cs;   // corner send
  T* b_cr = (T*)f_cr;   // corner recv
#endif
  
  
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
  pack_SXnode(src, gc_comm, b_ims, b_ips, nIDm, nIDp);
  if ( !IsendIrecv(b_ims, b_imr, b_ips, b_ipr, msz[0], nIDm, nIDp, &req[0]) ) return false;
  
  
  // Y direction
  nIDm = comm_tbl[J_minus];
  nIDp = comm_tbl[J_plus];
  pack_SYnode(src, gc_comm, b_jms, b_jps, nIDm, nIDp);
  if ( !IsendIrecv(b_jms, b_jmr, b_jps, b_jpr, msz[1], nIDm, nIDp, &req[4]) ) return false;
  
  
  // Z direction
  nIDm = comm_tbl[K_minus];
  nIDp = comm_tbl[K_plus];
  pack_SZnode(src, gc_comm, b_kms, b_kps, nIDm, nIDp);
  if ( !IsendIrecv(b_kms, b_kmr, b_kps, b_kpr, msz[2], nIDm, nIDp, &req[8]) ) return false;
  
  
#ifdef _DIAGONAL_COMM
  // edge
  if( !pack_SEnode(src, gc_comm, b_es, b_er, req) ) return false;
  
  // corner
  if( !pack_SCnode(src, gc_comm, b_cs, b_cr, req) ) return false;
#endif
  
  return true;
}



// #############################################################
template
bool BrickComm::Comm_S_cell(float* src, const int gc_comm, MPI_Request *req);

template
bool BrickComm::Comm_S_cell(double* src, const int gc_comm, MPI_Request *req);

template
bool BrickComm::Comm_S_cell(int* src, const int gc_comm, MPI_Request *req);

template
bool BrickComm::Comm_S_cell(unsigned* src, const int gc_comm, MPI_Request *req);

template
bool BrickComm::Comm_S_cell(long long* src, const int gc_comm, MPI_Request *req);


/* #########################################################
 * @brief スカラー変数 cell
 * @param [in,out]  src     スカラー変数
 * @param [in]      gc_comm 実際に通信する通信面数
 * @param [in,out]  req     MPI_Request
 * @retval true-success, false-fail
 */
template <class T>
bool BrickComm::Comm_S_cell(T* src,
                            const int gc_comm,
                            MPI_Request *req)
{
  T* b_ims = (T*)f_ims;  // I- direction send
  T* b_imr = (T*)f_imr;  // I- direction recv
  T* b_ips = (T*)f_ips;  // I+ direction send
  T* b_ipr = (T*)f_ipr;  // I+ direction recv
  T* b_jms = (T*)f_jms;  // J- direction send
  T* b_jmr = (T*)f_jmr;  // J- direction recv
  T* b_jps = (T*)f_jps;  // J+ direction send
  T* b_jpr = (T*)f_jpr;  // J+ direction recv
  T* b_kms = (T*)f_kms;  // K- direction send
  T* b_kmr = (T*)f_kmr;  // K- direction recv
  T* b_kps = (T*)f_kps;  // K+ direction send
  T* b_kpr = (T*)f_kpr;  // K+ direction recv
#ifdef _DIAGONAL_COMM
  T* b_es = (T*)f_es;   // edge send
  T* b_er = (T*)f_er;   // edge recv
  T* b_cs = (T*)f_cs;   // corner send
  T* b_cr = (T*)f_cr;   // corner recv
#endif
  
  
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
  pack_SXcell(src, gc_comm, b_ims, b_ips, nIDm, nIDp);
  if ( !IsendIrecv(b_ims, b_imr, b_ips, b_ipr, msz[0], nIDm, nIDp, &req[0]) ) return false;
  
  
  // Y direction
  nIDm = comm_tbl[J_minus];
  nIDp = comm_tbl[J_plus];
  pack_SYcell(src, gc_comm, b_jms, b_jps, nIDm, nIDp);
  if ( !IsendIrecv(b_jms, b_jmr, b_jps, b_jpr, msz[1], nIDm, nIDp, &req[4]) ) return false;
  
  
  // Z direction
  nIDm = comm_tbl[K_minus];
  nIDp = comm_tbl[K_plus];
  pack_SZcell(src, gc_comm, b_kms, b_kps, nIDm, nIDp);
  if ( !IsendIrecv(b_kms, b_kmr, b_kps, b_kpr, msz[2], nIDm, nIDp, &req[8]) ) return false;
  
  
#ifdef _DIAGONAL_COMM
  // edge
  if( !pack_SEcell(src, gc_comm, b_es, b_er, req) ) return false;
  
  // corner
  if( !pack_SCcell(src, gc_comm, b_cs, b_cr, req) ) return false;
#endif
  
  return true;
}


// #########################################################
template
bool BrickComm::Comm_S_wait_node(float* dest, const int gc_comm, MPI_Request *req);

template
bool BrickComm::Comm_S_wait_node(double* dest, const int gc_comm, MPI_Request *req);

template
bool BrickComm::Comm_S_wait_node(int* dest, const int gc_comm, MPI_Request *req);

template
bool BrickComm::Comm_S_wait_node(unsigned* dest, const int gc_comm, MPI_Request *req);

template
bool BrickComm::Comm_S_wait_node(long long* dest, const int gc_comm, MPI_Request *req);


/* #########################################################
 * @brief スカラー変数 node
 * @param [in,out]  dest    スカラー変数
 * @param [in]      gc_comm 実際に通信する通信面数
 * @param [out]     req     Array of MPI request
 * @retval true-success, false-fail
 */
template <class T>
bool BrickComm::Comm_S_wait_node(T* dest,
                                 const int gc_comm,
                                 MPI_Request *req)
{
#ifndef _DIAGONAL_COMM
  MPI_Status stat[4];
#else
  MPI_Status stat[26];
#endif
  
  T* b_ims = (T*)f_ims;  // I- direction send
  T* b_imr = (T*)f_imr;  // I- direction recv
  T* b_ips = (T*)f_ips;  // I+ direction send
  T* b_ipr = (T*)f_ipr;  // I+ direction recv
  T* b_jms = (T*)f_jms;  // J- direction send
  T* b_jmr = (T*)f_jmr;  // J- direction recv
  T* b_jps = (T*)f_jps;  // J+ direction send
  T* b_jpr = (T*)f_jpr;  // J+ direction recv
  T* b_kms = (T*)f_kms;  // K- direction send
  T* b_kmr = (T*)f_kmr;  // K- direction recv
  T* b_kps = (T*)f_kps;  // K+ direction send
  T* b_kpr = (T*)f_kpr;  // K+ direction recv
#ifdef _DIAGONAL_COMM
  T* b_es = (T*)f_es;   // edge send
  T* b_er = (T*)f_er;   // edge recv
  T* b_cs = (T*)f_cs;   // corner send
  T* b_cr = (T*)f_cr;   // corner recv
#endif
  
  
  //// X face ////
  int nIDm = comm_tbl[I_minus];
  int nIDp = comm_tbl[I_plus];
  if ( MPI_SUCCESS != MPI_Waitall( 4, &req[0], stat ) ) return false;
  unpack_SXnode(dest, gc_comm, b_imr, b_ipr, nIDm, nIDp);
  
  
  //// Y face ////
  nIDm = comm_tbl[J_minus];
  nIDp = comm_tbl[J_plus];
  if ( MPI_SUCCESS != MPI_Waitall( 4, &req[4], stat ) ) return false;
  unpack_SYnode(dest, gc_comm, b_jmr, b_jpr, nIDm, nIDp);
  
  
  //// Z face ////
  nIDm = comm_tbl[K_minus];
  nIDp = comm_tbl[K_plus];
  if ( MPI_SUCCESS != MPI_Waitall( 4, &req[8], stat ) ) return false;
  unpack_SZnode(dest, gc_comm, b_kmr, b_kpr, nIDm, nIDp);
  
  
#ifdef _DIAGONAL_COMM
  //// edge ////
  if ( MPI_SUCCESS != MPI_Waitall( 24, &req[12], stat ) ) return false;
  unpack_SEnode(dest, gc_comm, b_er);
  
  //// corner ////
  if ( MPI_SUCCESS != MPI_Waitall( 16, &req[36], stat ) ) return false;
  unpack_SCnode(dest, gc_comm, b_cr);
#endif
  
  return true;
}


// #########################################################
template
bool BrickComm::Comm_S_wait_cell(float* dest, const int gc_comm, MPI_Request *req);

template
bool BrickComm::Comm_S_wait_cell(double* dest, const int gc_comm, MPI_Request *req);

template
bool BrickComm::Comm_S_wait_cell(int* dest, const int gc_comm, MPI_Request *req);

template
bool BrickComm::Comm_S_wait_cell(unsigned* dest, const int gc_comm, MPI_Request *req);

template
bool BrickComm::Comm_S_wait_cell(long long* dest, const int gc_comm, MPI_Request *req);


/* #########################################################
 * @brief スカラー変数 cell
 * @param [in,out]  dest    スカラー変数
 * @param [in]      gc_comm 実際に通信する通信面数
 * @param [out]     req     Array of MPI request
 * @retval true-success, false-fail
 */
template <class T>
bool BrickComm::Comm_S_wait_cell(T* dest,
                                 const int gc_comm,
                                 MPI_Request *req)
{
#ifndef _DIAGONAL_COMM
  MPI_Status stat[4];
#else
  MPI_Status stat[26];
#endif
  
  T* b_ims = (T*)f_ims;  // I- direction send
  T* b_imr = (T*)f_imr;  // I- direction recv
  T* b_ips = (T*)f_ips;  // I+ direction send
  T* b_ipr = (T*)f_ipr;  // I+ direction recv
  T* b_jms = (T*)f_jms;  // J- direction send
  T* b_jmr = (T*)f_jmr;  // J- direction recv
  T* b_jps = (T*)f_jps;  // J+ direction send
  T* b_jpr = (T*)f_jpr;  // J+ direction recv
  T* b_kms = (T*)f_kms;  // K- direction send
  T* b_kmr = (T*)f_kmr;  // K- direction recv
  T* b_kps = (T*)f_kps;  // K+ direction send
  T* b_kpr = (T*)f_kpr;  // K+ direction recv
#ifdef _DIAGONAL_COMM
  T* b_es = (T*)f_es;   // edge send
  T* b_er = (T*)f_er;   // edge recv
  T* b_cs = (T*)f_cs;   // corner send
  T* b_cr = (T*)f_cr;   // corner recv
#endif
  
  
  //// X face ////
  int nIDm = comm_tbl[I_minus];
  int nIDp = comm_tbl[I_plus];
  if ( MPI_SUCCESS != MPI_Waitall( 4, &req[0], stat ) ) return false;
  unpack_SXcell(dest, gc_comm, b_imr, b_ipr, nIDm, nIDp);
  
  
  //// Y face ////
  nIDm = comm_tbl[J_minus];
  nIDp = comm_tbl[J_plus];
  if ( MPI_SUCCESS != MPI_Waitall( 4, &req[4], stat ) ) return false;
  unpack_SYcell(dest, gc_comm, b_jmr, b_jpr, nIDm, nIDp);
  
  
  //// Z face ////
  nIDm = comm_tbl[K_minus];
  nIDp = comm_tbl[K_plus];
  if ( MPI_SUCCESS != MPI_Waitall( 4, &req[8], stat ) ) return false;
  unpack_SZcell(dest, gc_comm, b_kmr, b_kpr, nIDm, nIDp);
  
  
#ifdef _DIAGONAL_COMM
  //// edge ////
  if ( MPI_SUCCESS != MPI_Waitall( 24, &req[12], stat ) ) return false;
  unpack_SEcell(dest, gc_comm, b_er);
  
  //// corner ////
  if ( MPI_SUCCESS != MPI_Waitall( 16, &req[36], stat ) ) return false;
  unpack_SCcell(dest, gc_comm, b_cr);
#endif
  
  return true;
}


// #########################################################
template
bool BrickComm::Comm_V_node(float* src, const int gc_comm, MPI_Request *req);

template
bool BrickComm::Comm_V_node(double* src, const int gc_comm, MPI_Request *req);

template
bool BrickComm::Comm_V_node(int* src, const int gc_comm, MPI_Request *req);

 
/* #########################################################
 * @brief ベクトル変数のノンブロッキング通信
 * @param [in,out]  src     ベクトル変数
 * @param [in]      gc_comm 実際に通信する通信面数
 * @param [in,out]  req     MPI_Request
 * @retval true-success, false-fail
 */
template <class T>
bool BrickComm::Comm_V_node(T* src,
                            const int gc_comm,
                            MPI_Request *req)
{
  T* b_ims = (T*)f_ims;  // I- direction send
  T* b_imr = (T*)f_imr;  // I- direction recv
  T* b_ips = (T*)f_ips;  // I+ direction send
  T* b_ipr = (T*)f_ipr;  // I+ direction recv
  T* b_jms = (T*)f_jms;  // J- direction send
  T* b_jmr = (T*)f_jmr;  // J- direction recv
  T* b_jps = (T*)f_jps;  // J+ direction send
  T* b_jpr = (T*)f_jpr;  // J+ direction recv
  T* b_kms = (T*)f_kms;  // K- direction send
  T* b_kmr = (T*)f_kmr;  // K- direction recv
  T* b_kps = (T*)f_kps;  // K+ direction send
  T* b_kpr = (T*)f_kpr;  // K+ direction recv
#ifdef _DIAGONAL_COMM
  T* b_es = (T*)f_es;   // edge send
  T* b_er = (T*)f_er;   // edge recv
  T* b_cs = (T*)f_cs;   // corner send
  T* b_cr = (T*)f_cr;   // corner recv
#endif
  
  
  // Communication identifier
  for (int i=0; i<NOFACE*2; i++) req[i] = MPI_REQUEST_NULL;
  
  // 実際に送受信するメッセージサイズ
  int msz[3];
  msz[0] = size[1] * size[2] * gc_comm * 3;
  msz[1] = size[0] * size[2] * gc_comm * 3;
  msz[2] = size[0] * size[1] * gc_comm * 3;
  
  
  // X direction
  int nIDm = comm_tbl[I_minus];
  int nIDp = comm_tbl[I_plus];
  pack_VXnode(src, gc_comm, b_ims, b_ips, nIDm, nIDp);
  if ( !IsendIrecv(b_ims, b_imr, b_ips, b_ipr, msz[0], nIDm, nIDp, &req[0]) ) return false;
  
  
  // Y direction
  nIDm = comm_tbl[J_minus];
  nIDp = comm_tbl[J_plus];
  pack_VYnode(src, gc_comm, b_jms, b_jps, nIDm, nIDp);
  if ( !IsendIrecv(b_jms, b_jmr, b_jps, b_jpr, msz[1], nIDm, nIDp, &req[4]) ) return false;
  
  
  // Z direction
  nIDm = comm_tbl[K_minus];
  nIDp = comm_tbl[K_plus];
  pack_VZnode(src, gc_comm, b_kms, b_kps, nIDm, nIDp);
  if ( !IsendIrecv(b_kms, b_kmr, b_kps, b_kpr, msz[2], nIDm, nIDp, &req[8]) ) return false;
  
  
#ifdef _DIAGONAL_COMM
  // edge
  if( !pack_VEnode(src, gc_comm, b_es, b_er, req) ) return false;
  
  // corner
  if( !pack_VCnode(src, gc_comm, b_cs, b_cr, req) ) return false;
#endif
  
  return true;
}


// #########################################################
template
bool BrickComm::Comm_V_cell(float* src, const int gc_comm, MPI_Request *req);

template
bool BrickComm::Comm_V_cell(double* src, const int gc_comm, MPI_Request *req);

template
bool BrickComm::Comm_V_cell(int* src, const int gc_comm, MPI_Request *req);


/* #########################################################
 * @brief ベクトル変数のノンブロッキング通信
 * @param [in,out]  src     ベクトル変数
 * @param [in]      gc_comm 実際に通信する通信面数
 * @param [in,out]  req     MPI_Request
 * @retval true-success, false-fail
 */
template <class T>
bool BrickComm::Comm_V_cell(T* src,
                            const int gc_comm,
                            MPI_Request *req)
{
  T* b_ims = (T*)f_ims;  // I- direction send
  T* b_imr = (T*)f_imr;  // I- direction recv
  T* b_ips = (T*)f_ips;  // I+ direction send
  T* b_ipr = (T*)f_ipr;  // I+ direction recv
  T* b_jms = (T*)f_jms;  // J- direction send
  T* b_jmr = (T*)f_jmr;  // J- direction recv
  T* b_jps = (T*)f_jps;  // J+ direction send
  T* b_jpr = (T*)f_jpr;  // J+ direction recv
  T* b_kms = (T*)f_kms;  // K- direction send
  T* b_kmr = (T*)f_kmr;  // K- direction recv
  T* b_kps = (T*)f_kps;  // K+ direction send
  T* b_kpr = (T*)f_kpr;  // K+ direction recv
#ifdef _DIAGONAL_COMM
  T* b_es = (T*)f_es;   // edge send
  T* b_er = (T*)f_er;   // edge recv
  T* b_cs = (T*)f_cs;   // corner send
  T* b_cr = (T*)f_cr;   // corner recv
#endif
  
  
  // Communication identifier
  for (int i=0; i<NOFACE*2; i++) req[i] = MPI_REQUEST_NULL;
  
  // 実際に送受信するメッセージサイズ
  int msz[3];
  msz[0] = size[1] * size[2] * gc_comm * 3;
  msz[1] = size[0] * size[2] * gc_comm * 3;
  msz[2] = size[0] * size[1] * gc_comm * 3;
  
  
  // X direction
  int nIDm = comm_tbl[I_minus];
  int nIDp = comm_tbl[I_plus];
  pack_VXcell(src, gc_comm, b_ims, b_ips, nIDm, nIDp);
  if ( !IsendIrecv(b_ims, b_imr, b_ips, b_ipr, msz[0], nIDm, nIDp, &req[0]) ) return false;
  
  
  // Y direction
  nIDm = comm_tbl[J_minus];
  nIDp = comm_tbl[J_plus];
  pack_VYcell(src, gc_comm, b_jms, b_jps, nIDm, nIDp);
  if ( !IsendIrecv(b_jms, b_jmr, b_jps, b_jpr, msz[1], nIDm, nIDp, &req[4]) ) return false;
  
  
  // Z direction
  nIDm = comm_tbl[K_minus];
  nIDp = comm_tbl[K_plus];
  pack_VZcell(src, gc_comm, b_kms, b_kps, nIDm, nIDp);
  if ( !IsendIrecv(b_kms, b_kmr, b_kps, b_kpr, msz[2], nIDm, nIDp, &req[8]) ) return false;
  
  
#ifdef _DIAGONAL_COMM
  // edge
  if( !pack_VEcell(src, gc_comm, b_es, b_er, req) ) return false;
  
  // corner
  if( !pack_VCcell(src, gc_comm, b_cs, b_cr, req) ) return false;
#endif
  
  return true;
}


// #########################################################
template
bool BrickComm::Comm_V_wait_node(float* dest, const int gc_comm, MPI_Request *req);

template
bool BrickComm::Comm_V_wait_node(double* dest, const int gc_comm, MPI_Request *req);

template
bool BrickComm::Comm_V_wait_node(int* dest, const int gc_comm, MPI_Request *req);


/* #########################################################
 * @brief ベクトル変数 node
 * @param [in,out]  dest    ベクトル変数
 * @param [in]      gc_comm 実際に通信する通信面数
 * @param [out]     req     Array of MPI request
 * @retval true-success, false-fail
 */
template <class T>
bool BrickComm::Comm_V_wait_node(T* dest,
                                 const int gc_comm,
                                 MPI_Request *req)
{
#ifndef _DIAGONAL_COMM
  MPI_Status stat[4];
#else
  MPI_Status stat[26];
#endif
  
  T* b_ims = (T*)f_ims;  // I- direction send
  T* b_imr = (T*)f_imr;  // I- direction recv
  T* b_ips = (T*)f_ips;  // I+ direction send
  T* b_ipr = (T*)f_ipr;  // I+ direction recv
  T* b_jms = (T*)f_jms;  // J- direction send
  T* b_jmr = (T*)f_jmr;  // J- direction recv
  T* b_jps = (T*)f_jps;  // J+ direction send
  T* b_jpr = (T*)f_jpr;  // J+ direction recv
  T* b_kms = (T*)f_kms;  // K- direction send
  T* b_kmr = (T*)f_kmr;  // K- direction recv
  T* b_kps = (T*)f_kps;  // K+ direction send
  T* b_kpr = (T*)f_kpr;  // K+ direction recv
#ifdef _DIAGONAL_COMM
  T* b_es = (T*)f_es;   // edge send
  T* b_er = (T*)f_er;   // edge recv
  T* b_cs = (T*)f_cs;   // corner send
  T* b_cr = (T*)f_cr;   // corner recv
#endif
  
  
  //// X face ////
  int nIDm = comm_tbl[I_minus];
  int nIDp = comm_tbl[I_plus];
  if ( MPI_SUCCESS != MPI_Waitall( 4, &req[0], stat ) ) return false;
  unpack_VXnode(dest, gc_comm, b_imr, b_ipr, nIDm, nIDp);
  
  
  //// Y face ////
  nIDm = comm_tbl[J_minus];
  nIDp = comm_tbl[J_plus];
  if ( MPI_SUCCESS != MPI_Waitall( 4, &req[4], stat ) ) return false;
  unpack_VYnode(dest, gc_comm, b_jmr, b_jpr, nIDm, nIDp);
  
  
  //// Z face ////
  nIDm = comm_tbl[K_minus];
  nIDp = comm_tbl[K_plus];
  if ( MPI_SUCCESS != MPI_Waitall( 4, &req[8], stat ) ) return false;
  unpack_VZnode(dest, gc_comm, b_kmr, b_kpr, nIDm, nIDp);
  
  
#ifdef _DIAGONAL_COMM
  //// edge ////
  if ( MPI_SUCCESS != MPI_Waitall( 24, &req[12], stat ) ) return false;
  unpack_VEnode(dest, gc_comm, b_er);
  
  //// corner ////
  if ( MPI_SUCCESS != MPI_Waitall( 16, &req[36], stat ) ) return false;
  unpack_VCnode(dest, gc_comm, b_cr);
#endif
  
  return true;
}


// #########################################################
template
bool BrickComm::Comm_V_wait_cell(float* dest, const int gc_comm, MPI_Request *req);

template
bool BrickComm::Comm_V_wait_cell(double* dest, const int gc_comm, MPI_Request *req);

template
bool BrickComm::Comm_V_wait_cell(int* dest, const int gc_comm, MPI_Request *req);


/* #########################################################
 * @brief ベクトル変数 node
 * @param [in,out]  dest    ベクトル変数
 * @param [in]      gc_comm 実際に通信する通信面数
 * @param [out]     req     Array of MPI request
 * @retval true-success, false-fail
 */
template <class T>
bool BrickComm::Comm_V_wait_cell(T* dest,
                                 const int gc_comm,
                                 MPI_Request *req)
{
#ifndef _DIAGONAL_COMM
  MPI_Status stat[4];
#else
  MPI_Status stat[26];
#endif
  
  T* b_ims = (T*)f_ims;  // I- direction send
  T* b_imr = (T*)f_imr;  // I- direction recv
  T* b_ips = (T*)f_ips;  // I+ direction send
  T* b_ipr = (T*)f_ipr;  // I+ direction recv
  T* b_jms = (T*)f_jms;  // J- direction send
  T* b_jmr = (T*)f_jmr;  // J- direction recv
  T* b_jps = (T*)f_jps;  // J+ direction send
  T* b_jpr = (T*)f_jpr;  // J+ direction recv
  T* b_kms = (T*)f_kms;  // K- direction send
  T* b_kmr = (T*)f_kmr;  // K- direction recv
  T* b_kps = (T*)f_kps;  // K+ direction send
  T* b_kpr = (T*)f_kpr;  // K+ direction recv
#ifdef _DIAGONAL_COMM
  T* b_es = (T*)f_es;   // edge send
  T* b_er = (T*)f_er;   // edge recv
  T* b_cs = (T*)f_cs;   // corner send
  T* b_cr = (T*)f_cr;   // corner recv
#endif
  
  
  //// X face ////
  int nIDm = comm_tbl[I_minus];
  int nIDp = comm_tbl[I_plus];
  if ( MPI_SUCCESS != MPI_Waitall( 4, &req[0], stat ) ) return false;
  unpack_VXcell(dest, gc_comm, b_imr, b_ipr, nIDm, nIDp);
  
  
  //// Y face ////
  nIDm = comm_tbl[J_minus];
  nIDp = comm_tbl[J_plus];
  if ( MPI_SUCCESS != MPI_Waitall( 4, &req[4], stat ) ) return false;
  unpack_VYcell(dest, gc_comm, b_jmr, b_jpr, nIDm, nIDp);
  
  
  //// Z face ////
  nIDm = comm_tbl[K_minus];
  nIDp = comm_tbl[K_plus];
  if ( MPI_SUCCESS != MPI_Waitall( 4, &req[8], stat ) ) return false;
  unpack_VZcell(dest, gc_comm, b_kmr, b_kpr, nIDm, nIDp);
  
  
#ifdef _DIAGONAL_COMM
  //// edge ////
  if ( MPI_SUCCESS != MPI_Waitall( 24, &req[12], stat ) ) return false;
  unpack_VEcell(dest, gc_comm, b_er);
  
  //// corner ////
  if ( MPI_SUCCESS != MPI_Waitall( 16, &req[36], stat ) ) return false;
  unpack_VCcell(dest, gc_comm, b_cr);
#endif
  
  return true;
}
