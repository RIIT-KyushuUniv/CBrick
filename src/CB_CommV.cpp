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
 * @file   CB_CommV.cpp
 * @brief  SubDomain class
 */

#include "CB_SubDomain.h"



/*
 * @brief ベクトル変数のノンブロッキング通信
 * @param [in,out]  src     ベクトル変数
 * @param [in]      gc_comm 実際に通信する通信面数
 * @param [in,out]  req     MPI_Request
 * @retval true-success, false-fail
 */
bool SubDomain::Comm_V_nonblocking(REAL_TYPE* src,
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
  int nIDm = comm_tbl[I_minus];
  int nIDp = comm_tbl[I_plus];

  pack_VI(src, gc_comm, f_ims, f_ips, nIDm, nIDp);
  if ( !IsendIrecv(f_ims, f_imr, f_ips, f_ipr, msz[0], nIDm, nIDp, &req[0]) ) return false;

  // Y direction
  nIDm = comm_tbl[J_minus];
  nIDp = comm_tbl[J_plus];

  pack_VJ(src, gc_comm, f_jms, f_jps, nIDm, nIDp);
  if ( !IsendIrecv(f_jms, f_jmr, f_jps, f_jpr, msz[1], nIDm, nIDp, &req[4]) ) return false;

  // Z direction
  nIDm = comm_tbl[K_minus];
  nIDp = comm_tbl[K_plus];

  pack_VK(src, gc_comm, f_kms, f_kps, nIDm, nIDp);
  if ( !IsendIrecv(f_kms, f_kmr, f_kps, f_kpr, msz[2], nIDm, nIDp, &req[8]) ) return false;

  return true;
}



/*
 * @brief ベクトル変数のノンブロッキング通信
 * @param [in,out]  dest    ベクトル変数
 * @param [in]      gc_comm 実際に通信する通信面数
 * @param [out]     req     Array of MPI request
 * @retval true-success, false-fail
 */
bool SubDomain::Comm_V_wait_nonblocking(REAL_TYPE* dest,
                                        const int gc_comm,
                                        MPI_Request *req)
{
  MPI_Status stat[4];

  //// X face ////
  int nIDm = comm_tbl[I_minus];
  int nIDp = comm_tbl[I_plus];
  if ( MPI_SUCCESS != MPI_Waitall( 4, &req[0], stat ) ) return false;
  unpack_VI(dest, gc_comm, f_imr, f_ipr, nIDm, nIDp);


  //// Y face ////
  nIDm = comm_tbl[J_minus];
  nIDp = comm_tbl[J_plus];
  if ( MPI_SUCCESS != MPI_Waitall( 4, &req[4], stat ) ) return false;
  unpack_VJ(dest, gc_comm, f_jmr, f_jpr, nIDm, nIDp);


  //// Z face ////
  nIDm = comm_tbl[K_minus];
  nIDp = comm_tbl[K_plus];
  if ( MPI_SUCCESS != MPI_Waitall( 4, &req[8], stat ) ) return false;
  unpack_VK(dest, gc_comm, f_kmr, f_kpr, nIDm, nIDp);

  return true;
}
