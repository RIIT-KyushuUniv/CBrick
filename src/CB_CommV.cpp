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
  int imax = size[0];
  int jmax = size[1];
  int kmax = size[2];

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

  pack_VX(src, gc_comm, f_xms, f_xps, nIDm, nIDp);
  if ( !IsendIrecv(f_xms, f_xmr, f_xps, f_xpr, msz[0], nIDm, nIDp, &req[0]) ) return false;

  // Y direction
  nIDm = comm_tbl[Y_minus];
  nIDp = comm_tbl[Y_plus];

  pack_VY(src, gc_comm, f_yms, f_yps, nIDm, nIDp);
  if ( !IsendIrecv(f_yms, f_ymr, f_yps, f_ypr, msz[1], nIDm, nIDp, &req[4]) ) return false;

  // Z direction
  nIDm = comm_tbl[Z_minus];
  nIDp = comm_tbl[Z_plus];

  pack_VZ(src, gc_comm, f_zms, f_zps, nIDm, nIDp);
  if ( !IsendIrecv(f_zms, f_zmr, f_zps, f_zpr, msz[2], nIDm, nIDp, &req[8]) ) return false;

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
  int nIDm = comm_tbl[X_minus];
  int nIDp = comm_tbl[X_plus];
  if ( MPI_SUCCESS != MPI_Waitall( 4, &req[0], stat ) ) return false;
  unpack_VX(dest, gc_comm, f_xmr, f_xpr, nIDm, nIDp);


  //// Y face ////
  nIDm = comm_tbl[Y_minus];
  nIDp = comm_tbl[Y_plus];
  if ( MPI_SUCCESS != MPI_Waitall( 4, &req[4], stat ) ) return false;
  unpack_VY(dest, gc_comm, f_ymr, f_ypr, nIDm, nIDp);


  //// Z face ////
  nIDm = comm_tbl[Z_minus];
  nIDp = comm_tbl[Z_plus];
  if ( MPI_SUCCESS != MPI_Waitall( 4, &req[8], stat ) ) return false;
  unpack_VZ(dest, gc_comm, f_zmr, f_zpr, nIDm, nIDp);

  return true;
}
