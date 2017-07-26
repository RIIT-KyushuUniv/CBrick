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
 * @fn Comm_V_blocking
 * @brief ベクトル変数（i,j,k,l）型のブロッキング通信
 * @param [in,out]  src     ベクトル変数
 * @param [in]      gc_comm 実際に通信する通信面数
 * @retval true-success, false-fail
 */
bool SubDomain::Comm_V_blocking(REAL_TYPE* src, const int gc_comm)
{
  int imax = size[0];
  int jmax = size[1];
  int kmax = size[2];
  int gc = halo_width;

  // 各成分の先頭アドレス
  size_t p_u = 0;
  size_t p_v = (imax+2*gc) * (jmax+2*gc) * (kmax+2*gc);
  size_t p_w = (imax+2*gc) * (jmax+2*gc) * (kmax+2*gc) * 2;

  // 実際に送受信するメッセージサイズ
  int msz[3];
  msz[0] = (size[1]+2*gc_comm) * (size[2]+2*gc_comm) * gc_comm;
  msz[1] = (size[0]+2*gc_comm) * (size[2]+2*gc_comm) * gc_comm;
  msz[2] = (size[0]+2*gc_comm) * (size[1]+2*gc_comm) * gc_comm;

  // X direction
  int nIDm = comm_tbl[X_minus];
  int nIDp = comm_tbl[X_plus];

  packX(&src[p_u], gc_comm, f_xms, f_xps, nIDm, nIDp);
  if ( !send_and_recv(f_xms, f_xmr, f_xps, f_xpr, msz[0], nIDm, nIDp) ) return false;
  //if ( !sendrecv(f_xms, f_xmr, f_xps, f_xpr, msz[0], nIDm, nIDp) ) return false;
  unpackX(&src[p_u], gc_comm, f_xmr, f_xpr, nIDm, nIDp);

  packX(&src[p_v], gc_comm, f_xms, f_xps, nIDm, nIDp);
  if ( !send_and_recv(f_xms, f_xmr, f_xps, f_xpr, msz[0], nIDm, nIDp) ) return false;
  //if ( !sendrecv(f_xms, f_xmr, f_xps, f_xpr, msz[0], nIDm, nIDp) ) return false;
  unpackX(&src[p_v], gc_comm, f_xmr, f_xpr, nIDm, nIDp);

  packX(&src[p_w], gc_comm, f_xms, f_xps, nIDm, nIDp);
  if ( !send_and_recv(f_xms, f_xmr, f_xps, f_xpr, msz[0], nIDm, nIDp) ) return false;
  //if ( !sendrecv(f_xms, f_xmr, f_xps, f_xpr, msz[0], nIDm, nIDp) ) return false;
  unpackX(&src[p_w], gc_comm, f_xmr, f_xpr, nIDm, nIDp);


  // Y direction
  nIDm = comm_tbl[Y_minus];
  nIDp = comm_tbl[Y_plus];

  packY(&src[p_u], gc_comm, f_yms, f_yps, nIDm, nIDp);
  if ( !send_and_recv(f_yms, f_ymr, f_yps, f_ypr, msz[1], nIDm, nIDp) ) return false;
  //if ( !sendrecv(f_yms, f_ymr, f_yps, f_ypr, msz[1], nIDm, nIDp) ) return false;
  unpackY(&src[p_u], gc_comm, f_ymr, f_ypr, nIDm, nIDp);

  packY(&src[p_v], gc_comm, f_yms, f_yps, nIDm, nIDp);
  if ( !send_and_recv(f_yms, f_ymr, f_yps, f_ypr, msz[1], nIDm, nIDp) ) return false;
  //if ( !sendrecv(f_yms, f_ymr, f_yps, f_ypr, msz[1], nIDm, nIDp) ) return false;
  unpackY(&src[p_v], gc_comm, f_ymr, f_ypr, nIDm, nIDp);

  packY(&src[p_w], gc_comm, f_yms, f_yps, nIDm, nIDp);
  if ( !send_and_recv(f_yms, f_ymr, f_yps, f_ypr, msz[1], nIDm, nIDp) ) return false;
  //if ( !sendrecv(f_yms, f_ymr, f_yps, f_ypr, msz[1], nIDm, nIDp) ) return false;
  unpackY(&src[p_w], gc_comm, f_ymr, f_ypr, nIDm, nIDp);


  // Z direction
  nIDm = comm_tbl[Z_minus];
  nIDp = comm_tbl[Z_plus];

  packZ(&src[p_u], gc_comm, f_zms, f_zps, nIDm, nIDp);
  if ( !send_and_recv(f_zms, f_zmr, f_zps, f_zpr, msz[2], nIDm, nIDp) ) return false;
  //if ( !sendrecv(f_zms, f_zmr, f_zps, f_zpr, msz[2], nIDm, nIDp) ) return false;
  unpackZ(&src[p_u], gc_comm, f_zmr, f_zpr, nIDm, nIDp);

  packZ(&src[p_v], gc_comm, f_zms, f_zps, nIDm, nIDp);
  if ( !send_and_recv(f_zms, f_zmr, f_zps, f_zpr, msz[2], nIDm, nIDp) ) return false;
  //if ( !sendrecv(f_zms, f_zmr, f_zps, f_zpr, msz[2], nIDm, nIDp) ) return false;
  unpackZ(&src[p_v], gc_comm, f_zmr, f_zpr, nIDm, nIDp);

  packZ(&src[p_w], gc_comm, f_zms, f_zps, nIDm, nIDp);
  if ( !send_and_recv(f_zms, f_zmr, f_zps, f_zpr, msz[2], nIDm, nIDp) ) return false;
  //if ( !sendrecv(f_zms, f_zmr, f_zps, f_zpr, msz[2], nIDm, nIDp) ) return false;
  unpackZ(&src[p_w], gc_comm, f_zmr, f_zpr, nIDm, nIDp);

  return true;
}
