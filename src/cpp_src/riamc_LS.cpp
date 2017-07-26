/*
###################################################################################
#
# RIAM-COMPACT
#
# Copyright (C) 2015-2017 Research Institute for Applied Mechanics(RIAM)
#                       / Research Institute for Information Technology(RIIT), Kyushu University.
# All rights reserved.
#
# Copyright (C) 2015-2016 Advanced Institute for Computational Science, RIKEN.
# All rights reserved.
#
###################################################################################
*/

/**
 * @file   riamc_LS.C
 * @brief  LS Class
 * @author aics
 */

#include "riamc.h"


// #################################################################
double RIAMC::Fdot1(REAL_TYPE* x)
{
  double flop_count=0.0;          /// 浮動小数点演算数
  double xy = 0.0;

  TIMING_start("Dot1");
  blas_dot1_(&xy, x, bcp, size, &flop_count);
  TIMING_stop("Dot1", flop_count);

  if ( numProc > 1 )
  {
    TIMING_start("A_R_Dot");
    double xy_tmp = xy;
    if  ( paraMngr->Allreduce(&xy_tmp, &xy, 1, MPI_SUM, procGrp) != CPM_SUCCESS ) Exit(0);
    TIMING_stop("A_R_Dot", 2.0*numProc*sizeof(double) );
  }

  return xy;
}


// #################################################################
double RIAMC::Fdot2(REAL_TYPE* x, REAL_TYPE* y)
{
  double flop_count=0.0;          /// 浮動小数点演算数
  double xy = 0.0;

  TIMING_start("Dot2");
  blas_dot2_(&xy, x, y, bcp, size, &flop_count);
  TIMING_stop("Dot2", flop_count);

  if ( numProc > 1 )
  {
    TIMING_start("A_R_Dot");
    double xy_tmp = xy;
    if  ( paraMngr->Allreduce(&xy_tmp, &xy, 1, MPI_SUM, procGrp) != CPM_SUCCESS ) Exit(0);
    TIMING_stop("A_R_Dot", 2.0*numProc*sizeof(double) );
  }

  return xy;
}


// #################################################################
void RIAMC::Preconditioner(REAL_TYPE* x, REAL_TYPE* b, const REAL_TYPE dt)
{

  double dummy = 1.0;

  cs = 0.0;

  // 前処理なし(コピー)
  if ( !isPreconditioned() )
  {
    TIMING_start("Blas_Copy");
    blas_copy_(x, b, size);
    TIMING_stop("Blas_Copy");
    return;
  }

  int lc_max = getInnerItr();

  // 前処理
  // 境界条件処理が実行される場合には、要注意
  if ( smoother == SOR2SMA )
  {
    SOR2_SMA(x, b, dt, lc_max, dummy, dummy, false);
  }
  else if ( smoother == SOR )
  {
    PointSOR(x, b, dt, lc_max, dummy, dummy, false);
  }

  //PointSSOR(x, b, dt, lc_max, dummy, dummy, false);
}



/* #################################################################
int RIAMC::SOR2_SMA()
{
  int ip;                         /// ローカルノードの基点(1,1,1)のカラーを示すインデクス
  /// ip=0 > R, ip=1 > B
  double flop_count=0.0;          /// 浮動小数点演算数
  REAL_TYPE omg = getOmega();     /// 加速係数
  double var[3];                  /// 誤差、残差、解
  int lc=0;                       /// ループカウント

  cs = 0.0;


  for (lc=1; lc<itrMax; lc++)
  {
    // 2色のマルチカラー(Red&Black)のセットアップ

    // ip = 0 基点(1,1,1)が Rからスタート
    //    = 1 基点(1,1,1)が Bからスタート
    if ( numProc > 1 )
    {
      ip = (head[0]+head[1]+head[2]+1) % 2;
    }
    else
    {
      ip = 0;
    }


    var[0] = 0.0; // 誤差
    var[1] = 0.0; // 残差
    var[2] = 0.0; // 解

    // 各カラー毎の間に同期, 残差は色間で積算する
    // R - color=0 / B - color=1
    for (int color=0; color<2; color++) {

      TIMING_start("Poisson_SOR2_SMA");
      flop_count = 0.0; // 色間で積算しない
      psor2sma_(x, size, pitch, &ip, &color, &omg, var, b, bcp, &cs, &flop_count);


      TIMING_stop("Poisson_SOR2_SMA", flop_count);


      // 境界条件
      TIMING_start("Poisson_BC");
      BC->OuterPBC(x, ensPeriodic);
      TIMING_stop("Poisson_BC", 0.0);


      // 同期処理
      if ( numProc > 1 )
      {
        TIMING_start("Sync_Poisson");
        int gd = GUIDE;
        if ( getSyncMode() == comm_sync )
        {
          if ( paraMngr->BndCommS3D(x, size[0], size[1], size[2], gd, 1, procGrp) != CPM_SUCCESS ) Exit(0); // 1 layer communication
        }
        else
        {
          int ireq[12];
          sma_comm_     (x, size, &gd, &color, &ip, cf_sz, cf_x, cf_y, cf_z, ireq, comm_tbl);
          sma_comm_wait_(x, size, &gd, &color, &ip, cf_sz, cf_x, cf_y, cf_z, ireq);
        }
        TIMING_stop("Sync_Poisson", face_comm_size*0.5*sizeof(REAL_TYPE));
      }
    }

    if ( converge_check )
    {
      // 収束判定 varは自乗量
      if ( Fcheck(var, b_l2, r0_l2) == true ) break;
    }

  }

  return lc;
}
*/


// #################################################################
// 反復変数の同期処理
void RIAMC::SyncScalar(REAL_TYPE* d_class, const int num_layer)
{
  if ( numProc > 1 )
  {
    TIMING_start("Sync_Poisson");
    int gd = GUIDE;
    if ( getSyncMode() == comm_sync )
    {
      if ( paraMngr->BndCommS3D(d_class, size[0], size[1], size[2], gd, num_layer, procGrp) != CPM_SUCCESS ) Exit(0);
    }
    else
    {
      MPI_Request req[12];
      for (int i=0; i<12; i++) req[i] = MPI_REQUEST_NULL;

      if ( paraMngr->BndCommS3D_nowait(d_class, size[0], size[1], size[2], gd, num_layer, req, procGrp) != CPM_SUCCESS ) Exit(0);
      if ( paraMngr->wait_BndCommS3D  (d_class, size[0], size[1], size[2], gd, num_layer, req, procGrp) != CPM_SUCCESS ) Exit(0);
    }
    TIMING_stop("Sync_Poisson", face_comm_size*(double)num_layer*sizeof(REAL_TYPE));
  }
}



// #################################################################
// PBiCBSTAB 収束判定は残差
// @note 反復回数が試行毎に異なる 内積のOpenMP並列のため
int LinearSolver::PBiCGstab(REAL_TYPE* x, REAL_TYPE* b, const REAL_TYPE dt, const double b_l2, const double r0_l2)
{
  double var[3];          /// 誤差, 残差, 解ベクトルのL2ノルム
  var[0] = var[1] = var[2] = 0.0;
  double flop = 0.0;

  REAL_TYPE cs = 0.0;

  TIMING_start("Blas_Clear");
  FBUtility::initS3D(pcg_q , size, 0.0);
  TIMING_stop("Blas_Clear", 0.0, 8);

  TIMING_start("Blas_Residual");
  flop = 0.0;
  blas_calc_rk_(pcg_r, x, b, bcp, size, pitch, &cs, &flop);
  TIMING_stop("Blas_Residual", flop);

  SyncScalar(pcg_r, 1);

  TIMING_start("Blas_Copy");
  blas_copy_(pcg_r0, pcg_r, size);
  TIMING_stop("Blas_Copy");

  double rho_old = 1.0;
  double alpha = 0.0;
  double omega  = 1.0;
  double r_omega = -omega;
  int lc=0;                      /// ループカウント

  for (lc=1; lc<getMaxIteration(); lc++)
  {
    double rho = Fdot2(pcg_r, pcg_r0);

    if( fabs(rho) < FLT_MIN )
    {
      lc = 0;
      break;
    }

    if( lc == 1 )
    {
      TIMING_start("Blas_Copy");
      blas_copy_(pcg_p, pcg_r, size);
      TIMING_stop("Blas_Copy");
    }
    else
    {
      double beta = rho / rho_old * alpha / omega;

      TIMING_start("Blas_BiCG_1");
      flop = 0.0;
      blas_bicg_1_(pcg_p, pcg_r, pcg_q, &beta, &omega, size, &flop);
      TIMING_stop("Blas_BiCG_1", flop);
    }
    SyncScalar(pcg_p, 1);

    TIMING_start("Blas_Clear");
    FBUtility::initS3D(pcg_p_, size, 0.0);
    TIMING_stop("Blas_Clear");

    Preconditioner(pcg_p_, pcg_p, dt);

    TIMING_start("Blas_AX");
    flop = 0.0;
    blas_calc_ax_(pcg_q, pcg_p_, bcp, size, pitch, &cs, &flop);
    TIMING_stop("Blas_AX", flop);

    alpha = rho / Fdot2(pcg_q, pcg_r0);

    double r_alpha = -alpha;
    TIMING_start("Blas_TRIAD");
    flop = 0.0;
    blas_triad_(pcg_s, pcg_q, pcg_r, &r_alpha, size, &flop);
    TIMING_stop("Blas_TRIAD", flop);

    SyncScalar(pcg_s, 1);

    TIMING_start("Blas_Clear");
    FBUtility::initS3D(pcg_s_, size, 0.0);
    TIMING_stop("Blas_Clear");

    Preconditioner(pcg_s_, pcg_s, dt);

    TIMING_start("Blas_AX");
    flop = 0.0;
    blas_calc_ax_(pcg_t_, pcg_s_, bcp, size, pitch, &cs, &flop);
    TIMING_stop("Blas_AX", flop);

    omega = Fdot2(pcg_t_, pcg_s) / Fdot1(pcg_t_);
    r_omega = -omega;

    TIMING_start("Blas_BiCG_2");
    flop = 0.0;
    blas_bicg_2_(x, pcg_p_, pcg_s_, &alpha , &omega, size, &flop);
    TIMING_stop("Blas_BiCG_2", flop);

    TIMING_start("Blas_TRIAD");
    flop = 0.0;
    blas_triad_  (pcg_r, pcg_t_, pcg_s, &r_omega, size, &flop);
    TIMING_stop("Blas_TRIAD", flop);

    var[1] = Fdot1(pcg_r);

    if ( Fcheck(var, b_l2, r0_l2) == true ) break;

    rho_old = rho;
  }


  TIMING_start("Poisson_BC");
  BC->OuterPBC(x, ensPeriodic);
  TIMING_stop("Poisson_BC");


  SyncScalar(x, 1);

  return lc;
}
