#ifndef _FFV_LS_H_
#define _FFV_LS_H_
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

#include "cpm_ParaManager.h"
#include "IterationControl.h"
#include "riamc.h"

#ifndef DISABLE_PMLIB
using namespace pm_lib;
#endif

  double face_comm_size; ///< 通信量

  SetBC3D* BC;       ///< BCクラス

#ifndef DISABLE_PMLIB
  PerfMonitor* PM;   ///< PerfMonitor class
#endif

  REAL_TYPE* pcg_p;  ///< work for BiCGstab
  REAL_TYPE* pcg_p_; ///< work for BiCGstab
  REAL_TYPE* pcg_r;  ///< work for BiCGstab
  REAL_TYPE* pcg_r0; ///< work for BiCGstab
  REAL_TYPE* pcg_q ; ///< work for BiCGstab
  REAL_TYPE* pcg_s;  ///< work for BiCGstab
  REAL_TYPE* pcg_s_; ///< work for BiCGstab
  REAL_TYPE* pcg_t ; ///< work for BiCGstab
  REAL_TYPE* pcg_t_; ///< work for BiCGstab

  int cf_sz[3];     ///< SOR2SMAの反復の場合のバッファサイズ
  REAL_TYPE *cf_x;  ///< i方向のバッファ
  REAL_TYPE *cf_y;  ///< j方向のバッファ
  REAL_TYPE *cf_z;  ///< k方向のバッファ

public:

  /** コンストラクタ */
  LinearSolver() : DomainInfo(), IterationCtl() {
    C   = NULL;
    BC  = NULL;
    bcp = NULL;
    bcd = NULL;
    pcg_p  = NULL;
    pcg_p_ = NULL;
    pcg_r  = NULL;
    pcg_r0 = NULL;
    pcg_q  = NULL;
    pcg_s  = NULL;
    pcg_s_ = NULL;
    pcg_t  = NULL;
    pcg_t_ = NULL;
    cf_x = NULL;
    cf_y = NULL;
    cf_z = NULL;

    face_comm_size = 0.0;

    for (int i=0; i<3; i++)
    {
      cf_sz[i] = 0;
    }

#ifndef DISABLE_PMLIB
    PM  = NULL;
#endif

  }

  /**　デストラクタ */
  ~LinearSolver() {}


  /**
   * @brief  Fcheck 非Div反復
   * @retval 収束したら true
   * @param [in]  var    誤差、残差、解ベクトルのL2ノルム
   * @param [in]  b_l2   右辺ベクトルのL2ノルム
   * @param [in]  r0_l2  初期残差ベクトルのL2ノルム
   */
  bool Fcheck(double* var, const double b_l2, const double r0_l2);


  /**
   * @brief Fdot for 1 array
   * @retval  内積値
   * @param [in]   x   vector1
   */
  double Fdot1(REAL_TYPE* x);


  /**
   * @brief Fdot for 2 arrays
   * @retval  内積値
   * @param [in]   x   vector1
   * @param [in]   y   vector2
   */
  double Fdot2(REAL_TYPE* x, REAL_TYPE* y);


  /**
   * @brief Preconditioner
   * @param [in,out] x   解ベクトル
   * @param [in]     b   RHS vector
   * @param [in]     dt  時間積分幅
   */
  void Preconditioner(REAL_TYPE* x, REAL_TYPE* b, const REAL_TYPE dt);



  /**
   * @brief 反復の同期処理
   * @param [in,out] d_class   対象データ
   * @param [in]     num_layer 通信の袖数
   */
  void SyncScalar(REAL_TYPE* d_class, const int num_layer);


  /**
   * @brief 初期化
   * @param [in]  C      Controlクラス
   * @param [in]  BC     SetBC3Dクラス
   * @param [in]  ModeT  タイミング測定モード
   * @param [in]  f_comm 通信量
   * @param [in]  PM     PerfMonitorクラス
   * @param [in]  bcp    BCindex P
   * @param [in]  bcd    BCindex ID
   * @param [in]  pcg_p  array for BiCGstab
   * @param [in]  pcg_p_ array for BiCGstab
   * @param [in]  pcg_r  array for BiCGstab
   * @param [in]  pcg_r0 array for BiCGstab
   * @param [in]  pcg_q  array for BiCGstab
   * @param [in]  pcg_s  array for BiCGstab
   * @param [in]  pcg_s_ array for BiCGstab
   * @param [in]  pcg_t  array for BiCGstab
   * @param [in]  pcg_t_ array for BiCGstab
   * @param [in]  ensP   周期境界の存在
   * @param [in]  cf_sz  バッファサイズ
   * @param [in]  cf_x   バッファ x方向
   * @param [in]  cf_y   バッファ y方向
   * @param [in]  cf_z   バッファ z方向
   */
#ifndef DISABLE_PMLIB
  void Initialize(Control* C,
                  SetBC3D* BC,
                  int ModeT,
                  double f_comm,
                  PerfMonitor* PM,
                  int* bcp,
                  int* bcd,
                  REAL_TYPE* pcg_p,
                  REAL_TYPE* pcg_p_,
                  REAL_TYPE* pcg_r,
                  REAL_TYPE* pcg_r0,
                  REAL_TYPE* pcg_q,
                  REAL_TYPE* pcg_s,
                  REAL_TYPE* pcg_s_,
                  REAL_TYPE* pcg_t,
                  REAL_TYPE* pcg_t_,
                  const int* ensP,
                  const int* cf_sz,
                  REAL_TYPE* cf_x,
                  REAL_TYPE* cf_y,
                  REAL_TYPE* cf_z);
#else
  void Initialize(Control* C,
                  SetBC3D* BC,
                  int ModeT,
                  double f_comm,
                  int* bcp,
                  int* bcd,
                  REAL_TYPE* pcg_p,
                  REAL_TYPE* pcg_p_,
                  REAL_TYPE* pcg_r,
                  REAL_TYPE* pcg_r0,
                  REAL_TYPE* pcg_q,
                  REAL_TYPE* pcg_s,
                  REAL_TYPE* pcg_s_,
                  REAL_TYPE* pcg_t,
                  REAL_TYPE* pcg_t_,
                  const int* ensP,
                  const int* cf_sz,
                  REAL_TYPE* cf_x,
                  REAL_TYPE* cf_y,
                  REAL_TYPE* cf_z);
#endif

  /**
   * @brief SOR法
   * @retval 反復数
   * @param [in,out] x              解ベクトル
   * @param [in]     b              RHS vector
   * @param [in]     dt             時間積分幅
   * @param [in]     itrMax         反復最大値
   * @param [in]     b_l2           L2 norm of b vector
   * @param [in]     r0_l2          初期残差ベクトルのL2ノルム
   * @param [in]     converge_check 収束判定を行う(true)
   */
  int PointSOR(REAL_TYPE* x, REAL_TYPE* b, const REAL_TYPE dt, const int itrMax, const double b_l2, const double r0_l2, bool converge_check=true);


  /**
   * @brief 2色オーダリングSORのストライドメモリアクセス版
   * @retval 反復数
   * @param [in,out] x              解ベクトル
   * @param [in]     b              RHS vector
   * @param [in]     dt             時間積分幅
   * @param [in]     itrMax         反復最大値
   * @param [in]     b_l2           L2 norm of b vector
   * @param [in]     r0_l2          初期残差ベクトルのL2ノルム
   * @param [in]     converge_check 収束判定を行う(true)

  int SOR2_SMA(REAL_TYPE* x, REAL_TYPE* b, const REAL_TYPE dt, const int itrMax, const double b_l2, const double r0_l2, bool converge_check=true);
*/

  /**
   * @brief 前処理つきBiCGstab
   * @retval 反復数
   * @param [in,out] x       解ベクトル
   * @param [in]     b       RHS vector
   * @param [in]     dt      時間積分幅
   * @param [in]     b_l2    L2 norm of b vector
   * @param [in]     r0_l2   初期残差ベクトルのL2ノルム
   */
  int PBiCGstab(REAL_TYPE* x, REAL_TYPE* b, const REAL_TYPE dt, const double b_l2, const double r0_l2);

};

#endif // _FFV_LS_H_
