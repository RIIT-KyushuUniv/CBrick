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

#ifndef _RIAMC_H_
#define _RIAMC_H_

#include <stdio.h>
#include <vector>
#include <string>

#include "riamc_Define.h"
#include "DomainInfo.h"
#include "riamc_Ffunc.h"
#include "riamc_util.h"
#include "rcVersion.h"

// Text parser
#include "TextParser.h"

// FX10 profiler
#if defined __K_FPCOLL
#include "fjcoll.h"
#elif defined __FX_FAPP
#include "fj_tool/fjcoll.h"
#endif

// Performance Monitor
#ifndef DISABLE_PMLIB
#include "PerfMonitor.h"
#endif

#ifdef _CDM_OUTPUT
// cdm_DFI
#include "cdm_DFI.h"
#endif


#ifndef _WIN32
#include <unistd.h>
#include <strings.h>
#else
#include "sph_win32_util.h"
#endif
#include <sys/types.h>

#if defined(IA32_LINUX) || defined(IA64_LINUX) || defined(SGI_ALTIX)
#include <sys/stat.h>
#endif

#ifdef MacOSX
#include <sys/uio.h>
#endif


using namespace std;

class RIAMC : public DomainInfo {

private:

  unsigned NSTEP;               ///< セッションの終了ステップ数

  REAL_TYPE RE;       ///< レイノルズ数
  REAL_TYPE ALPHA;    ///< 数値粘性の係数
  REAL_TYPE OMEGA;    ///< SORの加速係数
  REAL_TYPE EPS;      ///< 収束閾値
  REAL_TYPE CS;       ///< SMAGORINSKY CONSTANT
  int ITER;           ///< 最大反復回数

  int order_of_PM_key; ///< PMlib用の登録番号カウンタ < PM_NUM_MAX

  // 立方体の位置
  int G_cube_index[6];  ///< 立方体のインデクス　グローバル
  int bcf[6];           ///< 立方体のインデクス　ローカル, Fortran

  //TimeControl
  int ICON;             ///< 0 - initial / 1 - restart
  REAL_TYPE DT;         ///< 時間積分幅（無次元）
  int Interval_Log;     ///< logファイル出力間隔
  int Interval_Vis;     ///< visファイル出力間隔
  int Interval_Disp;    ///< 標準出力間隔

  // 制御
  int PM_Test;          ///< 性能テスト ONのとき、ファイルI/Oを禁止
  string Parallel_mode; ///< 並列モード文字列

  int LS_solver;  ///< 線形システムのソルバー
  int LS_norm;    ///< 線形システムのノルム

#ifndef DISABLE_PMLIB
  pm_lib::PerfMonitor PM;       ///< 性能モニタクラス
#endif

  double comm_size;

#ifdef _CDM_OUTPUT
  cdm_DFI *dfi_p;       ///< Pressure 出力用 cpm_DFI
  cdm_DFI *dfi_v;       ///< Velocity 出力用 cpm_DFI
#endif

public:

  REAL_TYPE TIME;

  REAL_TYPE* Z;
  REAL_TYPE* P;

  REAL_TYPE* V;
  REAL_TYPE* VD;
  REAL_TYPE* VV;
  REAL_TYPE* CV;

  REAL_TYPE* SGS;
  REAL_TYPE* RHS;

  REAL_TYPE* uvw;

public:
  /** コンストラクタ */
  RIAMC();

  /** デストラクタ */
  ~RIAMC() {};





private:

  /**
   * @brief タイミング測定開始
   * @param [in] key ラベル
   */
  inline void TIMING_start(const string key)
  {
#ifndef DISABLE_PMLIB
    // PMlib Intrinsic profiler
    PM.start(key);

    const char* s_label = key.c_str();

    // Venus FX profiler
#if defined __K_FPCOLL
    start_collection( s_label );
#elif defined __FX_FAPP
    fapp_start( s_label, 0, 0);
#endif

#endif // DISABLE_PMLIB
  }


  /**
   * @brief タイミング測定終了
   * @param [in] key             ラベル
   * @param [in] flopPerTask    「タスク」あたりの計算量/通信量(バイト) (ディフォルト0)
   * @param [in] iterationCount  実行「タスク」数 (ディフォルト1)
   */
  inline void TIMING_stop(const string key, double flopPerTask=0.0, int iterationCount=1)
  {
#ifndef DISABLE_PMLIB
    // Venus FX profiler
    const char* s_label = key.c_str();

#if defined __K_FPCOLL
    stop_collection( s_label );
#elif defined __FX_FAPP
    fapp_stop( s_label, 0, 0);
#endif

    // PMlib Intrinsic profiler
    PM.stop(key, flopPerTask, (unsigned)iterationCount);
#endif // DISABLE_PMLIB
  }


  // メモリ使用量の表示
  void displayMemoryInfo(FILE* fp, double G_mem, double L_mem, const char* str);


  // 履歴表示
  //void printHistory(FILE* fp, int step, REAL_TYPE time, int itr, REAL_TYPE rms);


  // パラメータの表示
  void printParameters(FILE* fp, int divtype);


  // 並列モード
  void setParallelism();


  // 時間積分幅や物理パラメータの設定
  int setParameters(TextParser* tp);

  // 配列のアロケート
  bool allocateArray(double &total);


#ifndef DISABLE_PMLIB
  // タイミング測定区間にラベルを与えるラッパー
  void set_label(const string label, pm_lib::PerfMonitor::Type type, bool exclusive=true);

  // タイミング測定区間にラベルを与える
  void set_timing_label();
#endif

  // 配列確保 S3D
  REAL_TYPE* Alloc_Real_S3D(const int* sz);

  // 配列確保 V3D
  REAL_TYPE* Alloc_Real_V3D(const int* sz);

  void printHistoryTitle(FILE* fp);

  void printHistory(FILE* fp, int step, REAL_TYPE tm, REAL_TYPE vmax, int ILAP,
                    REAL_TYPE RMSP, REAL_TYPE res, REAL_TYPE avp);

  bool CubeIndex();

  bool displayDomainInfo();


public:

  // 初期化格子生成など
  int Initialize(int argc, char **argv);


  // シミュレーションの1ステップの処理
  int MainLoop();


  // シミュレーションの終了時の処理
  int Post();

};


#endif // _RIAMC_H_
