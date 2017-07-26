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

#include "riamc.h"

// #################################################################
// constructor
RIAMC::RIAMC()
{

  NSTEP = 0;

  TIME = 0.0e0;

  ITER  = 0;
  OMEGA = 0.0;
  EPS   = 0.0;
  ALPHA = 0.0;

  DT = 0.001;
  Interval_Log  = 0;
  Interval_Vis  = 0;
  Interval_Disp = 0;

  RE = 0.0;
  CS = 0.0;

  LS_solver=0;
  LS_norm=0;

  order_of_PM_key = 0;

  ICON = 0;

  comm_size = 0.0;

  for (int i=0; i<6; i++) {
    G_cube_index[i] =0;
    bcf[i]=0;
  }

  // Performance test off
  PM_Test = OFF;

  Z = NULL;
  P = NULL;

  V  = NULL;
  VD = NULL;
  VV = NULL;
  CV = NULL;

  SGS = NULL;
  RHS = NULL;

  uvw = NULL;

}



// #################################################################
/* @brief メモリ消費情報を表示
 * @param [in]     fp    ファイルポインタ
 * @param [in]     G_mem グローバルメモリサイズ
 * @param [in]     L_mem ローカルメモリサイズ
 * @param [in]     str   表示用文字列
 */
void RIAMC::displayMemoryInfo(FILE* fp, double G_mem, double L_mem, const char* str)
{
  if ( numProc > 1 )
  {
    double tmp_memory = G_mem;
    if ( paraMngr->Allreduce(&tmp_memory, &G_mem, 1, MPI_SUM, procGrp) != CPM_SUCCESS ) Exit(0);
  }

  Hostonly_
  {
    FBUtility::MemoryRequirement(str, G_mem, L_mem, fp);
  }

  Hostonly_
  {
    printf("\n\n");
    //fprintf(fp, "\n\n");
  }
}



// #################################################################
/**
 * @brief シミュレーションの1ステップの処理
 */
int RIAMC::MainLoop()
{
  //TIMING_start("NS_Section");

  double flop_count = 0.0;  ///< flops計算用

  //自ランクの格子数
  int NX = size[0];
  int NY = size[1];
  int NZ = size[2];

  REAL_TYPE dh = (REAL_TYPE)pitch[0];

  // 圧力の平均計算に利用
  REAL_TYPE avp_normal = 1.0 / ( (REAL_TYPE)G_size[0]
                               * (REAL_TYPE)G_size[1]
                               * (REAL_TYPE)G_size[2] );

  // 1点あたりのrmsの計算に利用
  REAL_TYPE rms_normal = 1.0 / ( (REAL_TYPE)(G_size[0]-2)
                               * (REAL_TYPE)(G_size[1]-2)
                               * (REAL_TYPE)(G_size[2]-2) );

  int ILAP = 0;;
  REAL_TYPE res_p = 0.0;
  REAL_TYPE rms_p = 0.0;

  int ICOUNT = 0;
  int T_NSTEP = NSTEP;    //最終計算ステップをダミーにセット


  REAL_TYPE PP = 0.0;

  for (int ISTEP = 1; ISTEP <= NSTEP; ISTEP++)
  {
    /*
    if ( TerminateCntl::getTerminateFlag() )
    {
      fprintf(stdout, "\n\tForced terminate.\n\n");
      return 0; // forced terminate
    }*/

    TIME = TIME + DT;

    TIMING_start("Eddy_viscosity");
    flop_count = 0.0;
    eddy_viscosity_(size, &RE, &CS, Z, &dh, V, SGS, comm_tbl, &flop_count);
    TIMING_stop("Eddy_viscosity", flop_count);


    TIMING_start("BC_Eddy_viscosity");
    flop_count = 0.0;
    bc_eddy_(size, SGS, comm_tbl);
    TIMING_stop("BC_Eddy_viscosity", flop_count);


    // SGSの袖通信を行う
    if ( numProc > 1 )
    {
      TIMING_start("Comm_SGS");
      int gd = GUIDE;
      if( (paraMngr->BndCommS3D(SGS, NX, NY, NZ, gd, gd, 0)) != CPM_SUCCESS )
      {
          stamped_printf("\tError at SGS Comm Band Cell \n");
          return 0;
      }
      TIMING_stop("Comm_SGS", comm_size);
    }



    // STORE PREVIOUS VELOCITY
    TIMING_start("Copy_Vec");
    copy_vec_(size, V, VD);
    TIMING_stop("Copy_Vec", 0.0);


    if ( numProc > 1 )
    {
      TIMING_start("Comm_Vec_prev");
      int gd = GUIDE;
      if( (paraMngr->BndCommV3D(VD, NX, NY, NZ, gd, gd, 0)) != CPM_SUCCESS )
      {
          stamped_printf("\tError at VD Comm Band Cell \n");
          return 0;
      }
      TIMING_stop("Comm_Vec_prev", comm_size*3.0);
    }



    TIMING_start("Intermediate_v");
    flop_count = 0.0;
    intermediate_v_(size, &ALPHA, &DT, &RE, Z, &dh, V, VD, VV, SGS, comm_tbl, &flop_count);
    TIMING_stop("Intermediate_v", flop_count);


    // 計算結果 U,V,Wの袖通信を行う
    if ( numProc > 1 )
    {
      TIMING_start("Comm_InterV");
      int gd = GUIDE;
      if( (paraMngr->BndCommV3D(V, NX, NY, NZ, gd, gd, 0)) != CPM_SUCCESS )
      {
          stamped_printf("\tError at V Comm Band Cell \n");
          return 0;
      }
      TIMING_stop("Comm_InterV", comm_size*3.0);
    }



    TIMING_start("Contra_V_CC");
    flop_count = 0.0;
    contra_v_cc_(size, V, Z, &dh, CV, &flop_count);
    TIMING_stop("Contra_V_CC", flop_count);


    if ( numProc > 1 )
    {
      TIMING_start("Comm_CV");
      int gd = GUIDE;
      if( (paraMngr->BndCommV3D(CV, NX, NY, NZ, gd, gd, 0)) != CPM_SUCCESS )
      {
          stamped_printf("\tError at Comm Band Cell \n");
          return 0;
      }
      TIMING_stop("Comm_CV", comm_size*3.0);
    }



    TIMING_start("Contra_V_CF");
    flop_count = 0.0;
    contra_v_cf_(size, VV, CV, &flop_count);
    TIMING_stop("Contra_V_CF", flop_count);


    TIMING_start("RHS_Poisson");
    flop_count = 0.0;
    rhs_poisson_(size, &DT, VV, RHS, Z, &dh, &flop_count);
    TIMING_stop("RHS_Poisson", flop_count);



    // Poisson Iteration
    for (ILAP=1; ILAP<=ITER; ILAP++)
    {
      res_p = 0.0;
      rms_p = 0.0;

      // Switch SOLVER
      if ( LS_solver == LS_PSOR )
      {
        /*
        TIMING_start("Poisson_PSOR");
        flop_count = 0.0;
        psor_(size, &OMEGA, &rms_p, &res_p,
              P, C1, C2, C3, C5, C6, C7, C8, C9, RHS, comm_tbl, &flop_count);
        TIMING_stop("Poisson_PSOR", flop_count);
        */
      }
      else if ( LS_solver == LS_PSOR_MAF )
      {
        TIMING_start("Poisson_PSOR_MAF");
        flop_count = 0.0;
        psor_maf_(size, &OMEGA, &rms_p, &res_p,
              P, Z, &dh, RHS, comm_tbl, &flop_count);
        TIMING_stop("Poisson_PSOR_MAF", flop_count);
      }
      else if ( LS_solver == LS_SOR2SMA )
      {
        /*
        // 2色のマルチカラー(Red&Black)のセットアップ
        // ip = 0 基点(1,1,1)が Rからスタート
        //    = 1 基点(1,1,1)が Bからスタート
        int ip=0;
        if ( numProc > 1 )
        {
          ip = (head[0]+head[1]+head[2]+1) % 2;
        }
        else
        {
          ip = 0;
        }

        // 各カラー毎の間に同期, 残差は色間で積算する
        // R - color=0 / B - color=1
        TIMING_start("Poisson_2coloredSOR");
        flop_count = 0.0;
        for (int color=0; color<2; color++)
        {
          psor2sma_(size, &OMEGA, &rms_p, &res_p, &ip, &color,
                    P, C1, C2, C3, C5, C6, C7, C8, C9, RHS, comm_tbl, &flop_count);
        }
        TIMING_stop("Poisson_2coloredSOR", flop_count);
        */
      }
      else if ( LS_solver == LS_SOR2SMA_MAF )
      {
        // 2色のマルチカラー(Red&Black)のセットアップ
        // ip = 0 基点(1,1,1)が Rからスタート
        //    = 1 基点(1,1,1)が Bからスタート
        int ip=0;
        if ( numProc > 1 )
        {
          ip = (head[0]+head[1]+head[2]+1) % 2;
        }
        else
        {
          ip = 0;
        }

        // 各カラー毎の間に同期, 残差は色間で積算する
        // R - color=0 / B - color=1
        TIMING_start("Poisson_2_SOR_MAF");
        flop_count = 0.0;
        for (int color=0; color<2; color++)
        {
          psor2sma_maf_(size, &OMEGA, &rms_p, &res_p, &ip, &color,
                    P, Z, &dh, RHS, comm_tbl, &flop_count);
        }
        TIMING_stop("Poisson_2_SOR_MAF", flop_count);
      }
      else if ( LS_solver == LS_BiCGSTAB )
      {

      } // Switch SOLVER


      if ( numProc > 1 )
      {
        TIMING_start("Comm_Poisson");
        int gd = GUIDE;
        if( (paraMngr->BndCommS3D(P, NX, NY, NZ, gd, gd, 0)) != CPM_SUCCESS )
        {
          stamped_printf("\tError at P Comm Band Cell \n");
          return 0;
        }
        TIMING_stop("Comm_Poisson", comm_size);
      }


      if ( numProc > 1 )
      {
        TIMING_start("Comm_Res_Poisson");
        REAL_TYPE tmp[2] = {rms_p, res_p};
        REAL_TYPE buf[2] = {rms_p, res_p};
        if ( paraMngr->Allreduce(tmp, buf, 2, MPI_SUM, procGrp) != CPM_SUCCESS ) Exit(0);
        rms_p = buf[0];
        res_p = buf[1];
        TIMING_stop("Comm_Res_Poisson", 4.0*numProc*sizeof(REAL_TYPE));
      }


      rms_p *= rms_normal;
      res_p *= rms_normal;
      rms_p = sqrt(rms_p);
      res_p = sqrt(res_p);

      // convergence check
      if (PM_Test == OFF ) {

        if ( LS_norm == Norm_Rms )
        {
          if (rms_p <= EPS) break;
        }
        else { // residual
          if (res_p <= EPS) break;
        }
      }

    }  // Iteration



    REAL_TYPE p_av = 0.0;

    TIMING_start("Average_pressure");
    flop_count = 0.0;
    average_prs_(size, P, comm_tbl, &p_av, &flop_count);
    TIMING_stop("Average_pressure", flop_count);

    if ( numProc > 1 )
    {
      TIMING_start("Comm_Sum_avrP");
      REAL_TYPE tmp = p_av;
      if ( paraMngr->Allreduce(&tmp, &p_av, 1, MPI_SUM, procGrp) != CPM_SUCCESS ) Exit(0);
      TIMING_stop("Comm_Sum_avrP", 2.0*numProc*sizeof(REAL_TYPE));
    }


    TIMING_start("New_pressure");
    flop_count = 0.0;
    p_av *= avp_normal;
    new_pressure_(size, P, &p_av, &flop_count);
    TIMING_stop("New_pressure", flop_count);


    TIMING_start("BC_New_pressure");
    flop_count = 0.0;
    bc_new_prs_(size, &RE, P, Z, &dh, V, comm_tbl, &flop_count);
    TIMING_stop("BC_New_pressure", flop_count);


    // 計算結果 Pの袖通信を行う
    if ( numProc > 1 )
    {
      TIMING_start("Comm_New_P");
      int gd = GUIDE;
      if( (paraMngr->BndCommS3D( P, NX, NY, NZ, gd, gd, 0)) != CPM_SUCCESS )
      {
          stamped_printf("\tError at P Comm Band Cell \n");
          return 0;
      }
      TIMING_stop("Comm_New_P", comm_size);
    }


    TIMING_start("New_V_CC");
    flop_count = 0.0;
    new_vec_cc_(size, &DT, P, V, Z, &dh, comm_tbl, &flop_count);
    TIMING_stop("New_V_CC", flop_count);


    TIMING_start("New_V_CF");
    flop_count = 0.0;
    new_vec_cf_(size, &DT, P, VV, Z, &dh, &flop_count);
    TIMING_stop("New_V_CF", flop_count);


    TIMING_start("BC_New_velocity");
    flop_count = 0.0;
    bc_new_vel_(size, bcf, &DT, &dh, V, VD, comm_tbl, &flop_count);
    TIMING_stop("BC_New_velocity", flop_count);


    // 計算結果の袖通信を行う
    if ( numProc > 1 )
    {
      TIMING_start("Comm_New_V");
      int gd = GUIDE;
      if( (paraMngr->BndCommV3D(V, NX, NY, NZ, gd, gd, 0)) != CPM_SUCCESS )
      {
          stamped_printf("\tError at V Comm Band Cell \n");
          return 0;
      }
      if( (paraMngr->BndCommV3D(VV, NX, NY, NZ, gd, gd, 0)) != CPM_SUCCESS )
      {
          stamped_printf("\tError at VV Comm Band Cell \n");
          return 0;
      }
      TIMING_stop("Comm_New_V", comm_size*6.0);
    }



    if (ISTEP % Interval_Disp == 0 || ISTEP == NSTEP)
    {
      TIMING_start("Max_velocity");
      flop_count = 0.0;
      REAL_TYPE v_max;
      find_max_v_(size, &v_max, V, &flop_count);
      TIMING_stop("Max_velocity", flop_count);

      TIMING_start("Comm_Vmax");
      REAL_TYPE vtmp = v_max;
      if ( paraMngr->Allreduce(&vtmp, &v_max, 1, MPI_MAX, procGrp) != CPM_SUCCESS ) Exit(0);
      TIMING_stop("Comm_Vmax", 2.0*numProc*sizeof(REAL_TYPE));

      TIMING_start("History");
      Hostonly_ {
        printHistory(stdout, ISTEP, TIME, v_max, ILAP, rms_p, res_p, p_av);
        fflush(stdout);
      }
      TIMING_stop("History", 0.0);
    }

/*
    if (PM_Test == OFF ) {
      TIMING_start("Output_FieldData");
      display_(&MX, &MY, &MZ, &ISTEP, &T_NSTEP, &ILAP, &ICOUNT,
               &Interval_Disp, U, V, W, &Interval_Log, &Interval_Vis, &TIME, &RE, &rms_p, &myRank);
      TIMING_stop("Output_FieldData", 0.0);
    }
*/

    if (PM_Test == OFF ) {
#ifdef _CDM_OUTPUT
      if( ISTEP%Interval_Vis == 0 || ISTEP == NSTEP )
      {
        int gd = GUIDE;
        //// Pressure
        if( dfi_p )
        {
          // write
          dfi_p->WriteData( (unsigned)ISTEP
                         , TIME
                         , paraMngr->GetLocalNodeSize()
                         , 1
                         , gd
                         , P
                         , (REAL_TYPE*)NULL
                         , false
                         , 0
                         , 0.0
                         );
        }
        //// Velocity
        if( dfi_v )
        {
           // copy array
           copy_v2vex_(size, V, uvw);
         // write
           dfi_v->WriteData( (unsigned)ISTEP
                         , TIME
                         , paraMngr->GetLocalNodeSize()
                         , 3
                         , gd
                         , uvw
                         , (REAL_TYPE*)NULL
                         , false
                         , 0
                         , 0.0
                         );
        }
      }
#endif
    } // PM_Test
  }


  //TIMING_stop("NS_Section");

  return 1;
}


// #################################################################
/**
 * @brief シミュレーションの1ステップの処理
 */
int RIAMC::Post()
{

  double flop_count = 0.0;  ///< flops計算用

  if (PM_Test == OFF ) {
    TIMING_start("Output_FieldData");
    flop_count = 0.0;
    display2_(size, &TIME, &RE, V, P, VV, &myRank);
    TIMING_stop("Output_FieldData", flop_count);
  }

/*
  FILE *fp = NULL;

  Hostonly_
  {
    if ( !(fp=fopen("profiling.txt", "w")) )
    {
      stamped_printf("\tSorry, can't open 'profiling.txt' file. Write failed.\n");
      Exit(0);
    }
  }
*/

  char str[100];
  sprintf(str, "RIAM-Compact %s", RC_VERSION);

#ifndef DISABLE_PMLIB
  // 測定結果の集計(gathreメソッドは全ノードで呼ぶこと)
  PM.gather();

//  Hostonly_
//  {
    // 結果出力(排他測定のみ)
    string HostName = paraMngr->GetHostName();
    PM.print(stdout, HostName, str);
    //PM.print(fp, HostName, str);

    // 結果出力(非排他測定も)
    PM.printDetail(stdout);

//    if ( !fp ) fclose(fp);
//  }
#endif

  return 1;

}


// #################################################################
/* @brief 並列化と分割の方法を保持
 */
void RIAMC::setParallelism()
{
  // Serial or Parallel environment
  if( paraMngr->IsParallel() )
  {
    if ( numThreads > 1 )
    {
      Parallel_mode = "Hybrid";
    }
    else
    {
      Parallel_mode = "FlatMPI";
    }
  }
  else
  {
    numProc = 1;

    if ( numThreads > 1 )
    {
      Parallel_mode = "OpenMP";
    }
    else
    {
      Parallel_mode = "Serial";
    }
  }

}


#ifndef DISABLE_PMLIB
// #################################################################
/**
 * @brief タイミング測定区間にラベルを与えるラッパー
 * @param [in] label     ラベル
 * @param [in] type      測定対象タイプ(COMM or CALC)
 * @param [in] exclusive 排他測定フラグ(ディフォルトtrue)
 */
void RIAMC::set_label(const string label, pm_lib::PerfMonitor::Type type, bool exclusive)
{
  // 登録個数のチェック
  order_of_PM_key++;

  if ( order_of_PM_key > PM_NUM_MAX )
  {
    fprintf(stdout, "\tThe number of labels for Performance monitor goes over limit.\n");
    exit(0);
  }

  // 文字数がTM_LABEL_MAX-1を超えるものはカット
  if ( strlen(label.c_str()) > TM_LABEL_MAX-1 )
  {
    printf("\tWarning: Length of timing label must be less than %d\n", TM_LABEL_MAX-1);
  }

  // Performance Monitorへの登録
  PM.setProperties(label, type, exclusive);
}


// #################################################################
/**
 * @brief タイミング測定区間にラベルを与える
 */
void RIAMC::set_timing_label()
{
  using namespace pm_lib;
  // common
  set_label("Allocate_Arrays",         PerfMonitor::CALC);
  // common


  // Initialization_Section
  //set_label("Initialization_Section",  PerfMonitor::CALC, false);
  set_label("Domain_Decomposition",    PerfMonitor::CALC);
  set_label("Init_array",              PerfMonitor::CALC);
  set_label("Init_Vars",               PerfMonitor::CALC);
  set_label("Grid_gen",                PerfMonitor::CALC);
  set_label("Comm_Grid",               PerfMonitor::COMM);
  //set_label("Metrics",                 PerfMonitor::CALC);
  //set_label("Comm_Metrics",            PerfMonitor::COMM);

  // NS_Section
  //set_label("NS_Section",       PerfMonitor::CALC, false);
  set_label("Eddy_viscosity",      PerfMonitor::CALC);
  set_label("BC_Eddy_viscosity",   PerfMonitor::CALC);
  set_label("Comm_SGS",            PerfMonitor::COMM);
  set_label("Copy_Vec",            PerfMonitor::CALC);
  set_label("Comm_Vec_prev",       PerfMonitor::COMM);
  set_label("Intermediate_v",      PerfMonitor::CALC);
  set_label("Comm_InterV",         PerfMonitor::COMM);
  set_label("Contra_V_CC",         PerfMonitor::CALC);
  set_label("Comm_CV",             PerfMonitor::COMM);
  set_label("Contra_V_CF",         PerfMonitor::CALC);
  set_label("RHS_Poisson",         PerfMonitor::CALC);
  set_label("Poisson_PSOR",        PerfMonitor::CALC);
  set_label("Poisson_PSOR_MAF",    PerfMonitor::CALC);
  set_label("Poisson_2coloredSOR", PerfMonitor::CALC);
  set_label("Poisson_2_SOR_MAF",   PerfMonitor::CALC);
  set_label("Comm_Poisson",        PerfMonitor::COMM);
  set_label("Comm_Res_Poisson",    PerfMonitor::COMM);
  set_label("Average_pressure",    PerfMonitor::CALC);
  set_label("New_pressure",        PerfMonitor::CALC);
  set_label("Comm_Sum_avrP",       PerfMonitor::COMM);
  set_label("BC_New_pressure",     PerfMonitor::CALC);
  set_label("Comm_New_P",          PerfMonitor::COMM);
  set_label("New_V_CC",            PerfMonitor::CALC);
  set_label("New_V_CF",            PerfMonitor::CALC);
  set_label("BC_New_velocity",     PerfMonitor::CALC);
  set_label("Comm_New_V",          PerfMonitor::COMM);
  set_label("Max_velocity",        PerfMonitor::CALC);
  set_label("Comm_Vmax",           PerfMonitor::COMM);
  set_label("History",             PerfMonitor::CALC);

  // Post
  //set_label("Post_Section",     PerfMonitor::CALC, false);
  set_label("Output_FieldData",    PerfMonitor::CALC);

}
#endif

// #################################################################
/**
 * @brief 時刻歴出力
 * @param [in] fp   FILE pointer
 * @param [in] step ステップ
 * @param [in] tm   無次元時刻
 * @param [in] vmax 速度の最大値
 * @param [in] ILAP 反復回数
 * @param [in] rms_p 圧力反復のrms値
 * @param [in] res 圧力反復の残差
 * @param [in] avp 圧力平均値
 */
void RIAMC::printHistory(FILE* fp, int step, REAL_TYPE tm, REAL_TYPE vmax, int ILAP,
                         REAL_TYPE rms_p, REAL_TYPE res, REAL_TYPE avp)
{
  fprintf(fp, "%8d %14.6e %11.4e %5d    %12.5e   %12.5e   %11.4e\n",
          step, tm, vmax, ILAP, rms_p, res, avp);
}

// #################################################################
/**
 * @brief 時刻歴のタイトル出力
 * @param [in] fp   FILE pointer
 */
void RIAMC::printHistoryTitle(FILE* fp)
{
  fprintf(fp, "    step        time[-]    v_max[-]  ItrP            rmsP           resP          avrP\n");
}
