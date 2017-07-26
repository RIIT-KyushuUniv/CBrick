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
#include "rcVersion.h"
#include <climits>

int RIAMC::Initialize(int argc, char **argv)
{
  int ret = 0;                 // リターンコード
  int div_type = 0;            // 領域分割モード (1 - 分割数指定, 2 - 自動分割)
  double G_TotalMemory = 0.0;  ///< 計算に必要なメモリ量（グローバル）
  double TotalMemory   = 0.0;  ///< 計算に必要なメモリ量（ローカル）

  int NX, NY, NZ;


  // 並列管理クラスのインスタンスと初期化
  paraMngr = cpm_ParaManager::get_instance(argc, argv);
  if( !paraMngr ) return CPM_ERROR_PM_INSTANCE;

  myRank = paraMngr->GetMyRankID();

  setRankInfo(paraMngr, procGrp);

  MPI_Comm_size(MPI_COMM_WORLD, &numProc);

#ifdef _OPENMP
  numThreads = omp_get_max_threads();
  omp_set_num_threads(numThreads);
#endif

  setParallelism();



  double flop_count    = 0.0;  ///< flops計算用

  // パラメータローダのインスタンス生成
  TextParser tp;

  // パラメータのロードと保持
  if( argc>1 )
  {
    // 入力ファイルの指定
    string input_file = argv[1];
    int ierror=0;

    if ( (ierror = tp.read(input_file)) != TP_NO_ERROR )
    {
      Hostonly_ stamped_printf("\tError at reading '%s' file : %d\n", input_file.c_str(), ierror);
      Exit(0);
    }

    // 固定パラメータ
    if ( -1 == (div_type=setParameters(&tp)) ) return 0;

  }
  else {
    stamped_printf("\tError undefined input tp file\n");
    return 0;
  }

  // メッセージ表示
Hostonly_
{
  string ver_rc =RC_VERSION;
  string ver_tp =TP_VERSION;
  string ver_cpm=CPM_VERSION;

#ifndef DISABLE_PMLIB
  string ver_pm =PM_VERSION;
#endif

  fprintf(stdout, "\n\tWelcome to RIAM-COMPACT  \tVersion %s\n", ver_rc.c_str());
  fprintf(stdout, "\n\t      with\n");
  fprintf(stdout, "\t           TextParser    \tVersion %s\n", ver_tp.c_str());

#ifndef DISABLE_PMLIB
  fprintf(stdout, "\t           PMlib         \tVersion %s\n", ver_pm.c_str());
#endif

  fprintf(stdout, "\t           CPMlib        \tVersion %s\n", ver_cpm.c_str());

#if defined(CDM_VERSION)
  {
    string ver_cdm=CDM_VERSION;
    fprintf(stdout, "\t           CDMlib        \tVersion %s\n\n", ver_cdm.c_str());
  }
#endif
}

#ifndef DISABLE_PMLIB
  // タイミング測定の初期化
  PM.initialize( PM_NUM_MAX );
  PM.setRankInfo( myRank );
  PM.setParallelMode(Parallel_mode, numThreads, numProc);
  set_timing_label();
#endif


  // 領域分割 デフォルトで DIV_COMM_SIZE モード、プロセスグループ 0
  TIMING_start("Domain_Decomposition");
  size_t Nvc  = GUIDE; // 袖通信の最大数
  size_t Ncmp = 3;     // 最大はベクトル3成分

  if ( numProc > 1 )
  {

    if (div_type == 1)
    {
      // 分割数指定
      if( paraMngr->NodeInit(G_division, G_size, G_origin, G_region, Nvc, Ncmp) != CPM_SUCCESS )
      {
        Hostonly_ stamped_printf("\tError at CPM::NodeInit\n");
        Exit(0);
      }
    }
    else
    {
      // 自動分割
      if( paraMngr->NodeInit(G_size, G_origin, G_region, Nvc, Ncmp) != CPM_SUCCESS )
      {
        Hostonly_ stamped_printf("\tError at CPM::NodeInit\n");
        Exit(0);
      }
    }

    // 領域分割数取得
    const int* m_div = paraMngr->GetDivNum();
    G_division[0] = m_div[0];
    G_division[1] = m_div[1];
    G_division[2] = m_div[2];

    //自ランクのHeadIndexの取得
    const int* hdx = paraMngr->GetNodeHeadIndex();
    head[0] = hdx[0];
    head[1] = hdx[1];
    head[2] = hdx[2];

    //自ランクのTailIndexの取得
    const int* tal = paraMngr->GetNodeTailIndex();
    tail[0] = tal[0];
    tail[1] = tal[1];
    tail[2] = tal[2];

    //自ランクの格子数取得
    const int* lsize = paraMngr->GetLocalNodeSize();
    NX  = size[0] = lsize[0];
    NY  = size[1] = lsize[1];
    NZ  = size[2] = lsize[2];

    //自ランクの隣接ランク番号を取得
    const int* nbr = paraMngr->GetNeighborRankID();
    for (int i=0; i<6; i++) comm_tbl[i] = nbr[i];
  }
  else
  {
    G_division[0] = 1;
    G_division[1] = 1;
    G_division[2] = 1;

    NX  = size[0] = G_size[0];
    NY  = size[1] = G_size[1];
    NZ  = size[2] = G_size[2];

    head[0] = 0;
    head[1] = 0;
    head[2] = 0;

    tail[0] = NX-1;
    tail[1] = NY-1;
    tail[2] = NZ-1;
  }


  // 通信量　双方向 x ２面
   comm_size = (double)( (NX+2*GUIDE) * (NY+2*GUIDE)
                       + (NY+2*GUIDE) * (NZ+2*GUIDE)
                       + (NX+2*GUIDE) * (NZ+2*GUIDE) )
                       * 2.0 * 2.0 * sizeof(REAL_TYPE);


  // パラメータ表示
  Hostonly_ printParameters(stdout, div_type);

  if ( !displayDomainInfo() ) return 0;
  if ( !CubeIndex() ) return 0;

  TIMING_stop("Domain_Decomposition");



  // 配列のアロケート
  TIMING_start("Allocate_Arrays");
  if( !allocateArray(TotalMemory) )
  {
      stamped_printf("\tError at Allocate_Arrays\n");
      return 0;
  }
  TIMING_stop("Allocate_Arrays");


  // メモリ消費量の情報を表示

  Hostonly_
  {
    printf(    "\n----------\n\n");
  }
  G_TotalMemory = TotalMemory;

  displayMemoryInfo(stdout, G_TotalMemory, TotalMemory, "Solver");


  // 初期化
  TIMING_start("Init_Vars");
  flop_count = 0.0;
  var_init_(size, &ICON, V, P, VV);
  TIMING_stop("Init_Vars", flop_count);


  TIMING_start("Grid_gen");
  flop_count = 0.0;
  REAL_TYPE m_org[3] = {(REAL_TYPE)G_origin[0],
                        (REAL_TYPE)G_origin[1],
                        (REAL_TYPE)G_origin[2]
                       };
  REAL_TYPE m_pch[3] = {(REAL_TYPE)pitch[0],
                        (REAL_TYPE)pitch[1],
                        (REAL_TYPE)pitch[2]
                       };
  grid_gen_(size, Z, m_org, m_pch, head, comm_tbl, &flop_count);
  TIMING_stop("Grid_gen", flop_count);


  if ( numProc > 1 )
  {
    TIMING_start("Comm_Grid");
    int gd = GUIDE;

    if( (paraMngr->BndCommS3D( Z, NX, NY, NZ, gd, gd, 0)) != CPM_SUCCESS )
    {
      stamped_printf("\tError at Comm Band Cell \n");
      return 0;
    }
    TIMING_stop("Comm_Grid", comm_size*3.0);
  }




  #ifdef _CDM_OUTPUT
    //CDM initialize
    //CDM用headとtailの作成
    int cdm_head[3] = {head[0]+1, head[1]+1, head[2]+1};
    int cdm_tail[3] = {tail[0]+1, tail[1]+1, tail[2]+1};
    int gd = GUIDE;

    //Pressure
    dfi_p = cdm_DFI::WriteInit( MPI_COMM_WORLD
                              , cdm_DFI::Generate_DFI_Name("prs")
                              , "sph"
                              , "prs"
                              , CDM::E_CDM_FMT_SPH
                              , gd
                              , CDM::E_CDM_FLOAT64
                              , 1
                              , "proc.dfi"
                              , paraMngr->GetGlobalNodeSize()
                              , paraMngr->GetPitch()
                              , paraMngr->GetGlobalOrigin()
                              , paraMngr->GetDivNum()
                              , cdm_head
                              , cdm_tail
                              , "hogehoge"
                              , CDM::E_CDM_OFF );

    if( dfi_p )
    {
      dfi_p->setVariableName(0, "Pressure");
      dfi_p->WriteProcDfiFile(MPI_COMM_WORLD, false, 0, 0);
    }

    // Velocity
    dfi_v = cdm_DFI::WriteInit( MPI_COMM_WORLD
                              , cdm_DFI::Generate_DFI_Name("vel")
                              , "sph"
                              , "vel"
                              , CDM::E_CDM_FMT_SPH
                              , gd
                              , CDM::E_CDM_FLOAT64
                              , 3
                              , "proc.dfi"
                              , paraMngr->GetGlobalNodeSize()
                              , paraMngr->GetPitch()
                              , paraMngr->GetGlobalOrigin()
                              , paraMngr->GetDivNum()
                              , cdm_head
                              , cdm_tail
                              , "hogehoge"
                              , CDM::E_CDM_OFF );

    if( dfi_v )
    {
      dfi_v->setVariableName(0, "U");
      dfi_v->setVariableName(1, "V");
      dfi_v->setVariableName(2, "W");
    }
  #endif


  Hostonly_ printHistoryTitle(stdout);

  return 1;
}



// #################################################################
/**
 * @brief パラメータの表示
 * @param [in] fp       FILE pointer
 * @param [in] div_type 分割モード
 */
void RIAMC::printParameters(FILE* fp, int divtype)
{
  fprintf(fp, "\n------------\n");
  fprintf(fp, "PARAMETERS\n");
  fprintf(fp, "------------\n\n");

  fprintf(fp, "Domain\n");
  fprintf(fp, "\tGlobalOrigin   : (%16.6e, %16.6e, %16.6e)\n", G_origin[0], G_origin[1], G_origin[2]);
  fprintf(fp, "\tGlobalGrid     : (%16d, %16d, %16d)\n", G_size[0], G_size[1], G_size[2]);
  fprintf(fp, "\tGlobalRegion   : (%16.6e, %16.6e, %16.6e)\n", G_region[0], G_region[1], G_region[2]);
  fprintf(fp, "\tGlobalPitch    : (%16.6e, %16.6e, %16.6e)\n", pitch[0], pitch[1], pitch[2]);
  fprintf(fp, "\tDDM mode       : %s\n",(divtype==1)?"Specified":"AUTO");
  fprintf(fp, "\tGlobalDivision : (%16d, %16d, %16d)\n", G_division[0], G_division[1], G_division[2]);
  fprintf(fp, "\n\n");

  fprintf(fp, "Parameter\n");
  fprintf(fp, "\tBOX(global)\n");
  fprintf(fp, "\t\tMIN          : (%10d, %10d, %10d)\n",G_cube_index[0],G_cube_index[1],G_cube_index[2]);
  fprintf(fp, "\t\tMAX          : (%10d, %10d, %10d)\n",G_cube_index[3],G_cube_index[4],G_cube_index[5]);
  fprintf(fp, "\n");
  fprintf(fp, "\tReynolds                 : %16.6e\n", RE);
  fprintf(fp, "\tAlpha                    : %10.4f\n", ALPHA);
  fprintf(fp, "\tSmagorinsky Constant     : %16.6e\n", CS);
  fprintf(fp, "\n");

  fprintf(fp, "Linear System\n");
  fprintf(fp, "\tSolver                   : ");
  if      (LS_solver == LS_PSOR)        fprintf(fp, "Point SOR\n");
  else if (LS_solver == LS_PSOR_MAF)    fprintf(fp, "Point SOR metrics array free\n");
  else if (LS_solver == LS_SOR2SMA)     fprintf(fp, "2 color SOR\n");
  else if (LS_solver == LS_SOR2SMA_MAF) fprintf(fp, "2 color SOR metrics array free\n");
  else if (LS_solver == LS_BiCGSTAB)    fprintf(fp, "BiCGstab\n");
  fprintf(fp, "\tNorm                     : ");
  if (LS_norm == Norm_Res) fprintf(fp, "Residual\n");
  else if (LS_norm == Norm_Rms) fprintf(fp, "RMS\n");
  fprintf(fp, "\tIteration Max            : %10d\n", ITER);
  fprintf(fp, "\tOmega                    : %16.6e\n", OMEGA);
  fprintf(fp, "\tEpsilon                  : %16.6e\n", EPS);
  fprintf(fp, "\n");

  fprintf(fp, "Time Control\n");
  fprintf(fp, "\tRestart                  : %s\n", (ICON==0)?"Initial":"Restart");
  fprintf(fp, "\tLast step                : %16d\n", NSTEP);
  fprintf(fp, "\tDelta T                  : %16.6e\n", DT);
  fprintf(fp, "\tInterval for Log         : %10d\n", Interval_Log);
  fprintf(fp, "\tInterval for Vis         : %10d\n", Interval_Vis);
  fprintf(fp, "\tInterval for Display     : %10d\n", Interval_Disp);
  fprintf(fp, "\n");
  fprintf(fp, "Application Control\n");
  fprintf(fp, "\tPerformance test         : %s\n", (PM_Test==OFF)?"OFF":"ON");
  fprintf(fp, "\n\n");
}



// #################################################################
/**
 * @brief パラメータのロード
 * @param [in] tp      TextParser
 * @note 無次元パラメータ
 * @retval 分割モード
 */
int RIAMC::setParameters(TextParser* tp)
{
  string label, str;
  double ct;
  int type_div=1; // デフォルト、分割数指定

  // 全計算領域の基点
  label = "/DomainInfo/GlobalOrigin";
  if ( !tp->getInspectedVector(label, G_origin, 3) )
  {
    Hostonly_ stamped_printf("\tERROR : in parsing '%s'\n", label.c_str());
    return -1;
  }

  // 全計算格子数
  label = "/DomainInfo/GlobalGrid";
  if ( !tp->getInspectedVector(label, G_size, 3) )
  {
    Hostonly_ stamped_printf("\tERROR : in parsing '%s'\n", label.c_str());
    return -1;
  }

  // 全計算領域の大きさ
  label = "/DomainInfo/GlobalRegion";
  if ( !tp->getInspectedVector(label, G_region, 3) )
  {
    Hostonly_ stamped_printf("\tERROR : in parsing '%s'\n", label.c_str());
    return -1;
  }

  if ( (G_region[0]>0.0) && (G_region[1]>0.0) && (G_region[2]>0.0) )
  {
    ; // skip
  }
  else
  {
    Hostonly_ cout << "ERROR : in parsing [" << label << "]" << endl;
    Exit(0);
  }


  // 領域分割数
  label = "/DomainInfo/GlobalDivision";
  if ( !tp->getInspectedVector(label, G_division, 3) )
  {
    type_div = 2; // 自動分割
  }

  // プロセス分割数が指定されている場合のチェック
  if ( type_div == 1 )
  {
    if ( (G_division[0]>0) && (G_division[1]>0) && (G_division[2]>0) )
    {
      Hostonly_ printf("\tManual domain division is selected.\n");
    }
    else
    {
      Hostonly_ cout << "ERROR : in parsing [" << label << "]" << endl;
      Exit(0);
    }
  }


  //GlobalPitch
  pitch[0] = G_region[0] / (double)(G_size[0]-1);
  pitch[1] = G_region[1] / (double)(G_size[1]-1);
  pitch[2] = G_region[2] / (double)(G_size[2]-1);



  // 物理パラメータ
  label = "/Parameter/Reynolds";
  if ( !(tp->getInspectedValue(label, ct)) )
  {
    RE = 1.0e+4;
  } else {
    RE = (REAL_TYPE)ct;
  }

  // 風上化の数値粘性係数
  label = "/Parameter/Alpha";
  if ( !(tp->getInspectedValue(label, ct)) )
  {
    ALPHA = 0.5;
  } else {
    ALPHA = (REAL_TYPE)ct;
  }

  // SMAGORINSKY CONSTANT
  label = "/Parameter/Smagorinsky_Constant";
  if ( !(tp->getInspectedValue(label, ct)) )
  {
    CS = 0.1;
  } else {
    CS = (REAL_TYPE)ct;
  }


  // 線形システムのソルバー
  label = "/LinearSystem/Solver";
  if ( !(tp->getInspectedValue(label, str )) )
  {
    Hostonly_ printf("\tParsing error : No '%s'\n", label.c_str());
    Exit(0);
  }
  if ( !strcasecmp(str.c_str(), "pointsor") )
  {
    LS_solver = LS_PSOR;
  }
  else if ( !strcasecmp(str.c_str(), "pointsor_metrics_array_free") )
  {
    LS_solver = LS_PSOR_MAF;
  }
  else if ( !strcasecmp(str.c_str(), "sor2sma") )
  {
    LS_solver = LS_SOR2SMA;
  }
  else if ( !strcasecmp(str.c_str(), "sor2sma_metrics_array_free") )
  {
    LS_solver = LS_SOR2SMA_MAF;
  }
  else if ( !strcasecmp(str.c_str(), "bicgstab") )
  {
    LS_solver = LS_BiCGSTAB;
  }
  else{
    Hostonly_ printf("\tInvalid solver : No '%s'\n", str.c_str());
    Exit(0);
  }

  // 線形システムのノルム
  label = "/LinearSystem/Norm";
  if ( !(tp->getInspectedValue(label, str )) )
  {
    Hostonly_ printf("\tParsing error : No '%s'\n", label.c_str());
    Exit(0);
  }
  if ( !strcasecmp(str.c_str(), "residual") )
  {
    LS_norm = Norm_Res;
  }
  else if ( !strcasecmp(str.c_str(), "rms") )
  {
    LS_norm = Norm_Rms;
  }
  else{
    Hostonly_ printf("\tInvalid norm : No '%s'\n", str.c_str());
    Exit(0);
  }


  // 最大反復数
  label = "/LinearSystem/IterationMax";
  if ( !(tp->getInspectedValue(label, ct)) )
  {
    ITER = 100;
  }
  else {
    ITER = (int)ct;
  }

  // SORの加速係数
  label = "/LinearSystem/Omega";
  if ( !(tp->getInspectedValue(label, ct)) )
  {
    OMEGA = 1.0;
  }
  else {
    OMEGA = (REAL_TYPE)ct;
  }

  // 収束閾値
  label = "/LinearSystem/Epsilon";
  if ( !(tp->getInspectedValue(label, ct)) )
  {
    EPS = 0.001;
  }
  else {
    EPS = (REAL_TYPE)ct;
  }






  // リスタート指定
  label = "/TimeControl/Restart";
  if ( !(tp->getInspectedValue(label, str)) )
  {
    ICON = 0;
  } else {
    if ( !strcasecmp(str.c_str(), "initial") )
    {
      ICON = 0;
    }
    else if ( !strcasecmp(str.c_str(), "restart") )
    {
      ICON = 1;
    }
    else
    {
      fprintf(stdout, "\tParsing error : Invalid keyword for '%s'\n", label.c_str());
      return false;
    }
  }

  // 計算ステップ数
  label = "/TimeControl/LastStep";
  if ( !(tp->getInspectedValue(label, ct)) )
  {
    NSTEP = 50000;
  } else {
    NSTEP = (unsigned)ct;
  }

  // 時間積分幅
  label = "/TimeControl/DeltaT";
  if ( !(tp->getInspectedValue(label, ct)) )
  {
    DT = 0.001;
  } else {
    DT = (REAL_TYPE)ct;
  }

  // 可視化データの出力インターバル
  label = "/TimeControl/Interval_Vis_Output";
  if ( !(tp->getInspectedValue(label, ct)) )
  {
    Interval_Vis = 500;
  } else {
    Interval_Vis = (int)ct;
  }

  //  ログファイルの出力インターバル
  label = "/TimeControl/Interval_Log_Output";
  if ( !(tp->getInspectedValue(label, ct)) )
  {
    Interval_Log = 5000;
  } else {
    Interval_Log = (int)ct;
  }

  //  ログファイルの出力インターバル
  label = "/TimeControl/Interval_Log_Display";
  if ( !(tp->getInspectedValue(label, ct)) )
  {
    Interval_Disp = 1;
  } else {
    Interval_Disp = (int)ct;
  }


  label = "/ApplicationControl/PM_Test";

  if ( tp->chkLabel(label) )
  {
    if ( tp->getInspectedValue(label, str) )
    {
      if     ( !strcasecmp(str.c_str(), "on") )  PM_Test = ON;
      else if( !strcasecmp(str.c_str(), "off") ) PM_Test = OFF;
      else
      {
        Hostonly_ stamped_printf("\tInvalid keyword is described for '%s'\n", label.c_str());
        Exit(0);
      }
    }
    else
    {
      Exit(0);
    }
  }

  return type_div;
}

// #################################################################
/**
 * @brief 立方体インデクスの作成
 * @retval true/false
 */
bool RIAMC::CubeIndex()
{
  int cube_index[6];    ///< 立方体のインデクス　ローカル

  // 立方体のインデクス計算
  int box_sz[3];

  box_sz[0] = G_size[0] - 1;
  box_sz[1] = G_size[1] - 1;
  box_sz[2] = G_size[2] - 1;

  // 立方体のインデクス　Fortran表現
  G_cube_index[0] = box_sz[0]*2/10+1;
  G_cube_index[1] = box_sz[1]*2/5+1;
  G_cube_index[2] = 1;
  G_cube_index[3] = G_cube_index[0] + 50;
  G_cube_index[4] = G_cube_index[1] + 50;
  G_cube_index[5] = G_cube_index[2] + 50;

  if (G_cube_index[3] > G_size[0]) {
    stamped_printf("\tError : index range x %d\n", G_cube_index[3]);
    return false;
  }
  if (G_cube_index[4] > G_size[1]) {
    stamped_printf("\tError : index range y %d\n", G_cube_index[4]);
    return false;
  }
  if (G_cube_index[5] > G_size[2]) {
    stamped_printf("\tError : index range z %d\n", G_cube_index[5]);
    return false;
  }

  Hostonly_ {
    printf("Global  :   st_i    ed_i    st_j    ed_j    st_k    ed_k  // in Fortran index\n");
    printf("\t   %5d   %5d   %5d   %5d   %5d   %5d\n",
       G_cube_index[0], G_cube_index[3],
       G_cube_index[1], G_cube_index[4],
       G_cube_index[2], G_cube_index[5] );
  }


  // C表現へ変換
  for (int i=0; i<6; i++) G_cube_index[i] -= 1;

  // ローカルインデクスに変換
  paraMngr->Global2LocalIndex(G_cube_index[0],
                              G_cube_index[1],
                              G_cube_index[2],
                              cube_index[0],
                              cube_index[1],
                              cube_index[2]);

  paraMngr->Global2LocalIndex(G_cube_index[3],
                              G_cube_index[4],
                              G_cube_index[5],
                              cube_index[3],
                              cube_index[4],
                              cube_index[5]);


  /*
   PATTERN
                Rank region
            <================>
      (1)   |                |    (2)
     <--->  |      (4)       |   <--->
      (3) <--->   <--->    <---> (5)
            |                |
       <------------------------->
                   (6)


  */

  int st_x, st_y, st_z, ed_x, ed_y, ed_z;

  st_x = cube_index[0];
  st_y = cube_index[1];
  st_z = cube_index[2];
  ed_x = cube_index[3];
  ed_y = cube_index[4];
  ed_z = cube_index[5];


  // 各ランクで頂点の有無を調べる
  int idx_[6] = {INT_MAX, INT_MAX, INT_MAX, INT_MIN, INT_MIN, INT_MIN};

  // do loopで始点より終点の値が小さくしておき、ループを実行しないことでBCをスキップする

  // pattern 1
  if( st_x < 0 && ed_x < 0 ) { idx_[0] = -1; idx_[3] = -2; }
  if( st_y < 0 && ed_y < 0 ) { idx_[1] = -1; idx_[4] = -2; }
  if( st_z < 0 && ed_z < 0 ) { idx_[2] = -1; idx_[5] = -2; }

  // pattern 2
  if( st_x >= size[0] && ed_x >= size[0] ) { idx_[0] = -1; idx_[3] = -2; }
  if( st_y >= size[1] && ed_y >= size[1] ) { idx_[1] = -1; idx_[4] = -2; }
  if( st_z >= size[2] && ed_z >= size[2] ) { idx_[2] = -1; idx_[5] = -2; }

  // pattern 3 & 6
  if( st_x < 0 && ed_x >= 0 ) {
    idx_[0] = 0;
    idx_[3] = min(ed_x, size[0]-1);
  }
  if( st_y < 0 && ed_y >= 0 ) {
    idx_[1] = 0;
    idx_[4] = min(ed_y, size[1]-1);
  }
  if( st_z < 0 && ed_z >= 0 ) {
    idx_[2] = 0;
    idx_[5] = min(ed_z, size[2]-1);
  }

  // pattern 5 & 6
  if( ed_x >= size[0] && st_x < size[0] ) {
    idx_[0] = max(0, st_x);
    idx_[3] = size[0]-1;
  }
  if( ed_y >= size[1] && st_y < size[1] ) {
    idx_[1] = max(0, st_y);
    idx_[4] = size[1]-1;
  }
  if( ed_z >= size[2] && st_z < size[2] ) {
    idx_[2] = max(0, st_z);
    idx_[5] = size[2]-1;
  }

  // pattern 4
  if ( st_x >= 0 && ed_x < size[0] ) { idx_[0] = st_x; idx_[3] = ed_x; }
  if ( st_y >= 0 && ed_y < size[1] ) { idx_[1] = st_y; idx_[4] = ed_y; }
  if ( st_z >= 0 && ed_z < size[2] ) { idx_[2] = st_z; idx_[5] = ed_z; }

  // x, y, z の終点インデクスのどれかに-2の値があれば、領域外
  if (idx_[3] == -2 || idx_[4] == -2 || idx_[5] == -2) {
    idx_[0] = idx_[1] = idx_[2] = -1;
    idx_[3] = idx_[4] = idx_[5] = -2;
  }


  // Fortran index に変換 [1スタート]
  for (int i=0; i<6; i++) idx_[i] += 1;



  // ローカル情報をマスターに集約
  int* tmp = new int[6*numProc];
  for (int i=0; i<6*numProc; i++) tmp[i] = -1;

  // Fortran indexでを保持
  for (int i=0; i<6; i++) bcf[i] = idx_[i];

  // 集約
  for (int m=0; m<numProc; m++)
  {
    if ( paraMngr->Gather(idx_, 6, tmp, 6, 0, procGrp) != CPM_SUCCESS ) return false;
  }

  Hostonly_ {
    printf("\n\t Rank    st_i    ed_i    st_j    ed_j    st_k    ed_k  // in Fortran index\n");
    for (int i=0; i<numProc; i++)
    {
      printf("\t%5d   %5d   %5d   %5d   %5d   %5d   %5d\n",
             i,
             tmp[6*i+0], tmp[6*i+3],
             tmp[6*i+1], tmp[6*i+4],
             tmp[6*i+2], tmp[6*i+5] );
    }
    printf("\n");
  }

  if (tmp)
  {
    delete[] tmp;
    tmp = NULL;
  }

  return true;
}


// #################################################################
/**
 * @brief 領域情報の集約と表示
 * @retval true/false
 */
bool RIAMC::displayDomainInfo()
{
  int* m_xyz = new int[3*numProc];
  for (int i=0; i<3*numProc; i++) m_xyz[i] = -1;

  int* m_hd = new int[3*numProc];
  for (int i=0; i<3*numProc; i++) m_hd[i] = -1;

  // 集約
  for (int m=0; m<numProc; m++)
  {
    if ( paraMngr->Gather(size, 3, m_xyz, 3, 0, procGrp) != CPM_SUCCESS ) return false;
    if ( paraMngr->Gather(head, 3, m_hd,  3, 0, procGrp) != CPM_SUCCESS ) return false;
  }

  FILE* fp = NULL;

  Hostonly_
  {
    if ( !(fp=fopen("domaininfo.txt", "w")) )
    {
      stamped_printf("\tSorry, can't open 'domaininfo.txt' file. Write failed.\n");
      return false;
    }
  }

  // head[]は各ランクの開始格子点のグローバルインデクス（0スタート）
  Hostonly_ {
    fprintf(fp, "    Grid size  :  %5d   %5d   %5d\n", G_size[0], G_size[1], G_size[2]);
    fprintf(fp, "     Division  :  %5d   %5d   %5d\n", G_division[0], G_division[1], G_division[2]);

    fprintf(fp, "\n\t Rank  :   sz_x    sz_y    sz_z  :   hd_x    hd_y    hd_z\n");
    for (int i=0; i<numProc; i++)
    {
      fprintf(fp, "\t%5d  :  %5d   %5d   %5d  :  %5d   %5d   %5d\n",
             i,
             m_xyz[3*i+0],
             m_xyz[3*i+1],
             m_xyz[3*i+2],
             m_hd [3*i+0],
             m_hd [3*i+1],
             m_hd [3*i+2] );
    }
    printf("\n");

    if ( fp ) fclose(fp);
  }

  if (m_xyz)
  {
    delete[] m_xyz;
    m_xyz = NULL;
  }

  if (m_hd)
  {
    delete[] m_hd;
    m_hd = NULL;
  }

  return true;
}
