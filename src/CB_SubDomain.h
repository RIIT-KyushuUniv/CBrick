#ifndef _CB_SUBDOMAIN_H_
#define _CB_SUBDOMAIN_H_

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

/**
 * @file   CB_SubDomain.h
 * @brief  SubDomain class Header
 * @todo sort >> k>j スレッド化のためのオプション、オプションの与え方
 */

#include <mpi.h>
#include <string>
#include <stdlib.h>
#include "CB_Define.h"
#include "CB_Version.h"

// ワーク用の構造体
typedef struct {
  int sz[3];  ///< サブドメインのサイズ
  float srf;  ///< 通信量
  float sxy;  ///< 立方体への近さを表す自乗量
} score_tbl;


/****************************************************
 * 分割情報クラス （テンポラリ利用）
 */
class cntl_tbl {
public:
  int dsz[3];       ///< サブドメインの基準サイズ
  int mod[3];       ///< 基準サイズの個数
  int div[3];       ///< 全計算領域の各方向の分割数
  int  org_idx;     ///< セットした最初のインデクスを指す
  float sc_vol;     ///< 評価値：計算量のバランス
  float sc_com;     ///< 評価値：通信量
  float sc_len;     ///< 評価値：X方向長さ
  float sc_hex;     ///< 評価値：立方体度
  score_tbl* score;

  cntl_tbl() {
    for (int i=0; i<3; i++) {
      dsz[i] = 0;
      mod[i] = 0;
      div[i] = 0;
    }
    sc_vol = sc_com = sc_len = sc_hex = 0.0;
    org_idx = -1;
    score = NULL;
  }

  ~cntl_tbl() {
    if ( score ) delete [] score;
    score = NULL;
  }

  // copy constructor
  cntl_tbl(const cntl_tbl& src) {
    for (int i=0; i<3; i++) {
      dsz[i] = src.dsz[i];
      mod[i] = src.mod[i];
      div[i] = src.div[i];
    }
    sc_vol = src.sc_vol;
    sc_com = src.sc_com;
    sc_len = src.sc_len;
    sc_hex = src.sc_hex;
    score  = src.score;   // アドレスコピー、ポイント先の内容は変更しない
    org_idx= src.org_idx;
  }
};


/****************************************************
 * サブドメイン情報クラス （テンポラリ利用）
 */
class SubdomainInfo {
public:
  int sz[3];    ///< サブドメインの要素数
  int hd[3];    ///< サブドメインのヘッドインデクス（グローバル）
  int cm[6];    ///< 隣接ランクID

  SubdomainInfo() {
    for (int i=0; i<3; i++) {
      sz[i] = 0;
      hd[i] = -1;
    }
    for (int i=0; i<6; i++) cm[i] = -1;
  }

  ~SubdomainInfo() {}

  // copy constructor
  SubdomainInfo(const SubdomainInfo& src) {
    for (int i=0; i<3; i++) {
      sz[i] = src.sz[i];
      hd[i] = src.hd[i];
    }
    for (int i=0; i<6; i++) {
      cm[i] = src.cm[i];
    }
  }
};


/****************************************************
 * サブドメイン情報保持クラス
 */
class SubDomain {
public:
  int procGrp;          ///< プロセスグループ番号
  int myRank;           ///< 自ノードのランク番号
  int numProc;          ///< 全ランク数
  MPI_Comm mpi_comm;    ///< MPI コミュニケーター

  int G_div[3];         ///< 各軸方向の領域分割数
  int G_size[3];        ///< 全領域の要素数 (Global, Non-dimensional)
  int size[3];          ///< 各サブドメインの要素数 (Local, Non-dimensional
  int head[3];          ///< 開始インデクス（グローバルインデクス）
  int comm_tbl[6];      ///< 隣接ブロックのランク番号
  int halo_width;       ///< ガイドセル幅
  int ranking_opt;      ///< ランキングのオプション　（0=cubical, default, 1=vector）
  int div_mode;         ///< 分割モード (0=自動、1=指定)
  int f_index;          ///< Findex (0-OFF, 1-ON) @note 関連するところは head index
  int use_NB;           ///< ノンブロッキング通信の利用　（0-no_use, 1-use）


protected:
  std::string grid_type;///< "cell" or "node"
  SubdomainInfo* sd;    ///< 全サブドメインの分割要素数配列

private:
  // Buffer for asynchronous communication
  REAL_TYPE* f_xms;  // X- direction send
  REAL_TYPE* f_xmr;  // X- direction recv
  REAL_TYPE* f_xps;  // X+ direction send
  REAL_TYPE* f_xpr;  // X+ direction recv
  REAL_TYPE* f_yms;  // Y- direction send
  REAL_TYPE* f_ymr;  // Y- direction recv
  REAL_TYPE* f_yps;  // Y+ direction send
  REAL_TYPE* f_ypr;  // Y+ direction recv
  REAL_TYPE* f_zms;  // Z- direction send
  REAL_TYPE* f_zmr;  // Z- direction recv
  REAL_TYPE* f_zps;  // Z+ direction send
  REAL_TYPE* f_zpr;  // Z+ direction recv


public:
  // デフォルト　コンストラクタ
  SubDomain() {
    procGrp = -1;
    myRank  = -1;
    numProc = 1;
    halo_width = 0;
    ranking_opt = 0;
    div_mode = 0;
    f_index = 0;
    use_NB = 0;

    for (int i=0; i<NOFACE; i++) comm_tbl[i] = -1;

    for (int i=0; i<3; i++) {
      head[i]       = 0;
      size[i]       = 0;
      G_size[i]     = 0;
      G_div[i]      = 0;
    }
    sd = NULL;

    f_xms = NULL;  // X- direction send
    f_xmr = NULL;  // X- direction recv
    f_xps = NULL;  // X+ direction send
    f_xpr = NULL;  // X+ direction recv
    f_yms = NULL;  // Y- direction send
    f_ymr = NULL;  // Y- direction recv
    f_yps = NULL;  // Y+ direction send
    f_ypr = NULL;  // Y+ direction recv
    f_zms = NULL;  // Z- direction send
    f_zmr = NULL;  // Z- direction recv
    f_zps = NULL;  // Z+ direction send
    f_zpr = NULL;  // Z+ direction recv
  }


  // 領域サイズとプロセス数
  SubDomain(int m_gsz[],
            int m_halo,
            int m_np,
            int m_myrank,
            int m_procgrp,
            MPI_Comm m_comm,
            std::string m_type,
            std::string m_idxtyp,
            int m_useNB,
            int priority=0) {

    this->G_size[0]   = m_gsz[0];
    this->G_size[1]   = m_gsz[1];
    this->G_size[2]   = m_gsz[2];
    this->halo_width  = m_halo;
    this->numProc     = m_np;
    this->grid_type   = m_type;
    this->procGrp     = m_procgrp;
    this->myRank      = m_myrank;
    this->mpi_comm    = m_comm;
    this->ranking_opt = priority;
    this->use_NB      = m_useNB;

    if (m_type == "node" || m_type == "cell") {
      // ok
    }
    else {
      printf("Error : Invalid grid type [%s]\n", m_type.c_str());
      Exit(-1);
    }

    if (m_idxtyp == "Findex") {
      f_index = 1;
    }
    else if (m_idxtyp == "Cindex") {
      f_index = 0;
    }
    else {
      printf("Error : Invalid Index type [%s]\n", m_idxtyp.c_str());
      Exit(-1);
    }

    if (m_halo<0) Exit(-1);

    if ( !(sd = new SubdomainInfo[numProc]) ) {
      printf("\tFail to allocate memory\n");
      Exit(-1);
    }
  }

  /** デストラクタ */
  virtual ~SubDomain() {
    if ( sd ) delete [] sd;

    if ( f_xms ) delete [] f_xms;
    if ( f_xmr ) delete [] f_xmr;
    if ( f_xps ) delete [] f_xps;
    if ( f_xpr ) delete [] f_xpr;
    if ( f_yms ) delete [] f_yms;
    if ( f_ymr ) delete [] f_ymr;
    if ( f_yps ) delete [] f_yps;
    if ( f_ypr ) delete [] f_ypr;
    if ( f_zms ) delete [] f_zms;
    if ( f_zmr ) delete [] f_zmr;
    if ( f_zps ) delete [] f_zps;
    if ( f_zpr ) delete [] f_zpr;
  }

// CB_SubDomain.cpp
public:

  // @brief 通信テーブルを作成
  bool createRankTable();

  // @brief 最適な分割数を見つける
  bool findOptimalDivision();

  // @brief 　Global > Local　インデクス変換
  bool G2L_index(const int* Gi, int* Li);

  /*
   * @brief 分割数をセットする
   * @param [in] m_dv   分割数
   * @note G_div[0]*G_div[1]*G_div[2] != numProc を事前にチェックのこと
   */
  bool setDivision(const int* m_dv)
  {
    G_div[0] = m_dv[0];
    G_div[1] = m_dv[1];
    G_div[2] = m_dv[2];

    if ( G_div[0]>G_size[0] && G_div[1]>G_size[1] && G_div[2]>G_size[2] ) {
      printf("\tG_div[] size is our of range.\n");
      return false;
    }

    div_mode = 1;

    return true;
  }

  bool setSubDomain(int m_gsz[],
                    int m_halo,
                    int m_np,
                    int m_myrank,
                    int m_procgrp,
                    MPI_Comm m_comm,
                    std::string m_type,
                    std::string m_idxtyp,
                    int m_useNB,
                    int priority=0)
  {
    G_size[0]   = m_gsz[0];
    G_size[1]   = m_gsz[1];
    G_size[2]   = m_gsz[2];
    halo_width  = m_halo;
    numProc     = m_np;
    grid_type   = m_type;
    procGrp     = m_procgrp;
    myRank      = m_myrank;
    mpi_comm    = m_comm;
    ranking_opt = priority;
    use_NB      = m_useNB;

    if (m_type == "node" || m_type == "cell") {
      // ok
    }
    else {
      printf("Error : Invalid grid type [%s]\n", m_type.c_str());
      return false;
    }

    if (m_idxtyp == "Findex") {
      f_index = 1;
    }
    else if (m_idxtyp == "Cindex") {
      f_index = 0;
    }
    else {
      printf("Error : Invalid Index type [%s]\n", m_idxtyp.c_str());
      return false;
    }

    if (m_halo<0) Exit(-1);

    if ( !(sd = new SubdomainInfo[numProc]) ) {
      printf("\tFail to allocate memory\n");
      return false;
    }

    return true;
  }

protected:

  // @brief 開始インデクス=0の 3D=>1D インデクス変換、ガイドセルは考慮しない場合
  inline int rank_idx_0(const int _I, const int _J, const int _K, const int _NI, const int _NJ)
  {
    return _K * _NI * _NJ + _J * _NI + _I;
  }

  void Evaluation(cntl_tbl* t, const int tbl_sz, FILE* fp);

  bool findParameter();

  void getHeadIndex();

  int getNumCandidates();

  void getSizeCell(cntl_tbl* t, const int* in, const int m);

  void getSizeNode(score_tbl* t);

  void getSrf(score_tbl* t);

  void registerCandidates(cntl_tbl* t);

  int sortComm(cntl_tbl* t, const int c_sz, FILE* fp);

  int sortCube(cntl_tbl* t, const int c_sz, FILE* fp);

  int sortLenX(cntl_tbl* t, const int c_sz, FILE* fp);

  int sortVolume(cntl_tbl* t, const int tbl_sz, FILE* fp);



// CB_CommS.cpp
public:

  // @brief 通信バッファの確保
  bool initComm();

  // @brief スカラー変数のブロッキング通信
  // @param [in,out]  src     スカラー変数
  // @param [in]      gc_comm 実際に通信する通信面数
  // @retval true-success, false-fail
  bool Comm_S_blocking(REAL_TYPE* src, const int gc_comm);

  // @brief スカラー変数のノンブロッキング通信
  // @param [in,out]  src     スカラー変数
  // @param [in]      gc_comm 実際に通信する通信面数
  // @retval true-success, false-fail
  bool Comm_S_nonblocking(REAL_TYPE* src,
                          const int gc_comm,
                          MPI_Request *req);

  // @brief スカラー変数のノンブロッキング通信
  // @param [in,out]  dest    スカラー変数
  // @param [in]      gc_comm 実際に通信する通信面数
  // @param [out]     req     Array of MPI request
  // @retval true-success, false-fail
  bool Comm_S_wait_nonblocking(REAL_TYPE* dest,
                               const int gc_comm,
                               MPI_Request *req);


protected:

  bool send_and_recv(REAL_TYPE* ms,
                     REAL_TYPE* mr,
                     REAL_TYPE* ps,
                     REAL_TYPE* pr,
                     int msz,
                     int nIDm,
                     int nIDp);

  bool sendrecv(REAL_TYPE* ms,
                REAL_TYPE* mr,
                REAL_TYPE* ps,
                REAL_TYPE* pr,
                int msz,
                int nIDm,
                int nIDp);

  void packX(const REAL_TYPE *array,
             const int vc_comm,
             REAL_TYPE *sendm,
             REAL_TYPE *sendp,
             const int nIDm,
             const int nIDp);

  void unpackX(REAL_TYPE *array,
               const int vc_comm,
               const REAL_TYPE *recvm,
               const REAL_TYPE *recvp,
               const int nIDm,
               const int nIDp);

  void packY(const REAL_TYPE *array,
             const int vc_comm,
             REAL_TYPE *sendm,
             REAL_TYPE *sendp,
             const int nIDm,
             const int nIDp);

  void unpackY(REAL_TYPE *array,
               const int vc_comm,
               const REAL_TYPE *recvm,
               const REAL_TYPE *recvp,
               const int nIDm,
               const int nIDp);

  void packZ(const REAL_TYPE *array,
             const int vc_comm,
             REAL_TYPE *sendm,
             REAL_TYPE *sendp,
             const int nIDm,
             const int nIDp);

  void unpackZ(REAL_TYPE *array,
               const int vc_comm,
               const REAL_TYPE *recvm,
               const REAL_TYPE *recvp,
               const int nIDm,
               const int nIDp);

  bool IsendIrecv(REAL_TYPE* ms,
                  REAL_TYPE* mr,
                  REAL_TYPE* ps,
                  REAL_TYPE* pr,
                  int msz,
                  int nIDm,
                  int nIDp,
                  MPI_Request *req);


// CB_CommV.cpp

  /*
   * @fn Comm_V_blocking
   * @brief ベクトル変数（i,j,k,l）型のブロッキング通信
   * @param [in,out]  src     ベクトル変数
   * @param [in]      gc_comm 実際に通信する通信面数
   * @retval true-success, false-fail
   */
  bool Comm_V_blocking(REAL_TYPE* src, const int gc_comm);

};

#endif // _CB_SUBDOMAIN_H_
