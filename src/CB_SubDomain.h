#ifndef _CB_SUBDOMAIN_H_
#define _CB_SUBDOMAIN_H_

/*
###################################################################################
#
# CBrick
#
# Copyright (c) 2017-2018 Research Institute for Information Technology(RIIT),
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

// 分割モード auto_div
#define AUTO 0
#define SPEC 1

// ワーク用の構造体
typedef struct {
  int sz[3];  ///< サブドメインのサイズ
  float srf;  ///< 通信量
  float sxy;  ///< 立方体への近さを表す自乗量
} score_tbl;


/****************************************************
 * 分割情報クラス
 */
class cntl_tbl {
private:
  int mode;         ///< インスタンスのモード

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

  // デフォルトコンストラクタ
  cntl_tbl() {
    for (int i=0; i<3; i++) {
      dsz[i] = 0;
      mod[i] = 0;
      div[i] = 0;
    }
    sc_vol = sc_com = sc_len = sc_hex = 0.0;
    org_idx = -1;
    mode = 0;
    score = NULL;
  }

  // @param [in] num_array 評価するパラメータを保持する配列数
  cntl_tbl(const int num_array) {
    for (int i=0; i<3; i++) {
      dsz[i] = 0;
      mod[i] = 0;
      div[i] = 0;
    }
    sc_vol = sc_com = sc_len = sc_hex = 0.0;
    org_idx = -1;
    mode = 1;

    if (num_array < 1) exit(184);
    score = new score_tbl[num_array];
  }

  // コンストラクタ cntl_tbl p = src; からの自動変換に利用
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
    org_idx= src.org_idx;
    score  = src.score;   // アドレスのみコピー、ポイント先の内容は変更しない
    mode = 2;
  }

  virtual ~cntl_tbl() {
    if ( mode == 1 ) delete [] score;
  }
};


/****************************************************
 * サブドメイン情報クラス
 */
class SubdomainInfo {
public:
  int sz[3];      ///< サブドメインの要素数
  int hd[3];      ///< サブドメインのヘッドインデクス（グローバル）
  int cm[NOFACE]; ///< 隣接ランクID

  SubdomainInfo() {
    for (int i=0; i<3; i++) {
      sz[i] = 0;
      hd[i] = -1;
    }
    for (int i=0; i<NOFACE; i++) cm[i] = -1;
  }

  virtual ~SubdomainInfo() {}

  // copy constructor
  SubdomainInfo(const SubdomainInfo& src) {
    for (int i=0; i<3; i++) {
      sz[i] = src.sz[i];
      hd[i] = src.hd[i];
    }

    for (int i=0; i<NOFACE; i++) {
      cm[i] = src.cm[i];
    }
  }
};


/****************************************************
 * サブドメイン情報保持クラス
 */
class SubDomain {

protected:
  int G_div[3];         ///< 各軸方向の領域分割数
  int G_size[3];        ///< 全領域の要素数 (Global, Non-dimensional)
  int size[3];          ///< 各サブドメインの要素数 (Local, Non-dimensional
  int head[3];          ///< 開始インデクス（グローバルインデクス）
  int comm_tbl[NOFACE]; ///< 隣接ブロックのランク番号

  std::string grid_type;///< "cell" or "node"
  SubdomainInfo* sd;    ///< 全サブドメインの分割要素数配列

private:
  MPI_Comm mpi_comm;    ///< MPI コミュニケーター
  int procGrp;          ///< プロセスグループ番号
  int myRank;           ///< 自ノードのランク番号
  int halo_width;       ///< ガイドセル幅
  int auto_div;         ///< 分割モード (AUTO, SPEC)
  int f_index;          ///< Findex (0-OFF, 1-ON) @note 関連するところは head index
  int numProc;          ///< 全ランク数
  int ranking_opt;      ///< ランキングのオプション（0=cubical, default, 1=vector）

  // Buffer for asynchronous communication
  REAL_TYPE* f_ims;  // I- direction send
  REAL_TYPE* f_imr;  // I- direction recv
  REAL_TYPE* f_ips;  // I+ direction send
  REAL_TYPE* f_ipr;  // I+ direction recv
  REAL_TYPE* f_jms;  // J- direction send
  REAL_TYPE* f_jmr;  // J- direction recv
  REAL_TYPE* f_jps;  // J+ direction send
  REAL_TYPE* f_jpr;  // J+ direction recv
  REAL_TYPE* f_kms;  // K- direction send
  REAL_TYPE* f_kmr;  // K- direction recv
  REAL_TYPE* f_kps;  // K+ direction send
  REAL_TYPE* f_kpr;  // K+ direction recv
#ifdef _DIAGONAL_COMM
  REAL_TYPE* f_es;   // edge send
  REAL_TYPE* f_er;   // edge recv
  REAL_TYPE* f_cs;   // corner send
  REAL_TYPE* f_cr;   // corner recv
#endif

  int buf_flag;     // バッファを確保済みのときに1


public:
  // デフォルト コンストラクタ
  SubDomain() {
    procGrp = -1;
    myRank  = -1;
    numProc = 1;
    halo_width = 0;
    ranking_opt = 0;
    auto_div = AUTO;
    f_index = 0;

    for (int i=0; i<NOFACE; i++) comm_tbl[i] = -1;

    for (int i=0; i<3; i++) {
      head[i]       = 0;
      size[i]       = 0;
      G_size[i]     = 0;
      G_div[i]      = 0;
    }
    sd = NULL;

    f_ims = NULL;  // X- direction send
    f_imr = NULL;  // X- direction recv
    f_ips = NULL;  // X+ direction send
    f_ipr = NULL;  // X+ direction recv
    f_jms = NULL;  // Y- direction send
    f_jmr = NULL;  // Y- direction recv
    f_jps = NULL;  // Y+ direction send
    f_jpr = NULL;  // Y+ direction recv
    f_kms = NULL;  // Z- direction send
    f_kmr = NULL;  // Z- direction recv
    f_kps = NULL;  // Z+ direction send
    f_kpr = NULL;  // Z+ direction recv
#ifdef _DIAGONAL_COMM
    f_es  = NULL;  // edge send
    f_er  = NULL;  // edge recv
    f_cs  = NULL;  // corner send
    f_cr  = NULL;  // corner recv
#endif
    buf_flag = 0;
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

    if ( numProc < 1 ) Exit(-1);
    if ( !(sd=allocateSD()) ) Exit(-1);
  }


  /** デストラクタ */
  virtual ~SubDomain() {
    if ( buf_flag == 1 ) {
      delete [] f_ims;
      delete [] f_imr;
      delete [] f_ips;
      delete [] f_ipr;
      delete [] f_jms;
      delete [] f_jmr;
      delete [] f_jps;
      delete [] f_jpr;
      delete [] f_kms;
      delete [] f_kmr;
      delete [] f_kps;
      delete [] f_kpr;
#ifdef _DIAGONAL_COMM
      delete [] f_es;
      delete [] f_er;
      delete [] f_cs;
      delete [] f_cr;
#endif
    }
  }

// CB_SubDomain.cpp
public:

  // @brief 通信テーブルを作成
  bool createRankTable();

  // @brief 最適な分割数を見つける
  // @param [in] mode {0-IJK分割、デフォルト、1-JK分割}
  bool findOptimalDivision(int terrain_mode=0);

  // @brief Global > Localインデクス変換
  bool G2L_index(const int* Gi, int* Li);

  // @brief Global > Localインデクス変換 コンポーネントのみ
  bool G2L_index(const int Gi, int& Li, const int c);

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

    if ( G_div[0]>G_size[0] || G_div[1]>G_size[1] || G_div[2]>G_size[2] ) {
      Hostonly_ printf("\nERROR :  G_div[] value is greater than G_size[].\n\n");
      return false;
    }

    auto_div = SPEC;

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

    if (m_type == "node" || m_type == "cell") {
      // ok
    }
    else {
      Hostonly_ printf("Error : Invalid grid type [%s]\n", m_type.c_str());
      return false;
    }

    if (m_idxtyp == "Findex") {
      f_index = 1;
    }
    else if (m_idxtyp == "Cindex") {
      f_index = 0;
    }
    else {
      Hostonly_ printf("Error : Invalid Index type [%s]\n", m_idxtyp.c_str());
      return false;
    }

    if (m_halo<0) Exit(-1);

    if ( numProc < 1 ) return false;
    if ( !(sd=allocateSD()) ) return false;

    return true;
  }


  // @brief 領域の分割数を返す
  // @param [out] m_sz 分割数
  void getGlobalDivision(int* m_sz)
  {
    m_sz[0] = G_div[0];
    m_sz[1] = G_div[1];
    m_sz[2] = G_div[2];
  }

  // @brief 全計算領域の要素数を返す
  // @param [out] m_sz 要素数
  void getGlobalSize(int* m_sz)
  {
    m_sz[0] = G_size[0];
    m_sz[1] = G_size[1];
    m_sz[2] = G_size[2];
  }

  // @brief 部分領域の要素数を返す
  // @param [out] m_sz 要素数
  void getLocalSize(int* m_sz)
  {
    m_sz[0] = size[0];
    m_sz[1] = size[1];
    m_sz[2] = size[2];
  }

  // @brief 自領域の先頭インデクスを返す
  // @param [out] m_sz ランクmのhead[]
  void getLocalHead(int* m_sz)
  {
    m_sz[0] = head[0];
    m_sz[1] = head[1];
    m_sz[2] = head[2];
  }

  // @brief 通信テーブルを返す
  // @param [out] m_tbl 通信テーブル
  void getCommTable(int* m_tbl)
  {
    for (int i=0; i<NOFACE; i++) m_tbl[i] = comm_tbl[i];
  }

  // @brief subdomainの分割パラメータ領域サイズを返す
  // @param [in] m     sdのインデクス（ランク番号に相当）
  // @param [out] m_sz ランクmの分割数
  void getSubDomainSize(const int m, int m_sz[3])
  {
    m_sz[0] = sd[m].sz[0];
    m_sz[1] = sd[m].sz[1];
    m_sz[2] = sd[m].sz[2];
  }

  // @brief subdomainの分割パラメータ領域ヘッダを返す
  // @param [in] m     sdのインデクス（ランク番号に相当）
  // @param [out] m_sz ランクmのhead[]
  void getSubDomainHead(const int m, int m_sz[3])
  {
    m_sz[0] = sd[m].hd[0];
    m_sz[1] = sd[m].hd[1];
    m_sz[2] = sd[m].hd[2];
  }


private:

  // @brief SubDomainInfoクラスの確保
  SubdomainInfo* allocateSD()
  {
    return new SubdomainInfo[numProc];
  }

  void Evaluation(cntl_tbl* t, const int tbl_sz, FILE* fp);

  bool findParameter();

  void getHeadIndex();

  int getNumCandidates();

  int getNumCandidates4JK();

  void getSize(cntl_tbl* t, const int* in, const int m);

  void getSizeNode(score_tbl* t);

  void getSrf(score_tbl* t);

  void registerCandidates_Cell(cntl_tbl* t);

  void registerCandidates_Node(cntl_tbl* t);

  void registerCandidates4JK_Cell(cntl_tbl* tbl);

  void registerCandidates4JK_Node(cntl_tbl* tbl);

  int sortComm(cntl_tbl* t, const int c_sz, FILE* fp);

  int sortCube(cntl_tbl* t, const int c_sz, FILE* fp);

  int sortLenX(cntl_tbl* t, const int c_sz, FILE* fp);

  int sortVolume(cntl_tbl* t, const int tbl_sz, FILE* fp);



// CB_CommS.cpp
public:

  // @brief 通信バッファの確保
  // @param [in] gnum_compo バッファで利用する最大の層数（1-scalar, 3-vector, ?-others）
  bool initComm(const int num_compo);


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

private:

  bool IsendIrecv(REAL_TYPE* ms,
                  REAL_TYPE* mr,
                  REAL_TYPE* ps,
                  REAL_TYPE* pr,
                  int msz,
                  int nIDm,
                  int nIDp,
                  MPI_Request *req);


// CB_PackingScalarCell.cpp
private:

  void pack_SXcell(const REAL_TYPE *array,
               const int vc_comm,
               REAL_TYPE *sendm,
               REAL_TYPE *sendp,
               const int nIDm,
               const int nIDp);

  void unpack_SXcell(REAL_TYPE *array,
                 const int vc_comm,
                 const REAL_TYPE *recvm,
                 const REAL_TYPE *recvp,
                 const int nIDm,
                 const int nIDp);

  void pack_SYcell(const REAL_TYPE *array,
               const int vc_comm,
               REAL_TYPE *sendm,
               REAL_TYPE *sendp,
               const int nIDm,
               const int nIDp);

  void unpack_SYcell(REAL_TYPE *array,
                 const int vc_comm,
                 const REAL_TYPE *recvm,
                 const REAL_TYPE *recvp,
                 const int nIDm,
                 const int nIDp);

  void pack_SZcell(const REAL_TYPE *array,
               const int vc_comm,
               REAL_TYPE *sendm,
               REAL_TYPE *sendp,
               const int nIDm,
               const int nIDp);

  void unpack_SZcell(REAL_TYPE *array,
                 const int vc_comm,
                 const REAL_TYPE *recvm,
                 const REAL_TYPE *recvp,
                 const int nIDm,
                 const int nIDp);

  bool pack_SEcell(REAL_TYPE *array,
               const int vc_comm,
               REAL_TYPE *sendbuf,
               REAL_TYPE *recvbuf,
               MPI_Request *req);

  void unpack_SEcell(REAL_TYPE *array,
                 const int vc_comm,
                 const REAL_TYPE *recvbuf);

  bool pack_SCcell(REAL_TYPE *array,
               const int vc_comm,
               REAL_TYPE *sendbuf,
               REAL_TYPE *recvbuf,
               MPI_Request *req);

  void unpack_SCcell(REAL_TYPE *array,
                 const int vc_comm,
                 const REAL_TYPE *recvbuf);



  // CB_PackingScalarNode.cpp
  private:

    void pack_SXnode(const REAL_TYPE *array,
                 const int vc_comm,
                 REAL_TYPE *sendm,
                 REAL_TYPE *sendp,
                 const int nIDm,
                 const int nIDp);

    void unpack_SXnode(REAL_TYPE *array,
                   const int vc_comm,
                   const REAL_TYPE *recvm,
                   const REAL_TYPE *recvp,
                   const int nIDm,
                   const int nIDp);

    void pack_SYnode(const REAL_TYPE *array,
                 const int vc_comm,
                 REAL_TYPE *sendm,
                 REAL_TYPE *sendp,
                 const int nIDm,
                 const int nIDp);

    void unpack_SYnode(REAL_TYPE *array,
                   const int vc_comm,
                   const REAL_TYPE *recvm,
                   const REAL_TYPE *recvp,
                   const int nIDm,
                   const int nIDp);

    void pack_SZnode(const REAL_TYPE *array,
                 const int vc_comm,
                 REAL_TYPE *sendm,
                 REAL_TYPE *sendp,
                 const int nIDm,
                 const int nIDp);

    void unpack_SZnode(REAL_TYPE *array,
                   const int vc_comm,
                   const REAL_TYPE *recvm,
                   const REAL_TYPE *recvp,
                   const int nIDm,
                   const int nIDp);

  bool pack_SEnode(REAL_TYPE *array,
               const int vc_comm,
               REAL_TYPE *sendbuf,
               REAL_TYPE *recvbuf,
               MPI_Request *req);

  void unpack_SEnode(REAL_TYPE *array,
                 const int vc_comm,
                 const REAL_TYPE *recvbuf);

  bool pack_SCnode(REAL_TYPE *array,
               const int vc_comm,
               REAL_TYPE *sendbuf,
               REAL_TYPE *recvbuf,
               MPI_Request *req);

  void unpack_SCnode(REAL_TYPE *array,
                 const int vc_comm,
                 const REAL_TYPE *recvbuf);


// CB_CommV.cpp
public:

  // @brief ベクトル変数のノンブロッキング通信
  // @param [in,out]  src     ベクトル変数
  // @param [in]      gc_comm 実際に通信する通信面数
  // @retval true-success, false-fail
  bool Comm_V_nonblocking(REAL_TYPE* src,
                          const int gc_comm,
                          MPI_Request *req);

  // @brief ベクトル変数のノンブロッキング通信
  // @param [in,out]  dest    ベクトル変数
  // @param [in]      gc_comm 実際に通信する通信面数
  // @param [out]     req     Array of MPI request
  // @retval true-success, false-fail
  bool Comm_V_wait_nonblocking(REAL_TYPE* dest,
                               const int gc_comm,
                               MPI_Request *req);


// CB_PackingVectorCell.cpp
private:

  void pack_VXcell(const REAL_TYPE *array,
               const int vc_comm,
               REAL_TYPE *sendm,
               REAL_TYPE *sendp,
               const int nIDm,
               const int nIDp);

  void unpack_VXcell(REAL_TYPE *array,
                 const int vc_comm,
                 const REAL_TYPE *recvm,
                 const REAL_TYPE *recvp,
                 const int nIDm,
                 const int nIDp);

  void pack_VYcell(const REAL_TYPE *array,
               const int vc_comm,
               REAL_TYPE *sendm,
               REAL_TYPE *sendp,
               const int nIDm,
               const int nIDp);

  void unpack_VYcell(REAL_TYPE *array,
                 const int vc_comm,
                 const REAL_TYPE *recvm,
                 const REAL_TYPE *recvp,
                 const int nIDm,
                 const int nIDp);

  void pack_VZcell(const REAL_TYPE *array,
               const int vc_comm,
               REAL_TYPE *sendm,
               REAL_TYPE *sendp,
               const int nIDm,
               const int nIDp);

  void unpack_VZcell(REAL_TYPE *array,
                 const int vc_comm,
                 const REAL_TYPE *recvm,
                 const REAL_TYPE *recvp,
                 const int nIDm,
                 const int nIDp);

  bool pack_VEcell(REAL_TYPE *array,
               const int vc_comm,
               REAL_TYPE *sendbuf,
               REAL_TYPE *recvbuf,
               MPI_Request *req);

  void unpack_VEcell(REAL_TYPE *array,
                 const int vc_comm,
                 const REAL_TYPE *recvbuf);

  bool pack_VCcell(REAL_TYPE *array,
               const int vc_comm,
               REAL_TYPE *sendbuf,
               REAL_TYPE *recvbuf,
               MPI_Request *req);

  void unpack_VCcell(REAL_TYPE *array,
                 const int vc_comm,
                 const REAL_TYPE *recvbuf);

  // CB_PackingVectorNode.cpp
  private:

    void pack_VXnode(const REAL_TYPE *array,
                 const int vc_comm,
                 REAL_TYPE *sendm,
                 REAL_TYPE *sendp,
                 const int nIDm,
                 const int nIDp);

    void unpack_VXnode(REAL_TYPE *array,
                   const int vc_comm,
                   const REAL_TYPE *recvm,
                   const REAL_TYPE *recvp,
                   const int nIDm,
                   const int nIDp);

    void pack_VYnode(const REAL_TYPE *array,
                 const int vc_comm,
                 REAL_TYPE *sendm,
                 REAL_TYPE *sendp,
                 const int nIDm,
                 const int nIDp);

    void unpack_VYnode(REAL_TYPE *array,
                   const int vc_comm,
                   const REAL_TYPE *recvm,
                   const REAL_TYPE *recvp,
                   const int nIDm,
                   const int nIDp);

    void pack_VZnode(const REAL_TYPE *array,
                 const int vc_comm,
                 REAL_TYPE *sendm,
                 REAL_TYPE *sendp,
                 const int nIDm,
                 const int nIDp);

    void unpack_VZnode(REAL_TYPE *array,
                   const int vc_comm,
                   const REAL_TYPE *recvm,
                   const REAL_TYPE *recvp,
                   const int nIDm,
                   const int nIDp);

  bool pack_VEnode(REAL_TYPE *array,
               const int vc_comm,
               REAL_TYPE *sendbuf,
               REAL_TYPE *recvbuf,
               MPI_Request *req);

  void unpack_VEnode(REAL_TYPE *array,
                 const int vc_comm,
                 const REAL_TYPE *recvbuf);

  bool pack_VCnode(REAL_TYPE *array,
               const int vc_comm,
               REAL_TYPE *sendbuf,
               REAL_TYPE *recvbuf,
               MPI_Request *req);

  void unpack_VCnode(REAL_TYPE *array,
                 const int vc_comm,
                 const REAL_TYPE *recvbuf);

};

#endif // _CB_SUBDOMAIN_H_
