#ifndef _CB_COMM_H_
#define _CB_COMM_H_

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
 * @file   CB_Comm.h
 * @brief  BrickComm class Header
 */

#include <mpi.h>

#include <string>
#include <stdlib.h>
#include "CB_Define.h"
#include "CB_Pack.h"


class BrickComm {

private:
  int size[3];          ///< 各サブドメインの要素数 (Local, Non-dimensional
  int comm_tbl[NOFACE]; ///< 隣接ブロックのランク番号

  MPI_Comm mpi_comm;    ///< MPI コミュニケーター
  int halo_width;       ///< ガイドセル幅
  int buf_flag;         ///< バッファを確保済みのときに1
  std::string grid_type;///< "cell" or "node"

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



public:
  // デフォルト コンストラクタ
  BrickComm() {
    halo_width = 0;
    buf_flag = 0;

    for (int i=0; i<NOFACE; i++) comm_tbl[i] = -1;

    for (int i=0; i<3; i++) size[i] = 0;

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
  }

  /** デストラクタ */
  ~BrickComm() {
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



public:
  bool setBrickComm(const int m_sz[],
                    const int m_halo,
                    MPI_Comm m_comm,
                    const int m_tbl[],
                    std::string m_type) {

    this->size[0] = m_sz[0];
    this->size[1] = m_sz[1];
    this->size[2] = m_sz[2];
    this->halo_width  = m_halo;
    this->grid_type   = m_type;
    this->mpi_comm    = m_comm;

    for (int i=0; i<NOFACE; i++)
        this->comm_tbl[i] = m_tbl[i];

    if (m_type == "node" || m_type == "cell") {
      // ok
    }
    else {
      printf("Error : Invalid grid type [%s]\n", m_type.c_str());
      false;
    }

    if ( m_halo  < 0 ) false;

    return true;
  }


  /*
   * @fn init()
   * @brief 通信バッファの確保
   * @param [in] gnum_compo バッファで利用する最大の層数（1-scalar, 3-vector, ?-others）
   */
  bool init(const int num_compo)
  {
    int gc = halo_width;

    if (size[0]==0 || size[1]==0 || size[2]==0 || gc==0 || num_compo==0) {
      return false;
    }

    // バッファ領域としては、最大値で確保しておく
    int f_sz[3];
    f_sz[0] = size[1] * size[2] * gc * num_compo;
    f_sz[1] = size[0] * size[2] * gc * num_compo;
    f_sz[2] = size[0] * size[1] * gc * num_compo;

    if ( !(f_ims = new REAL_TYPE [f_sz[0]]) ) return false;
    if ( !(f_imr = new REAL_TYPE [f_sz[0]]) ) return false;
    if ( !(f_ips = new REAL_TYPE [f_sz[0]]) ) return false;
    if ( !(f_ipr = new REAL_TYPE [f_sz[0]]) ) return false;

    if ( !(f_jms = new REAL_TYPE [f_sz[1]]) ) return false;
    if ( !(f_jmr = new REAL_TYPE [f_sz[1]]) ) return false;
    if ( !(f_jps = new REAL_TYPE [f_sz[1]]) ) return false;
    if ( !(f_jpr = new REAL_TYPE [f_sz[1]]) ) return false;

    if ( !(f_kms = new REAL_TYPE [f_sz[2]]) ) return false;
    if ( !(f_kmr = new REAL_TYPE [f_sz[2]]) ) return false;
    if ( !(f_kps = new REAL_TYPE [f_sz[2]]) ) return false;
    if ( !(f_kpr = new REAL_TYPE [f_sz[2]]) ) return false;

  #ifdef _DIAGONAL_COMM
    // edge
    size_t lx = size[0] * gc * gc * num_compo;
    size_t ly = size[1] * gc * gc * num_compo;
    size_t lz = size[2] * gc * gc * num_compo;
    size_t le = lx*4 + ly*4 + lz*4;
    if ( !(f_es = new REAL_TYPE[le]) ) return false;
    if ( !(f_er = new REAL_TYPE[le]) ) return false;

    // corner
    size_t lc = gc * gc * gc * num_compo * 8;
    if ( !(f_cs = new REAL_TYPE[lc]) ) return false;
    if ( !(f_cr = new REAL_TYPE[lc]) ) return false;
  #endif

    buf_flag = 1; // バッファ確保ずみ

    return true;
  }

// CB_CommS.cpp
public:

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

#endif // _CB_COMM_H_
