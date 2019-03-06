#ifndef _CB_COMM_H_
#define _CB_COMM_H_

/*
###################################################################################
#
# CBrick
#
# Copyright (c) 2017-2019 Research Institute for Information Technology(RIIT),
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

  // バッファは8バイトで確保(double, long long)
  double* f_ims;  // I- direction send
  double* f_imr;  // I- direction recv
  double* f_ips;  // I+ direction send
  double* f_ipr;  // I+ direction recv
  double* f_jms;  // J- direction send
  double* f_jmr;  // J- direction recv
  double* f_jps;  // J+ direction send
  double* f_jpr;  // J+ direction recv
  double* f_kms;  // K- direction send
  double* f_kmr;  // K- direction recv
  double* f_kps;  // K+ direction send
  double* f_kpr;  // K+ direction recv
#ifdef _DIAGONAL_COMM
  double* f_es;   // edge send
  double* f_er;   // edge recv
  double* f_cs;   // corner send
  double* f_cr;   // corner recv
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
  
  // デストラクタ
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
  /* #########################################################
   * @brief パラメータセット
   * @param [in] m_sz   配列長
   * @param [in] m_halo ガイドセル長
   * @param [in] m_comm コミュニケータ
   * @param [in] m_tbl  隣接IDテーブル
   * @param [in] m_type "node" or "cell"
   */
  bool setBrickComm(const int m_sz[],
                    const int m_halo,
                    MPI_Comm m_comm,
                    const int m_tbl[],
                    std::string m_type)
  {
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

  
  /* #########################################################
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
    
    if ( !(f_ims = new double [f_sz[0]]) ) return false;
    if ( !(f_imr = new double [f_sz[0]]) ) return false;
    if ( !(f_ips = new double [f_sz[0]]) ) return false;
    if ( !(f_ipr = new double [f_sz[0]]) ) return false;
    
    if ( !(f_jms = new double [f_sz[1]]) ) return false;
    if ( !(f_jmr = new double [f_sz[1]]) ) return false;
    if ( !(f_jps = new double [f_sz[1]]) ) return false;
    if ( !(f_jpr = new double [f_sz[1]]) ) return false;
    
    if ( !(f_kms = new double [f_sz[2]]) ) return false;
    if ( !(f_kmr = new double [f_sz[2]]) ) return false;
    if ( !(f_kps = new double [f_sz[2]]) ) return false;
    if ( !(f_kpr = new double [f_sz[2]]) ) return false;
    
#ifdef _DIAGONAL_COMM
    // edge
    size_t lx = size[0] * gc * gc * num_compo;
    size_t ly = size[1] * gc * gc * num_compo;
    size_t lz = size[2] * gc * gc * num_compo;
    size_t le = lx*4 + ly*4 + lz*4;
    if ( !(f_es = new double [le]) ) return false;
    if ( !(f_er = new double [le]) ) return false;
    
    // corner
    size_t lc = gc * gc * gc * num_compo * 8;
    if ( !(f_cs = new double [lc]) ) return false;
    if ( !(f_cr = new double [lc]) ) return false;
#endif
    
    buf_flag = 1; // バッファ確保ずみ
    
    return true;
  }
  
  
  
// CB_Comm_inline.h
public:
  
  /** MPI_Datatypeを取得
   *  @param[in] ptr 取得したいデータのポインタ
   *  @return MPI_Datatype
   */
  template<class T> inline
  static MPI_Datatype GetMPI_Datatype(T *ptr);
  
  
  /*
   * @brief 隣接間通信 IsendIrecv Interface 1
   * @param [in]  ms   Send buffer to Minus direction
   * @param [out] mr   Recieve buffer from Mminus direction
   * @param [in]  ps   Send buffer to Plus direction
   * @param [out] pr   Recieve buffer from Plus direction
   * @param [in]  msz  send/recieve size
   * @param [in]  nIDm Neighbor ID for Minus direction
   * @param [in]  nIDp Neighbor ID for Plus direction
   * @param [out] req  Array of MPI request
   * @retval true-success, false-fail
   */
  template <class T> inline
  bool IsendIrecv(T* ms,
                  T* mr,
                  T* ps,
                  T* pr,
                  int msz,
                  int nIDm,
                  int nIDp,
                  MPI_Request *req);
  
  
  /*
   * @brief Irecv
   * @param [out]    ptr  Recieve pointer
   * @param [in]     sz   Recieve data size
   * @param [in]     nID  Neighbor ID
   * @param [in,out] req  Array of MPI request
   */
  template <class T> inline
  bool IrecvData(T* ptr,
                 int sz,
                 int nID,
                 MPI_Request *req);
  
  
  /*
   * @brief Isend
   * @param [in]     ptr  Send pointer
   * @param [in]     sz   Send data size
   * @param [in]     nID  Neighbor ID
   * @param [in,out] req  Array of MPI request
   */
  template <class T> inline
  bool IsendData(T* ptr,
                 int sz,
                 int nID,
                 MPI_Request *req);

  
  
  
// CB_Comm.cpp
public:
  
  /*
   * @brief 隣接間通信 IsendIrecv Interface 2
   * @param [in]  dtype  送受信データの型
   * @param [in]  ms     Send buffer to Minus direction
   * @param [out] mr     Recieve buffer from Mminus direction
   * @param [in]  ps     Send buffer to Plus direction
   * @param [out] pr     Recieve buffer from Plus direction
   * @param [in]  msz    send/recieve size
   * @param [in]  nIDm   Neighbor ID for Minus direction
   * @param [in]  nIDp   Neighbor ID for Plus direction
   * @param [out] req    Array of MPI request
   * @retval true-success, false-fail
   */
  bool IsendIrecv(MPI_Datatype dtype,
                  void* ms,
                  void* mr,
                  void* ps,
                  void* pr,
                  int msz,
                  int nIDm,
                  int nIDp,
                  MPI_Request* req);

  
  /*
   * @brief Irecv
   * @param [in]     dtype 受信データの型
   * @param [out]    ptr   Recieve pointer
   * @param [in]     sz    Recieve data size
   * @param [in]     nID   Neighbor ID
   * @param [in,out] req   Array of MPI request
   */
  bool IrecvData(MPI_Datatype dtype,
                 void* ptr,
                 int sz,
                 int nID,
                 MPI_Request* req);
  

  /*
   * @brief Isend
   * @param [in]     dtype 送信データの型
   * @param [in]     ptr   Send pointer
   * @param [in]     sz    Send data size
   * @param [in]     nID   Neighbor ID
   * @param [in,out] req   Array of MPI request
   */
  bool IsendData(MPI_Datatype dtype,
                 void* ptr,
                 int sz,
                 int nID,
                 MPI_Request* req);
  
  
  
  
public:
  
  /* #########################################################
   * @brief スカラー変数 node
   * @param [in,out]  src     スカラー変数
   * @param [in]      gc_comm 実際に通信する通信面数
   * @param [in,out]  req     MPI_Request
   * @retval true-success, false-fail
   */
  template <class T>
  bool Comm_S_node(T* src, const int gc_comm, MPI_Request *req);


  /* #########################################################
   * @brief スカラー変数 cell
   * @param [in,out]  src     スカラー変数
   * @param [in]      gc_comm 実際に通信する通信面数
   * @param [in,out]  req     MPI_Request
   * @retval true-success, false-fail
   */
  template <class T>
  bool Comm_S_cell(T* src, const int gc_comm, MPI_Request *req);
  
  
  /* #########################################################
   * @brief スカラー変数 node
   * @param [in,out]  dest    スカラー変数
   * @param [in]      gc_comm 実際に通信する通信面数
   * @param [out]     req     Array of MPI request
   * @retval true-success, false-fail
   */
  template <class T>
  bool Comm_S_wait_node(T* dest, const int gc_comm, MPI_Request *req);
  
  
  /* #########################################################
   * @brief スカラー変数 cell
   * @param [in,out]  dest    スカラー変数
   * @param [in]      gc_comm 実際に通信する通信面数
   * @param [out]     req     Array of MPI request
   * @retval true-success, false-fail
   */
  template <class T>
  bool Comm_S_wait_cell(T* dest, const int gc_comm, MPI_Request *req);

  
  
  
// CB_CommV.cpp
public:
  
  /* #########################################################
   * @brief ベクトル変数 node
   * @param [in,out]  src     ベクトル変数
   * @param [in]      gc_comm 実際に通信する通信面数
   * @param [in,out]  req     MPI_Request
   * @retval true-success, false-fail
   */
  template <class T>
  bool Comm_V_node(T* src, const int gc_comm, MPI_Request *req);

  
  /* #########################################################
   * @brief ベクトル変数 cell
   * @param [in,out]  src     ベクトル変数
   * @param [in]      gc_comm 実際に通信する通信面数
   * @param [in,out]  req     MPI_Request
   * @retval true-success, false-fail
   */
  template <class T>
  bool Comm_V_cell(T* src, const int gc_comm, MPI_Request *req);
  
  
  /* #########################################################
   * @brief ベクトル変数 node
   * @param [in,out]  dest    ベクトル変数
   * @param [in]      gc_comm 実際に通信する通信面数
   * @param [out]     req     Array of MPI request
   * @retval true-success, false-fail
   */
  template <class T>
  bool Comm_V_wait_node(T* dest, const int gc_comm, MPI_Request *req);

  
  /* #########################################################
   * @brief ベクトル変数 cell
   * @param [in,out]  dest    ベクトル変数
   * @param [in]      gc_comm 実際に通信する通信面数
   * @param [out]     req     Array of MPI request
   * @retval true-success, false-fail
   */
  template <class T>
  bool Comm_V_wait_cell(T* dest, const int gc_comm, MPI_Request *req);
  

  
  
// CB_PackingScalarCell.cpp
private:
  
  template <class T>
  void pack_SXcell(const T *array,
                   const int vc_comm,
                   T *sendm,
                   T *sendp,
                   const int nIDm,
                   const int nIDp);
  
  template <class T>
  void unpack_SXcell(T *array,
                     const int vc_comm,
                     const T *recvm,
                     const T *recvp,
                     const int nIDm,
                     const int nIDp);
  
  template <class T>
  void pack_SYcell(const T *array,
                   const int vc_comm,
                   T *sendm,
                   T *sendp,
                   const int nIDm,
                   const int nIDp);
  
  template <class T>
  void unpack_SYcell(T *array,
                     const int vc_comm,
                     const T *recvm,
                     const T *recvp,
                     const int nIDm,
                     const int nIDp);
  
  template <class T>
  void pack_SZcell(const T *array,
                   const int vc_comm,
                   T *sendm,
                   T *sendp,
                   const int nIDm,
                   const int nIDp);
  
  template <class T>
  void unpack_SZcell(T *array,
                     const int vc_comm,
                     const T *recvm,
                     const T *recvp,
                     const int nIDm,
                     const int nIDp);

  
#ifdef _DIAGONAL_COMM
  template <class T>
  bool pack_SEcell(T *array,
                   const int vc_comm,
                   T *sendbuf,
                   T *recvbuf,
                   MPI_Request *req);
  
  template <class T>
  void unpack_SEcell(T *array,
                     const int vc_comm,
                     const T *recvbuf);
  
  template <class T>
  bool pack_SCcell(T *array,
                   const int vc_comm,
                   T *sendbuf,
                   T *recvbuf,
                   MPI_Request *req);
  
  template <class T>
  void unpack_SCcell(T *array,
                     const int vc_comm,
                     const T *recvbuf);
#endif // _DIAGONAL_COMM

  
  
  // CB_PackingScalarNode.cpp
private:
  
  template <class T>
  void pack_SXnode(const T *array,
                   const int vc_comm,
                   T *sendm,
                   T *sendp,
                   const int nIDm,
                   const int nIDp);
  
  template <class T>
  void unpack_SXnode(T *array,
                     const int vc_comm,
                     const T *recvm,
                     const T *recvp,
                     const int nIDm,
                     const int nIDp);
  
  template <class T>
  void pack_SYnode(const T *array,
                   const int vc_comm,
                   T *sendm,
                   T *sendp,
                   const int nIDm,
                   const int nIDp);
  
  template <class T>
  void unpack_SYnode(T *array,
                     const int vc_comm,
                     const T *recvm,
                     const T *recvp,
                     const int nIDm,
                     const int nIDp);
  
  template <class T>
  void pack_SZnode(const T *array,
                   const int vc_comm,
                   T *sendm,
                   T *sendp,
                   const int nIDm,
                   const int nIDp);
  
  template <class T>
  void unpack_SZnode(T *array,
                     const int vc_comm,
                     const T *recvm,
                     const T *recvp,
                     const int nIDm,
                     const int nIDp);
  
  
#ifdef _DIAGONAL_COMM
  template <class T>
  bool pack_SEnode(T *array,
                   const int vc_comm,
                   T *sendbuf,
                   T *recvbuf,
                   MPI_Request *req);
  
  template <class T>
  void unpack_SEnode(T *array,
                     const int vc_comm,
                     const T *recvbuf);
  
  template <class T>
  bool pack_SCnode(T *array,
                   const int vc_comm,
                   T *sendbuf,
                   T *recvbuf,
                   MPI_Request *req);
  
  template <class T>
  void unpack_SCnode(T *array,
                     const int vc_comm,
                     const T *recvbuf);
#endif // _DIAGONAL_COMM
  
  
  
  // CB_PackingVectorCell.cpp
private:
  
  template <class T>
  void pack_VXcell(const T *array,
                   const int vc_comm,
                   T *sendm,
                   T *sendp,
                   const int nIDm,
                   const int nIDp);
  
  template <class T>
  void unpack_VXcell(T *array,
                     const int vc_comm,
                     const T *recvm,
                     const T *recvp,
                     const int nIDm,
                     const int nIDp);
  
  template <class T>
  void pack_VYcell(const T *array,
                   const int vc_comm,
                   T *sendm,
                   T *sendp,
                   const int nIDm,
                   const int nIDp);
  
  template <class T>
  void unpack_VYcell(T *array,
                     const int vc_comm,
                     const T *recvm,
                     const T *recvp,
                     const int nIDm,
                     const int nIDp);
  
  template <class T>
  void pack_VZcell(const T *array,
                   const int vc_comm,
                   T *sendm,
                   T *sendp,
                   const int nIDm,
                   const int nIDp);
  
  template <class T>
  void unpack_VZcell(T *array,
                     const int vc_comm,
                     const T *recvm,
                     const T *recvp,
                     const int nIDm,
                     const int nIDp);

  
#ifdef _DIAGONAL_COMM
  template <class T>
  bool pack_VEcell(T *array,
                   const int vc_comm,
                   T *sendbuf,
                   T *recvbuf,
                   MPI_Request *req);
  
  template <class T>
  void unpack_VEcell(T *array,
                     const int vc_comm,
                     const T *recvbuf);
  
  template <class T>
  bool pack_VCcell(T *array,
                   const int vc_comm,
                   T *sendbuf,
                   T *recvbuf,
                   MPI_Request *req);
  
  template <class T>
  void unpack_VCcell(T *array,
                     const int vc_comm,
                     const T *recvbuf);
#endif // _DIAGONAL_COMM
  
  
  
  // CB_PackingVectorNode.cpp
private:
  
  template <class T>
  void pack_VXnode(const T *array,
                   const int vc_comm,
                   T *sendm,
                   T *sendp,
                   const int nIDm,
                   const int nIDp);
  
  template <class T>
  void unpack_VXnode(T *array,
                     const int vc_comm,
                     const T *recvm,
                     const T *recvp,
                     const int nIDm,
                     const int nIDp);
  
  template <class T>
  void pack_VYnode(const T *array,
                   const int vc_comm,
                   T *sendm,
                   T *sendp,
                   const int nIDm,
                   const int nIDp);
  
  template <class T>
  void unpack_VYnode(T *array,
                     const int vc_comm,
                     const T *recvm,
                     const T *recvp,
                     const int nIDm,
                     const int nIDp);
  
  template <class T>
  void pack_VZnode(const T *array,
                   const int vc_comm,
                   T *sendm,
                   T *sendp,
                   const int nIDm,
                   const int nIDp);
  
  template <class T>
  void unpack_VZnode(T *array,
                     const int vc_comm,
                     const T *recvm,
                     const T *recvp,
                     const int nIDm,
                     const int nIDp);
  
  
#ifdef _DIAGONAL_COMM
  template <class T>
  bool pack_VEnode(T *array,
                   const int vc_comm,
                   T *sendbuf,
                   T *recvbuf,
                   MPI_Request *req);
  
  template <class T>
  void unpack_VEnode(T *array,
                     const int vc_comm,
                     const T *recvbuf);
  
  template <class T>
  bool pack_VCnode(T *array,
                   const int vc_comm,
                   T *sendbuf,
                   T *recvbuf,
                   MPI_Request *req);
  
  template <class T>
  void unpack_VCnode(T *array,
                     const int vc_comm,
                     const T *recvbuf);
#endif // _DIAGONAL_COMM
  
};


//インライン関数
#include "CB_Comm_inline.h"
#include "CB_PackingScalarCell.h"
#include "CB_PackingScalarNode.h"
#include "CB_PackingVectorCell.h"
#include "CB_PackingVectorNode.h"

#endif // _CB_COMM_H_
