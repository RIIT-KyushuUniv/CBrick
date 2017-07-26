//
//  comm.h
//  diff3dp
//
//  Created by keno on 2016/07/03.
//  Copyright © 2016年 keno. All rights reserved.
//

#ifndef comm_h
#define comm_h

#include <mpi.h>
#include <stdio.h>
#include "diff3dp.h"

/** 3次元インデクス(i,j,k) -> 1次元インデクス変換マクロ [C version]
 *  @param[in] _I  i方向インデクス
 *  @param[in] _J  j方向インデクス
 *  @param[in] _K  k方向インデクス
 *  @param[in] _NI i方向インデクスサイズ
 *  @param[in] _NJ j方向インデクスサイズ
 *  @param[in] _NK k方向インデクスサイズ
 *  @param[in] _VC 仮想セル数
 *  @return 1次元インデクス
 */
#define _IDX_S3D(_I,_J,_K,_NI,_NJ,_NK,_VC) \
( (size_t)(_K+_VC) * (size_t)(_NI+2*_VC) * (size_t)(_NJ+2*_VC) \
+ (size_t)(_J+_VC) * (size_t)(_NI+2*_VC) \
+ (size_t)(_I+_VC) \
)

#define _IDXFX(_I,_J,_K,_IS,_NJ,_NK,_VC) \
( (size_t)(_K+_VC) * (size_t)(_VC) * (size_t)(_NJ+2*_VC) \
+ (size_t)(_J+_VC) * (size_t)(_VC) \
+ (size_t)(_I-(_IS)) \
)

#define _IDXFY(_I,_J,_K,_NI,_JS,_NK,_VC) \
( (size_t)(_K+_VC)   * (size_t)(_NI+2*_VC) * (size_t)(_VC) \
+ (size_t)(_J-(_JS)) * (size_t)(_NI+2*_VC) \
+ (size_t)(_I+_VC) \
)

#define _IDXFZ(_I,_J,_K,_NI,_NJ,_KS,_VC) \
( (size_t)(_K-(_KS)) * (size_t)(_NI+2*_VC) * (size_t)(_NJ+2*_VC) \
+ (size_t)(_J+_VC)   * (size_t)(_NI+2*_VC) \
+ (size_t)(_I+_VC) \
)


void packX(const REAL_TYPE *array,
           const int imax, const int jmax, const int kmax, const int vc, const int vc_comm,
           REAL_TYPE *sendm,
           REAL_TYPE *sendp,
           const int nIDm, const int nIDp);

void unpackX(REAL_TYPE *array,
             const int imax, const int jmax, const int kmax, const int vc, const int vc_comm,
             const REAL_TYPE *recvm,
             const REAL_TYPE *recvp,
             const int nIDm, const int nIDp);

void packY(const REAL_TYPE *array,
           const int imax, const int jmax, const int kmax, const int vc, const int vc_comm,
           REAL_TYPE *sendm,
           REAL_TYPE *sendp,
           const int nIDm, const int nIDp);

void unpackY(REAL_TYPE *array,
             const int imax, const int jmax, const int kmax, const int vc, const int vc_comm,
             const REAL_TYPE *recvm,
             const REAL_TYPE *recvp,
             const int nIDm, const int nIDp);

void packZ(const REAL_TYPE *array,
           const int imax, const int jmax, const int kmax, const int vc, const int vc_comm,
           REAL_TYPE *sendm,
           REAL_TYPE *sendp,
           const int nIDm, const int nIDp);

void unpackZ(REAL_TYPE *array,
             const int imax, const int jmax, const int kmax, const int vc, const int vc_comm,
             const REAL_TYPE *recvm,
             const REAL_TYPE *recvp,
             const int nIDm, const int nIDp);

void send_and_recv(REAL_TYPE* ms, REAL_TYPE* mr, REAL_TYPE* ps, REAL_TYPE* pr, size_t msz, int nIDm, int nIDp);

void sendrecv(REAL_TYPE* ms, REAL_TYPE* mr, REAL_TYPE* ps, REAL_TYPE* pr, size_t msz, int nIDm, int nIDp);

void SyncGC_blocking(const int* tbl,
                     REAL_TYPE* xms, REAL_TYPE* xmr,
                     REAL_TYPE* xps, REAL_TYPE* xpr,
                     REAL_TYPE* yms, REAL_TYPE* ymr,
                     REAL_TYPE* yps, REAL_TYPE* ypr,
                     REAL_TYPE* zms, REAL_TYPE* zmr,
                     REAL_TYPE* zps, REAL_TYPE* zpr,
                     REAL_TYPE* src,
                     const int imax, const int jmax, const int kmax, const int gc, const int gc_comm);

void SyncGC_nonblocking(const int* tbl,
                        REAL_TYPE* xms, REAL_TYPE* xmr,
                        REAL_TYPE* xps, REAL_TYPE* xpr,
                        REAL_TYPE* yms, REAL_TYPE* ymr,
                        REAL_TYPE* yps, REAL_TYPE* ypr,
                        REAL_TYPE* zms, REAL_TYPE* zmr,
                        REAL_TYPE* zps, REAL_TYPE* zpr,
                        REAL_TYPE* src,
                        const int imax, const int jmax, const int kmax, const int gc, const int gc_comm, MPI_Request *req);

void IsendIrecv(REAL_TYPE* ms, REAL_TYPE* mr, REAL_TYPE* ps, REAL_TYPE* pr,
                size_t msz, int nIDm, int nIDp, MPI_Request *req);

void wait_nonblocking(const int* tbl,
                    REAL_TYPE* xms, REAL_TYPE* xmr,
                    REAL_TYPE* xps, REAL_TYPE* xpr,
                    REAL_TYPE* yms, REAL_TYPE* ymr,
                    REAL_TYPE* yps, REAL_TYPE* ypr,
                    REAL_TYPE* zms, REAL_TYPE* zmr,
                    REAL_TYPE* zps, REAL_TYPE* zpr,
                    REAL_TYPE* dest,
                    const int imax, const int jmax, const int kmax, const int gc, const int gc_comm, MPI_Request *req);

#endif /* comm_h */
