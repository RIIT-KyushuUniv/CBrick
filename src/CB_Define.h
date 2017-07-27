#ifndef _CB_DEFINE_H_
#define _CB_DEFINE_H_

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
 * @file   CB_Define.h
 * @brief  CBrick Definition Header
 */

#include <float.h>
#include <math.h>
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

// precision
#ifdef _REAL_IS_DOUBLE_
#define REAL_TYPE double
#else
/** 実数型の指定
 * - デフォルトでは、REAL_TYPE=float
 * - コンパイル時オプション-D_REAL_IS_DOUBLE_を付与することで
 *   REAL_TYPE=doubleになる
 */
#define REAL_TYPE float
#endif

#define _SIZE_DOUBLE_ 8

#define SINGLE_EPSILON 1.19e-7
#define DOUBLE_EPSILON 2.22e-16

#define NOFACE 6

#define ON 1
#define OFF 0


enum DIRection {
  X_minus=0,
  X_plus,
  Y_minus,
  Y_plus,
  Z_minus,
  Z_plus
};


/** 3次元インデクス(i,j,k) -> 1次元インデクス変換マクロ [C version]
 *  @param [in] _I  i方向インデクス
 *  @param [in] _J  j方向インデクス
 *  @param [in] _K  k方向インデクス
 *  @param [in] _NI i方向インデクスサイズ
 *  @param [in] _NJ j方向インデクスサイズ
 *  @param [in] _VC 仮想セル数
 *  @return 1次元インデクス
 */
#define _IDX_S3D(_I,_J,_K,_NI,_NJ,_VC) \
( (_K+_VC) * (_NI+2*_VC) * (_NJ+2*_VC) \
+ (_J+_VC) * (_NI+2*_VC) \
+ (_I+_VC) \
)

/* バッファへのインデクス変換 (I方向)
 *  @param [in] _I  i方向インデクス
 *  @param [in] _J  j方向インデクス
 *  @param [in] _K  k方向インデクス
 *  @param [in] _IS i方向の開始点インデクス
 *  @param [in] _NJ j方向インデクスサイズ
 *  @param [in] _VC 実際に送受信する仮想セル数
 *  @return 1次元インデクス
 */
#define _IDXFX(_I,_J,_K,_IS,_NJ,_VC) \
( (_K+_VC) * _VC * (_NJ+2*_VC) \
+ (_J+_VC) * _VC \
+ (_I-_IS) \
)

/* バッファへのインデクス変換 (J方向)
 *  @param [in] _I  i方向インデクス
 *  @param [in] _J  j方向インデクス
 *  @param [in] _K  k方向インデクス
 *  @param [in] _NI i方向インデクスサイズ
 *  @param [in] _JS j方向の開始点インデクス
 *  @param [in] _VC 実際に送受信する仮想セル数
 *  @return 1次元インデクス
 */
#define _IDXFY(_I,_J,_K,_NI,_JS,_VC) \
( (_K+_VC) * (_NI+2*_VC) * _VC \
+ (_J-_JS) * (_NI+2*_VC) \
+ (_I+_VC) \
)

/* バッファへのインデクス変換 (K方向)
 *  @param [in] _I  i方向インデクス
 *  @param [in] _J  j方向インデクス
 *  @param [in] _K  k方向インデクス
 *  @param [in] _NI i方向インデクスサイズ
 *  @param [in] _NJ j方向インデクスサイズ
 *  @param [in] _KS ｋ方向の開始点インデクス
 *  @param [in] _VC 実際に送受信する仮想セル数
 *  @return 1次元インデクス
 */
#define _IDXFZ(_I,_J,_K,_NI,_NJ,_KS,_VC) \
( (_K-_KS) * (_NI+2*_VC) * (_NJ+2*_VC) \
+ (_J+_VC) * (_NI+2*_VC) \
+ (_I+_VC) \
)


/** 3次元インデクス(i,j,k) -> 1次元インデクス変換マクロ
 *  @note i,j,kインデクスはF表記
 *  @param [in] _I  i方向インデクス
 *  @param [in] _J  j方向インデクス
 *  @param [in] _K  k方向インデクス
 *  @param [in] _NI i方向インデクスサイズ
 *  @param [in] _NJ j方向インデクスサイズ
 *  @param [in] _NK k方向インデクスサイズ
 *  @param [in] _VC 仮想セル数
 *  @return 1次元インデクス
 */
#define _F_IDX_S3D(_I,_J,_K,_NI,_NJ,_NK,_VC) \
( (_K+_VC-1) * (_NI+2*_VC) * (_NJ+2*_VC) \
+ (_J+_VC-1) * (_NI+2*_VC) \
+ (_I+_VC-1) \
)


#define stamped_printf printf("%s (%d):  ",__FILE__, __LINE__), printf
#define stamped_fprintf fprintf(fp, "%s (%d):  ",__FILE__, __LINE__), fprintf
#define mark() printf("%s (%d) [%d]:\n",__FILE__, __LINE__, myRank)
#define Hostonly_ if(myRank==0)

#define Exit(x) \
((void)printf("exit at %s:%u\n", __FILE__, __LINE__), exit((x)))

#endif // _CB_DEFINE_H_
