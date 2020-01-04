#ifndef _CB_PACK_H_
#define _CB_PACK_H_
/*
###################################################################################
#
# CBrick
#
# Copyright (c) 2017-2020 Research Institute for Information Technology(RIIT),
#                    Kyushu University.  All rights reserved.
#
####################################################################################
*/

/*
 * @file   CB_Pack.h
 * @brief  Packing header
 */



/* バッファへのインデクス変換 (I方向)
 *  @param [in] _I  i方向インデクス
 *  @param [in] _J  j方向インデクス
 *  @param [in] _K  k方向インデクス
 *  @param [in] _IS i方向の開始点インデクス
 *  @param [in] _NJ j方向インデクスサイズ
 *  @param [in] _VC 実際に送受信する仮想セル数
 *  @return 1次元インデクス
 */
// 横方向の袖まで含んだ形式
//#define _IDX_SI(_I,_J,_K,_IS,_NJ,_VC) \
//( (_K+(_VC)) * (_VC) * (_NJ+2*(_VC)) \
//+ (_J+(_VC)) * (_VC) \
//+ (_I-(_IS)) \
//)
#define _IDX_SI(_I,_J,_K,_NJ,_VC) \
    ( (_K) * (_VC) * (_NJ) \
    + (_J) * (_VC) \
    + (_I) )


/* バッファへのインデクス変換 (J方向)
 *  @param [in] _I  i方向インデクス
 *  @param [in] _J  j方向インデクス
 *  @param [in] _K  k方向インデクス
 *  @param [in] _NI i方向インデクスサイズ
 *  @param [in] _JS j方向の開始点インデクス
 *  @param [in] _VC 実際に送受信する仮想セル数
 *  @return 1次元インデクス
 */
//#define _IDX_SJ(_I,_J,_K,_NI,_JS,_VC) \
//( (_K+(_VC)) * (_NI+2*(_VC)) * (_VC) \
//+ (_J-(_JS)) * (_NI+2*(_VC)) \
//+ (_I+(_VC)) \
//)
#define _IDX_SJ(_I,_J,_K,_NI,_VC) \
    ( (_K) * (_NI) * (_VC) \
    + (_J) * (_NI) \
    + (_I) )

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
//#define _IDX_SK(_I,_J,_K,_NI,_NJ,_KS,_VC) \
//( (_K-(_KS)) * (_NI+2*(_VC)) * (_NJ+2*(_VC)) \
//+ (_J+(_VC)) * (_NI+2*(_VC)) \
//+ (_I+(_VC)) \
//)
#define _IDX_SK(_I,_J,_K,_NI,_NJ) \
    ( (_K) * (_NI) * (_NJ) \
    + (_J) * (_NI) \
    + (_I) )


/* バッファへのインデクス変換 (I方向)
 *  @param [in] _I  i方向インデクス
 *  @param [in] _J  j方向インデクス
 *  @param [in] _K  k方向インデクス
 *  @param [in] _L  ベクトル成分インデクス {0,1,2}
 *  @param [in] _IS i方向の開始点インデクス
 *  @param [in] _NJ j方向インデクスサイズ
 *  @param [in] _NK k方向インデクスサイズ
 *  @param [in] _VC 実際に送受信する仮想セル数
 *  @return 1次元インデクス
 */
//#define _IDX_VI(_I,_J,_K,_L,_IS,_NJ,_NK,_VC) \
//( (_L) * (_VC) * (_NJ+2*(_VC)) * (_NK+2*(_VC)) \
//+ (_K+(_VC)) * (_VC) * (_NJ+2*(_VC)) \
//+ (_J+(_VC)) * (_VC) \
//+ (_I-(_IS)) \
//)
#define _IDX_VI(_I,_J,_K,_L,_NJ,_NK,_VC) \
( (_L) * (_VC) * (_NJ) * (_NK) \
+ _IDX_SI(_I,_J,_K,_NJ,_VC) \
)

/* バッファへのインデクス変換 (J方向)
 *  @param [in] _I  i方向インデクス
 *  @param [in] _J  j方向インデクス
 *  @param [in] _K  k方向インデクス
 *  @param [in] _L  ベクトル成分インデクス {0,1,2}
 *  @param [in] _NI i方向インデクスサイズ
 *  @param [in] _JS j方向の開始点インデクス
 *  @param [in] _NK k方向インデクスサイズ
 *  @param [in] _VC 実際に送受信する仮想セル数
 *  @return 1次元インデクス
 */
//#define _IDX_VJ(_I,_J,_K,_L,_NI,_JS,_NK,_VC) \
//( (_L) * (_NI+2*(_VC)) * (_VC) * (_NK+2*(_VC)) \
//+ (_K+(_VC)) * (_NI+2*(_VC)) * (_VC) \
//+ (_J-(_JS)) * (_NI+2*(_VC)) \
//+ (_I+(_VC)) \
//)
#define _IDX_VJ(_I,_J,_K,_L,_NI,_NK,_VC) \
( (_L) * (_NI) * (_VC) * (_NK) \
+ _IDX_SJ(_I,_J,_K,_NI,_VC) \
)

/* バッファへのインデクス変換 (K方向)
 *  @param [in] _I  i方向インデクス
 *  @param [in] _J  j方向インデクス
 *  @param [in] _K  k方向インデクス
 *  @param [in] _L  ベクトル成分インデクス {0,1,2}
 *  @param [in] _NI i方向インデクスサイズ
 *  @param [in] _NJ j方向インデクスサイズ
 *  @param [in] _KS ｋ方向の開始点インデクス
 *  @param [in] _VC 実際に送受信する仮想セル数
 *  @return 1次元インデクス
 */
//#define _IDX_VK(_I,_J,_K,_L,_NI,_NJ,_KS,_VC) \
//( (_L) * (_NI+2*(_VC)) * (_NJ+2*(_VC)) * (_VC) \
//+ (_K-(_KS)) * (_NI+2*(_VC)) * (_NJ+2*(_VC)) \
//+ (_J+(_VC)) * (_NI+2*(_VC)) \
//+ (_I+(_VC)) \
//)
#define _IDX_VK(_I,_J,_K,_L,_NI,_NJ,_VC) \
( (_L) * (_NI) * (_NJ) * (_VC) \
+ _IDX_SK(_I,_J,_K,_NI,_NJ) \
)

#endif // _CB_PACK_H_
