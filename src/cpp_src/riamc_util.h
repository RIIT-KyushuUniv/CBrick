#ifndef _RIAM_UTY_H_
#define _RIAM_UTY_H_

//##################################################################################
//
// Flow Base class
//
// Copyright (c) 2007-2011 VCAD System Research Program, RIKEN.
// All rights reserved.
//
// Copyright (c) 2011-2015 Institute of Industrial Science, The University of Tokyo.
// All rights reserved.
//
// Copyright (c) 2012-2016 Advanced Institute for Computational Science, RIKEN.
// All rights reserved.
//
//##################################################################################

/**
 * @file   RIAM_util.h
 * @brief  FlowBase FBUtility class Header for FFV-C
 * @author aics
 */

#include "cpm_Define.h"
#include <math.h>
#include <string>
#include <string.h>
#include <stdio.h>
#include <algorithm>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdlib.h>

#include "riamc_Define.h"

#include "omp.h"

using namespace std;

// #################################################################
class FBUtility {

public:
  /** コンストラクタ */
  FBUtility() {}

  /** デストラクタ */
  ~FBUtility() {}


public:

  /** 文字列を小文字にして比較
   * @param [in] str1  比較string
   * @param [in] str2  比較string
   * @return true-同じ / false-異なる
   */
  static bool compare(const string str1, const string str2)
  {
    if ( !strcasecmp(str1.c_str(), str2.c_str()) )
    {
      return true;
    }
    return false;
  }


  // ディレクトリがなければ作成、既存なら何もしない（単一ディレクトリ）
  static int c_mkdir(const char* path);



  /**
   * @brief dirの方向ラベルを返す
   * @param [in] dir 方向コード
   * @return 方向ラベル
   */
  static string getDirection(const int dir)
  {
    string face;
    if      (dir == X_minus) face = "X-";
    else if (dir == X_plus)  face = "X+";
    else if (dir == Y_minus) face = "Y-";
    else if (dir == Y_plus)  face = "Y+";
    else if (dir == Z_minus) face = "Z-";
    else if (dir == Z_plus)  face = "Z+";
    return face;
  }

  /**
   * @brief dirの方向ラベルを返す
   * @param [in] dir 方向コード
   * @return 方向ラベル
   */
  static string getDirStr(const int dir)
  {
    string face;
    if      (dir == X_minus) face = "Xminus";
    else if (dir == X_plus)  face = "Xplus";
    else if (dir == Y_minus) face = "Yminus";
    else if (dir == Y_plus)  face = "Yplus";
    else if (dir == Z_minus) face = "Zminus";
    else if (dir == Z_plus)  face = "Zplus";
    return face;
  }


  // 階層ディレクトリの作成
  static int mkdirs(string path);


  /**
   * @brief メモリ使用量を表示する
   * @param [in] mode     処理モード
   * @param [in] Memory   必要メモリ量
   * @param [in] l_memory local
   * @param [in] fp       ファイルポインタ
   */
  static void MemoryRequirement(const char* mode, const double Memory, const double l_memory, FILE* fp);

};

#endif // _RIAM_UTY_H_
