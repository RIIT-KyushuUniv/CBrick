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
 * @file   FBUtility.C
 * @brief  FlowBase FBUtility class
 * @author aics
 */


#include "riamc_util.h"


// #################################################################
/**
 * @brief ディレクトリがなければ作成、既存なら何もしない（単一ディレクトリ）
 * @param [in] path ディレクトリパス
 */
int FBUtility::c_mkdir(const char* path)
{
  // 標準ライブラリ呼び出し
  // パーミッションはumaskで指定
  umask(022);

  int ret = mkdir(path, 0777); // rwx

  if ( 0 != ret )
  {
    // 既存以外のエラー
    if ( EEXIST != errno )
    {
      printf( "\tError(errno)=[%s]\n", strerror(errno) );
      return 0;
    }
  }

  return 1;
}


// #################################################################
// メモリ使用量を表示する
void FBUtility::MemoryRequirement(const char* mode, const double Memory, const double l_memory, FILE* fp)
{
  const double mem = Memory;
  const double lmem= l_memory;
  const double KB = 1024.0;
  const double MB = 1024.0*KB;
  const double GB = 1024.0*MB;
  const double TB = 1024.0*GB;
  const double PB = 1024.0*TB;
  const double factor = 1.05; // estimate 5% for addtional

  fprintf (fp,"\t>> Memory required for %s : ", mode);

  // Global memory
  fprintf (fp," Global=");
  if ( mem > PB ) {
    fprintf (fp,"%6.2f (PB)", mem / PB *factor);
  }
  else if ( mem > TB ) {
    fprintf (fp,"%6.2f (TB)", mem / TB *factor);
  }
  else if ( mem > GB ) {
    fprintf (fp,"%6.2f (GB)", mem / GB *factor);
  }
  else if ( mem > MB ) {
    fprintf (fp,"%6.2f (MB)", mem / MB *factor);
  }
  else if ( mem > KB ) {
    fprintf (fp,"%6.2f (KB)", mem / KB *factor);
  }
  else if ( mem <= KB ){
    fprintf (fp,"%6.2f (B)", mem *factor);
  }
  else {
    fprintf (fp,"Caution! Memory required : %d (Byte)", (int)(mem *factor) );
  }

  // Local memory
  fprintf (fp," : Local=");
  if ( lmem > PB ) {
    fprintf (fp,"%6.2f (PB)\n", lmem / PB *factor);
  }
  else if ( lmem > TB ) {
    fprintf (fp,"%6.2f (TB)\n", lmem / TB *factor);
  }
  else if ( lmem > GB ) {
    fprintf (fp,"%6.2f (GB)\n", lmem / GB *factor);
  }
  else if ( lmem > MB ) {
    fprintf (fp,"%6.2f (MB)\n", lmem / MB *factor);
  }
  else if ( lmem > KB ) {
    fprintf (fp,"%6.2f (KB)\n", lmem / KB *factor);
  }
  else if ( lmem <= KB ){
    fprintf (fp,"%6.2f (B)\n", lmem *factor);
  }
  else {
    fprintf (fp,"Caution! Memory required : %d (Byte)\n", (int)(lmem *factor) );
  }

  fflush(fp);
}


// #################################################################
/**
 * @brief 階層ディレクトリの作成
 * @param [in] path ディレクトリパス
 */
int FBUtility::mkdirs(string path)
{
  int len = path.size() + 4;
  char* buf = new char[len];

  if ( !buf )
  {
    printf("Error: create buffer(%d) %s\n", errno, strerror(errno));
    return(-1);
  }
  strcpy(buf, path.c_str());

  // 階層的にディレクトリを作成
  char *p = NULL;
  int ret = 0;

  for(p=strchr(buf+1, '/'); p; p=strchr(p+1, '/'))
  {
    *p = '\0';
    ret = c_mkdir(buf);
    if (ret != 1)
    {
      delete [] buf;
      return(-1);
    }
    *p = '/';
  }

  if (buf)
  {
    delete [] buf;
  }
  return(1);
}
