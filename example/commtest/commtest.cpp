//
//  commtest.cpp
//
//  Created by keno on 2017/07/01.
//  Copyright © 2017 keno. All rights reserved.
//

// Compile
// $ gcc commtest.cpp -o commtest

// Execution
// $ mpirun -np X commtest nx ny nz gc grid index
// (ex)
// $ mpirun -np 4 commtest 64 64 64 2 node F
// $ mpirun -np 4 commtest 65 65 65 2 cell C

// 袖通信テスト

#define BASE 1000

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

#include <CB_SubDomain.h>
#include <CB_Comm.h>
#include <stdlib.h>
#include <string.h>

////////////////////////////////////////////////////////////////////////////////
REAL_TYPE* alloc_real(const int* sz, const int gc, const int sv=1)
{
  REAL_TYPE* p=NULL;
  size_t len = (size_t)( (sz[0]+2*gc) * (sz[1]+2*gc) * (sz[2]+2*gc) * sv );
  if ( !(p = new REAL_TYPE[len]) ) {
    printf("fail to allocate memory\n");
  }
  for( size_t i=0; i<len; i++ ) p[i] = -1.0;

  return p;
}



////////////////////////////////////////////////////////////////////////////////
bool checkTransferX(const int* rsize,
                    const int gc,
                    const int myRank,
                    const REAL_TYPE* X,
                    const REAL_TYPE* Y,
                    const REAL_TYPE* Z,
                    const int* nID,
                    FILE* fp)
{
  int NI = rsize[3*myRank+0];
  int NJ = rsize[3*myRank+1];
  int NK = rsize[3*myRank+2];

  int dm = nID[I_minus];
  int dp = nID[I_plus];

/*
  NodeのC表記
                        gc = 0     1
                             <--gc->
         NI-3  NI-2  NI-1   NI   NI+1
      -----+-----+-----|-----+-----+-----
          -2    -1     0     1     2
                 <--gc->
            gc = 0     1

*/
  // MINUS DIR
  if (dm >= 0) {
    fprintf(fp,"\nI_MINUS dir from rank[%d] ========================================\n", dm);
    fprintf(fp,"\nX  -----------------------------------\n\n");
    for( int k=0; k<NK; k++ ){
    for( int j=0; j<NJ; j++ ){
    for (int i=0; i<gc; i++) {
      REAL_TYPE bf = X[_IDX_S3D(i-1,j,k,NI,NJ,gc)];
      int rk = (int)bf/(int)BASE;
      int val= (int)(bf-rk*BASE);
      if ( (rk != dm) || (val != rsize[3*dm+0]-2+i) ) {
        fprintf(fp,"(%3d %3d %3d) rk= %d val= %d gc= %d\n", i-1,j,k,rk,val,i);
      }
    }}}

    fprintf(fp,"\nY  -----------------------------------\n\n");
    for( int k=0; k<NK; k++ ){
    for( int j=0; j<NJ; j++ ){
    for (int i=0; i<gc; i++) {
      REAL_TYPE bf = Y[_IDX_S3D(i-1,j,k,NI,NJ,gc)];
      int rk = (int)bf/(int)BASE;
      int val= (int)(bf-rk*BASE);
      if ( (rk != dm) || (val != j) ) {
        fprintf(fp,"(%3d %3d %3d) rk= %d val= %d gc= %d\n", i-1,j,k,rk,val,i);
      }
    }}}

    fprintf(fp,"\nZ  -----------------------------------\n\n");
    for( int k=0; k<NK; k++ ){
    for( int j=0; j<NJ; j++ ){
    for (int i=0; i<gc; i++) {
      REAL_TYPE bf = Z[_IDX_S3D(i-1,j,k,NI,NJ,gc)];
      int rk = (int)bf/(int)BASE;
      int val= (int)(bf-rk*BASE);
      if ( (rk != dm) || (val != k) ) {
        fprintf(fp,"(%3d %3d %3d) rk= %d val= %d gc= %d\n", i-1,j,k,rk,val,i);
      }
    }}}
  }


  // PLUS DIR
  if (dp >= 0) {
    fprintf(fp,"\nI_PLUS dir from rank[%d] ========================================\n", dp);
    fprintf(fp,"\nX  -----------------------------------\n\n");
    for( int k=0; k<NK; k++ ){
    for( int j=0; j<NJ; j++ ){
    for (int i=0; i<gc; i++) {
      REAL_TYPE bf = X[_IDX_S3D(NI+i,j,k,NI,NJ,gc)];
      int rk = (int)bf/(int)BASE;
      int val= (int)(bf-rk*BASE);
      if ( (rk != dp) || (val != 1+i) ) {
        fprintf(fp,"(%3d %3d %3d) rk= %d val= %d gc= %d\n", NI+i,j,k,rk,val,i);
      }
    }}}

    fprintf(fp,"\nY  -----------------------------------\n\n");
    for( int k=0; k<NK; k++ ){
    for( int j=0; j<NJ; j++ ){
    for (int i=0; i<gc; i++) {
      REAL_TYPE bf = Y[_IDX_S3D(NI+i,j,k,NI,NJ,gc)];
      int rk = (int)bf/(int)BASE;
      int val= (int)(bf-rk*BASE);
      if ( (rk != dp) || (val != j) ) {
        fprintf(fp,"(%3d %3d %3d) rk= %d val= %d gc= %d\n", NI+i,j,k,rk,val,i);
      }
    }}}

    fprintf(fp,"\nZ  -----------------------------------\n\n");
    for( int k=0; k<NK; k++ ){
    for( int j=0; j<NJ; j++ ){
    for (int i=0; i<gc; i++) {
      REAL_TYPE bf = Z[_IDX_S3D(NI+i,j,k,NI,NJ,gc)];
      int rk = (int)bf/(int)BASE;
      int val= (int)(bf-rk*BASE);
      if ( (rk != dp) || (val != k) ) {
        fprintf(fp,"(%3d %3d %3d) rk= %d val= %d gc= %d\n", NI+i,j,k,rk,val,i);
      }
    }}}
  }

  return true;
}




////////////////////////////////////////////////////////////////////////////////
bool checkTransferY(const int* rsize,
                    const int gc,
                    const int myRank,
                    const REAL_TYPE* X,
                    const REAL_TYPE* Y,
                    const REAL_TYPE* Z,
                    const int* nID,
                    FILE* fp)
{
  int NI = rsize[3*myRank+0];
  int NJ = rsize[3*myRank+1];
  int NK = rsize[3*myRank+2];

  int dm = nID[J_minus];
  int dp = nID[J_plus];

  // MINUS DIR
  if (dm >= 0) {
    fprintf(fp,"\nJ_MINUS from rank[%d] ========================================\n", dm);
    fprintf(fp,"\nX  -----------------------------------\n\n");
    for( int k=0; k<NK; k++ ){
    for( int j=0; j<gc; j++ ){
    for (int i=0; i<NI; i++) {
      REAL_TYPE bf = X[_IDX_S3D(i,j-1,k,NI,NJ,gc)];
      int rk = (int)bf/(int)BASE;
      int val= (int)(bf-rk*BASE);
      if ( (rk != dm) || (val != i) ) {
        fprintf(fp,"(%3d %3d %3d) rk = %d val = %dn", i, j-1, k, rk, val );
      }
    }}}

    fprintf(fp,"\nY  -----------------------------------\n\n");
    for( int k=0; k<NK; k++ ){
    for( int j=0; j<gc; j++ ){
    for (int i=0; i<NI; i++) {
      REAL_TYPE bf = Y[_IDX_S3D(i,j-1,k,NI,NJ,gc)];
      int rk = (int)bf/(int)BASE;
      int val= (int)(bf-rk*BASE);
      if ( (rk != dm) || (val != rsize[3*dm+1]-2+j) ) {
        fprintf(fp,"(%3d %3d %3d) rk = %d val = %d\n", i, 1-j, k, rk, val);
      }
    }}}

    fprintf(fp,"\nZ  -----------------------------------\n\n");
    for( int k=0; k<NK; k++ ){
    for( int j=0; j<gc; j++ ){
    for (int i=0; i<NI; i++) {
      REAL_TYPE bf = Z[_IDX_S3D(i,j-1,k,NI,NJ,gc)];
      int rk = (int)bf/(int)BASE;
      int val= (int)(bf-rk*BASE);
      if ( (rk != dm) || (val != k) ) {
        fprintf(fp,"(%3d %3d %3d) rk = %d val = %d\n",i, 1-j, k, rk, val);
      }
    }}}
  }


  // PLUS DIR
  if (dp >= 0) {
    fprintf(fp,"\nJ_PLUS from rank[%d] ========================================\n",dp);
    fprintf(fp,"\nX  -----------------------------------\n\n");
    for( int k=0; k<NK; k++ ){
    for( int j=0; j<gc; j++ ){
    for (int i=0; i<NI; i++) {
      REAL_TYPE bf = X[_IDX_S3D(i,NJ+j,k,NI,NJ,gc)];
      int rk = (int)bf/(int)BASE;
      int val= (int)(bf-rk*BASE);
      if ( (rk != dp) || (val != i) ) {
        fprintf(fp,"(%3d %3d %3d) rk = %d val = %d\n", i, NJ+j, k, rk, val);
      }
    }}}

    fprintf(fp,"\nY  -----------------------------------\n\n");
    for( int k=0; k<NK; k++ ){
    for( int j=0; j<gc; j++ ){
    for (int i=0; i<NI; i++) {
      REAL_TYPE bf = Y[_IDX_S3D(i,NJ+j,k,NI,NJ,gc)];
      int rk = (int)bf/(int)BASE;
      int val= (int)(bf-rk*BASE);
      if ( (rk != dp) || (val != 1+j) ) {
        fprintf(fp,"(%3d %3d %3d) rk = %d val = %d\n", i, NJ+j, k, rk, val);
      }
    }}}

    fprintf(fp,"\nZ  -----------------------------------\n\n");
    for( int k=0; k<NK; k++ ){
    for( int j=0; j<gc; j++ ){
    for (int i=0; i<NI; i++) {
      REAL_TYPE bf = Z[_IDX_S3D(i,NJ+j,k,NI,NJ,gc)];
      int rk = (int)bf/(int)BASE;
      int val= (int)(bf-rk*BASE);
      if ( (rk != dp) || (val != k) ) {
        fprintf(fp,"(%3d %3d %3d) rk = %d val = %d\n", i, NJ+j, k, rk, val);
      }
    }}}
  }

  return true;
}



////////////////////////////////////////////////////////////////////////////////
bool checkTransferZ(const int* rsize,
                    const int gc,
                    const int myRank,
                    const REAL_TYPE* X,
                    const REAL_TYPE* Y,
                    const REAL_TYPE* Z,
                    const int* nID,
                    FILE* fp)
{
  int NI = rsize[3*myRank+0];
  int NJ = rsize[3*myRank+1];
  int NK = rsize[3*myRank+2];

  int dm = nID[K_minus];
  int dp = nID[K_plus];

  // MINUS DIR
  if (dm >= 0) {
    fprintf(fp,"\nK_MINUS from rank[%d] ========================================\n",dm);
    fprintf(fp,"\nX  -----------------------------------\n\n");
    for( int k=0; k<gc; k++ ){
    for( int j=0; j<NJ; j++ ){
    for (int i=0; i<NI; i++) {
      REAL_TYPE bf = X[_IDX_S3D(i,j,k-1,NI,NJ,gc)];
      int rk = (int)bf/(int)BASE;
      int val= (int)(bf-rk*BASE);
      if ( (rk != dm) || (val != i) ) {
        fprintf(fp,"(%3d %3d %3d) rk = %d val = %d\n", i, j, k-1, rk, val);
      }
    }}}

    fprintf(fp,"\nY  -----------------------------------\n\n");
    for( int k=0; k<gc; k++ ){
    for( int j=0; j<NJ; j++ ){
    for (int i=0; i<NI; i++) {
      REAL_TYPE bf = Y[_IDX_S3D(i,j,k-1,NI,NJ,gc)];
      int rk = (int)bf/(int)BASE;
      int val= (int)(bf-rk*BASE);
      if ( (rk != dm) || (val != j) ) {
        fprintf(fp,"(%3d %3d %3d) rk = %d val = %d\n", i, j, k-1, rk, val);
      }
    }}}

    fprintf(fp,"\nZ  -----------------------------------\n\n");
    for( int k=0; k<gc; k++ ){
    for( int j=0; j<NJ; j++ ){
    for (int i=0; i<NI; i++) {
      REAL_TYPE bf = Z[_IDX_S3D(i,j,k-1,NI,NJ,gc)];
      int rk = (int)bf/(int)BASE;
      int val= (int)(bf-rk*BASE);
      if ( (rk != dm) || (val != rsize[3*dm+2]-2+k) ) {
        fprintf(fp,"(%3d %3d %3d) rk = %d val = %d\n", i, j, k-1, rk, val);
      }
    }}}
  }


  // PLUS DIR
  if (dp >= 0) {
    fprintf(fp,"\nK_PLUS from rank[%d] ========================================\n",dp);
    fprintf(fp,"\nX  -----------------------------------\n\n");
    for( int k=0; k<gc; k++ ){
    for( int j=0; j<NJ; j++ ){
    for (int i=0; i<NI; i++) {
      REAL_TYPE bf = X[_IDX_S3D(i,j,NK+k,NI,NJ,gc)];
      int rk = (int)bf/(int)BASE;
      int val= (int)(bf-rk*BASE);
      if ( (rk != dp) || (val != i) ) {
        fprintf(fp,"(%3d %3d %3d) rk = %d val = %d\n",i, j, NK+k, rk, val);
      }
    }}}

    fprintf(fp,"\nY  -----------------------------------\n\n");
    for( int k=0; k<gc; k++ ){
    for( int j=0; j<NJ; j++ ){
    for (int i=0; i<NI; i++) {
      REAL_TYPE bf = Y[_IDX_S3D(i,j,NK+k,NI,NJ,gc)];
      int rk = (int)bf/(int)BASE;
      int val= (int)(bf-rk*BASE);
      if ( (rk != dp) || (val != j) ) {
        fprintf(fp,"(%3d %3d %3d) rk = %d val = %d\n",i, j, NK+k, rk, val);
      }
    }}}

    fprintf(fp,"\nZ  -----------------------------------\n\n");
    for( int k=0; k<gc; k++ ){
    for( int j=0; j<NJ; j++ ){
    for (int i=0; i<NI; i++) {
      REAL_TYPE bf = Z[_IDX_S3D(i,j,NK+k,NI,NJ,gc)];
      int rk = (int)bf/(int)BASE;
      int val= (int)(bf-rk*BASE);
      if ( (rk != dp) || (val != 1+k) ) {
        fprintf(fp,"(%3d %3d %3d) rk = %d val = %d\n",i, j, NK+k, rk, val);
      }
    }}}
  }

  return true;
}



////////////////////////////////////////////////////////////////////////////////
// @program commtest
// @brief communication of guide cell test
// @param [in] argc  number of aruguents on command line
// @param [in] argv  each string of argument
int main(int argc, char * argv[]) {

  // Check command line arguments
  if ( argc != 8 && argc != 11) {
    printf("Usage:\n");
    printf("\t$ mpirun -np N %s nx ny nz gc grid index d_mode\n", argv[0]);
    printf("\t       or \n");
    printf("\t$ mpirun -np N %s nx ny nz gc grid index d_mode dvx dvy dvz\n", argv[0]);
    printf("\tnx  ; Number of cells for x-dir.\n");
    printf("\tny  ; Number of cells for y-dir.\n");
    printf("\tnz  ; Number of cells for z-dir.\n");
    printf("\tgc  ; Number of guide cell.\n");
    printf("\tgrid; grid type (node, cell).\n");
    printf("\tindex; index type (F, C).\n");
    printf("\td_mode; division mode (IJK, IJ, JK).\n");
    printf("\tdvx, dvy, dvz : Number of division for each axis.\n");
    return -1;
  }


  int np, myRank, proc_grp = 0;
  int nID[NOFACE];           ///< 隣接ランクテーブル
  int head[3];               ///< 各領域の先頭インデクス
  int lsz[3]={0,0,0};        ///< ローカルサイズ
  int G_div[3]={0, 0, 0};    ///< 領域全体の分割数
  int G_size[3]={0,0,0};     ///< 領域全体のサイズ
  int d_mode=-1;             ///< 分割方法の指定
  std::string idx_str;       ///< インデクス指定の文字列
  std::string grd_str;       ///< node or cell指定の文字列
  int div_type=-1;           ///< 自動分割 or 分割指定

  if ( argc == 8) div_type = AUTO;
  else div_type = SPEC;

  for (int i=0; i<NOFACE; i++) nID[i]=-1;

  // conversion from ASCII char to digit
  G_size[0] = atoi(argv[1]);
  G_size[1] = atoi(argv[2]);
  G_size[2] = atoi(argv[3]);
  int gc    = atoi(argv[4]);
  char *grid  = argv[5];
  char *index = argv[6];
  char *div   = argv[7];


  if (!strcasecmp(grid, "node")) {
    grd_str = "node";
  }
  else if (!strcasecmp(grid, "cell")) {
    grd_str = "cell";
  }
  else {
    printf("Grid error >> %s\n", grid);
    exit(1);
  }

  if (!strcasecmp(index, "f")) {
    idx_str = "Findex";
  }
  else if (!strcasecmp(index, "c")) {
    idx_str = "Cindex";
  }
  else {
    printf("Index error >> %s\n", index);
    exit(1);
  }

  if (!strcasecmp(div, "ijk")) {
    d_mode=0;
  }
  else if (!strcasecmp(div, "ij")) {
    d_mode=1;
  }
  else if (!strcasecmp(div, "jk")) {
    d_mode=2;
  }
  else {
    printf("DIV error >> %s\n", div);
    exit(1);
  }


  // Initialize MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  if ( argc == 11) {
    G_div[0]    = atoi(argv[8]);
    G_div[1]    = atoi(argv[9]);
    G_div[2]    = atoi(argv[10]);
    if (G_div[0]*G_div[1]*G_div[2] != np) {
      printf("division != np\n");
      exit(0);
    }
  }

  Hostonly_ {
    printf("%d x %d x %d / %d\n", G_size[0], G_size[1], G_size[2], np);
    printf("guide         : %d\n", gc);

    printf("Grid type     : ");
    if (!strcasecmp(grid, "node")) {
      printf("NODE\n");
    }
    else if (!strcasecmp(grid, "cell")) {
      printf("CELL\n");
    }

    printf("Index         : ");
    if (!strcasecmp(index, "f")) {
      printf("F index\n");
    }
    else if (!strcasecmp(index, "c")) {
      printf("C index\n");
    }

    printf("Division mode : ");
    if (d_mode==0) {
      printf("IJK\n");
    }
    else if (d_mode==1) {
      printf("IJ\n");
    }
    else {
      printf("JK\n");
    }

    if ( argc == 11) {
      printf("Domain division  : %d %d %d\n", G_div[0], G_div[1],G_div[2]);
    }
  }



  SubDomain D;
  D.setSubDomain(G_size, gc, np, myRank, proc_grp, MPI_COMM_WORLD,
                 grd_str, idx_str);


  if (div_type == SPEC) {
    if ( !D.setDivision(G_div) ) {
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD, -1);
    }
  }

  if ( !D.findOptimalDivision(d_mode) )
  {
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD, -1);
  }


  if ( !D.createRankTable() ) {
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  // 領域分割数取得 （自動分割の場合）
  D.getGlobalDivision(G_div);

  //自ランクのHeadIndexの取得
  D.getLocalHead(head);

  //自ランクの格子数取得
  D.getLocalSize(lsz);

  //自ランクの隣接ランク番号を取得
  D.getCommTable(nID);



  // 各ランクの情報を共有
  int* asize = new int [np*3];
  MPI_Allgather(lsz, 3, MPI_INT,
                asize, 3, MPI_INT, MPI_COMM_WORLD);


  if (myRank==0) {
    printf("\n");
    printf("Rank :    NI    NJ    NK\n");
    for (int i=0; i<np; i++) {
      printf("[%2d] : %5d %5d %5d\n", i,
         asize[3*i+0],
         asize[3*i+1],
         asize[3*i+2]);
    }
    printf("\n");
  }


  int* rankID = new int [np*6];
  MPI_Allgather(nID, 6, MPI_INT,
                rankID, 6, MPI_INT, MPI_COMM_WORLD);

  if (myRank==0) {
    printf("Rank :  I-  I+  J-  J+  K-  K+\n");
    for (int i=0; i<np; i++) {
      printf("[%2d] : %3d %3d %3d %3d %3d %3d\n", i,
            rankID[6*i+0],
            rankID[6*i+1],
            rankID[6*i+2],
            rankID[6*i+3],
            rankID[6*i+4],
            rankID[6*i+5]);
    }
    printf("\n");
  }


  int NI = lsz[0];
  int NJ = lsz[1];
  int NK = lsz[2];

  if (NI>BASE || NJ>BASE || NK>BASE) {
    printf("Local size exceeds limit %d : Rank %d (%d, %d, %d)\n",
         BASE, myRank, NI, NJ, NK);
    MPI_Abort(MPI_COMM_WORLD, -1);
  }

  //printf("NIJK = %d %d %d\n", NI, NJ, NK);


  // Definition of Array
  REAL_TYPE* X = NULL;  ///< scalar work
  REAL_TYPE* Y = NULL;  ///< scalar work
  REAL_TYPE* Z = NULL;  ///< scalar work
  REAL_TYPE* V = NULL;  ///< vector work
  if ( !(X=alloc_real(lsz, gc)) ) MPI_Abort(MPI_COMM_WORLD, -1);
  if ( !(Y=alloc_real(lsz, gc)) ) MPI_Abort(MPI_COMM_WORLD, -1);
  if ( !(Z=alloc_real(lsz, gc)) ) MPI_Abort(MPI_COMM_WORLD, -1);
  if ( !(V=alloc_real(lsz, gc, 3)) ) MPI_Abort(MPI_COMM_WORLD, -1);


  // 初期値
  /* X, Y, Zの各配列はBASEをシフト値として値をセット
     X配列にはインデクスiとして、myRank*BASE+i, Y, Zも同様
     アドレス計算はC表記、配列の値もC表記のインデクスをいれておく

Cell image
             Rank 0   <--|-->   Rank 1
               0*BASE+(NI-1)_0
                     v
rank0       [NI-2] [NI-1]  [NI]  [NI+1]
      -----+------+------|------+------+---------> i
rank1        [-2]   [-1]    [0]    [1]    [2]
                             ^
                          1*BASE+1

セルの場合、Rank1の要素0にはrank0*BASE+lsz[0] (Rnak0の値)
 BASEの桁の値で送信された値のランクを確認し、下位の数値で
 i方向の位置を確認する

====
Node image アドレスはC表記

            Rank A   <--|-->   Rank B
    rank0           0*BASE+(NI-1)_A
                        v
rankA  [NI-3] [NI-2] [NI-1]  [NI]  [NI+1]
     -----+------+------|------+------+-------> i
rankB   [-2]   [-1]    [0]    [1]    [2]
                               ^
    rank1                   1*BASE+1

ノードの場合、Rank1の要素0にはrank0*BASE+lsz[0]-1 (Rnak0の値)

===
     隣接間通信によりガイドセルを交換すると
     Fortran IndexでRank1の要素0は rank*BASE+lsz[0] (Rnak0の値)となるので
     想定の値になっているかをチェック

     スカラーのチェックはX, Y, Zの各方向をみる
     ベクトルはVにX, Y, Zの値を入れて、隣接間通信を行い、同様にチェックする
  */


  for( int k=0; k<NK; k++ ){
  for( int j=0; j<NJ; j++ ){
  for( int i=0; i<NI; i++ ){
    X[_IDX_S3D(i,j,k,NI,NJ,gc)] = myRank*BASE+i;
  }}}

  for( int k=0; k<NK; k++ ){
  for( int j=0; j<NJ; j++ ){
  for( int i=0; i<NI; i++ ){
    Y[_IDX_S3D(i,j,k,NI,NJ,gc)] = myRank*BASE+j;
  }}}

  for( int k=0; k<NK; k++ ){
  for( int j=0; j<NJ; j++ ){
  for( int i=0; i<NI; i++ ){
    Z[_IDX_S3D(i,j,k,NI,NJ,gc)] = myRank*BASE+k;
  }}}

  for( int k=0; k<NK; k++ ){
  for( int j=0; j<NJ; j++ ){
  for( int i=0; i<NI; i++ ){
    V[_IDX_V3D(i,j,k,0,NI,NJ,NK,gc)] = X[_IDX_S3D(i,j,k,NI,NJ,gc)];
    V[_IDX_V3D(i,j,k,1,NI,NJ,NK,gc)] = Y[_IDX_S3D(i,j,k,NI,NJ,gc)];
    V[_IDX_V3D(i,j,k,2,NI,NJ,NK,gc)] = Z[_IDX_S3D(i,j,k,NI,NJ,gc)];
  }}}


  // Communication identifier for nonblocking
  // 斜め袖通信有/無の両方に対応するにはNOFACE*2で確保する
  MPI_Request req[NOFACE*2];
  for (int i=0; i<NOFACE*2; i++) req[i] = MPI_REQUEST_NULL;


  BrickComm CM;      ///< 通信クラス


  // 通信クラス設定
  if ( !CM.setBrickComm(lsz, gc, MPI_COMM_WORLD, nID, grd_str) ) {
    stamped_printf("\tBrickComm settng error.\n");
    return 0;
  }


  // 通信バッファ確保  # of component = 3
  if  ( !CM.init(3) ) {
    stamped_printf("\tBrickComm initialize error.\n");
    return 0;
  }

  // 袖通信(scalar)
  if (!strcasecmp(grid, "node")) {
    CM.Comm_S_node(X, gc, req);
    CM.Comm_S_wait_node(X, gc, req);
  }
  else {
    CM.Comm_S_cell(X, gc, req);
    CM.Comm_S_wait_cell(X, gc, req);
  }

  if (!strcasecmp(grid, "node")) {
    CM.Comm_S_node(Y, gc, req);
    CM.Comm_S_wait_node(Y, gc, req);
  }
  else {
    CM.Comm_S_cell(Y, gc, req);
    CM.Comm_S_wait_cell(Y, gc, req);
  }

  if (!strcasecmp(grid, "node")) {
    CM.Comm_S_node(Z, gc, req);
    CM.Comm_S_wait_node(Z, gc, req);
  }
  else {
    CM.Comm_S_cell(Z, gc, req);
    CM.Comm_S_wait_cell(Z, gc, req);
  }

  if (!strcasecmp(grid, "node")) {
    CM.Comm_V_node(V, gc, req);
    CM.Comm_V_wait_node(V, gc, req);
  }
  else {
    CM.Comm_V_cell(V, gc, req);
    CM.Comm_V_wait_cell(V, gc, req);
  }



  // scalar
  char fname[30];
  sprintf( fname, "log_S_%03d.txt", myRank );
  FILE* fp=fopen(fname, "w");
  if ( !checkTransferX(asize,
                      gc,
                      myRank,
                      X, Y, Z,
                      nID, fp) ) MPI_Abort(MPI_COMM_WORLD, -1);

  if ( !checkTransferY(asize,
                      gc,
                      myRank,
                      X, Y, Z,
                      nID, fp) ) MPI_Abort(MPI_COMM_WORLD, -1);

  if ( !checkTransferZ(asize,
                      gc,
                      myRank,
                      X, Y, Z,
                      nID, fp) ) MPI_Abort(MPI_COMM_WORLD, -1);
  fclose(fp);



  // Vector
  sprintf( fname, "log_V_%03d.txt", myRank );
  fp=fopen(fname, "w");

  size_t n = (NI+2*gc) * (NJ+2*gc) * (NK+2*gc) ;

  if ( !checkTransferX(asize,
                      gc,
                      myRank,
                      &V[0],
                      &V[n],
                      &V[2*n],
                      nID, fp) ) MPI_Abort(MPI_COMM_WORLD, -1);

  if ( !checkTransferY(asize,
                      gc,
                      myRank,
                      &V[0],
                      &V[n],
                      &V[2*n],
                      nID, fp) ) MPI_Abort(MPI_COMM_WORLD, -1);

  if ( !checkTransferZ(asize,
                      gc,
                      myRank,
                      &V[0],
                      &V[n],
                      &V[2*n],
                      nID, fp) ) MPI_Abort(MPI_COMM_WORLD, -1);
  fclose(fp);

  // deallocate
  delete [] X;
  delete [] Y;
  delete [] Z;
  delete [] V;
  delete [] asize;
  delete [] rankID;


  // finalize MPI
  MPI_Finalize();
  Hostonly_ printf("Successfully terminated.\n\n");

  return 0;
}
