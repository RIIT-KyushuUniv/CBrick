//
//  commtest.cpp
//
//  Created by keno on 2017/07/01.
//  Copyright © 2017 keno. All rights reserved.
//

// Compile
// $ gcc commtest.cpp -o commtest

// Execution
// $ mpirun -np X commtest nx ny nz gc grid
// (ex)
// $ mpirun -np 4 commtest 64 64 64 2 node
// $ mpirun -np 4 commtest 65 65 65 2 cell

// 袖通信テスト

#include <CB_SubDomain.h>
#include <CB_Comm.h>
#include <stdlib.h>

// @fn alloc_real
// @brief allocatin of real array
// @param [in]     sz  Size of array
// @param [in]     val initial value
// @retval pointer
// @note Type of real is defined by REAL_TYPE
REAL_TYPE* alloc_real(const size_t sz, REAL_TYPE val=REAL_TYPE(0))
{
  REAL_TYPE* p=NULL;
  if ( !(p = new REAL_TYPE[sz]) ) {
    printf("fail to allocate memory\n");
  }
  for( size_t i=0;i<sz;i++ )
  {
    p[i] = val;
  }
  return p;
}

void printArray(FILE *fp, const char *title, REAL_TYPE *a, int imax, int jmax, int kmax, int lmax, int gc)
{
  fprintf(fp, "\n");
  fprintf(fp, "***** %s *****\n", title);
  fprintf(fp, "imax:%d\n", imax);
  fprintf(fp, "jmax:%d\n", jmax);
  fprintf(fp, "kmax:%d\n", kmax);
  fprintf(fp, "lmax:%d\n", lmax);
  fprintf(fp, "gc  :%d\n", gc);

  for( int l=0;l<lmax;l++ )
  {
    fprintf(fp, "** l=%d\n", l);
    for( int k=kmax+gc-1;k>=0-gc;k-- )
    {
      fprintf(fp, "\n* k=%d\n", k);
      for( int j=jmax+gc-1;j>=0-gc;j-- )
      {
        fprintf(fp, "%4d, ", j);
        for( int i=0-gc;i<imax+gc;i++ )
        {
          fprintf(fp, "%10d,", int(a[_IDX_V3D(i,j,k,l,imax,jmax,kmax,gc)]));
        }
        fprintf(fp, "\n");
      }
    }
  }
}

// @program commtest
// @brief communication of guide cell test
// @param [in] argc  number of aruguents on command line
// @param [in] argv  each string of argument
int main(int argc, char * argv[]) {

  // Check command line arguments
  if ( argc != 6) {
    printf("Usage:\n");
    printf("\t$ mpirun -np N %s nx ny nz gc grid\n", argv[0]);
    printf("\tnx  ; Number of cells for x-dir.\n");
    printf("\tny  ; Number of cells for y-dir.\n");
    printf("\tnz  ; Number of cells for z-dir.\n");
    printf("\tgc  ; Number of guide cell.\n");
    printf("\tgrid; grid type (node, cell).\n");
    return -1;
  }

  int m_sz[3], np, myRank;
  int proc_grp = 0;
  int div_mode=0;
  int nID[6]={0,0,0,0,0,0};  ///< 隣接ランクテーブル

  // conversion from ASCII char to digit
  m_sz[0]    = atoi(argv[1]);
  m_sz[1]    = atoi(argv[2]);
  m_sz[2]    = atoi(argv[3]);
  int gc     = atoi(argv[4]);
  char *grid = argv[5];

  // Initialize MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  Hostonly_ printf("%d x %d x %d / %d\n", m_sz[0], m_sz[1], m_sz[2], np);

  // priorityはデフォルト
  SubDomain D(m_sz, gc, np, myRank, proc_grp, MPI_COMM_WORLD, grid, "Cindex");

  // 最適な分割方法を検索
  if ( !D.findOptimalDivision(div_mode) )
  {
    MPI_Barrier(MPI_COMM_WORLD);
    exit(-1);
  }

  // create rank table
  if ( !D.createRankTable() )
  {
    MPI_Barrier(MPI_COMM_WORLD);
    exit(-1);
  }


  // subdomain size
  int lsz[3]={0,0,0}; // サブドメインサイズ
  D.getLocalSize(lsz);

  int imax = lsz[0];
  int jmax = lsz[1];
  int kmax = lsz[2];
  size_t len = (size_t)( (imax+2*gc) * (jmax+2*gc) * (kmax+2*gc) );

  D.getCommTable(nID);

  // Definition of Array
  REAL_TYPE* s = NULL;  ///< scalar work
  REAL_TYPE* v = NULL;  ///< vector work
  if ( !(s=alloc_real(len  ,REAL_TYPE(-1))  ) ) MPI_Abort(MPI_COMM_WORLD, -1);
  if ( !(v=alloc_real(len*3,REAL_TYPE(-1))) ) MPI_Abort(MPI_COMM_WORLD, -1);

  // デバッグ用初期値をセット
  int sidx=0;
  if( !strcmp(grid,"node") ) sidx=1; //nodeのときは下側開始インデクスを1にする
  int cnt_s = (myRank+1)*10000;
  for( int k=sidx;k<kmax;k++ ){
  for( int j=sidx;j<jmax;j++ ){
  for( int i=sidx;i<imax;i++ ){
    s[_IDX_S3D(i,j,k,imax,jmax,gc)] = cnt_s++;
  }}}
  for( int l=0;l<3;l++ ){
    int cnt_v = (myRank+1)*100000 + l*10000;
    for( int k=sidx;k<kmax;k++ ){
    for( int j=sidx;j<jmax;j++ ){
    for( int i=sidx;i<imax;i++ ){
      v[_IDX_V3D(i,j,k,l,imax,jmax,kmax,gc)] = cnt_v++;
    }}}
  }

  // Communication identifier for nonblocking
  // 斜め袖通信有/無の両方に対応するにはNOFACE*2で確保する
  MPI_Request req[NOFACE*2];
  for (int i=0; i<NOFACE*2; i++) req[i] = MPI_REQUEST_NULL;


  BrickComm CM;      ///< 通信クラス

  // 通信クラス設定
  if ( !CM.setBrickComm(lsz, gc, MPI_COMM_WORLD, nID, "cell") ) {
    stamped_printf("\tBrickComm settng error.\n");
    return 0;
  }

  // 通信バッファ確保  # of component = 3
  if  ( !CM.init(3) ) {
    stamped_printf("\tBrickComm initialize error.\n");
    return 0;
  }

  // 袖通信(scalar)
  CM.Comm_S_nonblocking(s, gc, req);
  CM.Comm_S_wait_nonblocking(s, gc, req);

  // 袖通信(vector)
  CM.Comm_V_nonblocking(v, gc, req);
  CM.Comm_V_wait_nonblocking(v, gc, req);

  // log
  char fname1[512], fname2[512];
  sprintf(fname1, "commS_%04d.log", myRank);
  sprintf(fname2, "commV_%04d.log", myRank);
  FILE *fp1 = fopen(fname1, "wt");
  FILE *fp2 = fopen(fname2, "wt");
  printArray(fp1, "scalar", s, imax, jmax, kmax, 1, gc);
  printArray(fp2, "vector", v, imax, jmax, kmax, 3, gc);
  fclose(fp1);
  fclose(fp2);

  // deallocate
  delete [] s;
  delete [] v;

  // finalize MPI
  MPI_Finalize();
  Hostonly_ printf("Successfully terminated.\n\n");

  return 0;
}
