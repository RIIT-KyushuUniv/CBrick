// Diff3D

#include "Diff3D.h"
#include <CB_SubDomain.h>

// @fn alloc_int
// @brief allocatin of integer array
// @param [in]     sz Size of array
// @retval pointer
int* alloc_int(const size_t sz)
{
  int* p=NULL;
  if ( !(p = new int[sz]) ) {
    printf("fail to allocate memory\n");
  }
  return p;
}


// @fn alloc_real
// @brief allocatin of real array
// @param [in]     sz Size of array
// @retval pointer
// @note Type of real is defined by REAL_TYPE
REAL_TYPE* alloc_real(const size_t sz)
{
  REAL_TYPE* p=NULL;
  if ( !(p = new REAL_TYPE[sz]) ) {
    printf("fail to allocate memory\n");
  }
  return p;
}


// @fn usage
// @brief print usage
void usage()
{
  printf("Usage:\n");

  printf("\t$ mpirun -np n diff3d config_file\n\n");
  printf("\tContents of config_file\n");
  printf("\tnx;       Number of cells for x-dir.\n");
  printf("\tny;       Number of cells for y-dir.\n");
  printf("\tnz;       Number of cells for z-dir.\n");
  printf("\tDivMode;  Division mode (0-IJK, 1-JK).\n");
  printf("\tXdiv;     Number of division for x-dir.\n");
  printf("\tYdiv;     Number of division for y-dir.\n");
  printf("\tZdiv;     Number of division for z-dir.\n");
  printf("\tdt;       Time increment of time marching (non-dimensional)\n");
  printf("\talpha;    Coefficient of diffusion (non-dimensional)\n");
  printf("\tlaststep; Time step to calculate\n");
  printf("\tfileout;  Interval for writing a file\n");
  printf("\tcomm.;    Communication mode (0-blocking, 1-nonblocking)\n\n");
}


// @fn read_config
// @brief read parameters in config file
// @param [in,out] p      Struct of physical parameters
// @param [in,out] c      Struct of control parameters
// @param [out]    m_sz   全計算領域の要素数
// @param [out]    div    領域分割数
// @param [in]     fname  File name of config
// @param [in]     myRank Rank ID of my own
// @retval 0-fail, 1-success
bool read_config(Phys_Param* p,
                 Cntl_Param* c,
                 int* m_sz,
                 int* div,
                 const char* fname,
                 const int myRank)
{
  // config fileのオープン
  Hostonly_
  {
    FILE* fp = NULL;

    if ( !(fp=fopen(fname, "r")) )
    {
      stamped_printf("\tSorry, can't open 'config.txt' file.\n");
      return false;
    }

    // for string
    char buf[1024];
    float m_dt=0.0, m_alp=0.0;

    fscanf(fp, "%d %s", &m_sz[0], buf);
    fscanf(fp, "%d %s", &m_sz[1], buf);
    fscanf(fp, "%d %s", &m_sz[2], buf);

    fscanf(fp, "%d %s", &div[0], buf);
    fscanf(fp, "%d %s", &div[1], buf);
    fscanf(fp, "%d %s", &div[2], buf);

    fscanf(fp, "%e %s", &m_dt, buf);
    fscanf(fp, "%e %s", &m_alp, buf);

    fscanf(fp, "%d %s", &c->laststep, buf);
    fscanf(fp, "%d %s", &c->fileout, buf);
    fscanf(fp, "%d %s", &c->blocking, buf);
    fscanf(fp, "%d %s", &c->div_mode, buf);

    fclose(fp);
    fp = NULL;

    p->dh = 1.0 / (REAL_TYPE)m_sz[0];
    p->dt = (REAL_TYPE)m_dt;
    p->alpha = (REAL_TYPE)m_alp;

    printf("\nPARAMETERS\n");
    printf("\tnx   = %4d ;         Number of cells for x-dir.\n", m_sz[0]);
    printf("\tny   = %4d ;         Number of cells for y-dir.\n", m_sz[1]);
    printf("\tnz   = %4d ;         Number of cells for z-dir.\n", m_sz[2]);

    printf("\tDivMode =      %s\n", (c->div_mode==0)?"IJK":"JK");

    if (div[0]*div[1]*div[2] == 0) {
      printf("\tNumber of division is not specified\n");
    }
    else {
      printf("\tXdiv = %4d ;         Number of division for x-dir.\n", div[0]);
      printf("\tYdiv = %4d ;         Number of division for y-dir.\n", div[1]);
      printf("\tZdiv = %4d ;         Number of division for z-dir.\n", div[2]);
    }

    printf("\tdh   = %10.6e ; Mesh width (isotropic, non-dimensional)\n", p->dh);
    printf("\tdt   = %10.6e ; Time increment of time marching (non-dimensional)\n", p->dt);
    printf("\talpha= %10.6e ; Coefficient of diffusion (non-dimensional)\n", p->alpha);
    printf("\tlaststep = %8d ; Time step to calculate\n", c->laststep);
    printf("\tfileout  = %8d ; Interval for writing a file\n\n", c->fileout);
    printf("\tComm.    = ");
    if ( c->blocking == 0) {
      printf("blocking\n\n");
    }
    else {
      printf("non-blocking\n\n");
    }
  } // Hostonly


  // buffer to pack variables
  REAL_TYPE* param=NULL;
  if ( !(param= new REAL_TYPE[3]) ) {
    printf("fail to allocate memory\n");
    return false;
  }

  int* cntl=NULL;
  if ( !(cntl= new int[10]) ) {
    printf("fail to allocate memory\n");
    return false;
  }

  // packing to send buffer
  Hostonly_ {
    param[0] = p->dh;
    param[1] = p->dt;
    param[2] = p->alpha;

    cntl[0] = m_sz[0];
    cntl[1] = m_sz[1];
    cntl[2] = m_sz[2];
    cntl[3] = div[0];
    cntl[4] = div[1];
    cntl[5] = div[2];
    cntl[6] = c->laststep;
    cntl[7] = c->fileout;
    cntl[8] = c->blocking;
    cntl[9] = c->div_mode;

  }

  // Broadcat to all processes
  if ( sizeof(REAL_TYPE) == _SIZE_DOUBLE_ ) {
    MPI_Bcast(param, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  else {
    MPI_Bcast(param, 3, MPI_FLOAT, 0, MPI_COMM_WORLD);
  }

  p->dh    = param[0];
  p->dt    = param[1];
  p->alpha = param[2];



  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Bcast(cntl, 10, MPI_INT, 0, MPI_COMM_WORLD);

  m_sz[0]     = cntl[0];
  m_sz[1]     = cntl[1];
  m_sz[2]     = cntl[2];
  div[0]      = cntl[3];
  div[1]      = cntl[4];
  div[2]      = cntl[5];
  c->laststep = cntl[6];
  c->fileout  = cntl[7];
  c->blocking = cntl[8];
  c->div_mode = cntl[9];


  // release buffer
  if ( param ) delete [] param;
  if ( cntl  ) delete [] cntl;

  return true;
}


// @program diff3d
// @brief 3D unsteady diffusion equation
// @param [in] argc  number of aruguents on command line
// @param [in] argv  each string of argument
int main(int argc, char * argv[])
{
  Phys_Param P_phys;
  Cntl_Param P_cntl;

  int np=0;         ///< Number of processes
  int myRank=-1;    ///< Rank number
  int proc_grp = 0; ///< プロセスグループ番号
  int mode=0;       ///< 0-serial, 1-parallel
  int gc = 1;       ///< ガイドセル幅=1
  int div_type=0;   ///< 分割指定 (0-自動、1-指定)

  int dsz[3]={0,0,0};    ///< 全計算領域の要素数
  REAL_TYPE origin[3]={0.0, 0.0, 0.0};  ///< global origin
  int m_dv[3]={0,0,0};   ///< 領域分割数


  // Initialize MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  if ( np > 1 ) mode = 1;

  // Check command line arguments
  if ( argc != 2) {
    usage();
    return 0;
  }

  Hostonly_  {
    if ( sizeof(REAL_TYPE) == _SIZE_DOUBLE_ ) {
      printf("Real is DOUBLE\n");
    }
    else {
      printf("Real is FLOAT\n");
    }
  }


  // config file
  if ( !read_config(&P_phys, &P_cntl, dsz, m_dv, argv[1], myRank) )
  {
    MPI_Finalize();
    return 1;
  }

  // コンフィギュレーションファイルの分割数がゼロでなければ、指定分割
  if (m_dv[0]*m_dv[1]*m_dv[2] != 0)  div_type = 1;


  if ( div_type == 1 && m_dv[0]*m_dv[1]*m_dv[2] != np) {
    printf("\tThe number of proceees does not agree with the div[] size.\n");
    return 1;
  }

  SubDomain D(dsz, gc, np, myRank, proc_grp, MPI_COMM_WORLD, "cell", "Findex");

  // 分割数指定
  if (div_type == 1) D.setDivision(m_dv);

  if ( !D.findOptimalDivision(P_cntl.div_mode) ) MPI_Abort(MPI_COMM_WORLD, -1);
  if ( !D.createRankTable() ) MPI_Abort(MPI_COMM_WORLD, -1);


  // origin of each SubDomain, Fortran Index
  P_phys.org[0] = origin[0] + (D.head[0]-1)*P_phys.dh;
  P_phys.org[1] = origin[1] + (D.head[1]-1)*P_phys.dh;
  P_phys.org[2] = origin[2] + (D.head[2]-1)*P_phys.dh;


  int lsz[3]={0,0,0}; //サブドメインサイズ
  lsz[0] = D.size[0];
  lsz[1] = D.size[1];
  lsz[2] = D.size[2];

  // Definition of Array
  REAL_TYPE* q = NULL;  ///< unknown variable
  REAL_TYPE* w = NULL;  ///< scalar work


  // Array allocation
  size_t len = (size_t)( (lsz[0]+2*gc) * (lsz[1]+2*gc) * (lsz[2]+2*gc) );

  if ( !(q=alloc_real(len)) ) MPI_Abort(MPI_COMM_WORLD, -1);
  if ( !(w=alloc_real(len)) ) MPI_Abort(MPI_COMM_WORLD, -1);


  // 通信バッファ確保 scalar
  D.initComm(1);

  // Communication identifier for nonblocking
  MPI_Request req[12];
  for (int i=0; i<12; i++) req[i] = MPI_REQUEST_NULL;

  // initialization
  initialize_(lsz, &gc, q, w);


  D.Comm_S_nonblocking(q, gc, req);
  D.Comm_S_wait_nonblocking(q, gc, req);

  D.Comm_S_nonblocking(w, gc, req);
  D.Comm_S_wait_nonblocking(w, gc, req);



  REAL_TYPE time = 0.0;
  REAL_TYPE res;

  char fname[20];
  sprintf( fname, "result_%04d.sph", myRank );

  // time step loop
  Hostonly_ printf("    step         time     residual\n");

  for (int step=1; step<=P_cntl.laststep; step++)
  {
    time += P_phys.dt;

    // boundary condition
    bc_(lsz, &gc, q, &P_phys.dh, P_phys.org, D.comm_tbl);


    D.Comm_S_nonblocking(q, gc, req);
    D.Comm_S_wait_nonblocking(q, gc, req);


    // time marching
    res = 0.0;
    euler_explicit_(lsz, &gc, q, w, &P_phys.dh, &P_phys.dt, &P_phys.alpha, &res);

    REAL_TYPE tmp=res;

    if ( sizeof(REAL_TYPE) == _SIZE_DOUBLE_ ) {
      MPI_Allreduce(&tmp, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    else {
      MPI_Allreduce(&tmp, &res, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    }

    res = sqrt(res);


    D.Comm_S_nonblocking(q, gc, req);
    D.Comm_S_wait_nonblocking(q, gc, req);


    // display history
    Hostonly_ printf("%8d %12.6e %12.6e\n", step, time, res);

    // file out
    if (0 == step-(step/P_cntl.fileout)*P_cntl.fileout)
    {
      write_sph_(lsz, &gc, &step, &time, &P_phys.dh, P_phys.org, fname, q);
    }

  }

  // post


  // release
  if ( q ) delete [] q;
  if ( w ) delete [] w;


  MPI_Finalize();

  return 0;
}
