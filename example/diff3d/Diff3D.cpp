//
//  Diff3D.cpp
//
//  Created by keno on 2016/06/26.
//  Copyright © 2016年 keno. All rights reserved.
//

#include "mpi.h"
#include "Diff3D.h"
#include "utility.h"
#include "mydebug.h"
#include "table.h"
#include "comm.h"

/* C header */
int allocateArray(const int* sz, const int gc, REAL_TYPE* q, REAL_TYPE* w);

int read_config(Phys_Param* p, Cntl_Param* c, const char* fname, const int myrank);




// Fortran functions
extern "C" {
extern void initialize_(int* sz,
                        int* gc,
                        REAL_TYPE* q,
                        REAL_TYPE* w);

extern void bc_(int* sz,
                int* gc,
                REAL_TYPE* q,
                REAL_TYPE* dh,
                REAL_TYPE* org,
                int* tbl);

extern void euler_explicit_(int* sz,
                            int* gc,
                            REAL_TYPE* q,
                            REAL_TYPE* w,
                            REAL_TYPE* dh,
                            REAL_TYPE* dt,
                            REAL_TYPE* alpha,
                            REAL_TYPE* res);

extern void write_sph_(int* sz,
                       int* gc,
                       int* step,
                       REAL_TYPE* time,
                       REAL_TYPE* dh,
                       REAL_TYPE* org,
                       char* fname,
                       REAL_TYPE* q);

};


// @program diff3d
// @brief 3D unsteady diffusion equation
// @param [in] argc  number of aruguents on command line
// @param [in] argv  each string of argument
int main(int argc, char * argv[])
{
  Phys_Param P_phys;
  Cntl_Param P_cntl;

  // Check command line arguments
  if ( argc != 2) {
    usage();
    return 0;
  }

  int np=0;        ///< Number of processes
  int myrank=-1;   ///< Rank number
  int mode=0;      ///< 0-serial, 1-parallel

  int dsz[3]={0,0,0};    ///< Number of cell sizes of entire domain
  int comm_tbl[6]={-1,-1,-1,-1,-1,-1};  ///< Communication table
  REAL_TYPE origin[3]={0.0, 0.0, 0.0};  ///< global origin

  Hostonly_  {
    if ( sizeof(REAL_TYPE) == 8 ) {
      printf("Real is DOUBLE\n");
    }
    else {
      printf("Real is FLOAT\n");
    }
  }


  // Initialize MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  if ( np > 1 ) mode = 1;

  // config file
  if ( !read_config(&P_phys, &P_cntl, argv[1], myrank) )
  {
    MPI_Finalize();
    return 1;
  }

  // compute communication table
  if ( !mktable(np, myrank, mode, &P_phys, &P_cntl, comm_tbl, origin) )
  {
    MPI_Finalize();
    return 1;
  }




  // Definition of Array
  REAL_TYPE* q = NULL;  ///< unknown variable
  REAL_TYPE* w = NULL;  ///< scalar work
  int gc = 1;           ///< guide cell


  // Array allocation
  int m_sz[3] = {P_cntl.sz[0], P_cntl.sz[1], P_cntl.sz[2]};
  size_t len = (size_t)( (m_sz[0]+2*gc) * (m_sz[1]+2*gc) * (m_sz[2]+2*gc) );

  if ( !(q=alloc_real(len)) ) return 0;
  if ( !(w=alloc_real(len)) ) return 0;


  // Buffer allocation
  REAL_TYPE* face_xms = NULL;  // X- direction send
  REAL_TYPE* face_xmr = NULL;  // X- direction recv
  REAL_TYPE* face_xps = NULL;  // X+ direction send
  REAL_TYPE* face_xpr = NULL;  // X+ direction recv
  REAL_TYPE* face_yms = NULL;  // Y- direction send
  REAL_TYPE* face_ymr = NULL;  // Y- direction recv
  REAL_TYPE* face_yps = NULL;  // Y+ direction send
  REAL_TYPE* face_ypr = NULL;  // Y+ direction recv
  REAL_TYPE* face_zms = NULL;  // Z- direction send
  REAL_TYPE* face_zmr = NULL;  // Z- direction recv
  REAL_TYPE* face_zps = NULL;  // Z+ direction send
  REAL_TYPE* face_zpr = NULL;  // Z+ direction recv

  size_t face_sz[3];
  face_sz[0] = (size_t)(m_sz[1]+2*gc) * (size_t)(m_sz[2]+2*gc);
  face_sz[1] = (size_t)(m_sz[0]+2*gc) * (size_t)(m_sz[2]+2*gc);
  face_sz[2] = (size_t)(m_sz[0]+2*gc) * (size_t)(m_sz[1]+2*gc);

  if ( !(face_xms=(REAL_TYPE*)malloc(sizeof(REAL_TYPE)*face_sz[0])) ) {
    printf("fail to allocate memory\n");
    exit(0);
  }
  if ( !(face_xmr=(REAL_TYPE*)malloc(sizeof(REAL_TYPE)*face_sz[0])) ) {
    printf("fail to allocate memory\n");
    exit(0);
  }
  if ( !(face_xps=(REAL_TYPE*)malloc(sizeof(REAL_TYPE)*face_sz[0])) ) {
    printf("fail to allocate memory\n");
    exit(0);
  }
  if ( !(face_xpr=(REAL_TYPE*)malloc(sizeof(REAL_TYPE)*face_sz[0])) ) {
    printf("fail to allocate memory\n");
    exit(0);
  }

  if ( !(face_yms=(REAL_TYPE*)malloc(sizeof(REAL_TYPE)*face_sz[1])) ) {
    printf("fail to allocate memory\n");
    exit(0);
  }
  if ( !(face_ymr=(REAL_TYPE*)malloc(sizeof(REAL_TYPE)*face_sz[1])) ) {
    printf("fail to allocate memory\n");
    exit(0);
  }
  if ( !(face_yps=(REAL_TYPE*)malloc(sizeof(REAL_TYPE)*face_sz[1])) ) {
    printf("fail to allocate memory\n");
    exit(0);
  }
  if ( !(face_ypr=(REAL_TYPE*)malloc(sizeof(REAL_TYPE)*face_sz[1])) ) {
    printf("fail to allocate memory\n");
    exit(0);
  }

  if ( !(face_zms=(REAL_TYPE*)malloc(sizeof(REAL_TYPE)*face_sz[2])) ) {
    printf("fail to allocate memory\n");
    exit(0);
  }
  if ( !(face_zmr=(REAL_TYPE*)malloc(sizeof(REAL_TYPE)*face_sz[2])) ) {
    printf("fail to allocate memory\n");
    exit(0);
  }
  if ( !(face_zps=(REAL_TYPE*)malloc(sizeof(REAL_TYPE)*face_sz[2])) ) {
    printf("fail to allocate memory\n");
    exit(0);
  }
  if ( !(face_zpr=(REAL_TYPE*)malloc(sizeof(REAL_TYPE)*face_sz[2])) ) {
    printf("fail to allocate memory\n");
    exit(0);
  }

  // Communication identifier for nonblocking
  MPI_Request req[12];
  for (int i=0; i<12; i++) req[i] = MPI_REQUEST_NULL;

  // initialization
  initialize_(m_sz, &gc, q, w);

  if ( P_cntl.blocking == 1) {
    SyncGC_nonblocking(comm_tbl,
                  face_xms, face_xmr, face_xps, face_xpr,
                  face_yms, face_ymr, face_yps, face_ypr,
                  face_zms, face_zmr, face_zps, face_zpr,
                  q, m_sz[0], m_sz[1], m_sz[2], gc, gc, req);
    wait_nonblocking(comm_tbl,
                  face_xms, face_xmr, face_xps, face_xpr,
                  face_yms, face_ymr, face_yps, face_ypr,
                  face_zms, face_zmr, face_zps, face_zpr,
                  q, m_sz[0], m_sz[1], m_sz[2], gc, gc, req);

    SyncGC_nonblocking(comm_tbl,
                  face_xms, face_xmr, face_xps, face_xpr,
                  face_yms, face_ymr, face_yps, face_ypr,
                  face_zms, face_zmr, face_zps, face_zpr,
                  w, m_sz[0], m_sz[1], m_sz[2], gc, gc, req);
    wait_nonblocking(comm_tbl,
                  face_xms, face_xmr, face_xps, face_xpr,
                  face_yms, face_ymr, face_yps, face_ypr,
                  face_zms, face_zmr, face_zps, face_zpr,
                  w, m_sz[0], m_sz[1], m_sz[2], gc, gc, req);

  }
  else {
    SyncGC_blocking(comm_tbl,
                  face_xms, face_xmr, face_xps, face_xpr,
                  face_yms, face_ymr, face_yps, face_ypr,
                  face_zms, face_zmr, face_zps, face_zpr,
                  q, m_sz[0], m_sz[1], m_sz[2], gc, gc);

    SyncGC_blocking(comm_tbl,
                  face_xms, face_xmr, face_xps, face_xpr,
                  face_yms, face_ymr, face_yps, face_ypr,
                  face_zms, face_zmr, face_zps, face_zpr,
                  w, m_sz[0], m_sz[1], m_sz[2], gc, gc);
  }

  REAL_TYPE time = 0.0;
  REAL_TYPE res;

  char fname[20];
  sprintf( fname, "result_%04d.sph", myrank );

  // time step loop
  Hostonly_ printf("    step         time     residual\n");

  for (int step=1; step<=P_cntl.laststep; step++)
  {
    time += P_phys.dt;

    // boundary condition
    bc_(m_sz, &gc, q, &P_phys.dh, P_phys.org, comm_tbl);

    if ( P_cntl.blocking == 1) {
      SyncGC_nonblocking(comm_tbl,
                    face_xms, face_xmr, face_xps, face_xpr,
                    face_yms, face_ymr, face_yps, face_ypr,
                    face_zms, face_zmr, face_zps, face_zpr,
                    q, m_sz[0], m_sz[1], m_sz[2], gc, gc, req);
      wait_nonblocking(comm_tbl,
                    face_xms, face_xmr, face_xps, face_xpr,
                    face_yms, face_ymr, face_yps, face_ypr,
                    face_zms, face_zmr, face_zps, face_zpr,
                    q, m_sz[0], m_sz[1], m_sz[2], gc, gc, req);
    }
    else {
      SyncGC_blocking(comm_tbl,
                    face_xms, face_xmr, face_xps, face_xpr,
                    face_yms, face_ymr, face_yps, face_ypr,
                    face_zms, face_zmr, face_zps, face_zpr,
                    q, m_sz[0], m_sz[1], m_sz[2], gc, gc);
    }

    // time marching
    res = 0.0;
    euler_explicit_(m_sz, &gc, q, w, &P_phys.dh, &P_phys.dt, &P_phys.alpha, &res);

    REAL_TYPE tmp=res;

    if ( sizeof(REAL_TYPE) == 8 ) {
      MPI_Allreduce(&tmp, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    }
    else {
      MPI_Allreduce(&tmp, &res, 1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
    }

    res = sqrt(res);

    if ( P_cntl.blocking == 1) {
      SyncGC_nonblocking(comm_tbl,
                      face_xms, face_xmr, face_xps, face_xpr,
                      face_yms, face_ymr, face_yps, face_ypr,
                      face_zms, face_zmr, face_zps, face_zpr,
                      q, m_sz[0], m_sz[1], m_sz[2], gc, gc, req);
      wait_nonblocking(comm_tbl,
                      face_xms, face_xmr, face_xps, face_xpr,
                      face_yms, face_ymr, face_yps, face_ypr,
                      face_zms, face_zmr, face_zps, face_zpr,
                      q, m_sz[0], m_sz[1], m_sz[2], gc, gc, req);
    }
    else {
      SyncGC_blocking(comm_tbl,
                      face_xms, face_xmr, face_xps, face_xpr,
                      face_yms, face_ymr, face_yps, face_ypr,
                      face_zms, face_zmr, face_zps, face_zpr,
                      q, m_sz[0], m_sz[1], m_sz[2], gc, gc);
    }



    // display history
    Hostonly_ printf("%8d %12.6e %12.6e\n", step, time, res);

    // file out
    if (0 == step-(step/P_cntl.fileout)*P_cntl.fileout)
    {
      write_sph_(m_sz, &gc, &step, &time, &P_phys.dh, P_phys.org, fname, q);
    }

  }

  // post




  // release
  if ( !q ) free(q);
  if ( !w ) free(w);

  if ( !face_xms ) free(face_xms);
  if ( !face_xmr ) free(face_xmr);
  if ( !face_xps ) free(face_xps);
  if ( !face_xpr ) free(face_xpr);
  if ( !face_yms ) free(face_yms);
  if ( !face_ymr ) free(face_ymr);
  if ( !face_yps ) free(face_yps);
  if ( !face_ypr ) free(face_ypr);
  if ( !face_zms ) free(face_zms);
  if ( !face_zmr ) free(face_zmr);
  if ( !face_zps ) free(face_zps);
  if ( !face_zpr ) free(face_zpr);


  MPI_Finalize();

  return 0;
}






// @fn read_config
// @brief read parameters in config file
// @param [in,out] p      Struct of physical parameters
// @param [in,out] c      Struct of control parameters
// @param [in]     fname  File name of config
// @param [in]     myRank Rank ID of my own
// @retval 0-fail, 1-success
int read_config(Phys_Param* p, Cntl_Param* c, const char* fname, const int myrank)
{
  FILE* fp = NULL;

  // config fileのオープン
  Hostonly_
  {
    if ( !(fp=fopen(fname, "r")) )
    {
      stamped_printf("\tSorry, can't open 'config.txt' file.\n");
      return 0; // < ここの動作には並列処理時のバグが残っています。何でしょう？
    }

    // for string
    char buf[1024];

    fscanf(fp, "%d %s", &c->sz[0], buf);
    fscanf(fp, "%d %s", &c->sz[1], buf);
    fscanf(fp, "%d %s", &c->sz[2], buf);

    fscanf(fp, "%d %s", &c->div[0], buf);
    fscanf(fp, "%d %s", &c->div[1], buf);
    fscanf(fp, "%d %s", &c->div[2], buf);

    //fscanf(fp, "%e %s", &p->dh, buf);
    fscanf(fp, "%e %s", &p->dt, buf);
    fscanf(fp, "%e %s", &p->alpha, buf);

    fscanf(fp, "%d %s", &c->laststep, buf);
    fscanf(fp, "%d %s", &c->fileout, buf);
    fscanf(fp, "%d %s", &c->blocking, buf);

    fclose(fp);
    fp = NULL;

    p->dh = 1.0 / (REAL_TYPE)c->sz[0];

    printf("\nPARAMETERS\n");
    printf("\tnx   = %4d ;         Number of cells for x-dir.\n", c->sz[0]);
    printf("\tny   = %4d ;         Number of cells for y-dir.\n", c->sz[1]);
    printf("\tnz   = %4d ;         Number of cells for z-dir.\n", c->sz[2]);
    printf("\tXdiv = %4d ;         Number of division for x-dir.\n", c->div[0]);
    printf("\tYdiv = %4d ;         Number of division for y-dir.\n", c->div[1]);
    printf("\tZdiv = %4d ;         Number of division for z-dir.\n", c->div[2]);
    printf("\tdh   = %10.6e ; Mesh width (isotropic, non-dimensional)\n", p->dh);
    printf("\tdt   = %10.6e ; Time increment of time marching (non-dimensional)\n", p->dt);
    printf("\talpha= %10.6e ; Coefficient of diffusion (non-dimensional)\n", p->alpha);
    printf("\tlaststep = %8d ; Time step to calculate\n", c->laststep);
    printf("\tfileout  = %8d ; Interval for writing a file\n\n", c->fileout);
    printf("\tblocking = ");
    if ( c->blocking == 0) {
      printf("blocking\n\n");
    }
    else {
      printf("nonblocking\n\n");
    }
  }

  // buffer to pack variables
  REAL_TYPE* param=NULL;
  if ( !(param=(REAL_TYPE*)malloc(sizeof(REAL_TYPE)*3)) ) {
    printf("fail to allocate memory\n");
    exit(0);
  }

  int* cntl=NULL;
  if ( !(cntl=(int*)malloc(sizeof(int)*9)) ) {
    printf("fail to allocate memory\n");
    exit(0);
  }

  // packing to send buffer
  Hostonly_ {
    param[0] = p->dh;
    param[1] = p->dt;
    param[2] = p->alpha;

    cntl[0] = c->sz[0];
    cntl[1] = c->sz[1];
    cntl[2] = c->sz[2];
    cntl[3] = c->div[0];
    cntl[4] = c->div[1];
    cntl[5] = c->div[2];
    cntl[6] = c->laststep;
    cntl[7] = c->fileout;
    cntl[8] = c->blocking;
  }

  // Broadcat to all processes
  if ( sizeof(REAL_TYPE) == 8 ) {
    MPI_Bcast(param, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  else {
    MPI_Bcast(param, 3, MPI_FLOAT, 0, MPI_COMM_WORLD);
  }

  p->dh    = param[0];
  p->dt    = param[1];
  p->alpha = param[2];

#if 0
  printf("rank=%3d %10.6e %10.6e %10.6e\n", myrank, p->dh, p->dt, p->alpha);
#endif


  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Bcast(cntl, 9, MPI_INT, 0, MPI_COMM_WORLD);

  c->sz[0]    = cntl[0];
  c->sz[1]    = cntl[1];
  c->sz[2]    = cntl[2];
  c->div[0]   = cntl[3];
  c->div[1]   = cntl[4];
  c->div[2]   = cntl[5];
  c->laststep = cntl[6];
  c->fileout  = cntl[7];
  c->blocking = cntl[8];


#if 0
  printf("rank=%3d %4d %4d %4d %4d %4d %4d %8d %8d\n",
         myrank,
         c->sz[0], c->sz[1], c->sz[2],
         c->div[0], c->div[1], c->div[2],
         c->laststep, c->fileout);
#endif

  // release buffer
  if ( !param ) free(param);
  if ( !cntl  ) free(cntl);

  return 1;
}
