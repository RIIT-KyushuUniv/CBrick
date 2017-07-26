//
//  table.c
//  sample
//
//  Created by keno on 2016/06/18.
//  Copyright © 2016 keno. All rights reserved.
//

#include "mpi.h"
#include <stdlib.h>
#include "table.h"

// @fn mktable
// @brief prepare communication table
// @param [in]      np      Number of processes
// @param [in]      myrank  Rank number
// @param [in]      mode    0-serial, 1-parallel
// @param [in,out]  m_phys  Physical parameter
// @param [in,out]  m_cntl  Control parameter
// @param [out]     my_tbl  Communication table
// @param [in]      g_org   Global origin
// @retval 0-fail, 1-success
int mktable(const int np, const int myrank, const int mode, Phys_Param* m_phys, Cntl_Param* m_cntl, int my_tbl[6], REAL_TYPE g_org[3])
{
  const int dsz[3] = {m_cntl->sz[0], m_cntl->sz[1], m_cntl->sz[2]};
  const int div[3] = {m_cntl->div[0], m_cntl->div[1], m_cntl->div[2]};


  // check parallel mode
  if ( mode == 1) {

    // check arguments
    if ( div[0]*div[1]*div[2] != np || np == 1 ) {
      if ( myrank == 0 ) printf("\t  Xdiv*Ydiv*Zdiv must be n > 1.\n");
      return 0; // <= if return value is -1, what's happen?
    }
  }

  // display mode
  if ( mode == 1 ) {
    if ( myrank == 0 ) printf("\t## Parallel mode  %d processes ##\n\n", np);
  }
  else {
    printf("\t## Serial mode ##\n\n");
    return 0;
  }


  // domain parameters
  int my_size[3] = {0,0,0};
  int my_head[3] = {0,0,0};

  // create comm table
  setTable(div, np, myrank, dsz, my_size, my_head, my_tbl);

  m_cntl->sz[0] = my_size[0];
  m_cntl->sz[1] = my_size[1];
  m_cntl->sz[2] = my_size[2];

  m_phys->org[0] = g_org[0] + (REAL_TYPE)(my_head[0]-1) * m_phys->dh;
  m_phys->org[1] = g_org[1] + (REAL_TYPE)(my_head[1]-1) * m_phys->dh;
  m_phys->org[2] = g_org[2] + (REAL_TYPE)(my_head[2]-1) * m_phys->dh;

  // gather info and printout
  display(np, myrank, my_size, my_head, my_tbl, m_phys->org);

  return 1;
}



/*
 * @fn setTable
 * @brief calculate parameters and preserve in tbl[]
 * @param [in]  div    number of division
 * @param [in]  np     number of processes
 * @param [in]  myrank rank numbber
 * @param [in]  dsz    number of cell size of a whole domain
 * @param [out] my_sz  size of own subdomain
 * @param [out] my_hd  head of own subdomain
 * @param [out] my_tbl communication table of own subdomain
 */
void setTable(const int div[3],
              const int np,
              const int myrank,
              const int dsz[3],
              int my_sz[3],
              int my_hd[3],
              int my_tbl[6])
{
  // region parameters
  // div[3] ; Number of division for each axis
  // np     ; Number of whole processes
  // dsz[3] ; Number of cells for each axis
  // myrank ; Rank number

  // sz[3] ; default domain size for each dir.
  int sz[3];
  sz[0] = dsz[0]/div[0];  if ( dsz[0] != sz[0]*div[0] ) sz[0] +=1;
  sz[1] = dsz[1]/div[1];  if ( dsz[1] != sz[1]*div[1] ) sz[1] +=1;
  sz[2] = dsz[2]/div[2];  if ( dsz[2] != sz[2]*div[2] ) sz[2] +=1;

  // compute head position
  int head[3];
  int comm[6];
  getHead(sz, div, dsz, head, myrank, comm);

  // registration
  my_sz[0] = sz[0];
  my_sz[1] = sz[1];
  my_sz[2] = sz[2];

  my_hd[0] = head[0];
  my_hd[1] = head[1];
  my_hd[2] = head[2];

  my_tbl[0] = comm[0];
  my_tbl[1] = comm[1];
  my_tbl[2] = comm[2];
  my_tbl[3] = comm[3];
  my_tbl[4] = comm[4];
  my_tbl[5] = comm[5];
}



/*
 * @fn getHead
 * @brief Calculate head position, subdomain size, and comm table
 * @param [in,out] sz     subdomain size
 * @param [in]     dv     number of division
 * @param [in]     dsz    number of cell size of a whole domain
 * @param [out]    head   head position
 * @param [in]     myrank rank numbber
 * @param [out]    comm   communication table
 */
void getHead(int* sz, const int* dv, const int* dsz, int* head, const int myrank, int* comm)
{
  // head index and domain size for each direction
  int* head_x = NULL;
  int* head_y = NULL;
  int* head_z = NULL;
  int* sz_x = NULL;
  int* sz_y = NULL;
  int* sz_z = NULL;
  int* rt=NULL;

  // sizes which consider guide cells
  if ( !(head_x=(int*)malloc(sizeof(int)*(dv[0]+2))) ) {
    printf("fail to allocate memory\n");
    exit(0);
  }
  if ( !(head_y=(int*)malloc(sizeof(int)*(dv[1]+2))) ) {
    printf("fail to allocate memory\n");
    exit(0);
  }
  if ( !(head_z=(int*)malloc(sizeof(int)*(dv[2]+2))) ) {
    printf("fail to allocate memory\n");
    exit(0);
  }

  if ( !(sz_x=(int*)malloc(sizeof(int)*(dv[0]+2))) ) {
    printf("fail to allocate memory\n");
    exit(0);
  }
  if ( !(sz_y=(int*)malloc(sizeof(int)*(dv[1]+2))) ) {
    printf("fail to allocate memory\n");
    exit(0);
  }
  if ( !(sz_z=(int*)malloc(sizeof(int)*(dv[2]+2))) ) {
    printf("fail to allocate memory\n");
    exit(0);
  }

  // array size => (0:dv[0]+1, 0:dv[1]+1, 0:dv[2]+1)
  if ( !(rt=(int*)malloc(sizeof(int)*(dv[0]+2)*(dv[1]+2)*(dv[2]+2))) ) {
    printf("fail to allocate memory\n");
    exit(0);
  }

  // initialization
  for (int i=0; i<dv[0]+2; i++) {
    head_x[i] = 0;
    sz_x[i] = 0;
  }
  for (int i=0; i<dv[1]+2; i++) {
    head_y[i] = 0;
    sz_y[i] = 0;
  }
  for (int i=0; i<dv[2]+2; i++) {
    head_z[i] = 0;
    sz_z[i] = 0;
  }

#pragma omp parallel for schedule(static) collapse(2)
  for (int k=0; k<dv[2]+2; k++) {
    for (int j=0; j<dv[1]+2; j++) {
      for (int i=0; i<dv[0]+2; i++) {
        size_t m = _F_IDX_S3D(i, j, k, dv[0], dv[1], dv[2], 1);
        rt[m] = -1;
      }
    }
  }

  // rank number
  int c = 0;
  for (int k=1; k<=dv[2]; k++) {
    for (int j=1; j<=dv[1]; j++) {
      for (int i=1; i<=dv[0]; i++) {
        size_t m = _F_IDX_S3D(i, j, k, dv[0], dv[1], dv[2], 1);
        rt[m] = c++;
      }
    }
  }

  // crate table : head position and size
  head_x[1] = 1;
  sz_x[1] = sz[0];
  for (int i=2; i<=dv[0]; i++) {
    int tmp = head_x[i-1]+sz_x[i-1];
    if ( tmp <= dsz[0] ) {
      head_x[i] = tmp;
      if (head_x[i]+sz[0] <= dsz[0]) {
        sz_x[i] = sz[0];
      }
      else {
        sz_x[i] = dsz[0] - head_x[i] + 1;
      }
    }
  }

  head_y[1] = 1;
  sz_y[1] = sz[1];
  for (int j=2; j<=dv[1]; j++) {
    int tmp = head_y[j-1]+sz_y[j-1];
    if ( tmp <= dsz[1] ) {
      head_y[j] = tmp;
      if (head_y[j]+sz[1] <= dsz[1]) {
        sz_y[j] = sz[1];
      }
      else {
        sz_y[j] = dsz[1] - head_y[j] + 1;
      }
    }
  }

  head_z[1] = 1;
  sz_z[1] = sz[2];
  for (int k=2; k<=dv[2]; k++) {
    int tmp = head_z[k-1]+sz_z[k-1];
    if ( tmp <= dsz[2] ) {
      head_z[k] = tmp;
      if (head_z[k]+sz[2] <= dsz[2]) {
        sz_z[k] = sz[2];
      }
      else {
        sz_z[k] = dsz[2] - head_z[k] + 1;
      }
    }
  }



  // determine head and size from table
  for (int k=1; k<=dv[2]; k++) {
    for (int j=1; j<=dv[1]; j++) {
      for (int i=1; i<=dv[0]; i++) {
        size_t m = _F_IDX_S3D(i, j, k, dv[0], dv[1], dv[2], 1);

        if (rt[m] == myrank) {
          head[0] = head_x[i];
          head[1] = head_y[j];
          head[2] = head_z[k];

          sz[0] = sz_x[i];
          sz[1] = sz_y[j];
          sz[2] = sz_z[k];

          comm[X_minus] = rt[_F_IDX_S3D(i-1, j  , k  , dv[0], dv[1], dv[2], 1)];
          comm[X_plus]  = rt[_F_IDX_S3D(i+1, j  , k  , dv[0], dv[1], dv[2], 1)];
          comm[Y_minus] = rt[_F_IDX_S3D(i  , j-1, k  , dv[0], dv[1], dv[2], 1)];
          comm[Y_plus]  = rt[_F_IDX_S3D(i  , j+1, k  , dv[0], dv[1], dv[2], 1)];
          comm[Z_minus] = rt[_F_IDX_S3D(i  , j  , k-1, dv[0], dv[1], dv[2], 1)];
          comm[Z_plus]  = rt[_F_IDX_S3D(i  , j  , k+1, dv[0], dv[1], dv[2], 1)];
        }
      }
    }
  }


  // release working array
  if ( !head_x ) free(head_x);
  if ( !head_y ) free(head_y);
  if ( !head_z ) free(head_z);
  if ( !sz_x ) free(sz_x);
  if ( !sz_y ) free(sz_y);
  if ( !sz_z ) free(sz_z);
  if ( !rt   ) free(rt);
}



/*
 * @fn display
 * @brief gather table and display
 * @param [in] np      number of processes
 * @param [in] myrank  rank numbber
 * @param [in] m_size  size of own subdomain
 * @param [in] m_head  head of own subdomain
 * @param [in] m_tbl   communication table of own subdomain
 * @param [in] m_org   Local origin
 */
void display(const int np,
             const int myrank,
             const int* m_size,
             const int* m_head,
             const int* m_tbl,
             const REAL_TYPE* m_org)
{
  int* g_size=NULL;  ///< buffer to gather size
  int* g_head=NULL;  ///< buffer to gather head
  int* g_tbl =NULL;  ///< buffer to gather table
  REAL_TYPE* g_org=NULL;

  if ( !(g_size=(int*)malloc(sizeof(int)*np*3)) ) {
    printf("fail to allocate memory\n");
    exit(0);
  }

  if ( !(g_head=(int*)malloc(sizeof(int)*np*3)) ) {
    printf("fail to allocate memory\n");
    exit(0);
  }

  if ( !(g_tbl=(int*)malloc(sizeof(int)*np*6)) ) {
    printf("fail to allocate memory\n");
    exit(0);
  }

  if ( !(g_org=(REAL_TYPE*)malloc(sizeof(REAL_TYPE)*np*3)) ) {
    printf("fail to allocate memory\n");
    exit(0);
  }

  int my_sz[3];
  int my_hd[3];
  int my_cm[6];
  REAL_TYPE my_og[3];

  for (int i=0; i<3; i++) {
    my_sz[i]=m_size[i];
    my_hd[i]=m_head[i];
    my_og[i]=m_org[i];
  }

  for (int i=0; i<6; i++) {
    my_cm[i] = m_tbl[i];
  }

  // gather size from all processes to root process
  MPI_Gather(my_sz, 3, MPI_INT, g_size, 3, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(my_hd, 3, MPI_INT, g_head, 3, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(my_cm, 6, MPI_INT, g_tbl,  6, MPI_INT, 0, MPI_COMM_WORLD);
  if ( sizeof(REAL_TYPE) == 8 ) {
    MPI_Gather(my_og, 3, MPI_DOUBLE, g_org, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
  else {
    MPI_Gather(my_og, 3, MPI_FLOAT, g_org, 3, MPI_FLOAT, 0, MPI_COMM_WORLD);
  }

  // display, root only
  if ( myrank == 0 ) {
    printf("Rank : size (  X,   Y,   Z)  : head (  X,   Y,   Z)  : connection [ X-  X+  Y-  Y+  Z-  Z+] : Local origin (X, Y, Z)\n");

    for (int m=0; m<np; m++) {
      printf("%4d :      (%3d, %3d, %3d)  :      (%3d, %3d, %3d)  :            [%3d %3d %3d %3d %3d %3d] : (%10.5e, %10.5e, %10.5e) \n",
             m,
             g_size[3*m],
             g_size[3*m+1],
             g_size[3*m+2],
             g_head[3*m],
             g_head[3*m+1],
             g_head[3*m+2],
             g_tbl[6*m],
             g_tbl[6*m+1],
             g_tbl[6*m+2],
             g_tbl[6*m+3],
             g_tbl[6*m+4],
             g_tbl[6*m+5],
             g_org[3*m],
             g_org[3*m+1],
             g_org[3*m+2]
             );
    }
    printf("\n\n");
  }


  // release array
  if ( !g_size ) free(g_size);
  if ( !g_head ) free(g_head);
  if ( !g_tbl  ) free(g_tbl);
  if ( !g_org  ) free(g_org);
}
