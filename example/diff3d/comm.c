//
//  comm.c
//  diff3dp
//
//  Created by keno on 2016/07/03.
//  Copyright © 2016年 keno. All rights reserved.
//

#include <mpi.h>
#include "comm.h"

// @brief NonBlocking commuication
void SyncGC_nonblocking(const int* tbl,
                     REAL_TYPE* xms, REAL_TYPE* xmr,
                     REAL_TYPE* xps, REAL_TYPE* xpr,
                     REAL_TYPE* yms, REAL_TYPE* ymr,
                     REAL_TYPE* yps, REAL_TYPE* ypr,
                     REAL_TYPE* zms, REAL_TYPE* zmr,
                     REAL_TYPE* zps, REAL_TYPE* zpr,
                     REAL_TYPE* src,
                     const int imax, const int jmax, const int kmax, const int gc, const int gc_comm, MPI_Request *req)
{
  // Communication identifier
  for (int i=0; i<12; i++) req[i] = MPI_REQUEST_NULL;

  // X direction

  int nIDm = tbl[X_minus];
  int nIDp = tbl[X_plus];
  size_t msz = gc * (jmax+2*gc) * (kmax+2*gc);

  packX(src, imax, jmax, kmax, gc, gc_comm, xms, xps, nIDm, nIDp);
  IsendIrecv(xms, xmr, xps, xpr, msz, nIDm, nIDp, &req[0]);

  // Y direction
  msz = gc * (imax+2*gc) * (kmax+2*gc);
  nIDm = tbl[Y_minus];
  nIDp = tbl[Y_plus];

  packY(src, imax, jmax, kmax, gc, gc_comm, yms, yps, nIDm, nIDp);
  IsendIrecv(yms, ymr, yps, ypr, msz, nIDm, nIDp, &req[4]);

  // Z direction
  msz = gc * (imax+2*gc) * (jmax+2*gc);
  nIDm = tbl[Z_minus];
  nIDp = tbl[Z_plus];

  packZ(src, imax, jmax, kmax, gc, gc_comm, zms, zps, nIDm, nIDp);
  IsendIrecv(zms, zmr, zps, zpr, msz, nIDm, nIDp, &req[8]);
}

void IsendIrecv(REAL_TYPE* ms, REAL_TYPE* mr, REAL_TYPE* ps, REAL_TYPE* pr, size_t msz, int nIDm, int nIDp, MPI_Request *req)
{
  // Identifier
  MPI_Request r0 = MPI_REQUEST_NULL;
  MPI_Request r1 = MPI_REQUEST_NULL;
  MPI_Request r2 = MPI_REQUEST_NULL;
  MPI_Request r3 = MPI_REQUEST_NULL;

  int tag_p=0;
  int tag_m=0;

  // Recieve Minus
  if ( sizeof(REAL_TYPE) == 8 ) {
    if ( nIDm >= 0 ) MPI_Irecv(mr, msz, MPI_DOUBLE, nIDm, tag_m, MPI_COMM_WORLD, &r1);
  }
  else {
    if ( nIDm >= 0 ) MPI_Irecv(mr, msz, MPI_FLOAT, nIDm, tag_m, MPI_COMM_WORLD, &r1);
  }

  // Recieve Plus
  if ( sizeof(REAL_TYPE) == 8 ) {
    if ( nIDp >= 0 ) MPI_Irecv(pr, msz, MPI_DOUBLE, nIDp, tag_p, MPI_COMM_WORLD, &r3);
  }
  else {
    if ( nIDp >= 0 ) MPI_Irecv(pr, msz, MPI_FLOAT, nIDp, tag_p, MPI_COMM_WORLD, &r3);
  }

  // Send Plus
  if ( sizeof(REAL_TYPE) == 8 ) {
    if ( nIDp >= 0 ) MPI_Isend(ps, msz, MPI_DOUBLE, nIDp, tag_p, MPI_COMM_WORLD, &r2);
  }
  else {
    if ( nIDp >= 0 ) MPI_Isend(ps, msz, MPI_FLOAT, nIDp, tag_p, MPI_COMM_WORLD, &r2);
  }

  // Send Minus
  if ( sizeof(REAL_TYPE) == 8 ) {
    if ( nIDm >= 0 ) MPI_Isend(ms, msz, MPI_DOUBLE, nIDm, tag_m, MPI_COMM_WORLD, &r0);
  }
  else {
    if ( nIDm >= 0 ) MPI_Isend(ms, msz, MPI_FLOAT, nIDm, tag_m, MPI_COMM_WORLD, &r0);
  }


  req[0] = r0;
  req[1] = r1;
  req[2] = r2;
  req[3] = r3;
}

// @brief wait for IsendIrecv()
void wait_nonblocking(const int* tbl,
                    REAL_TYPE* xms, REAL_TYPE* xmr,
                    REAL_TYPE* xps, REAL_TYPE* xpr,
                    REAL_TYPE* yms, REAL_TYPE* ymr,
                    REAL_TYPE* yps, REAL_TYPE* ypr,
                    REAL_TYPE* zms, REAL_TYPE* zmr,
                    REAL_TYPE* zps, REAL_TYPE* zpr,
                    REAL_TYPE* dest,
                    const int imax, const int jmax, const int kmax, const int gc, const int gc_comm, MPI_Request *req)
{
  MPI_Status stat[4];

  //// X face ////
  int nIDm = tbl[X_minus];
  int nIDp = tbl[X_plus];
  MPI_Waitall( 4, &req[0], stat );
  unpackX(dest, imax, jmax, kmax, gc, gc_comm, xmr, xpr, nIDm, nIDp);


  //// Y face ////
  nIDm = tbl[Y_minus];
  nIDp = tbl[Y_plus];
  MPI_Waitall( 4, &req[4], stat );
  unpackY(dest, imax, jmax, kmax, gc, gc_comm, ymr, ypr, nIDm, nIDp);


  //// Z face ////
  nIDm = tbl[Z_minus];
  nIDp = tbl[Z_plus];
  MPI_Waitall( 4, &req[8], stat );
  unpackZ(dest, imax, jmax, kmax, gc, gc_comm, zmr, zpr, nIDm, nIDp);
}


// @brief Blocking commuication
void SyncGC_blocking(const int* tbl,
                     REAL_TYPE* xms, REAL_TYPE* xmr,
                     REAL_TYPE* xps, REAL_TYPE* xpr,
                     REAL_TYPE* yms, REAL_TYPE* ymr,
                     REAL_TYPE* yps, REAL_TYPE* ypr,
                     REAL_TYPE* zms, REAL_TYPE* zmr,
                     REAL_TYPE* zps, REAL_TYPE* zpr,
                     REAL_TYPE* src,
                     const int imax, const int jmax, const int kmax, const int gc, const int gc_comm)
{
  // X direction

  int nIDm = tbl[X_minus];
  int nIDp = tbl[X_plus];
  size_t msz = gc * (jmax+2*gc) * (kmax+2*gc);

  packX(src, imax, jmax, kmax, gc, gc_comm, xms, xps, nIDm, nIDp);
  send_and_recv(xms, xmr, xps, xpr, msz, nIDm, nIDp);
  //sendrecv(xms, xmr, xps, xpr, msz, nIDm, nIDp);
  unpackX(src, imax, jmax, kmax, gc, gc_comm, xmr, xpr, nIDm, nIDp);


  // Y direction
  msz = gc * (imax+2*gc) * (kmax+2*gc);
  nIDm = tbl[Y_minus];
  nIDp = tbl[Y_plus];

  packY(src, imax, jmax, kmax, gc, gc_comm, yms, yps, nIDm, nIDp);
  send_and_recv(yms, ymr, yps, ypr, msz, nIDm, nIDp);
  //sendrecv(yms, ymr, yps, ypr, msz, nIDm, nIDp);
  unpackY(src, imax, jmax, kmax, gc, gc_comm, ymr, ypr, nIDm, nIDp);

  // Z direction
  msz = gc * (imax+2*gc) * (jmax+2*gc);
  nIDm = tbl[Z_minus];
  nIDp = tbl[Z_plus];

  packZ(src, imax, jmax, kmax, gc, gc_comm, zms, zps, nIDm, nIDp);
  send_and_recv(zms, zmr, zps, zpr, msz, nIDm, nIDp);
  //sendrecv(zms, zmr, zps, zpr, msz, nIDm, nIDp);
  unpackZ(src, imax, jmax, kmax, gc, gc_comm, zmr, zpr, nIDm, nIDp);

}

void send_and_recv(REAL_TYPE* ms, REAL_TYPE* mr, REAL_TYPE* ps, REAL_TYPE* pr, size_t msz, int nIDm, int nIDp)
{
  // Plus side of subdomain
  int tag_p=0;
  MPI_Status *stat_p;

  if ( sizeof(REAL_TYPE) == 8 ) {
    if ( nIDp >= 0 ) MPI_Send(ps, msz, MPI_DOUBLE, nIDp, tag_p, MPI_COMM_WORLD);
    if ( nIDp >= 0 ) MPI_Recv(pr, msz, MPI_DOUBLE, nIDp, tag_p, MPI_COMM_WORLD, stat_p);
  }
  else {
    if ( nIDp >= 0 ) MPI_Send(ps, msz, MPI_FLOAT, nIDp, tag_p, MPI_COMM_WORLD);
    if ( nIDp >= 0 ) MPI_Recv(pr, msz, MPI_FLOAT, nIDp, tag_p, MPI_COMM_WORLD, stat_p);
  }


  // Minus side of subdomain
  int tag_m=0;
  MPI_Status *stat_m;

  if ( sizeof(REAL_TYPE) == 8 ) {
    if ( nIDm >= 0 ) MPI_Recv(mr, msz, MPI_DOUBLE, nIDm, tag_m, MPI_COMM_WORLD, stat_m);
    if ( nIDm >= 0 ) MPI_Send(ms, msz, MPI_DOUBLE, nIDm, tag_m, MPI_COMM_WORLD);
  }
  else {
    if ( nIDm >= 0 ) MPI_Recv(mr, msz, MPI_FLOAT, nIDm, tag_m, MPI_COMM_WORLD, stat_m);
    if ( nIDm >= 0 ) MPI_Send(ms, msz, MPI_FLOAT, nIDm, tag_m, MPI_COMM_WORLD);
  }

}

void sendrecv(REAL_TYPE* ms, REAL_TYPE* mr, REAL_TYPE* ps, REAL_TYPE* pr, size_t msz, int nIDm, int nIDp)
{
  // Plus side of subdomain
  int tag_p=0;
  MPI_Status *stat_p;

  if ( sizeof(REAL_TYPE) == 8 ) {
    if ( nIDp >= 0 ) MPI_Sendrecv(ps, msz, MPI_DOUBLE, nIDp, tag_p,
                                  pr, msz, MPI_DOUBLE, nIDp, tag_p,
                                  MPI_COMM_WORLD, stat_p);
  }
  else {
    if ( nIDp >= 0 ) MPI_Sendrecv(ps, msz, MPI_FLOAT, nIDp, tag_p,
                                  pr, msz, MPI_FLOAT, nIDp, tag_p,
                                  MPI_COMM_WORLD, stat_p);
  }


  // Minus side of subdomain
  int tag_m=0;
  MPI_Status *stat_m;

  if ( sizeof(REAL_TYPE) == 8 ) {
    if ( nIDm >= 0 ) MPI_Sendrecv(ms, msz, MPI_DOUBLE, nIDm, tag_m,
                                  mr, msz, MPI_DOUBLE, nIDm, tag_m,
                                  MPI_COMM_WORLD, stat_m);
  }
  else {
    if ( nIDm >= 0 ) MPI_Sendrecv(ms, msz, MPI_FLOAT, nIDm, tag_m,
                                  mr, msz, MPI_FLOAT, nIDm, tag_m,
                                  MPI_COMM_WORLD, stat_m);
  }
}


/*
 * @brief pack send data for X direction
 * @param [in]  array   source array
 * @param [in]  imax    size x
 * @param [in]  jmax    size y
 * @param [in]  kmax    size z
 * @param [in]  vc      guide cell size
 * @param [in]  vc_comm number of guide cell layer to be sent
 * @param [out] sendm   send buffer of X- direction
 * @param [out] sendp   send buffer of X+ direction
 * @param [in]  nIDm    Rank number of X- direction
 * @param [in]  nIDp    Rank number of X+ direction
 */
void packX(const REAL_TYPE *array, const int imax, const int jmax, const int kmax, const int vc, const int vc_comm,
           REAL_TYPE *sendm, REAL_TYPE *sendp, const int nIDm, const int nIDp)
{
  int ix = imax;
  int jx = jmax;
  int kx = kmax;
  int gc = vc;
  int gc_comm = vc_comm;

  if( nIDm >= 0 )
  {
#pragma omp parallel for firstprivate(ix, jx, kx, gc, gc_comm) schedule(static) collapse(2)
    for( int k=0-gc_comm;k<kx+gc_comm;k++ ){
      for( int j=0-gc_comm;j<jx+gc_comm;j++ ){
        for( int i=0;i<gc_comm;i++ ){
          sendm[_IDXFX(i,j,k,0,jx,kx,gc_comm)] = array[_IDX_S3D(i,j,k,ix,jx,kx,gc)];
        }
      }
    }
  }

  if( nIDp >= 0 )
  {
#pragma omp parallel for firstprivate(ix, jx, kx, gc, gc_comm) schedule(static) collapse(2)
    for( int k=0-gc_comm;k<kx+gc_comm;k++ ){
      for( int j=0-gc_comm;j<jx+gc_comm;j++ ){
        for( int i=ix-gc_comm;i<ix;i++ ){
          sendp[_IDXFX(i,j,k,ix-gc_comm,jx,kx,gc_comm)] = array[_IDX_S3D(i,j,k,ix,jx,kx,gc)];
        }
      }
    }
  }
}


/*
 * @brief unpack send data for X direction
 * @param [in,out]  array   dest array
 * @param [in]  imax    size x
 * @param [in]  jmax    size y
 * @param [in]  kmax    size z
 * @param [in]  vc      guide cell size
 * @param [in]  vc_comm number of guide cell layer to be sent
 * @param [in]  recvm   recv buffer of X- direction
 * @param [in]  recvp   recv buffer of X+ direction
 * @param [in]  nIDm    Rank number of X- direction
 * @param [in]  nIDp    Rank number of X+ direction
 */
void unpackX(REAL_TYPE *array, const int imax, const int jmax, const int kmax, const int vc, const int vc_comm,
             const REAL_TYPE *recvm, const REAL_TYPE *recvp, const int nIDm, const int nIDp)
{
  int ix = imax;
  int jx = jmax;
  int kx = kmax;
  int gc = vc;
  int gc_comm = vc_comm;

  if( nIDm >= 0 )
  {
#pragma omp parallel for firstprivate(ix, jx, kx, gc, gc_comm) schedule(static) collapse(2)
    for( int k=0-gc_comm;k<kx+gc_comm;k++ ){
      for( int j=0-gc_comm;j<jx+gc_comm;j++ ){
        for( int i=0-gc_comm;i<0;i++ ){
          array[_IDX_S3D(i,j,k,ix,jx,kx,gc)] = recvm[_IDXFX(i,j,k,0-gc_comm,jx,kx,gc_comm)];
        }
      }
    }
  }

  if( nIDp >= 0 )
  {
#pragma omp parallel for firstprivate(ix, jx, kx, gc, gc_comm) schedule(static) collapse(2)
    for( int k=0-gc_comm;k<kx+gc_comm;k++ ){
      for( int j=0-gc_comm;j<jx+gc_comm;j++ ){
        for( int i=ix;i<ix+gc_comm;i++ ){
          array[_IDX_S3D(i,j,k,ix,jx,kx,gc)] = recvp[_IDXFX(i,j,k,ix,jx,kx,gc_comm)];
        }
      }
    }
  }
}

/*
 * @brief pack send data for Y direction
 * @param [in]  array   source array
 * @param [in]  imax    size x
 * @param [in]  jmax    size y
 * @param [in]  kmax    size z
 * @param [in]  vc      guide cell size
 * @param [in]  vc_comm number of guide cell layer to be sent
 * @param [out] sendm   send buffer of Y- direction
 * @param [out] sendp   send buffer of Y+ direction
 * @param [in]  nIDm    Rank number of Y- direction
 * @param [in]  nIDp    Rank number of Y+ direction
 */
void packY(const REAL_TYPE *array, const int imax, const int jmax, const int kmax, const int vc, const int vc_comm,
           REAL_TYPE *sendm, REAL_TYPE *sendp, const int nIDm, const int nIDp)
{
  int ix = imax;
  int jx = jmax;
  int kx = kmax;
  int gc = vc;
  int gc_comm = vc_comm;

  if( nIDm >= 0 )
  {
#pragma omp parallel for firstprivate(ix, jx, kx, gc, gc_comm) schedule(static) collapse(2)
    for( int k=0-gc_comm;k<kx+gc_comm;k++ ){
      for( int j=0;j<gc_comm;j++ ){
        for( int i=0-gc_comm;i<ix+gc_comm;i++ ){
          sendm[_IDXFY(i,j,k,ix,0,kx,gc_comm)] = array[_IDX_S3D(i,j,k,ix,jx,kx,gc)];
        }
      }
    }
  }

  if( nIDp >= 0 )
  {
#pragma omp parallel for firstprivate(ix, jx, kx, gc, gc_comm) schedule(static) collapse(2)
    for( int k=0-gc_comm;k<kx+gc_comm;k++ ){
      for( int j=jx-gc_comm;j<jx;j++ ){
        for( int i=0-gc_comm;i<ix+gc_comm;i++ ){
          sendp[_IDXFY(i,j,k,ix,jx-gc_comm,kx,gc_comm)] = array[_IDX_S3D(i,j,k,ix,jx,kx,gc)];
        }
      }
    }
  }
}


/*
 * @brief unpack send data for Y direction
 * @param [in,out]  array   dest array
 * @param [in]  imax    size x
 * @param [in]  jmax    size y
 * @param [in]  kmax    size z
 * @param [in]  vc      guide cell size
 * @param [in]  vc_comm number of guide cell layer to be sent
 * @param [in]  recvm   recv buffer of Y- direction
 * @param [in]  recvp   recv buffer of Y+ direction
 * @param [in]  nIDm    Rank number of Y- direction
 * @param [in]  nIDp    Rank number of Y+ direction
 */
void unpackY(REAL_TYPE *array, const int imax, const int jmax, const int kmax, const int vc, const int vc_comm,
             const REAL_TYPE *recvm, const REAL_TYPE *recvp, const int nIDm, const int nIDp)
{
  int ix = imax;
  int jx = jmax;
  int kx = kmax;
  int gc = vc;
  int gc_comm = vc_comm;

  if( nIDm >= 0 )
  {
#pragma omp parallel for firstprivate(ix, jx, kx, gc, gc_comm) schedule(static) collapse(2)
    for( int k=0-gc_comm;k<kx+gc_comm;k++ ){
      for( int j=0-gc_comm;j<0;j++ ){
        for( int i=0-gc_comm;i<ix+gc_comm;i++ ){
          array[_IDX_S3D(i,j,k,ix,jx,kx,gc)] = recvm[_IDXFY(i,j,k,ix,0-gc_comm,kx,gc_comm)];
        }
      }
    }
  }

  if( nIDp >= 0 )
  {
#pragma omp parallel for firstprivate(ix, jx, kx, gc, gc_comm) schedule(static) collapse(2)
    for( int k=0-gc_comm;k<kx+gc_comm;k++ ){
      for( int j=jx;j<jx+gc_comm;j++ ){
        for( int i=0-gc_comm;i<ix+gc_comm;i++ ){
          array[_IDX_S3D(i,j,k,ix,jx,kx,gc)] = recvp[_IDXFY(i,j,k,ix,jx,kx,gc_comm)];
        }
      }
    }
  }
}

/*
 * @brief pack send data for Z direction
 * @param [in]  array   source array
 * @param [in]  imax    size x
 * @param [in]  jmax    size y
 * @param [in]  kmax    size z
 * @param [in]  vc      guide cell size
 * @param [in]  vc_comm number of guide cell layer to be sent
 * @param [out] sendm   send buffer of Z- direction
 * @param [out] sendp   send buffer of Z+ direction
 * @param [in]  nIDm    Rank number of Z- direction
 * @param [in]  nIDp    Rank number of Z+ direction
 */
void packZ(const REAL_TYPE *array, const int imax, const int jmax, const int kmax, const int vc, const int vc_comm,
           REAL_TYPE *sendm, REAL_TYPE *sendp, const int nIDm, const int nIDp)
{
  int ix = imax;
  int jx = jmax;
  int kx = kmax;
  int gc = vc;
  int gc_comm = vc_comm;

  if( nIDm >= 0 )
  {
#pragma omp parallel for firstprivate(ix, jx, kx, gc, gc_comm) schedule(static) collapse(2)
    for( int k=0;k<gc_comm;k++ ){
      for( int j=0-gc_comm;j<jx+gc_comm;j++ ){
        for( int i=0-gc_comm;i<ix+gc_comm;i++ ){
          sendm[_IDXFZ(i,j,k,ix,jx,0,gc_comm)] = array[_IDX_S3D(i,j,k,ix,jx,kx,gc)];
        }
      }
    }
  }

  if( nIDp >= 0 )
  {
#pragma omp parallel for firstprivate(ix, jx, kx, gc, gc_comm) schedule(static) collapse(2)
    for( int k=kx-gc_comm;k<kx;k++ ){
      for( int j=0-gc_comm;j<jx+gc_comm;j++ ){
        for( int i=0-gc_comm;i<ix+gc_comm;i++ ){
          sendp[_IDXFZ(i,j,k,ix,jx,kx-gc_comm,gc_comm)] = array[_IDX_S3D(i,j,k,ix,jx,kx,gc)];
        }
      }
    }
  }
}


/*
 * @brief unpack send data for Z direction
 * @param [in,out]  array   dest array
 * @param [in]  imax    size x
 * @param [in]  jmax    size y
 * @param [in]  kmax    size z
 * @param [in]  vc      guide cell size
 * @param [in]  vc_comm number of guide cell layer to be sent
 * @param [in]  recvm   recv buffer of Z- direction
 * @param [in]  recvp   recv buffer of Z+ direction
 * @param [in]  nIDm    Rank number of Z- direction
 * @param [in]  nIDp    Rank number of Z+ direction
 */
void unpackZ(REAL_TYPE *array, const int imax, const int jmax, const int kmax, const int vc, const int vc_comm,
             const REAL_TYPE *recvm, const REAL_TYPE *recvp, const int nIDm, const int nIDp)
{
  int ix = imax;
  int jx = jmax;
  int kx = kmax;
  int gc = vc;
  int gc_comm = vc_comm;

  if( nIDm >= 0 )
  {
#pragma omp parallel for firstprivate(ix, jx, kx, gc, gc_comm) schedule(static) collapse(2)
    for( int k=0-gc_comm;k<0;k++ ){
      for( int j=0-gc_comm;j<jx+gc_comm;j++ ){
        for( int i=0-gc_comm;i<ix+gc_comm;i++ ){
          array[_IDX_S3D(i,j,k,ix,jx,kx,gc)] = recvm[_IDXFZ(i,j,k,ix,jx,0-gc_comm,gc_comm)];
        }
      }
    }
  }

  if( nIDp >= 0 )
  {
#pragma omp parallel for firstprivate(ix, jx, kx, gc, gc_comm) schedule(static) collapse(2)
    for( int k=kx;k<kx+gc_comm;k++ ){
      for( int j=0-gc_comm;j<jx+gc_comm;j++ ){
        for( int i=0-gc_comm;i<ix+gc_comm;i++ ){
          array[_IDX_S3D(i,j,k,ix,jx,kx,gc)] = recvp[_IDXFZ(i,j,k,ix,jx,kx,gc_comm)];
        }
      }
    }
  }
}
