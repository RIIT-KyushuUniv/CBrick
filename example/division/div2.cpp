//
//  div2.cpp
//
//  Created by keno on 2017/07/01.
//  Copyright © 2017 keno. All rights reserved.
//

// Compile
// $ gcc div2.cpp -o div2

// Execution
// $ mpirun -np  4 div2 15  1  1 2 2 1
// $ mpirun -np 12 div2 31 12 20 2 2 3
// $ mpirun -np 16 div2 30 24 24 2 4 2

// プロセス数、領域サイズと領域分割数を与えて、最適な分割パターンを得る（セルベース FVM）

#include <CB_SubDomain.h>
#include <stdlib.h>


int main(int argc, char * argv[]) {

  // Check command line arguments
  if ( argc != 7) {
    printf("Usage:\n");
    printf("\t$ mpirun -np N div2 nx ny nz divx divy divz\n");
    return -1;
  }

  int m_sz[3], np, myrank, m_dv[3];
  int proc_grp = 0;
  int gc = 1;      // ガイドセル幅=1

  // conversion from ASCII char to digit
  m_sz[0] = atoi(argv[1]);
  m_sz[1] = atoi(argv[2]);
  m_sz[2] = atoi(argv[3]);
  m_dv[0] = atoi(argv[4]);
  m_dv[1] = atoi(argv[5]);
  m_dv[2] = atoi(argv[6]);

  // Initialize MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  if (myrank==0) printf("%d x %d x %d / %d\n", m_sz[0], m_sz[1], m_sz[2], np);

  if (m_dv[0]*m_dv[1]*m_dv[2] != np) {
    printf("\tThe number of proceees does not agree with the div[] size.\n");
    exit(-1);
  }

  // priorityはデフォルト
  SubDomain D(m_sz, gc, m_dv, np, myrank, proc_grp, MPI_COMM_WORLD, "cell");

  if ( !D.findParameter() ) MPI_Abort(MPI_COMM_WORLD, -1);

  if ( !D.createRankTable() ) MPI_Abort(MPI_COMM_WORLD, -1);

  MPI_Finalize();

  return 0;
}
