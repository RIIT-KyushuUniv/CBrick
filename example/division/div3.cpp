//
//  div3.cpp
//
//  Created by keno on 2017/07/01.
//  Copyright © 2017 keno. All rights reserved.
//

// Compile
// $ gcc div3.cpp -o div3

// Execution
// $ mpirun -np  4 div3 15  1  1 0
// $ mpirun -np 12 div3 31 12 20 0
// $ mpirun -np 16 div3 30 24 24 0

// 領域サイズと領域分割数を与えて、最適な分割パターンを得る（ノードベース FDM）

#include <CB_SubDomain.h>
#include <stdlib.h>


int main(int argc, char * argv[]) {

  // Check command line arguments
  if ( argc != 5) {
    printf("Usage:\n");
    printf("\t$ mpirun -np N div1 nx ny nz div_mode\n");
    return -1;
  }

  int m_sz[3], np, myRank;
  int proc_grp = 0;
  int gc = 1;      // ガイドセル幅=1
  int div_mode=0;

  // conversion from ASCII char to digit
  m_sz[0] = atoi(argv[1]);
  m_sz[1] = atoi(argv[2]);
  m_sz[2] = atoi(argv[3]);
  div_mode= atoi(argv[4]);

  // Initialize MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);

  Hostonly_ printf("%d x %d x %d / %d\n", m_sz[0], m_sz[1], m_sz[2], np);

  // priorityはデフォルト
  SubDomain D(m_sz, gc, np, myRank, proc_grp, MPI_COMM_WORLD, "node", "Cindex");

  if ( !D.findOptimalDivision(div_mode) )
  {
    MPI_Barrier(MPI_COMM_WORLD);
    exit(-1);
  }

  if ( !D.createRankTable() )
  {
    MPI_Barrier(MPI_COMM_WORLD);
    exit(-1);
  }

  MPI_Finalize();
  Hostonly_ printf("Successfully terminated.\n\n");

  return 0;
}
