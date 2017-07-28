//
//  div1.cpp
//
//  Created by keno on 2017/07/01.
//  Copyright © 2017 keno. All rights reserved.
//

// Compile
// $ gcc div1.cpp -o div1

// Execution
// $ mpirun -np  4 div1 15  1  1
// $ mpirun -np 12 div1 31 12 20
// $ mpirun -np 16 div1 30 24 24

// プロセス数と領域サイズを与えて、最適な分割パターンを得る（セルベース FVM）

#include <CB_SubDomain.h>
#include <stdlib.h>


int main(int argc, char * argv[]) {

  // Check command line arguments
  if ( argc != 4) {
    printf("Usage:\n");
    printf("\t$ mpirun -np N div1 nx ny nz\n");
    return -1;
  }

  int m_sz[3], np, myrank;
  int proc_grp = 0;
  int gc = 1;      // ガイドセル幅=1

  // conversion from ASCII char to digit
  m_sz[0] = atoi(argv[1]);
  m_sz[1] = atoi(argv[2]);
  m_sz[2] = atoi(argv[3]);

  // Initialize MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &np);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  if (myrank==0) printf("%d x %d x %d / %d\n", m_sz[0], m_sz[1], m_sz[2], np);

  // priorityはデフォルト
  SubDomain D(m_sz, gc, np, myrank, proc_grp, MPI_COMM_WORLD, "cell", "Cindex", 0);

  if ( !D.findOptimalDivision() ) MPI_Abort(MPI_COMM_WORLD, -1);

  if ( !D.createRankTable() ) MPI_Abort(MPI_COMM_WORLD, -1);

  MPI_Finalize();

  return 0;
}
