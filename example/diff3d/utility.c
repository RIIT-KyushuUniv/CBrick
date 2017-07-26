s
//
//  utility.c
//  diff3d
//
//  Created by keno on 2016/06/26.
//  Copyright © 2016年 keno. All rights reserved.
//

#include "utility.h"


// @fn usage
// @brief print usage
void usage()
{
  printf("Usage:\n");
  //printf("\t$ ./diff3d config_file\n\n");
  printf("\t$ mpirun -np n diff3d config_file\n\n");
  printf("\tContents of config_file\n");
  printf("\tnx;       Number of cells for x-dir.\n");
  printf("\tny;       Number of cells for y-dir.\n");
  printf("\tnz;       Number of cells for z-dir.\n");
  printf("\tXdiv;     Number of division for x-dir.\n");
  printf("\tYdiv;     Number of division for y-dir.\n");
  printf("\tZdiv;     Number of division for z-dir.\n");
  //printf("\tdh;       Mesh width (isotropic, non-dimensional)\n");
  printf("\tdt;       Time increment of time marching (non-dimensional)\n");
  printf("\talpha;    Coefficient of diffusion (non-dimensional)\n");
  printf("\tlaststep; Time step to calculate\n");
  printf("\tfileout;  Interval for writing a file\n\n");
}
