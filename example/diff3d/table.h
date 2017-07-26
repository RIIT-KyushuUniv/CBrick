//
//  table.h
//  diff3dp
//
//  Created by keno on 2016/07/01.
//  Copyright © 2016年 keno. All rights reserved.
//

#ifndef table_h
#define table_h

#include "diff3dp.h"

void getHead(int* sz,
             const int* dv,
             const int* dsz,
             int head[3],
             const int myrank,
             int* comm);

void setTable(const int div[3],
              const int np,
              const int myrank,
              const int dsz[3],
              int my_sz[3],
              int my_hd[3],
              int my_tbl[6]);

void display(const int np,
             const int myrank,
             const int* my_size,
             const int* my_head,
             const int* my_tbl,
             const REAL_TYPE* org);

int mktable(const int np,
            const int myrank,
            const int mode,
            Phys_Param* m_phys,
            Cntl_Param* m_cntl,
            int comm_tbl[6],
            REAL_TYPE g_org[3]);

#endif /* table_h */
