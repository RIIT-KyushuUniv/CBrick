/*
###################################################################################
#
# RIAM-COMPACT
#
# Copyright (C) 2015-2017 Research Institute for Applied Mechanics(RIAM)
#                       / Research Institute for Information Technology(RIIT), Kyushu University.
# All rights reserved.
#
# Copyright (C) 2015-2016 Advanced Institute for Computational Science, RIKEN.
# All rights reserved.
#
###################################################################################
*/

#ifndef _RIAMC_F_FUNC_H_
#define _RIAMC_F_FUNC_H_


//RIAM_COMPACT;
extern "C" {

extern void grid_gen_(
                  int* sz,
                  REAL_TYPE* Z,
                  REAL_TYPE* org,
                  REAL_TYPE* pch,
                  int* head,
                  int* comm_tbl,
                  double* flop);

extern void var_init_(
                  int* sz,
                  int* ICON,
                  REAL_TYPE* V,
                  REAL_TYPE* P,
                  REAL_TYPE* VV);

extern void eddy_viscosity_(
                  int* sz,
                  REAL_TYPE* RE,
                  REAL_TYPE* CS2,
                  REAL_TYPE* Z,
                  REAL_TYPE* pitch,
                  REAL_TYPE* V,
                  REAL_TYPE* SGS,
                  int* comm_tbl,
                  double* fcount);

extern void bc_eddy_(
                  int* sz,
                  REAL_TYPE* SGS,
                  int* comm_tbl);

extern void copy_vec_(
                  int* sz,
                  REAL_TYPE* V,
                  REAL_TYPE* VD);

extern void intermediate_v_(
                  int* sz,
                  REAL_TYPE* Alpha,
                  REAL_TYPE* DeltaT,
                  REAL_TYPE* RE,
                  REAL_TYPE* Z,
                  REAL_TYPE* pitch,
                  REAL_TYPE* V,
                  REAL_TYPE* VD,
                  REAL_TYPE* VV,
                  REAL_TYPE* SGS,
                  int* comm_tbl,
                  double* fcount);

extern void contra_v_cc_(
                  int* sz,
                  REAL_TYPE* V,
                  REAL_TYPE* Z,
                  REAL_TYPE* pitch,
                  REAL_TYPE* CV,
                  double* fcount);

extern void contra_v_cf_(
                  int* sz,
                  REAL_TYPE* VV,
                  REAL_TYPE* CV,
                  double* fcount);

extern void rhs_poisson_(
                  int* sz,
                  REAL_TYPE* DT,
                  REAL_TYPE* VV,
                  REAL_TYPE* RHS,
                  REAL_TYPE* Z,
                  REAL_TYPE* pitch,
                  double* fcount);

extern void average_prs_(
                  int* sz,
                  REAL_TYPE* P,
                  int* comm_tbl,
                  REAL_TYPE* PP,
                  double* fcount);

extern void new_pressure_(
                  int* sz,
                  REAL_TYPE* P,
                  REAL_TYPE* PP,
                  double* fcount);

extern void bc_new_prs_(
                  int* sz,
                  REAL_TYPE* RE,
                  REAL_TYPE* P,
                  REAL_TYPE* Z,
                  REAL_TYPE* pitch,
                  REAL_TYPE* V,
                  int* comm_tbl,
                  double* fcount);

extern void new_vec_cc_(
                  int* sz,
                  REAL_TYPE* DeltaT,
                  REAL_TYPE* P,
                  REAL_TYPE* V,
                  REAL_TYPE* Z,
                  REAL_TYPE* pitch,
                  int* comm_tbl,
                  double* fcount);

extern void new_vec_cf_(
                  int* sz,
                  REAL_TYPE* DeltaT,
                  REAL_TYPE* P,
                  REAL_TYPE* VV,
                  REAL_TYPE* Z,
                  REAL_TYPE* pitch,
                  double* fcount);

extern void bc_new_vel_(
                  int* sz,
                  int* index,
                  REAL_TYPE* DeltaT,
                  REAL_TYPE* pitch,
                  REAL_TYPE* V,
                  REAL_TYPE* VD,
                  int* comm_tbl,
                  double* fcount);

extern void display_(
                  int* sz,
                  int* ISTEP, int* T_NSTEP, int* ILAP,
                  int* ICOUNT, int* Interval_Disp,
                  REAL_TYPE *V,
                  int* Interval_Log, int* Interval_Vis,
                  REAL_TYPE* TIME, REAL_TYPE* RE, REAL_TYPE* RMSP, int* myRank);

extern void display2_(
                  int* sz,
                  REAL_TYPE* TIME, REAL_TYPE* RE,
                  REAL_TYPE* V,
                  REAL_TYPE* P,
                  REAL_TYPE* VV,
                  int* myRank);

extern void copy_v2vex_(
                  int* sz,
                  REAL_TYPE* V,
                  REAL_TYPE* VEX);

extern void find_max_v_(
                  int* sz,
                  REAL_TYPE* v_max,
                  REAL_TYPE* V,
                  double* flop);


// Linear solver

extern void psor2sma_maf_(
                  int* sz,
                  REAL_TYPE* OMG,
                  REAL_TYPE* RSMP,
                  REAL_TYPE* residual,
                  int* offset,
                  int* color,
                  REAL_TYPE* P,
                  REAL_TYPE* Z,
                  REAL_TYPE* pitch,
                  REAL_TYPE* RHS,
                  int* comm_tbl,
                  double* fcount);

extern void psor_maf_(
                  int* sz,
                  REAL_TYPE* OMG,
                  REAL_TYPE* RSMP,
                  REAL_TYPE* residual,
                  REAL_TYPE* P,
                  REAL_TYPE* Z,
                  REAL_TYPE* pitch,
                  REAL_TYPE* RHS,
                  int* comm_tbl,
                  double* fcount);
}



#endif // _RIAMC_F_FUNC_H_
