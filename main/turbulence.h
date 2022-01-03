/*
 * File: turbulence.h
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Header file for turbulence.c.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#ifndef _turbulence_h
#define _turbulence_h

// Von Karman's Constant
#define KAPPA_VK 0.42

// Background length scale (as q^2 l)
#define LBACKGROUND (1e-10)
// Background velocity scale (as q^2)
#define QBACKGROUND (1e-4*1e-4)

#define TURB_KKL 1
#define TURB_GEN 2
#define TURB_KEPS 3
#define TURB_KOMEGA 4
#define TURB_KKL_GLS 5
#define TURB_PARABOLIC 10

void my25(gridT *grid, physT *phys, propT *prop, REAL **wnew, REAL **q, REAL **l,
          REAL **Cn_q, REAL **Cn_l, REAL **nuT, REAL **kappaT, MPI_Comm comm, int myproc);

void gls(gridT *grid, physT *phys, propT *prop, REAL **wnew, REAL **q, REAL **l,
         REAL **Cn_q, REAL **Cn_l, REAL **nuT, REAL **kappaT, MPI_Comm comm, int myproc);

void parabolic_viscosity(gridT *grid, physT *phys, propT *prop, REAL **wnew, REAL **q, REAL **l,
                         REAL **Cn_q, REAL **Cn_l, REAL **nuT, REAL **kappaT, MPI_Comm comm, int myproc);

void cell_centered_bed_stress_interp(physT *phys, gridT *grid, REAL *taux, REAL *tauy);
void edge_centered_bed_stress(physT *phys, gridT *grid, REAL *tau);

#endif
