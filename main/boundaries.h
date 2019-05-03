/*
 * File: boundaries.h
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Header file for boundaries.c.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#ifndef _boundaries_h
#define _boundaries_h

#include "suntans.h"
#include "phys.h"
#include "grid.h"
#include "met.h"
#include "util.h"

#define NT 3

// Enumerated type for open/specified bc specification
enum {
  specified, open
};


// struct to hold netcdf boundary condition data for a single
// scalar.  Pulling this out into its own struct to allow for
// more generic handling across T,S and sediment.
typedef struct _scalar_boundT {
  REAL ***boundary_scal_t;
  REAL **boundary_scal_f; // forward timestep
  REAL **boundary_scal_b; // backward timestep
  REAL **boundary_scal;   // current value
  REAL *point_scal; // point sources
  REAL ***scal_t;  // type 3, cell centered boundaries
  REAL **scal_f;
  REAL **scal_b;
  REAL **scal;
  // string used to form input/output variable names for netcdf
  // T,S, sed1, sed2, etc.
  char varname[256]; 
} scalar_boundT;

// temperature, salinity, sediments, etc.
#define MAXSCALARS 10

// Structure array to store netcdf boundary condition data
typedef struct _boundT{
  
  // Dimension sizes
  size_t Ntype2;
  size_t Ntype3;
  size_t Npoint_source;
  size_t Nseg;
  size_t Nt;
  size_t Nk;
  
  // boolean operators
  int hasType2;
  int hasType3;
  int hasSeg;

  // Grid cell indices
  int *edgep;
  int *localedgep;
  int *cellp;
  int *segedgep;
  int *segp;
  int *point_cell;
  int *point_layer;

  // Indices that point the grid to the cell in the file
  int *ind2;
  int *ind3;
  int *ind3edge;
  int *ind_point;

  // Boundary coordinates
  REAL *xe;
  REAL *ye;
  REAL *xv;
  REAL *yv;
  REAL *z;
  REAL *time;	
  REAL *segarea;
  REAL *localsegarea;

  // Time record locators
  int t0;
  int t1;
  int t2; 

  // Data arrays at forward (_f) and backward (_b) timestep
  // Type-2 (edge centred) boundaries

  REAL ***boundary_u_t;
  REAL ***boundary_v_t;
  REAL ***boundary_w_t;
  REAL ***boundary_T_t;
  REAL ***boundary_S_t;
  REAL **boundary_Q_t;

  REAL **boundary_u_f;
  REAL **boundary_v_f;
  REAL **boundary_w_f;
  REAL **boundary_T_f;
  REAL **boundary_S_f;
  REAL *boundary_Q_f;

  REAL **boundary_u_b;
  REAL **boundary_v_b;
  REAL **boundary_w_b;
  REAL **boundary_T_b;
  REAL **boundary_S_b;
  REAL *boundary_Q_b;

  REAL **boundary_u;
  REAL **boundary_v;
  REAL **boundary_w;
  REAL **boundary_T;
  REAL **boundary_S;
  REAL *boundary_Q;

  // Point sources
  REAL *point_Q;
  REAL *point_T;
  REAL *point_S;
  
  // Type-3 (cell centred) boundaries
  REAL ***uc_t;
  REAL ***vc_t;
  REAL ***wc_t;
  REAL ***T_t;
  REAL ***S_t;
  REAL ** h_t;

  REAL **uc_f;
  REAL **vc_f;
  REAL **wc_f;
  REAL **T_f;
  REAL **S_f;
  REAL * h_f;

  REAL **uc_b;
  REAL **vc_b;
  REAL **wc_b;
  REAL **T_b;
  REAL **S_b;
  REAL *h_b;

  REAL **uc;
  REAL **vc;
  REAL **wc;
  REAL **T;
  REAL **S;
  REAL *h;

  int num_scalars;
  scalar_boundT scalars[MAXSCALARS];
  // references to specific elements of scalars
  scalar_boundT *T_scal,*S_scal;
} boundT;

// Declare the boundary structure global
boundT *bound;

FILE *windFID;

void OpenBoundaryFluxes(REAL **q, REAL **ub, REAL **ubn, gridT *grid, physT *phys, propT *prop,
                        int myproc);
void BoundaryVelocities(gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm);
void BoundaryScalars(gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm);
void BoundarySediment(gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm);
void PointSources(gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm);
void PointSourcesContinuity(REAL **w, gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm);
void PointSourceTemp(REAL **A, REAL **B, gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm);
void PointSourceSalt(REAL **A, REAL **B, gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm);
void WindStress(gridT *grid, physT *phys, propT *prop, metT *met, int myproc);
void InitBoundaryData(propT *prop, gridT *grid, int myproc, MPI_Comm comm);
scalar_boundT *AllocateBoundaryScalar(const char *name, boundT *bound, int myproc);
void AllocateBoundaryData(propT *prop, gridT *grid, boundT **bound, int myproc, MPI_Comm comm);
void UpdateBdyNC(propT *prop, gridT *grid, int myproc, MPI_Comm comm);

// Copies netcdf-derived data to the arrays ready for phys, sediments
void BoundaryScalarGeneric(REAL **boundary_scal, // [jptr-edgedist[2], k]
                           REAL **scal, // [cellp,k]
                           scalar_boundT *scalar,
                           gridT *grid, physT *phys, propT *prop,int myproc, MPI_Comm comm);

#endif
