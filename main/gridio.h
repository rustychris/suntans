/*
 * File: physio.h
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Header file for gridio.c.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#ifndef _gridio_h
#define _gridio_h

#include "grid.h"
#include "mympi.h"

void ReadGrid(gridT **grid, int myproc, int numprocs, MPI_Comm comm);
void OutputGridData(gridT *maingrid, gridT *grid, int myproc, int numprocs);
void ReadGridFileNames(int myproc);
void ReadDepth(gridT *grid, int myproc);
void InitializeEdgeDepths(gridT *grid, int myproc, MPI_Comm comm, grid_type_t grid_type);
void WriteEdgeDepths(gridT *grid,int myproc, grid_type_t grid_type);
void ReadMainGrid(gridT *grid, int myproc);

#endif
