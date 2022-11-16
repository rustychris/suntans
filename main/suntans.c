/*
 * File: suntans.c
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * SUNTANS - Stanford Unstructured Nonhydrostatic Terrain-following Adaptive
 * Navier-Stokes Simulator
 *
 * http://suntans.stanford.edu
 * 
 * Written by Oliver B. Fringer
 * Dept. of Civil and Environmental Engineering
 * Stanford University
 * 
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 * 
 */
#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>

#include "suntans.h"
#include "memory.h"
#include "mympi.h"
#include "grid.h"
#include "gridio.h"
#include "phys.h"
#include "sediments.h"
#include "physio.h"
#include "report.h"


/* Handler to print a stack trace on segfaults
 */
void handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}

// define global variables for filenames
char DATADIR[BUFFERLENGTH],
  DATAFILE[2*BUFFERLENGTH], // *2 for sprintf target
  PSLGFILE[BUFFERLENGTH], 
  POINTSFILE[BUFFERLENGTH], 
  EDGEFILE[BUFFERLENGTH], 
  CELLSFILE[BUFFERLENGTH], 
  NODEFILE[BUFFERLENGTH], 
  INPUTDEPTHFILE[BUFFERLENGTH],
  CELLCENTEREDFILE[BUFFERLENGTH], 
  EDGECENTEREDFILE[BUFFERLENGTH], 
  VERTSPACEFILE[BUFFERLENGTH], 
  TOPOLOGYFILE[BUFFERLENGTH];
// define global variables
int TRIANGULATE, GRID, SOLVE, VERBOSE, WARNING, ASCII, RESTART, NUMPROCS, STEPSPERFILE;

int main(int argc, char *argv[])
{
  int myproc, numprocs, j;
  MPI_Comm comm;
  gridT *grid;
  physT *phys;
  propT *prop;

  signal(SIGSEGV, handler); // install SIGSEGV handler
  
  StartMpi(&argc,&argv,&comm,&myproc,&numprocs);

  ParseFlags(argc,argv,myproc);
  
  if(GRID)
    GetGrid(&grid,myproc,numprocs,comm);
  else
    ReadGrid(&grid,myproc,numprocs,comm);
  

  if(SOLVE) {
    //read parameters in suntans.dat into the solver
    ReadProperties(&prop,grid,myproc);
    // give space and initialize dzf(edge) dzz(center) dzzold(center)
    InitializeVerticalGrid(&grid,myproc);
    AllocatePhysicalVariables(grid,&phys,prop);
    AllocateTransferArrays(&grid,myproc,numprocs,comm);
    
    InitializeEdgeDepths(grid,myproc,comm,LOCAL_GRID);
    OpenFiles(prop,myproc);
    if(RESTART)
      ReadPhysicalVariables(grid,phys,prop,myproc,comm);
    else
      InitializePhysicalVariables(grid,phys,prop,myproc,comm);
    
    Solve(grid,phys,prop,myproc,numprocs,comm);
    //    FreePhysicalVariables(grid,phys,prop);
    //    FreeTransferArrays(grid,myproc,numprocs,comm);
  }

  EndMpi(&comm);
}





