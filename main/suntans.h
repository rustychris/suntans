/*
 * File: suntans.h
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Main header file
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#ifndef _suntans_h
#define _suntans_h

#include "math.h"
#include "time.h"

// number of faces (for triangle)
#define REAL double
#define BUFFERLENGTH 256
#define NUMEDGECOLUMNS 3
#define BREAK printf("%d\n",*((int *)0));
#define PI 3.141592654
#define RHO0 1000.0
#define GRAV 9.81
#define INFTY 1e20
#define CONSERVED 1e-5
#define EMPTY 999999
#define SMALL 1e-15
#define CHECKCONSISTENCY 0
// Cell with depth less than this is considered inactive
// and dry.
// possible that this should be significantly larger than
// the CG epsilon
#define DRYCELLHEIGHT 5e-4

// merge top layer with next below when it is thinner 
// than this
#define DZMIN_SURFACE (2*DRYCELLHEIGHT)

// Cell with depth less than this will be propped up to
// to this deep.  should be less than DRYCELLHEIGHT
#define CLAMPHEIGHT (0.2*DRYCELLHEIGHT)

// Edge (or is it cell?) with depth less than this get a very high drag
// coefficient (100)
#define BUFFERHEIGHT 1e-2
#define DEFAULT_NFACES 3

// enable low-level output for a particular model element
// #define DBG_PROC 0
// #define DBG_EDGE 13333
// #define DBG_CELL 4130


// Error/Exit codes
#define EXIT_WRITING 1

#define DEFAULTDATAFILE "suntans.dat"

// define global variables for filenames
char DATADIR[BUFFERLENGTH],
  DATAFILE[BUFFERLENGTH],
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

#endif
