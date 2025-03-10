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
#define INFTY 1e20
#define CONSERVED 1e-5
#define EMPTY 999999
#define SMALL 1e-15
#define CHECKCONSISTENCY 0
// Cell with depth less than this is considered inactive
// and dry.
// possible that this should be significantly larger than
// the CG epsilon
#define DRYCELLHEIGHT 1e-3

// merge top layer with next below when it is thinner
// dzmin_surface moved to runtime configuration

// Evaporation (negative rain) will not decrease a cell's depth
// below this number.
#define DZMIN_EVAP 0.2

// Cell with depth less than this will be propped up to
// to this deep.  should be less than DRYCELLHEIGHT
#define CLAMPHEIGHT (0.2*DRYCELLHEIGHT)

// Edge with depth less than this get a very high drag
// coefficient (100)
#define BUFFERHEIGHT 1e-2
// alternatively, employ a gross weir when d eta/dx exceeds this
// #define WEIR_GRADIENT 1e-2

#define DEFAULT_NFACES 3

// enable low-level output for a particular model element
// #define DBG_PROC nnn
// #define DBG_EDGE nnn
// #define DBG_CELL nnn

// Error/Exit codes
#define EXIT_WRITING 1

#define DEFAULTDATAFILE "suntans.dat"

// 0: disable horizontal advection of turbulence
// 1: enable
#define HOR_ADV_TURBULENCE 1

// for debugging. global scaling of nut and kappaT.
// #define TURB_SCALE 0.2

// define global variables for filenames
extern char DATADIR[BUFFERLENGTH],
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
extern int TRIANGULATE, GRID, SOLVE, VERBOSE, WARNING, ASCII, RESTART, NUMPROCS, STEPSPERFILE;

#endif
