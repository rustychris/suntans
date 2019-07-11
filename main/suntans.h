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
#define DRYCELLHEIGHT 1e-3

// merge top layer with next below when it is thinner 
// than this
// Had been (20*DRYCELLHEIGHT)
// with SFB grid, that was getting down to 2.4s time step limitation.
// Changing this from 2cm to 5cm gives a similar increase in the min_time_step.
// 2019-07-01: for breach flows, testing 0.10m to see if the timesteps are
// are limited by vertical or horizontal fluxes
#define DZMIN_SURFACE (0.25)

// Evaporation (negative rain) will not decrease a cell's depth
// below this number.
#define DZMIN_EVAP 0.2

// Cell with depth less than this will be propped up to
// to this deep.  should be less than DRYCELLHEIGHT
#define CLAMPHEIGHT (0.2*DRYCELLHEIGHT)

// Edge (or is it cell?) with depth less than this get a very high drag
// coefficient (100)
// #define BUFFERHEIGHT 1e-2
// alternatively, employ a gross weir when d eta/dx exceeds this
#define WEIR_GRADIENT 1e-2

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
