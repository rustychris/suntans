/*
 * File: profiles.h
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Header file containing global variable definitions for profiles.c.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#ifndef _profiles_h
#define _profiles_h

#include "suntans.h"

#define MAXOUTPUTVARIABLES 8
#define ALLPROFILEVARIABLES "husbTqnkC"
#define DEFAULTPROFILEVARIABLES "husb"

extern FILE *FreeSurfaceProfFID, *HorizontalVelocityProfFID, *VerticalVelocityProfFID,
  *SalinityProfFID, *BGSalinityProfFID, *TemperatureProfFID, *PressureProfFID, 
  *EddyViscosityProfFID, *ScalarDiffusivityProfFID, *ProfileDataFID,
  **SediProfFID;

extern int existProfs, numInterpPoints, ntoutProfs, NkmaxProfs, numTotalDataPoints, numLocalDataPoints, *dataIndices, *interpIndices;
extern int *total2d, all2d, *total3d, *allIndices;
extern REAL *dataXY, *merge_tmp, *merge_tmp2;
extern char ProfileVariables[BUFFERLENGTH];

void InterpData(gridT *grid, physT *phys, propT *prop, MPI_Comm, int numprocs, int myproc);

#endif
