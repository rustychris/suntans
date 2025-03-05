/*
 * File: initialization.c
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Contains the functions that are used
 * to initialize the depth, free-surface, and salinity.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#include "math.h"
#include "fileio.h"
#include "suntans.h"
#include "initialization.h"
#include "mynetcdf.h"

/* this mimics grid->dz_lump_bed_frac - it is not strictly necessary that the value
   defined here be the same as what is used elsewhere - it is used just for determining
   computational cost of each water column.
*/
#define DZ_LUMP_BED_FRAC 0.2


#define sech 1/cosh
/*
 * Function: GetDZ
 * Usage: GetDZ(dz,depth,Nkmax,myproc);
 * ------------------------------------
 * Returns the vertical grid spacing in the array dz.
 *
 */
int GetDZ(REAL *dz, REAL depth, REAL localdepth, int Nkmax, int myproc) {
  int k, status;
  REAL z=0, dz0, r = GetValue(DATAFILE,"rstretch",&status);
  /* Added for reading vertspace.at */
  char fname[BUFFERLENGTH];
  FILE *vertspace_fp=NULL;
  static REAL *vertspace_cached=NULL;

  /* RH: special case, if r==0, then use a pre-existing vertspace.dat */
  REAL dz_k;
  static int warned=0;

  if ( r == 0.0) {
    if(!vertspace_cached ) {
      vertspace_cached = (REAL*)SunMalloc(sizeof(REAL)*Nkmax,"GetDZ");

      sprintf(fname,"%s/vertspace.dat.in",DATADIR);

      if ( !warned ) {
        printf("RH: using pre-existing vertspace.dat.in in GetDZ\n");
        printf("     reading from %s\n",fname);
        warned = 1;
      }

      vertspace_fp = MPI_FOpen(fname,"r","GetDZ",myproc);

      for(k=0;k<Nkmax;k++) {
        vertspace_cached[k] = getfield(vertspace_fp,fname);
      }
      fclose(vertspace_fp);
    }
    if ( dz != NULL ) {
      for(k=0;k<Nkmax;k++) {
        dz[k] = vertspace_cached[k];
      }
      return -999;
    } else {
      z=0;
      for (k=0;k<Nkmax;k++) {
        dz_k = vertspace_cached[k];

        z += dz_k; // z is now the depth of the bottom of layer k

        if ( z >= localdepth ) {
          // Check to see if this cell should be lumped
          if ( localdepth - (z-dz_k) < DZ_LUMP_BED_FRAC * dz_k) {
            return k;
          }
          return k+1; // return the number of layers
        }
      }
      printf("\nWARNING: Assigning Nkmax, localdepth = %f\n\n", localdepth);
      printf("Depth = %f\n", depth);
      printf("z = %f\n", z);
      return Nkmax;
    }
  } else {
    /* Existing handling for r != 0 left unchanged */
    if(dz!=NULL) {
      if(r==1) 
        for(k=0;k<Nkmax;k++)
  	  dz[k]=depth/Nkmax;
      else if(r>1 && r<=1.3) {    
        dz[0] = depth*(r-1)/(pow(r,Nkmax)-1);
        if(VERBOSE>2) printf("Minimum vertical grid spacing is %.2f\n",dz[0]);
        for(k=1;k<Nkmax;k++) 
	  dz[k]=r*dz[k-1];
      } else if(r>-1.3 && r<-1) {    
        r=fabs(r);
        dz[Nkmax-1] = depth*(r-1)/(pow(r,Nkmax)-1);
        if(VERBOSE>2) printf("Minimum vertical grid spacing is %.2f\n",dz[Nkmax-1]);
        for(k=Nkmax-2;k>=0;k--) 
	  dz[k]=r*dz[k+1];
      } else {
        printf("Error in GetDZ when trying to create vertical grid:\n");
        printf("Absolute value of stretching parameter rstretch must  be in the range (1,1.1).\n");
        exit(1);
      }
    } else {
      r=fabs(r);
      if(r!=1)
        dz0 = depth*(r-1)/(pow(r,Nkmax)-1);
      else
        dz0 = depth/Nkmax;
      z = dz0;
      for(k=1;k<Nkmax;k++) {
        dz0*=r;
        z+=dz0;
        if(z>=localdepth) {
	  return k;
        }
      }
    }
  }
  return -1;
}
  
/*
 * Function: ReturnDepth
 * Usage: grid->dv[n]=ReturnDepth(grid->xv[n],grid->yv[n]);
 * --------------------------------------------------------
 * Helper function to create a bottom bathymetry.  Used in
 * grid.c in the GetDepth function when IntDepth is 0.
 *
 */
REAL ReturnDepth(REAL x, REAL y) {
  REAL length, xmid, shelfdepth, depth;

//  length = 10000;
//  xmid = 65000;
//  shelfdepth = 500;
//  depth = 3000;
//  if(x<=xmid-length/2)
//    return depth;
//  else if(x>xmid-length/2 && x<=xmid+length/2 && length>0)
//    return depth-(depth-shelfdepth)*(x-xmid+length/2)/length;
//  else
//    return shelfdepth;
    return 1;
}

 /*
  * Function: ReturnFreeSurface
  * Usage: grid->h[n]=ReturnFreeSurface(grid->xv[n],grid->yv[n]);
  * -------------------------------------------------------------
  * Helper function to create an initial free-surface. Used
  * in phys.c in the InitializePhysicalVariables function.
  *
  */
REAL ReturnFreeSurface(REAL x, REAL y, REAL d) {
  return 0;
}

/*
 * Function: ReturnSalinity
 * Usage: grid->s[n]=ReturnSalinity(grid->xv[n],grid->yv[n],z);
 * ------------------------------------------------------------
 * Helper function to create an initial salinity field.  Used
 * in phys.c in the InitializePhysicalVariables function.
 *
 */
REAL ReturnSalinity(REAL x, REAL y, REAL z) {
  REAL thermocline_depth=20;

  if(z>-thermocline_depth)
    return 3.4286*pow(fabs(thermocline_depth),0.0187)-3.6;
  return 3.4286*pow(fabs(z),0.0187)-3.6;
}

/*
 * Function: ReturnTemperature
 * Usage: grid->T[n]=ReturnTemperaturegrid->xv[n],grid->yv[n],z);
 * ------------------------------------------------------------
 * Helper function to create an initial temperature field.  Used
 * in phys.c in the InitializePhysicalVariables function.
 *
 */
REAL ReturnTemperature(REAL x, REAL y, REAL z, REAL depth) {
  if(x<500)
    return 1;
  return 0;
}

/*
 * Function: ReturnHorizontalVelocity
 * Usage: grid->u[n]=ReturnHorizontalVelocity(grid->xv[n],grid->yv[n],
 *                                            grid->n1[n],grid->n2[n],z);
 * ------------------------------------------------------------
 * Helper function to create an initial velocity field.  Used
 * in phys.c in the InitializePhysicalVariables function.
 *
 */
REAL ReturnHorizontalVelocity(REAL x, REAL y, REAL n1, REAL n2, REAL z) {
  return 0;
}

/*
 * Function: ReturnSediment
 * Usage: SediC[Nsize][n][Nk]=ReturnSediment(grid->xv[n],grid->yv[n],z);
 * ------------------------------------------------------------
 * Helper function to create an initial sediment concentration field.  Used
 * in sediment.c IntitalizeSediment function
 *
 */
REAL ReturnSediment(REAL x, REAL y, REAL z, int sizeno) {
  if(z>-2)
    return 1;
  return 0;
}

/*
 * Function: ReturnBedSedimentRatio
 * Usage: SediC[Nsize][n][Nk]=ReturnBedSedimentRatio(grid->xv[n],grid->yv[n],z);
 * ------------------------------------------------------------
 * Helper function to create an initial bed sediment concentration field.  Used
 * in sediment.c IntitalizeSediment function
 * the sum of ratio should be 1
 */
REAL ReturnBedSedimentRatio(REAL x, REAL y, int layer, int sizeno,int nsize) {
  REAL a;
  a=1.0/nsize;
  return a;
}
