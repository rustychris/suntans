/*
* NetCDF IO functions
* -------------------
* All functions either generic or specific involved in netcdf io should go here
* 
*/ 

#include "mynetcdf.h"
#include "merge.h"
#include "sediments.h"
#include "turbulence.h" 

/***********************************************
* Private functions
***********************************************/

const void* FillValue(int empty);
static void ravel(REAL **tmparray, REAL *tmpvec, gridT *grid);
static void ravelW(REAL **tmparray, REAL *tmpvec,gridT *grid);
static void ravelEdge(REAL **tmparray, REAL *tmpvec, gridT *grid);
static void nc_addattr(int ncid, int varid, char *attname, char *attvalue);
static void nc_addattr_int(int ncid, int varid, char *attname, int *attvalue);
static void nc_addattr_real(int ncid, int varid, char *attname, REAL *attvalue);

void nc_read_3D(int ncid, char *vname, size_t *start, size_t *count, REAL ***tmparray);
void nc_read_2D(int ncid, char *vname, size_t *start, size_t *count, REAL **tmparray, int myproc);
void nc_write_double(int ncid, char *vname, REAL *tmparray, int myproc);
void nc_write_int(int ncid, char *vname, int *tmparray, int myproc);
void nc_write_intvar(int ncid, char *vname, gridT *grid, int *tmparray, int myproc);
void nc_write_doublevar(int ncid, char *vname, gridT *grid, REAL *tmparray, int myproc);

static void nc_write_2D_merge(int ncid, int tstep, REAL *array, propT *prop, gridT *grid,
                              char *varname, int numprocs, int myproc, MPI_Comm comm);
static void nc_write_2D_merge_int(int ncid, int tstep, int *array, propT *prop, gridT *grid,
                                  char *varname, int numprocs, int myproc, MPI_Comm comm);
static void nc_write_3D_merge(int ncid, int tstep, REAL **array, propT *prop, gridT *grid,
                              char *varname, int isw, int numprocs, int myproc, MPI_Comm comm);
static void nc_write_3Dedge_merge(int ncid, int tstep, REAL **array, propT *prop, gridT *grid,
                                  char *varname,int isw, int numprocs, int myproc, MPI_Comm comm);
static void nc_write_2Dedge_merge(int ncid, int tstep, REAL *array, propT *prop, gridT *grid,
                                  char *varname, int numprocs, int myproc, MPI_Comm comm);

static void InitialiseOutputNCugridMerge(propT *prop, physT *phys, gridT *grid, metT *met, int myproc);

static void cell_centered_bed_stress_center(physT *phys, gridT *grid,REAL *taux,REAL *tauy);

/* Compute cell-centered components of bed stress in the same manner as sediments.c,
   using cell-centered velocity.
*/
static void cell_centered_bed_stress_center(physT *phys, gridT *grid, REAL *taux, REAL *tauy) {
  /*
   * taux: destination for x component of cell-centered be stress [Nc]
   * tauy: destination for y component of cell-centered be stress [Nc]
   */
  int k, n, j, nf, i, iptr;
  REAL u,v,u2,CdB,tauB,umag;

  for(i=0;i<grid->Nc;i++) { taux[i]=tauy[i]=0.0; }
  
  // for each computational cell (non-stage defined)
  for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    n=grid->cellp[iptr];
    
    u=phys->uc[n][grid->Nk[n]-1];
    v=phys->vc[n][grid->Nk[n]-1];
    
    u2=pow(u,2)+pow(v,2);

    // cell-averaged drag coefficient
    CdB=0.0;
    for(nf=0;nf<grid->nfaces[n];nf++) {
      j = grid->face[n*grid->maxfaces+nf];
      CdB+=phys->CdB[j];
    }
    CdB /= grid->nfaces[n];
    
    tauB=CdB*u2;
    umag=sqrt(u2);
    taux[n]=tauB * u/umag;
    tauy[n]=tauB * v/umag;
  }
}

/*########################################################
*
* General functions
*
*#########################################################*/
/*
 * Function: MPI_NCOpen
 * Usage: fid = MPI_NCOpen(string,NC_NOWRITE,"GetValue",myproc);
 * -----------------------------------------------------
 * Exits if the requested file does not exist and closes
 * MPI cleanly. The third string is useful for determining which
 * function the function was called from.  When two processes
 * are trying to read the same file at the same time, an error
 * code of EAGAIN results.  This is ommitted as a possible error
 * code.
 *
 */
int MPI_NCOpen(char *file, int perms, char *caller, int myproc) {
  extern int errno;
  int ncid;
  int retval;
  char str[BUFFERLENGTH];

  if (perms==NC_NOWRITE){
    // Just open the file for read access
    if ( (VERBOSE>1) && (myproc==0) ) printf("Opening netcdf file: %s\n",file) ;
    if ((retval = nc_open(file,perms, &ncid)))
      ERR(retval);
  } else {
    // Create a new netcdf dataset
    if (VERBOSE>1) printf("Creating netcdf file: %s\n",file) ;
    // RH: add share flag on the off-chance it permits access to data
    // during the run
    if ((retval = nc_create(file,perms|NC_SHARE, &ncid)))
      ERR(retval);
  }
  if (retval){
    printf("Error in Function %s while trying to open %s ",caller,file);
    printf("Error: %s\n", nc_strerror(retval));
    MPI_Finalize();
    exit(EXIT_FAILURE);
  } else {
    //if ( (VERBOSE>2) ) printf("Successfully opened file: %s on processor %d \n",file,myproc) ;
    return ncid;
  }
}

/*
 * Function: MPI_NCClose(int ncid)
 * -------------------------------
 * Wrapper function for nc_close()
 */
int MPI_NCClose(int ncid){
    int retval;
    if ((retval = nc_close(ncid)))
	ERR(retval);
    return retval;
}


/*
* Function: nc_read_3D()
* ----------------------
* Reads a 3D array from a netcdf file and returns the output in an array (not a vector)
*
* Warning: there are no dimension checks performed here so be careful.
* The size of the array dimension should equal 'count' ie [n, k ,j]
*/

void nc_read_3D(int ncid, char *vname, size_t start[3], size_t count[3], REAL ***tmparray){

    int j, k, n, ii;
    int varid, retval;
    //REAL tmpvec[ (int)count[0] * (int)count[1] * (int)count[2] ];
    REAL outdata[(int)count[0]][(int)count[1]][(int)count[2]];

    //Read the data
    if ((retval = nc_inq_varid(ncid, vname, &varid))) {
      ERRM(retval,vname);
    }
    //    if ((retval = nc_get_vara_double(ncid, varid, start, count, &tmparray[0][0][0])))
    //	ERR(retval);
    if ((retval = nc_get_vara_double(ncid, varid, start, count, &outdata[0][0][0])))  {
      // RH: gcc complains that size_t is effectively long unsigned, and %d
      // is for int.
      printf("start: %ld %ld %ld\n",start[0],start[1],start[2]);
      printf("count: %ld %ld %ld\n",count[0],count[1],count[2]);
      printf("varname: %s\n",vname);
      ERRM(retval,"nc_get_vara_double");
    }

    // Loop through and insert the vector values into an array
    for(n=0;n<(int)count[0];n++){
	for(k=0;k<(int)count[1];k++){
	    for(j=0;j<(int)count[2];j++){
		//Linear index
		//ii = n*(int)count[1]*(int)count[2]+k*(int)count[2]+j;
		//tmparray[n][k][j]=tmpvec[ii];
		tmparray[n][k][j]=outdata[n][k][j];
	    }
	}
    }

}// End function


/*
* Function: nc_read_2D()
* ----------------------
* Reads a 2D array from a netcdf file and returns the output in an array (not a vector)
*
* Warning: there are no dimension checks performed here so be careful.
* The size of the array dimension should equal 'count' ie [n,,j]
*/

void nc_read_2D(int ncid, char *vname, size_t start[2], size_t count[2], REAL **tmparray, int myproc){

    int j, n, ii;
    int varid, retval;
    REAL outdata[(int)count[0]][(int)count[1]];

    // printf("nc_read_2d -- Proc: %d, vname: %s, start[%d][%d], count[%d][%d]",myproc,vname,(int)start[0],(int)start[1],(int)count[0],(int)count[1]);
  
    //Read the data
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    
    if ((retval = nc_get_vara_double(ncid, varid, start, count, &outdata[0][0]))) 
	ERR(retval);
    
    // Loop through and insert the vector values into an array
    for(n=0;n<(int)count[0];n++){
	for(j=0;j<(int)count[1];j++){
	    //Linear index
	    //ii = n*(int)count[1]+j;
	    //tmparray[n][j]=tmpvec[ii];
	    tmparray[n][j]=outdata[n][j];
	    //printf("myproc: %d, start[0]: %d, n: %d of %d, j: %d of %d, outdata[n][j]: %f, tmparray[n][j]: %f\n",myproc, (int)start[0], n,(int)count[0],j,(int)count[1],outdata[n][j], tmparray[n][j]);
	}
    }
    //printf(" Done\n");
}// End function

/*
 * Function: nc_write_double()
 * --------------------------
 *
 * Wrapper function for writing a double variable
 *
 */
void nc_write_double(int ncid, char *vname, REAL *tmparray, int myproc){
    int varid, retval;

    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_put_var_double(ncid, varid, tmparray)))
      ERR(retval);

} //end function

/*
 * Function: nc_write_int()
 * --------------------------
 *
 * Wrapper function for writing an integer variable
 *
 */
void nc_write_int(int ncid, char *vname, int *tmparray, int myproc){
    int varid, retval;

    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_put_var_int(ncid, varid, tmparray)))
      ERR(retval);

} //end function

/*
 * Function: nc_write_intvar()
 * --------------------------
 *
 * Wrapper function for writing a variable length integer variable
 *
 */
void nc_write_intvar(int ncid, char *vname, gridT *grid, int *tmparray, int myproc){
    int varid, retval,i,nf;
    int *tmpvar;
    tmpvar = (int *)SunMalloc((grid->Nc*grid->maxfaces)*sizeof(int),"nc_write_intvar");

    for(i=0;i<grid->Nc;i++){
	for(nf=0;nf<grid->maxfaces;nf++){
	    if(nf < grid->nfaces[i])
	    	tmpvar[i*grid->maxfaces+nf]=tmparray[i*grid->maxfaces+nf];
	    else
	    	tmpvar[i*grid->maxfaces+nf]=(int)EMPTY;
	}
    }
    
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_put_var_int(ncid, varid, tmpvar)))
      ERR(retval);

} //end function

/*
 * Function: nc_write_doublevar()
 * --------------------------
 *
 * Wrapper function for writing a variable length double variable
 *
 */
void nc_write_doublevar(int ncid, char *vname, gridT *grid, REAL *tmparray, int myproc){
    int varid, retval,j,nf;
    REAL *tmpvar;
    tmpvar = (REAL *)SunMalloc(grid->Nc*grid->maxfaces*sizeof(REAL),"nc_write_intvar");

    for(j=0;j<grid->Nc;j++){
	for(nf=0;nf<grid->maxfaces;nf++){
	    if(nf < grid->nfaces[j])
	    	tmpvar[j*grid->maxfaces+nf]=tmparray[j*grid->maxfaces+nf];
	    else
	    	tmpvar[j*grid->maxfaces+nf]=(REAL)EMPTY;
	}
    }
    
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_put_var_double(ncid, varid, tmpvar)))
      ERR(retval);

} //end function



/*
* Function: getTimeRec()
* -----------------------------
*  Retuns the index of the first preceding time step in the vector time
*  May return out-of-bounds values, -1 or nt.
*/
int getTimeRec(REAL nctime, REAL *time, int nt){
   int j;

   for(j=0;j<nt;j++){
     if (time[j]>=nctime) {
       return j-1;
     }
   }
   return nt;
}

/*
 * Function: nc_write_2D_merge()
 * -------------------------------
 * Merges a 2D (time-varying) variable and writes to netcdf
 */
static void nc_write_2D_merge(int ncid, int tstep, REAL *array, propT *prop, gridT *grid, char *varname, int numprocs, int myproc, MPI_Comm comm){

   int varid, retval;
   //size_t starttwo[] = {prop->nctimectr,0};
   size_t starttwo[] = {tstep,0};
   size_t counttwo[] = {1,grid->Nc};

    MergeCellCentered2DArray(array,grid,numprocs,myproc,comm);

    if(myproc==0){
    	counttwo[1] = mergedGrid->Nc;
	if ((retval = nc_inq_varid(ncid, varname, &varid)))
	    ERR(retval);
	if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, merged2DArray )))
	    ERR(retval);
    }
}


/*
 * Function: nc_write_2D_merge_int()
 * -------------------------------
 * Merges a 2D (time-varying) integer variable and writes to netcdf
 */
static void nc_write_2D_merge_int(int ncid, int tstep, int *array, propT *prop, gridT *grid,
                                  char *varname, int numprocs, int myproc, MPI_Comm comm){
   int varid, retval;
   size_t starttwo[] = {tstep,0};
   size_t counttwo[] = {1,grid->Nc};

   MergeCellCentered2DArray_int(array,grid,numprocs,myproc,comm);

   if(myproc==0){
     counttwo[1] = mergedGrid->Nc;
     if ((retval = nc_inq_varid(ncid, varname, &varid)))
       ERR(retval);
     if ((retval = nc_put_vara_int(ncid, varid, starttwo, counttwo, merged2DArray_int )))
       ERR(retval);
   }
}



/*
 * Function: nc_write_3D_merge()
 * -------------------------------
 * Merges a 3D (time-varying) variable and writes to netcdf
 */
static void nc_write_3D_merge(int ncid, int tstep, REAL **array, propT *prop, gridT *grid, char *varname,int isw, int numprocs, int myproc, MPI_Comm comm){

   int varid, retval,i,k;
   //size_t startthree[] = {prop->nctimectr,0,0};
   size_t startthree[] = {tstep,0,0};
   size_t countthree[] = {1,grid->Nkmax,grid->Nc};

    MergeCellCentered3DArray(array,grid,numprocs,myproc,comm);   

    if(myproc==0){
    	countthree[2]=mergedGrid->Nc;
	if ((retval = nc_inq_varid(ncid, varname, &varid)))
	    ERR(retval);

	//Roll the array out into a vector
        for(i=0;i<mergedGrid->Nc;i++){
	    for(k=0;k<mergedGrid->Nkmax;k++){
	      if(k<mergedGrid->Nk[i]){
		merged3DVector[k*mergedGrid->Nc+i] = merged3DArray[i][k];
	      }else{
		merged3DVector[k*mergedGrid->Nc+i] = (REAL)EMPTY;
	      }
	    }
        }
	if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, &merged3DVector[0])))
	    ERR(retval);
    }
}


/*
 * Function: nc_write_2Dedge_merge()
 * -------------------------------
 * Merges a 2D (time-varying) edge variable and writes to netcdf
 */
static void nc_write_2Dedge_merge(int ncid, int tstep, REAL *array, propT *prop, gridT *grid, char *varname, int numprocs, int myproc, MPI_Comm comm){

   int varid, retval,i,k;
   size_t starttwo[] = {tstep,0};
   size_t counttwo[] = {1,-1 /* Ne */};

   // Places result directly in merged3DVector, using only the first Ne
   // elements
   MergeEdgeCentered2DArray(array,grid,numprocs,myproc,comm);   

   if(myproc==0){
     counttwo[1]=mergedGrid->Ne;
     if ((retval = nc_inq_varid(ncid, varname, &varid)))
       ERR(retval);

     if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, &merged3DVector[0])))
       ERR(retval);
   }
}

/*
 * Function: nc_write_3Dedge_merge()
 * -------------------------------
 * Merges a 3D (time-varying) edge variable and writes to netcdf
 */
static void nc_write_3Dedge_merge(int ncid, int tstep, REAL **array, propT *prop, gridT *grid, char *varname,int isw, int numprocs, int myproc, MPI_Comm comm){

   int varid, retval,i,k;
   //size_t startthree[] = {prop->nctimectr,0,0};
   size_t startthree[] = {tstep,0,0};
   size_t countthree[] = {1,grid->Nkmax,grid->Ne};

    MergeEdgeCentered3DArray(array,grid,numprocs,myproc,comm);   

    if(myproc==0){
    	countthree[2]=mergedGrid->Ne;
	if ((retval = nc_inq_varid(ncid, varname, &varid)))
	    ERR(retval);

	//Roll the array out into a vector
        for(i=0;i<mergedGrid->Ne;i++){
	    for(k=0;k<mergedGrid->Nkmax;k++){
	      if(k<mergedGrid->Nke[i]){
		merged3DVector[k*mergedGrid->Ne+i] = merged3DEArray[i][k];
	      }else{
		merged3DVector[k*mergedGrid->Ne+i] = (REAL)EMPTY;
	      }
	    }
        }
	if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, &merged3DVector[0])))
	    ERR(retval);
    }
}

/*###############################################################
*
* SUNTANS output file functions
*
#################################################################*/
/*
* Function: WriteOutputNCmerge()
* -----------------------------
* Write SUNTANS output to netcdf file/s merged onto one core.
* 
*/
void WriteOutputNCmerge(propT *prop, gridT *grid, physT *phys, metT *met, int blowup, int numprocs, int myproc, MPI_Comm comm){
  int iptr,i; // RH for debugging

  int ncid;
  int varid, retval, k;
  // Start and count vectors for one, two and three dimensional arrays
  size_t startone[] = {prop->nctimectr};
  size_t countone[] = {1};
  size_t starttwo[] = {prop->nctimectr,0};
  size_t counttwo[] = {1,grid->Nc};
  size_t startthree[] = {prop->nctimectr,0,0};
  size_t countthree[] = {1,grid->Nkmax,grid->Nc};
  const size_t countthreew[] = {1,grid->Nkmax+1,grid->Nc};
  const REAL time[] = {prop->nctime};
  char str[2*BUFFERLENGTH], filename[BUFFERLENGTH];

  nc_set_log_level(3); // This helps with debugging errors
   
  // RH: this used to be prop->n==1+prop->nstart, but that misses
  // the first output step.
  if(!(prop->n%prop->ntout) || prop->n==prop->nstart || blowup) {
    
    // RH: same, used to be prop->n==1+prop->nstart
    if( (prop->nctimectr>=prop->nstepsperncfile) || prop->n==prop->nstart){
      if(prop->n > 1+prop->nstart){
	// Close the old netcdf file
	if(myproc==0){
	  printf("Closing opened output netcdf file...\n");
	  MPI_NCClose(prop->outputNetcdfFileID);
	}
      }
       
      // Open the new netcdf file
      MPI_GetFile(filename,DATAFILE,"outputNetcdfFile","WriteOutputNCmerge",myproc);

      sprintf(str,"%s_%04d.nc",filename,prop->ncfilectr);
      if(myproc==0){
	prop->outputNetcdfFileID = MPI_NCOpen(str,NC_CLASSIC_MODEL|NC_NETCDF4,"WriteOutputNCmerge",myproc);
	// Initialise a new output file
	InitialiseOutputNCugridMerge(prop, phys, grid, met, myproc);
      }else{
	prop->outputNetcdfFileID=-1;
      }
	
      // Reset the time counter
      prop->nctimectr = 0;
      prop->ncfilectr += 1;
      startone[0] = prop->nctimectr;
    }
    ncid = prop->outputNetcdfFileID;

    if(myproc==0 && VERBOSE>1){ 
      if(!blowup) {
	printf("Outputting data to netcdf at step %d of %d\n",prop->n,prop->nsteps+prop->nstart);
      } else {
	printf("Outputting blowup data to netcdf at step %d of %d\n",prop->n,prop->nsteps+prop->nstart);
      }
    }
    if(myproc==0){ 
      /* Write the time data*/
      if ((retval = nc_inq_varid(ncid, "time", &varid))) {
	printf("ncid: %d\n",ncid);
	ERRM(retval,"nc_inq_varid('time'), proc 0");
      }
      if ((retval = nc_put_vara_double(ncid, varid, startone, countone, time ))) {
	ERRM(retval,"proc 0, save time");
      }
       
      countthree[2] = mergedGrid->Nc;
      counttwo[1] = mergedGrid->Nc;
    }
     
    /* Write to the physical variables*/
     
    // 2D cell-centered variables
    nc_write_2D_merge(ncid, prop->nctimectr, phys->h, prop, grid, "eta", numprocs, myproc, comm);

    // Need a place to work up cell centered data
    //  phys->tmpvar is 3D cell centered, Nkmax
    //  phys->tmpvarW is 3D cell centered, Nkmax+1
    double *tmpx=phys->tmpvar; // [Nc*Nk], but we use only need to be [Nc]
    double *tmpy=phys->tmpvarW; // [Nc*(Nk+1)], but here we use only [Nc]
    cell_centered_bed_stress_interp(phys,grid,tmpx,tmpy); // also updates phys->tau_B
    nc_write_2D_merge(ncid, prop->nctimectr, tmpx, prop, grid, "tauB_x", numprocs, myproc, comm);
    nc_write_2D_merge(ncid, prop->nctimectr, tmpy, prop, grid, "tauB_y", numprocs, myproc, comm);
    nc_write_2Dedge_merge(ncid, prop->nctimectr, phys->tau_B, prop, grid, "tauB", numprocs, myproc, comm); // no var defined yet

    cell_centered_bed_stress_center(phys,grid,tmpx,tmpy);
    nc_write_2D_merge(ncid, prop->nctimectr, tmpx, prop, grid, "tauBc_x", numprocs, myproc, comm);
    nc_write_2D_merge(ncid, prop->nctimectr, tmpy, prop, grid, "tauBc_y", numprocs, myproc, comm);

    // RH: variables to diagnose time limiting cells.
    nc_write_2D_merge_int(ncid, prop->nctimectr, phys->limiting_cell, prop, grid, "limiting_cell", numprocs, myproc, comm);
    nc_write_2D_merge(ncid, prop->nctimectr, phys->min_time_step, prop, grid, "min_time_step", numprocs, myproc, comm);
     
    if(prop->metmodel>0) {
      // Atmospheric flux variables
      nc_write_2D_merge(ncid,prop->nctimectr, met->Uwind, prop, grid, "Uwind", numprocs, myproc, comm);
      nc_write_2D_merge(ncid,prop->nctimectr,  met->Vwind, prop, grid, "Vwind", numprocs, myproc, comm);
      nc_write_2D_merge(ncid,prop->nctimectr,  met->Tair, prop, grid, "Tair", numprocs, myproc, comm);
      nc_write_2D_merge(ncid,prop->nctimectr,  met->Pair, prop, grid, "Pair", numprocs, myproc, comm);
      nc_write_2D_merge(ncid,prop->nctimectr,  met->rain, prop, grid, "rain", numprocs, myproc, comm);
      nc_write_2D_merge(ncid,prop->nctimectr,  met->RH, prop, grid, "RH", numprocs, myproc, comm);
      nc_write_2D_merge(ncid,prop->nctimectr,  met->cloud, prop, grid, "cloud", numprocs, myproc, comm);
      nc_write_2D_merge(ncid,prop->nctimectr,  met->Hs, prop, grid, "Hs", numprocs, myproc, comm);
      nc_write_2D_merge(ncid,prop->nctimectr,  met->Hl, prop, grid, "Hl", numprocs, myproc, comm);
      nc_write_2D_merge(ncid,prop->nctimectr,  met->Hlw, prop, grid, "Hlw", numprocs, myproc, comm);
      nc_write_2D_merge(ncid,prop->nctimectr,  met->Hsw, prop, grid, "Hsw", numprocs, myproc, comm);
      nc_write_2D_merge(ncid,prop->nctimectr,  met->tau_x, prop, grid, "tau_x", numprocs, myproc, comm);
      nc_write_2D_merge(ncid,prop->nctimectr,  met->tau_y, prop, grid, "tau_y", numprocs, myproc, comm);
      if(prop->beta > 0) {
	nc_write_2D_merge(ncid,prop->nctimectr,  met->EP, prop, grid, "EP", numprocs, myproc, comm);
      }
    }
    // 3D cell-centered variables
    nc_write_3D_merge(ncid,prop->nctimectr,  phys->uc, prop, grid, "uc",0, numprocs, myproc, comm);
    nc_write_3D_merge(ncid,prop->nctimectr,  phys->vc, prop, grid, "vc",0, numprocs, myproc, comm);
    nc_write_3D_merge(ncid,prop->nctimectr,  phys->nu_tv, prop, grid, "nu_v",0, numprocs, myproc, comm);
     
    nc_write_3D_merge(ncid,prop->nctimectr,  phys->qT, prop, grid, "turb_q",0, numprocs, myproc, comm);
    nc_write_3D_merge(ncid,prop->nctimectr,  phys->lT, prop, grid, "turb_l",0, numprocs, myproc, comm);
     
    if(prop->beta>0)
      nc_write_3D_merge(ncid,prop->nctimectr,  phys->s, prop, grid, "salt",0, numprocs, myproc, comm);
     
    if(prop->gamma>0)
      nc_write_3D_merge(ncid,prop->nctimectr,  phys->T, prop, grid, "temp",0, numprocs, myproc, comm);
     
    if( (prop->gamma>0) || (prop->beta>0) ) 
      nc_write_3D_merge(ncid,prop->nctimectr,  phys->rho, prop, grid, "rho",0, numprocs, myproc, comm);
     
    if(prop->calcage){
      nc_write_3D_merge(ncid,prop->nctimectr,  age->agec, prop, grid, "agec",0, numprocs, myproc, comm);
      nc_write_3D_merge(ncid,prop->nctimectr,  age->agealpha, prop, grid, "agealpha",0, numprocs, myproc, comm);
    }
     
    // Vertical velocity 
    nc_write_3D_merge(ncid,prop->nctimectr,  phys->w, prop, grid, "w",1, numprocs, myproc, comm);
     
    // 3D edge-based variables 
    nc_write_3Dedge_merge(ncid,prop->nctimectr,  phys->u, prop, grid, "U",0, numprocs, myproc, comm);
     
     
    /* Update the time counter*/
    prop->nctimectr += 1;  
  }
   
} // End of function


/*
* Function: WriteOutputNC()
* -----------------------------
* Main function for writing SUNTANS output to netcdf file/s
* 
*/
void WriteOutputNC(propT *prop, gridT *grid, physT *phys, metT *met, int blowup, int myproc){
   int ncid = prop->outputNetcdfFileID;
   int varid, retval, k;
   char varname[256];
   // Start and count vectors for one, two and three dimensional arrays
   const size_t startone[] = {prop->nctimectr};
   const size_t countone[] = {1};
   const size_t starttwo[] = {prop->nctimectr,0};
   const size_t counttwo[] = {1,grid->Nc};
   const size_t count_t_Ne[] = {1,grid->Ne};
   const size_t startthree[] = {prop->nctimectr,0,0};
   const size_t countthree[] = {1,grid->Nkmax,grid->Nc};
   const size_t countthree_e[] = {1,grid->Nkmax,grid->Ne};
   const size_t countthreew[] = {1,grid->Nkmax+1,grid->Nc};
   const REAL time[] = {prop->nctime};

   nc_set_log_level(3); // This helps with debugging errors
   
   if(!(prop->n%prop->ntout) || prop->n==1+prop->nstart || blowup) {

    if(myproc==0 && VERBOSE>1){ 
      if(!blowup) 
        printf("Outputting data to netcdf at step %d of %d\n",prop->n,prop->nsteps+prop->nstart);
      else
        printf("Outputting blowup data to netcdf at step %d of %d\n",prop->n,prop->nsteps+prop->nstart);

      // debugging a netcdf error:
      printf("time start: %d\n",prop->nctimectr);
    }
    
    /* Write the time data*/
    if ((retval = nc_inq_varid(ncid, "time", &varid)))
      ERRM(retval,"time id");
    if ((retval = nc_put_vara_double(ncid, varid, startone, countone, time )))
      ERRM(retval,"time");
    
    /* Write to the physical variables*/
    if ((retval = nc_inq_varid(ncid, "eta", &varid)))
      ERRM(retval,"eta id");
    if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, phys->h )))
      ERRM(retval,"eta");

    if ((retval = nc_inq_varid(ncid, "limiting_cell", &varid)))
      ERRM(retval,"limit_cell id");
    if ((retval = nc_put_vara_int(ncid, varid, starttwo, counttwo, phys->limiting_cell )))
      ERRM(retval,"limiting_cell");

    if ((retval = nc_inq_varid(ncid, "min_time_step", &varid)))
      ERRM(retval,"min_time_step");
    if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, phys->min_time_step )))
      ERRM(retval,"min_time_step");
    
    if ((retval = nc_inq_varid(ncid, "uc", &varid)))
      ERRM(retval,"uc id");
    ravel(phys->uc, phys->tmpvar, grid);
    if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, phys->tmpvar )))
      ERRM(retval,"uc");
    
    if ((retval = nc_inq_varid(ncid, "vc", &varid)))
      ERRM(retval,"vc id");
    ravel(phys->vc, phys->tmpvar, grid);
    if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, phys->tmpvar )))
      ERRM(retval,"vc");
      
    // write w at cell top and bottom
    if ((retval = nc_inq_varid(ncid, "w", &varid)))
      ERRM(retval,"w id");
    ravelW(phys->w, phys->tmpvarW, grid);
    if ((retval = nc_put_vara_double(ncid, varid, startthree, countthreew, phys->tmpvarW )))
      ERRM(retval,"w");

    if ((retval = nc_inq_varid(ncid, "nu_v", &varid)))
      ERRM(retval,"nu_v id");
    ravel(phys->nu_tv, phys->tmpvar, grid);
    if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, phys->tmpvar )))
      ERRM(retval,"nu_v");

    if ((retval = nc_inq_varid(ncid, "turb_q", &varid)))
      ERRM(retval,"turb_q id");
    ravel(phys->qT, phys->tmpvar, grid);
    if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, phys->tmpvar )))
      ERRM(retval,"turb_q");
    
    if ((retval = nc_inq_varid(ncid, "turb_l", &varid)))
      ERRM(retval,"turb_l id");
    ravel(phys->lT, phys->tmpvar, grid);
    if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, phys->tmpvar )))
      ERRM(retval,"turb_l");
    
    // Tracers
    if(prop->beta>0){
      if ((retval = nc_inq_varid(ncid, "salt", &varid)))
        ERRM(retval,"salt id");
      ravel(phys->s, phys->tmpvar, grid);
      if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, phys->tmpvar )))
        ERRM(retval,"salt");
    }
     
    if(prop->gamma>0){
	if ((retval = nc_inq_varid(ncid, "temp", &varid)))
	  ERRM(retval,"temp id");
	ravel(phys->T, phys->tmpvar, grid);
	if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, phys->tmpvar )))
	  ERRM(retval,"temp");
     }
      
     if( (prop->gamma>0) || (prop->beta>0) ){ 
	if ((retval = nc_inq_varid(ncid, "rho", &varid)))
	  ERRM(retval,"rho id");
	ravel(phys->rho, phys->tmpvar, grid);
	if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, phys->tmpvar )))
	  ERRM(retval,"rho");
     }

     if(prop->calcage>0){ 
       if ((retval = nc_inq_varid(ncid, "agec", &varid)))
         ERRM(retval,"agec id");
       ravel(age->agec, phys->tmpvar, grid);
       if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, phys->tmpvar )))
         ERRM(retval,"agec");
       
       if ((retval = nc_inq_varid(ncid, "agealpha", &varid)))
         ERRM(retval,"agealpha id");
       ravel(age->agealpha, phys->tmpvar, grid);
       if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, phys->tmpvar )))
         ERRM(retval,"agealpha");
     }

     // Sediment output:
     if(prop->computeSediments){
       for(k=0;k<sediments->Nsize;k++){
         sprintf(varname,"sed%d",k+1); // 1-based to match suntans.dat
         if ((retval = nc_inq_varid(ncid, varname, &varid)))
           ERRM(retval,"sed var id");
         ravel(sediments->SediC[k], phys->tmpvar, grid);
         if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, phys->tmpvar )))
           ERRM(retval,"sed conc");
       }
     }
     
     // Vertical grid spacing
     if ((retval = nc_inq_varid(ncid, "dzz", &varid)))
       ERR(retval);
     ravel(grid->dzz, phys->tmpvar, grid);
     if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, phys->tmpvar )))
       ERR(retval);
     
     if ((retval = nc_inq_varid(ncid, "dzf", &varid)))
       ERR(retval);
     ravelEdge(grid->dzf, phys->tmpvarE, grid);
     if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree_e, phys->tmpvarE )))
       ERR(retval);

     if ((retval = nc_inq_varid(ncid, "etop", &varid)))
       ERR(retval);
     if ((retval = nc_put_vara_int(ncid, varid, starttwo, count_t_Ne, grid->etop ) ))
       ERR(retval);

     // Edge normal velocity
     if ((retval = nc_inq_varid(ncid, "U", &varid)))
       ERR(retval);
     ravelEdge(phys->u, phys->tmpvarE, grid);
     if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree_e, phys->tmpvarE )))
       ERR(retval);

     // Wind variables
     if(prop->metmodel>0){
       if ((retval = nc_inq_varid(ncid, "Uwind", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, met->Uwind )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "Vwind", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, met->Vwind )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "Tair", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, met->Tair )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "Pair", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, met->Pair )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "rain", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, met->rain )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "RH", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, met->RH )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "cloud", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, met->cloud )))
	 ERR(retval);
       
       // Heat flux variables
       if ((retval = nc_inq_varid(ncid, "Hs", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, met->Hs )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "Hl", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, met->Hl )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "Hlw", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, met->Hlw )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "Hsw", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, met->Hsw )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "tau_x", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, met->tau_x )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "tau_y", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, met->tau_y )))
	 ERR(retval);
       
       if(prop->beta > 0.0){
	  if ((retval = nc_inq_varid(ncid, "EP", &varid)))
	     ERR(retval);
	  if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, met->EP )))
	     ERR(retval);
       }
     }
     
    /* Update the time counter*/
    prop->nctimectr += 1;
   }
   // Free the temporary vector
   //SunFree(tmpvar,grid->Nc*grid->Nkmax*sizeof(REAL),"WriteOuputNC");
   //SunFree(tmpvarE,grid->Ne*grid->Nkmax*sizeof(REAL),"WriteOuputNC");

} // End of function


/*
* Function: InitialiseOutputNCugridMerge()
* ------------------------------------
*
* Initialises the output netcdf file/s
* Files conform to the UGRID-0.9 CF conventions
* 
* Merges all data onto one processor
* The pointer to each file is stored in prop->outputNetcdfFileID
*
* The properties of the grid are stored in the global mergedGrid structure.
* 
*/
static void InitialiseOutputNCugridMerge(propT *prop, physT *phys, gridT *grid, metT *met, int myproc){
   int ncid = prop->outputNetcdfFileID;
   int retval, k, n, j, i;
   int varid;
   int dimid_Nc, dimid_Ne , dimid_Np, dimid_time, dimid_numsides, dimid_Two, dimid_Nkw, dimid_Nk; 
   int dimidone[1];
   int dimidtwo[2];
   int dimidthree[3];
   int nofill=0;
   const size_t starttwo[] = {0,0};
   const size_t counttwo[] = {mergedGrid->Nkmax,mergedGrid->Nc};
   REAL *z_r;
   REAL *z_w;
   int *edges;
   const int DEFLATE=1;
   const int DEFLATELEVEL=2;
   const REAL FILLVALUE = (REAL)EMPTY;
   int num_unlimdims;


   //REAL *tmpvar;
   // Need to write the 3-D arrays as vectors
   //tmpvar = (REAL *)SunMalloc(grid->Nc*grid->Nkmax*sizeof(REAL),"InitialiseOutputNC");  
   
   ///* Initialise the depth arrays */
   z_r = (REAL *)SunMalloc((grid->Nkmax)*sizeof(REAL),"InitialiseOutputNCugrid");
   z_w = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"InitialiseOutputNCugrid");

   ///* Initialize an edge array */
   //edges = (int *)SunMalloc(2*(grid->Ne)*sizeof(int),"InitialiseOutputNCugrid");
   
   /**************
    *
    * Start writing...
    *
    **************/
   if(VERBOSE>1 && myproc==0) printf("Initialising output netcdf files...");

   // Set the netcdf time ctr to 0
   prop->nctimectr=0;
      
   /********************************************************************** 
    *
    * Define the dimensions
    *
    **********************************************************************/

   if ((retval = nc_def_dim(ncid, "Nc", mergedGrid->Nc, &dimid_Nc)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "Np", grid->Np, &dimid_Np)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "Ne", mergedGrid->Ne, &dimid_Ne)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "Nk", grid->Nkmax, &dimid_Nk)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "Nkw", grid->Nkmax+1, &dimid_Nkw)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "numsides", grid->maxfaces, &dimid_numsides)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "Two", 2, &dimid_Two)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "time", NC_UNLIMITED, &dimid_time)))
	ERR(retval);


   /*
    * Define the global attributes - this should be expanded to include model input parameters 
    *
    */
   nc_addattr(ncid, NC_GLOBAL,"title","SUNTANS NetCDF output file");
    /********************************************************************** 
    *
    * Define the grid topology variables and attributes
    *
    **********************************************************************/

    //suntans_mesh
    if ((retval = nc_def_var(ncid,"suntans_mesh",NC_INT,0,0,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"cf_role","mesh_topology");
    nc_addattr(ncid, varid,"long_name","Topology data of 2D unstructured mesh");
    nc_addattr(ncid, varid,"topology_dimension","2");
    nc_addattr(ncid, varid,"node_coordinates","xp yp");
    nc_addattr(ncid, varid,"face_node_connectivity","cells");
    nc_addattr(ncid, varid,"edge_node_connectivity","edges");
    nc_addattr(ncid, varid,"face_coordinates","xv yv");
    nc_addattr(ncid, varid,"edge_coordinates","xe ye");
    nc_addattr(ncid, varid,"face_edge_connectivity","face");
    nc_addattr(ncid, varid,"edge_face_connectivity","grad");

    // cells
    dimidtwo[0] = dimid_Nc;
    dimidtwo[1] = dimid_numsides;
    if ((retval = nc_def_var(ncid,"cells",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"cf_role","face_node_connectivity");
    nc_addattr(ncid, varid,"long_name","Maps every face to its corner nodes");
    //if ((retval = nc_put_var_int(ncid,varid, grid->cells)))
    //  ERR(retval);

        //nfaces
    dimidone[0] = dimid_Nc;
    if ((retval = nc_def_var(ncid,"nfaces",NC_INT,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Number of cell faces");
    nc_addattr(ncid, varid,"units","");
    nc_addattr(ncid, varid,"coordinates","xv yv");

    //face
    dimidtwo[0] = dimid_Nc;
    dimidtwo[1] = dimid_numsides;
    if ((retval = nc_def_var(ncid,"face",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"cf_role","face_edge_connectivity");
    nc_addattr(ncid, varid,"long_name","Maps every face to its edges");
    //if ((retval = nc_put_var_int(ncid,varid, grid->face)))
    //  ERR(retval);

    //edges
    dimidtwo[0] = dimid_Ne;
    dimidtwo[1] = dimid_Two;
    if ((retval = nc_def_var(ncid,"edges",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"cf_role","edge_node_connectivity");
    nc_addattr(ncid, varid,"long_name","Maps every edge to the two nodes it connects");
    //if ((retval = nc_put_var_int(ncid,varid, grid->edges)))
    //  ERR(retval);

    //neigh
    dimidtwo[0] = dimid_Nc;
    dimidtwo[1] = dimid_numsides;
    if ((retval = nc_def_var(ncid,"neigh",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"cf_role","face_face_connectivity");
    nc_addattr(ncid, varid,"long_name","Maps every face to its neighbouring faces");
    //if ((retval = nc_put_var_int(ncid,varid, grid->neigh)))
    //  ERR(retval);

    //grad
    dimidtwo[0] = dimid_Ne;
    dimidtwo[1] = dimid_Two;
    if ((retval = nc_def_var(ncid,"grad",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"cf_role","edge_face_connectivity");
    nc_addattr(ncid, varid,"long_name","Maps every edge to the two faces it connects ");
    //if ((retval = nc_put_var_int(ncid,varid, grid->grad)))
    //  ERR(retval);


    /********************************************************************** 
    *
    * Define the grid coordinate variables and attributes 
    *
    **********************************************************************/

    //xv
    dimidone[0] = dimid_Nc;
    if ((retval = nc_def_var(ncid,"xv",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"standard_name","Easting");
    nc_addattr(ncid, varid,"long_name","Easting of 2D mesh face");
    //if ((retval = nc_put_var_double(ncid,varid, grid->xv)))
    //   ERR(retval);
      
    //yv
    dimidone[0] = dimid_Nc;
    if ((retval = nc_def_var(ncid,"yv",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"standard_name","Northing");
    nc_addattr(ncid, varid,"long_name","Northing of 2D mesh face");
    //if ((retval = nc_put_var_double(ncid,varid, grid->yv)))
    //   ERR(retval);
       
    //xp
    dimidone[0] = dimid_Np;
    if ((retval = nc_def_var(ncid,"xp",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"standard_name","Easting");
    nc_addattr(ncid, varid,"long_name","Easting of 2D mesh node");
    //if ((retval = nc_put_var_double(ncid,varid, grid->xp)))
    //   ERR(retval);
        
    //yp
    dimidone[0] = dimid_Np;
    if ((retval = nc_def_var(ncid,"yp",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"standard_name","Northing");
    nc_addattr(ncid, varid,"long_name","Northing of 2D mesh node");
    //if ((retval = nc_put_var_double(ncid,varid, grid->yp)))
    //   ERR(retval);
         
    //xe
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"xe",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"standard_name","Easting");
    nc_addattr(ncid, varid,"long_name","Easting of 2D mesh edge");
    //if ((retval = nc_put_var_double(ncid,varid, grid->xe)))
    //   ERR(retval);
          
    //ye
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"ye",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"standard_name","Northing");
    nc_addattr(ncid, varid,"long_name","Northing of 2D mesh edge");
    //if ((retval = nc_put_var_double(ncid,varid, grid->ye)))
    //   ERR(retval);
        
    /********************************************************************** 
    *
    * Define the grid metric variables 
    *
    **********************************************************************/
    //normal
    dimidtwo[0] = dimid_Nc;
    dimidtwo[1] = dimid_numsides;
    if ((retval = nc_def_var(ncid,"normal",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Dot product of unique normal with outward normal of each edge");

    //n1
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"n1",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","x-component of the edge normal");
    nc_addattr(ncid, varid,"coordinates","xe ye");

    //n2
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"n2",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","y-component of the edge normal");
    nc_addattr(ncid, varid,"coordinates","xe ye");


    //df
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"df",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","edge length");
    nc_addattr(ncid, varid,"units","m");
    nc_addattr(ncid, varid,"coordinates","xe ye");
    //if ((retval = nc_put_var_double(ncid,varid, grid->df)))
    //  ERR(retval);

    //dg
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"dg",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","distance between faces on either side of edge");
    nc_addattr(ncid, varid,"units","m");
    nc_addattr(ncid, varid,"coordinates","xe ye");
    //if ((retval = nc_put_var_double(ncid,varid, grid->dg)))
    //  ERR(retval);

    //def
    dimidtwo[0] = dimid_Nc;
    dimidtwo[1] = dimid_numsides;
    if ((retval = nc_def_var(ncid,"def",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Distance between faces and edges");
    //if ((retval = nc_put_var_double(ncid,varid, grid->def)))
    //  ERR(retval);
    nc_addattr(ncid, varid,"coordinates","xe ye");
    nc_addattr(ncid, varid,"units","m");
    
    //mark
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"mark",NC_INT,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Edge marker type");
    nc_addattr(ncid, varid,"units","0 - computational; 1 - closed; 2 flux BC; 3 - stage BC; 4 - other BC; 5 - interproc; 6 - ghost cell.");
    nc_addattr(ncid, varid,"coordinates","xe ye");

    //Ac
    dimidone[0] = dimid_Nc;
    if ((retval = nc_def_var(ncid,"Ac",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Horizontal area of 2D mesh");
    nc_addattr(ncid, varid,"units","m2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","xv yv");
    //if ((retval = nc_put_var_double(ncid,varid, grid->Ac)))
    //  ERR(retval);
    /********************************************************************** 
    *
    * Define the vertical grid variables and attributes 
    *
    **********************************************************************/
   //dz
   dimidone[0] = dimid_Nk;   
   if ((retval = nc_def_var(ncid,"dz",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","z layer spacing");
   nc_addattr(ncid, varid,"units","m");
   //if ((retval = nc_put_var_double(ncid,varid, grid->dz)))
   //  ERR(retval);

   // Calculate and write the vertical coordinate z levels 
   z_w[0]=0.0;
   for(k=0;k<grid->Nkmax;k++){
      z_w[k+1] = z_w[k] + grid->dz[k];
      if(k==0){
	 z_r[k] = grid->dz[k]*0.5;
      }else{
	 z_r[k] = z_r[k-1]+grid->dz[k];
      }
   }
   //z_r
   dimidone[0] = dimid_Nk;
   if ((retval = nc_def_var(ncid,"z_r",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"standard_name","ocean_z_coordinate");
   nc_addattr(ncid, varid,"long_name","depth at layer mid points");
   nc_addattr(ncid, varid,"units","m");
   nc_addattr(ncid, varid,"positive","down");
   //if ((retval = nc_put_var_double(ncid,varid, z_r)))
   //  ERR(retval);

   //z_w
   dimidone[0] = dimid_Nkw;
   if ((retval = nc_def_var(ncid,"z_w",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"standard_name","ocean_z_coordinate");
   nc_addattr(ncid, varid,"long_name","depth at layer edges");
   nc_addattr(ncid, varid,"units","m");
   nc_addattr(ncid, varid,"positive","down");
   //if ((retval = nc_put_var_double(ncid,varid, z_w)))
   //  ERR(retval);
   //Nk
   dimidone[0] = dimid_Nc;   
   if ((retval = nc_def_var(ncid,"Nk",NC_INT,1,dimidone,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Number of layers at face"); 
   //if ((retval = nc_put_var_int(ncid,varid, grid->Nk)))
   //  ERR(retval);

   //Nke
   dimidone[0] = dimid_Ne;   
   if ((retval = nc_def_var(ncid,"Nke",NC_INT,1,dimidone,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Number of layers at edge");
   //if ((retval = nc_put_var_int(ncid,varid, grid->Nke)))
   //  ERR(retval);

    //dv
    dimidone[0] = dimid_Nc;
    if ((retval = nc_def_var(ncid,"dv",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"standard_name","sea_floor_depth_below_geoid");
    nc_addattr(ncid, varid,"long_name","seafloor depth");
    nc_addattr(ncid, varid,"units","m");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","xv yv");
    nc_addattr(ncid, varid,"positive","down");
    //if ((retval = nc_put_var_double(ncid,varid, grid->dv)))
    //  ERR(retval);

    //time
    dimidone[0] = dimid_time;
    if ((retval = nc_def_var(ncid,"time",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"long_name","time");
    nc_addattr(ncid, varid,"units","seconds since 1990-01-01 00:00:00");  

   /********************************************************************** 
    * 
    * Define the physical variables and attributes 
    *
    **********************************************************************/
   
   dimidtwo[0] = dimid_time;
   dimidtwo[1] = dimid_Nc;
   
   dimidthree[0] = dimid_time;
   dimidthree[1] = dimid_Nk;
   dimidthree[2] = dimid_Nc;
   
   // eta
   if ((retval = nc_def_var(ncid,"eta",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Sea surface elevation");
   nc_addattr(ncid, varid,"units","m");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"positive","up");
   nc_addattr(ncid, varid,"coordinates","time yv xv");

   // cfl_count
   if ((retval = nc_def_var(ncid,"limiting_cell",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Cell frequency of CFL limit");
   nc_addattr(ncid, varid,"units","-");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time yv xv");

   // minimum time step per cell
   if ((retval = nc_def_var(ncid,"min_time_step",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Minimum time step limitation");
   nc_addattr(ncid, varid,"units","s");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time yv xv");

   
   //u
   if ((retval = nc_def_var(ncid,"uc",NC_DOUBLE,3,dimidthree,&varid)))
      ERR(retval);
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Eastward water velocity component");
   nc_addattr(ncid, varid,"units","m s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time z_r yv xv");

   //v
   if ((retval = nc_def_var(ncid,"vc",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval);   
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Northward water velocity component");
   nc_addattr(ncid, varid,"units","m s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time z_r yv xv");
   
   //w
   dimidthree[1] = dimid_Nkw;
   if ((retval = nc_def_var(ncid,"w",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval); 
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Vertical water velocity component");
   nc_addattr(ncid, varid,"units","m s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time z_w yv xv");  

   dimidthree[1] = dimid_Nk;
   
   //nu_v
   if ((retval = nc_def_var(ncid,"nu_v",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval); 
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Vertical eddy viscosity");
   nc_addattr(ncid, varid,"units","m2 s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time z_r yv xv");

   // turb_q
   if ((retval = nc_def_var(ncid,"turb_q",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval); 
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","turbulent velocity scale");
   nc_addattr(ncid, varid,"units","m s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time z_r yv xv");

   // turb_l
   if ((retval = nc_def_var(ncid,"turb_l",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval); 
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","turbulent length scale");
   nc_addattr(ncid, varid,"units","m");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time z_r yv xv");

   //salinity
   if(prop->beta>0){
     if ((retval = nc_def_var(ncid,"salt",NC_DOUBLE,3,dimidthree,&varid)))
      ERR(retval);
     if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
     if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Salinity");
    nc_addattr(ncid, varid,"units","ppt");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time z_r yv xv");
   }
   
   //temperature
   if(prop->gamma>0){
     if ((retval = nc_def_var(ncid,"temp",NC_DOUBLE,3,dimidthree,&varid)))
      ERR(retval); 
     if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
     if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Water temperature");
    nc_addattr(ncid, varid,"units","degrees C");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time z_r yv xv");
   }
   
   //rho
   if( (prop->gamma>0) || (prop->beta>0) ){
     if ((retval = nc_def_var(ncid,"rho",NC_DOUBLE,3,dimidthree,&varid)))
      ERR(retval);
     if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
     if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Water density");
    nc_addattr(ncid, varid,"units","kg m-3");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time z_r yv xv");
   }
   
   //age
   if(prop->calcage>0){
     if ((retval = nc_def_var(ncid,"agec",NC_DOUBLE,3,dimidthree,&varid)))
       ERR(retval);
     if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
       ERR(retval);
     if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
       ERR(retval);
     nc_addattr(ncid, varid,"long_name","Age concentration");
     nc_addattr(ncid, varid,"units","");
     nc_addattr(ncid, varid,"mesh","suntans_mesh");
     nc_addattr(ncid, varid,"location","face");
     nc_addattr(ncid, varid,"coordinates","time z_r yv xv");
     
     if ((retval = nc_def_var(ncid,"agealpha",NC_DOUBLE,3,dimidthree,&varid)))
       ERR(retval);
     if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
       ERR(retval);
     if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
       ERR(retval);
     nc_addattr(ncid, varid,"long_name","Age alpha parameter");
     nc_addattr(ncid, varid,"units","seconds");
     nc_addattr(ncid, varid,"mesh","suntans_mesh");
     nc_addattr(ncid, varid,"location","face");
     nc_addattr(ncid, varid,"coordinates","time z_r yv xv");
   }

   { // Bed stress - evaluated at edges, interpolated to center
     if ((retval = nc_def_var(ncid,"tauB_x",NC_DOUBLE,2,dimidtwo,&varid)))
       ERR(retval); 
     nc_addattr(ncid, varid,"long_name","Bed stress x-component, cell-centered from u_edge");
     nc_addattr(ncid, varid,"units","m2 s-2");
     nc_addattr(ncid, varid,"mesh","suntans_mesh");
     nc_addattr(ncid, varid,"location","face");
     nc_addattr(ncid, varid,"coordinates","time yv xv");
     
     if ((retval = nc_def_var(ncid,"tauB_y",NC_DOUBLE,2,dimidtwo,&varid)))
       ERR(retval); 
     nc_addattr(ncid, varid,"long_name","Bed stress y-component, cell-centered from u_edge");
     nc_addattr(ncid, varid,"units","m2 s-2");
     nc_addattr(ncid, varid,"mesh","suntans_mesh");
     nc_addattr(ncid, varid,"location","face");
     nc_addattr(ncid, varid,"coordinates","time yv xv");   
     
     if ((retval = nc_def_var(ncid,"tauB",NC_DOUBLE,2,(int[]){dimid_time,dimid_Ne},&varid)))
       ERR(retval); 
     nc_addattr(ncid, varid,"long_name","Bed stress edge-centered");
     nc_addattr(ncid, varid,"units","m2 s-2");
     nc_addattr(ncid, varid,"mesh","suntans_mesh");
     nc_addattr(ncid, varid,"location","edge");
     nc_addattr(ncid, varid,"coordinates","time ye xe");   
   }
   
   { // Bed stress evaluated at cell center
     if ((retval = nc_def_var(ncid,"tauBc_x",NC_DOUBLE,2,dimidtwo,&varid)))
       ERR(retval); 
     nc_addattr(ncid, varid,"long_name","Bed stress x-component, cell-centered from uc,vc");
     nc_addattr(ncid, varid,"units","N m-2");
     nc_addattr(ncid, varid,"mesh","suntans_mesh");
     nc_addattr(ncid, varid,"location","face");
     nc_addattr(ncid, varid,"coordinates","time yv xv");
     
     if ((retval = nc_def_var(ncid,"tauBc_y",NC_DOUBLE,2,dimidtwo,&varid)))
       ERR(retval); 
     nc_addattr(ncid, varid,"long_name","Bed stress y-component, cell-centered from uc,vc");
     nc_addattr(ncid, varid,"units","N m-2");
     nc_addattr(ncid, varid,"mesh","suntans_mesh");
     nc_addattr(ncid, varid,"location","face");
     nc_addattr(ncid, varid,"coordinates","time yv xv");   
   }

   
   //U
   dimidthree[2] = dimid_Ne;
   if ((retval = nc_def_var(ncid,"U",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval); 
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Edge normal velocity");
   nc_addattr(ncid, varid,"units","m s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","edge");
   nc_addattr(ncid, varid,"coordinates","time z_r ye xe");  

   // Meteorological variables (2-D) //
   
   if(prop->metmodel>0){
    // Uwind
    if ((retval = nc_def_var(ncid,"Uwind",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Eastward wind velocity component");
    nc_addattr(ncid, varid,"units","m s-1");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");

    // Vwind
    if ((retval = nc_def_var(ncid,"Vwind",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Eastward wind velocity component");
    nc_addattr(ncid, varid,"units","m s-1");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");

    // Tair
    if ((retval = nc_def_var(ncid,"Tair",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval);
    nc_addattr(ncid, varid,"long_name","Air temperature");
    nc_addattr(ncid, varid,"units","degrees C");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");

    // Pair
    if ((retval = nc_def_var(ncid,"Pair",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Air pressure");
    nc_addattr(ncid, varid,"units","millibar");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");

    // rain
    if ((retval = nc_def_var(ncid,"rain",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Rain fall rate");
    nc_addattr(ncid, varid,"units","kg m2 s-1");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");

    //RH
    if ((retval = nc_def_var(ncid,"RH",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Relative humidity");
    nc_addattr(ncid, varid,"units","Percent (%)");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");

    //cloud
    if ((retval = nc_def_var(ncid,"cloud",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Cloud cover fraction");
    nc_addattr(ncid, varid,"units","");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");
    
    // Surface flux variables //
    // Hs
    if ((retval = nc_def_var(ncid,"Hs",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Sensible heat flux");
    nc_addattr(ncid, varid,"units","W m-2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");   
    nc_addattr(ncid, varid,"positive","down");   

    // Hl
    if ((retval = nc_def_var(ncid,"Hl",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
     nc_addattr(ncid, varid,"long_name","Latent heat flux");
    nc_addattr(ncid, varid,"units","W m-2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");   
    nc_addattr(ncid, varid,"positive","down");   

    // Hlw
    if ((retval = nc_def_var(ncid,"Hlw",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
     nc_addattr(ncid, varid,"long_name","Net longwave radiation flux");
    nc_addattr(ncid, varid,"units","W m-2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");      
    nc_addattr(ncid, varid,"positive","down");   

    // Hsw
    if ((retval = nc_def_var(ncid,"Hsw",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Net shortwave radiation flux");
    nc_addattr(ncid, varid,"units","W m-2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");   
    nc_addattr(ncid, varid,"positive","down");   

    // tau_x
    if ((retval = nc_def_var(ncid,"tau_x",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Eastward component surface wind stress");
    nc_addattr(ncid, varid,"units","N m-2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");   

    // tau_y
    if ((retval = nc_def_var(ncid,"tau_y",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
     nc_addattr(ncid, varid,"long_name","Northward component surface wind stress");
    nc_addattr(ncid, varid,"units","N m-2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");   

    //EP
    if(prop->beta > 0.0){
	if ((retval = nc_def_var(ncid,"EP",NC_DOUBLE,2,dimidtwo,&varid)))
	    ERR(retval); 
	nc_addattr(ncid, varid,"long_name","Evaporation minus precipiaton");
	nc_addattr(ncid, varid,"units","m s-1");
	nc_addattr(ncid, varid,"mesh","suntans_mesh");
	nc_addattr(ncid, varid,"location","face");
	nc_addattr(ncid, varid,"coordinates","time yv xv");   
    }

   }

   //End file definition mode
   if ((retval = nc_enddef(ncid)))
	ERR(retval);

   
   /**********************************************************
   *
   * Write data (needs to be done out of definition mode for classic model)
   *
   ****************************************************************/
   
   nc_write_intvar(ncid,"cells",mergedGrid,mergedGrid->cells,myproc);
   nc_write_intvar(ncid,"face",mergedGrid,mergedGrid->face,myproc);
   nc_write_int(ncid,"nfaces",mergedGrid->nfaces,myproc);
   nc_write_int(ncid,"edges",mergedGrid->edges,myproc);
   //nc_write_intvar(ncid,"neigh",grid,grid->neigh,myproc);
   nc_write_int(ncid,"grad",mergedGrid->grad,myproc);
   //nc_write_int(ncid,"gradf",grid->gradf,myproc);
   nc_write_int(ncid,"mark",mergedGrid->mark,myproc);
   //nc_write_int(ncid,"mnptr",grid->mnptr,myproc);
   //nc_write_int(ncid,"eptr",grid->eptr,myproc);
   //
   nc_write_double(ncid,"xv",mergedGrid->xv,myproc);
   nc_write_double(ncid,"yv",mergedGrid->yv,myproc);
   nc_write_double(ncid,"xe",mergedGrid->xe,myproc);
   nc_write_double(ncid,"ye",mergedGrid->ye,myproc);
   nc_write_double(ncid,"xp",grid->xp,myproc);
   nc_write_double(ncid,"yp",grid->yp,myproc);

   nc_write_intvar(ncid,"normal",mergedGrid,mergedGrid->normal,myproc);
   nc_write_double(ncid,"n1",mergedGrid->n1,myproc);
   nc_write_double(ncid,"n2",mergedGrid->n2,myproc);
   nc_write_double(ncid,"df",mergedGrid->df,myproc);
   nc_write_double(ncid,"dg",mergedGrid->dg,myproc);
   //nc_write_doublevar(ncid,"def",grid,grid->def,myproc);
   nc_write_double(ncid,"Ac",mergedGrid->Ac,myproc);

   nc_write_double(ncid,"dz",grid->dz,myproc);
   nc_write_double(ncid,"z_r",z_r,myproc);
   nc_write_double(ncid,"z_w",z_w,myproc);
   nc_write_int(ncid,"Nk",mergedGrid->Nk,myproc);
   nc_write_int(ncid,"Nke",mergedGrid->Nke,myproc);
   nc_write_double(ncid,"dv",mergedGrid->dv,myproc);


      // Free the temporary vectors
   //SunFree(z_r,grid->Nkmax*sizeof(REAL),"InitialiseOutputNCugrid");
   //SunFree(z_w,(grid->Nkmax+1)*sizeof(REAL),"InitialiseOutputNCugrid");
   //SunFree(tmpvar,grid->Nc*grid->Nkmax,"InitialiseOutputNC");
   if(VERBOSE>1 && myproc==0) printf("Done.\n");

}// End function


/*
* Function: InitialiseOutputNCugrid()
* ------------------------------------
*
* Initialises the output netcdf file/s
* Files conform to the UGRID-0.9 CF conventions
* 
* One file per processor
* The pointer to each file is stored in prop->outputNetcdfFileID
* 
*/
void InitialiseOutputNCugrid(propT *prop, gridT *grid, physT *phys, metT *met, int myproc){
   int ncid = prop->outputNetcdfFileID;
   int retval, k, n, j, i;
   int varid;
   int dimid_Nc, dimid_Ne , dimid_Np, dimid_time, dimid_numsides, dimid_Two, dimid_Nkw, dimid_Nk; 
   int dimidone[1];
   int dimidtwo[2];
   int dimidthree[3];
   int nofill=0;
   const size_t starttwo[] = {0,0};
   const size_t counttwo[] = {grid->Nkmax,grid->Nc};
   REAL *z_r;
   REAL *z_w;
   int *edges;
   const int DEFLATE=1;
   const int DEFLATELEVEL=2;
   const REAL FILLVALUE = (REAL)EMPTY;
   char varname[256];

   //REAL *tmpvar;
   // Need to write the 3-D arrays as vectors
   //tmpvar = (REAL *)SunMalloc(grid->Nc*grid->Nkmax*sizeof(REAL),"InitialiseOutputNC");  
   
   /* Initialise the depth arrays */
   z_r = (REAL *)SunMalloc((grid->Nkmax)*sizeof(REAL),"InitialiseOutputNCugrid");
   z_w = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"InitialiseOutputNCugrid");

   /* Initialize an edge array */
   edges = (int *)SunMalloc(2*(grid->Ne)*sizeof(int),"InitialiseOutputNCugrid");
   
   /**************
    *
    * Start writing...
    *
    **************/
   if(VERBOSE>1 && myproc==0) printf("Initialising output netcdf files...");
   
   // Set the netcdf time ctr to 0
   prop->nctimectr=0;
   
   /* Define the global attributes - this should be expanded to include model input parameters*/
   nc_addattr(ncid, NC_GLOBAL,"title","SUNTANS NetCDF output file");
   
   /********************************************************************** 
    *
    * Define the dimensions
    *
    **********************************************************************/
   if ((retval = nc_def_dim(ncid, "Nc", grid->Nc, &dimid_Nc)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "Np", grid->Np, &dimid_Np)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "Ne", grid->Ne, &dimid_Ne)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "Nk", grid->Nkmax, &dimid_Nk)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "Nkw", grid->Nkmax+1, &dimid_Nkw)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "numsides", grid->maxfaces, &dimid_numsides)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "Two", 2, &dimid_Two)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "time", NC_UNLIMITED, &dimid_time)))
	ERR(retval);
   
    /********************************************************************** 
    *
    * Define the grid topology variables and attributes
    *
    **********************************************************************/

    //suntans_mesh
    if ((retval = nc_def_var(ncid,"suntans_mesh",NC_INT,0,0,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"cf_role","mesh_topology");
    nc_addattr(ncid, varid,"long_name","Topology data of 2D unstructured mesh");
    nc_addattr(ncid, varid,"topology_dimension","2");
    nc_addattr(ncid, varid,"node_coordinates","xp yp");
    nc_addattr(ncid, varid,"face_node_connectivity","cells");
    nc_addattr(ncid, varid,"edge_node_connectivity","edges");
    nc_addattr(ncid, varid,"face_coordinates","xv yv");
    nc_addattr(ncid, varid,"edge_coordinates","xe ye");
    nc_addattr(ncid, varid,"face_edge_connectivity","face");
    nc_addattr(ncid, varid,"edge_face_connectivity","grad");

    // cells
    dimidtwo[0] = dimid_Nc;
    dimidtwo[1] = dimid_numsides;
    if ((retval = nc_def_var(ncid,"cells",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"cf_role","face_node_connectivity");
    nc_addattr(ncid, varid,"long_name","Maps every face to its corner nodes");
    //if ((retval = nc_put_var_int(ncid,varid, grid->cells)))
    //  ERR(retval);

    //face
    dimidtwo[0] = dimid_Nc;
    dimidtwo[1] = dimid_numsides;
    if ((retval = nc_def_var(ncid,"face",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"cf_role","face_edge_connectivity");
    nc_addattr(ncid, varid,"long_name","Maps every face to its edges");
    //if ((retval = nc_put_var_int(ncid,varid, grid->face)))
    //  ERR(retval);

    //nfaces
    dimidone[0] = dimid_Nc;
    if ((retval = nc_def_var(ncid,"nfaces",NC_INT,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Number of cell faces");
    nc_addattr(ncid, varid,"units","");
    nc_addattr(ncid, varid,"coordinates","xv yv");

    //edges
    dimidtwo[0] = dimid_Ne;
    dimidtwo[1] = dimid_Two;
    if ((retval = nc_def_var(ncid,"edges",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"cf_role","edge_node_connectivity");
    nc_addattr(ncid, varid,"long_name","Maps every edge to the two nodes it connects");

    //neigh
    dimidtwo[0] = dimid_Nc;
    dimidtwo[1] = dimid_numsides;
    if ((retval = nc_def_var(ncid,"neigh",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"cf_role","face_face_connectivity");
    nc_addattr(ncid, varid,"long_name","Maps every face to its neighbouring faces");
    //if ((retval = nc_put_var_int(ncid,varid, grid->neigh)))
    //  ERR(retval);

    //grad
    dimidtwo[0] = dimid_Ne;
    dimidtwo[1] = dimid_Two;
    if ((retval = nc_def_var(ncid,"grad",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"cf_role","edge_face_connectivity");
    nc_addattr(ncid, varid,"long_name","Maps every edge to the two faces it connects ");
    //if ((retval = nc_put_var_int(ncid,varid, grid->grad)))
    //  ERR(retval);

    //gradf
    dimidtwo[0] = dimid_Ne;
    dimidtwo[1] = dimid_Two;
    if ((retval = nc_def_var(ncid,"gradf",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Index of face number for a particular cell ");
    //if ((retval = nc_put_var_int(ncid,varid, grid->grad)))
    //  ERR(retval);

    //mark
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"mark",NC_INT,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Edge marker type");
    nc_addattr(ncid, varid,"units","0 - computational; 1 - closed; 2 flux BC; 3 - stage BC; 4 - other BC; 5 - interproc; 6 - ghost cell.");
    nc_addattr(ncid, varid,"coordinates","xe ye");


    //mnptr
    dimidone[0] = dimid_Nc;
    if ((retval = nc_def_var(ncid,"mnptr",NC_INT,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Maps face indices between partitioned and unpartioned grid");
    //if ((retval = nc_put_var_int(ncid,varid, grid->mnptr)))
    //  ERR(retval);

    //eptr
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"eptr",NC_INT,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Maps edge indices between partitioned and unpartioned grid");
    //if ((retval = nc_put_var_int(ncid,varid, grid->eptr)))
    //  ERR(retval);

   /********************************************************************** 
    *
    * Define the grid coordinate variables and attributes 
    *
    **********************************************************************/

    //xv
    dimidone[0] = dimid_Nc;
    if ((retval = nc_def_var(ncid,"xv",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"standard_name","Easting");
    nc_addattr(ncid, varid,"long_name","Easting of 2D mesh face");
    //if ((retval = nc_put_var_double(ncid,varid, grid->xv)))
    //   ERR(retval);
      
    //yv
    dimidone[0] = dimid_Nc;
    if ((retval = nc_def_var(ncid,"yv",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"standard_name","Northing");
    nc_addattr(ncid, varid,"long_name","Northing of 2D mesh face");
    //if ((retval = nc_put_var_double(ncid,varid, grid->yv)))
    //   ERR(retval);
       
    //xp
    dimidone[0] = dimid_Np;
    if ((retval = nc_def_var(ncid,"xp",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"standard_name","Easting");
    nc_addattr(ncid, varid,"long_name","Easting of 2D mesh node");
    //if ((retval = nc_put_var_double(ncid,varid, grid->xp)))
    //   ERR(retval);
        
    //yp
    dimidone[0] = dimid_Np;
    if ((retval = nc_def_var(ncid,"yp",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"standard_name","Northing");
    nc_addattr(ncid, varid,"long_name","Northing of 2D mesh node");
    //if ((retval = nc_put_var_double(ncid,varid, grid->yp)))
    //   ERR(retval);
         
    //xe
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"xe",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"standard_name","Easting");
    nc_addattr(ncid, varid,"long_name","Easting of 2D mesh edge");
    //if ((retval = nc_put_var_double(ncid,varid, grid->xe)))
    //   ERR(retval);
          
    //ye
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"ye",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"standard_name","Northing");
    nc_addattr(ncid, varid,"long_name","Northing of 2D mesh edge");
    //if ((retval = nc_put_var_double(ncid,varid, grid->ye)))
    //   ERR(retval);
        
    /********************************************************************** 
    *
    * Define the grid metric variables 
    *
    **********************************************************************/

    //normal
    dimidtwo[0] = dimid_Nc;
    dimidtwo[1] = dimid_numsides;
    if ((retval = nc_def_var(ncid,"normal",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Dot product of unique normal with outward normal of each edge");
    //if ((retval = nc_put_var_int(ncid,varid, grid->normal)))
    //  ERR(retval);

    //n1
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"n1",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","x-component of the edge normal");
    nc_addattr(ncid, varid,"coordinates","xe ye");
    //if ((retval = nc_put_var_double(ncid,varid, grid->n1)))
    //  ERR(retval);

    //n2
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"n2",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","y-component of the edge normal");
    nc_addattr(ncid, varid,"coordinates","xe ye");
    //if ((retval = nc_put_var_double(ncid,varid, grid->n2)))
    //  ERR(retval);

    //df
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"df",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","edge length");
    nc_addattr(ncid, varid,"units","m");
    nc_addattr(ncid, varid,"coordinates","xe ye");
    //if ((retval = nc_put_var_double(ncid,varid, grid->df)))
    //  ERR(retval);

    //dg
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"dg",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","distance between faces on either side of edge");
    nc_addattr(ncid, varid,"units","m");
    nc_addattr(ncid, varid,"coordinates","xe ye");
    //if ((retval = nc_put_var_double(ncid,varid, grid->dg)))
    //  ERR(retval);

    //def
    dimidtwo[0] = dimid_Nc;
    dimidtwo[1] = dimid_numsides;
    if ((retval = nc_def_var(ncid,"def",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Distance between faces and edges");
    nc_addattr(ncid, varid,"units","m");
    nc_addattr(ncid, varid,"coordinates","xe ye");
    //if ((retval = nc_put_var_double(ncid,varid, grid->def)))
    //  ERR(retval);

    //Ac
    dimidone[0] = dimid_Nc;
    if ((retval = nc_def_var(ncid,"Ac",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Horizontal area of 2D mesh");
    nc_addattr(ncid, varid,"units","m2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","xv yv");
    //if ((retval = nc_put_var_double(ncid,varid, grid->Ac)))
    //  ERR(retval);

    /********************************************************************** 
    *
    * Define the vertical grid variables and attributes 
    *
    **********************************************************************/
   //dz
   dimidone[0] = dimid_Nk;   
   if ((retval = nc_def_var(ncid,"dz",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","z layer spacing");
   nc_addattr(ncid, varid,"units","m");
   //if ((retval = nc_put_var_double(ncid,varid, grid->dz)))
   //  ERR(retval);

   // Calculate and write the vertical coordinate z levels 
   z_w[0]=0.0;
   for(k=0;k<grid->Nkmax;k++){
      z_w[k+1] = z_w[k] + grid->dz[k];
      if(k==0){
	 z_r[k] = grid->dz[k]*0.5;
      }else{
	 z_r[k] = z_r[k-1]+grid->dz[k];
      }
   }
   //z_r
   dimidone[0] = dimid_Nk;   
   if ((retval = nc_def_var(ncid,"z_r",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"standard_name","ocean_z_coordinate");
   nc_addattr(ncid, varid,"long_name","depth at layer mid points");
   nc_addattr(ncid, varid,"units","m");
   nc_addattr(ncid, varid,"positive","down");
   //if ((retval = nc_put_var_double(ncid,varid, z_r)))
   //  ERR(retval);

   //z_w
   dimidone[0] = dimid_Nkw;
   if ((retval = nc_def_var(ncid,"z_w",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"standard_name","ocean_z_coordinate");
   nc_addattr(ncid, varid,"long_name","depth at layer edges");
   nc_addattr(ncid, varid,"units","m");
   nc_addattr(ncid, varid,"positive","down");
   //if ((retval = nc_put_var_double(ncid,varid, z_w)))
   //  ERR(retval);

   //Nk
   dimidone[0] = dimid_Nc;   
   if ((retval = nc_def_var(ncid,"Nk",NC_INT,1,dimidone,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Number of layers at face"); 
   //if ((retval = nc_put_var_int(ncid,varid, grid->Nk)))
   //  ERR(retval);

   //Nke
   dimidone[0] = dimid_Ne;   
   if ((retval = nc_def_var(ncid,"Nke",NC_INT,1,dimidone,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Number of layers at edge");
   //if ((retval = nc_put_var_int(ncid,varid, grid->Nke)))
   //  ERR(retval);

    //dv
    dimidone[0] = dimid_Nc;
    if ((retval = nc_def_var(ncid,"dv",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"standard_name","sea_floor_depth_below_geoid");
    nc_addattr(ncid, varid,"long_name","seafloor depth");
    nc_addattr(ncid, varid,"units","m");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","xv yv");
    nc_addattr(ncid, varid,"positive","down");
    //if ((retval = nc_put_var_double(ncid,varid, grid->dv)))
    //  ERR(retval);

    //dzz
    dimidthree[0] = dimid_time;
    dimidthree[1] = dimid_Nk;
    dimidthree[2] = dimid_Nc;
    if ((retval = nc_def_var(ncid,"dzz",NC_DOUBLE,3,dimidthree,&varid)))
       ERR(retval);
    if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
       ERR(retval);
    if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
       ERR(retval);
    nc_addattr(ncid, varid,"long_name","z layer spacing at faces");
    nc_addattr(ncid, varid,"units","m");

    //dzf
    dimidthree[0] = dimid_time;
    dimidthree[1] = dimid_Nk;
    dimidthree[2] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"dzf",NC_DOUBLE,3,dimidthree,&varid)))
       ERR(retval);
    if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
       ERR(retval);
    if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
       ERR(retval);
    nc_addattr(ncid, varid,"long_name","z layer spacing at edges");
    nc_addattr(ncid, varid,"units","m");

#ifdef NC_OUT_ETOP
    // etop
    dimidtwo[0]=dimid_time;
    dimidtwo[1]=dimid_Ne;
    if ((retval = nc_def_var(ncid,"etop",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","top active layer at edges");
#endif // NC_OUT_ETOP

    //time
    dimidone[0] = dimid_time;
    if ((retval = nc_def_var(ncid,"time",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"long_name","time");
    nc_addattr(ncid, varid,"units","seconds since 1990-01-01 00:00:00");

   /**********************************************************************
    *
    * Define the physical variables and attributes
    *
    **********************************************************************/

   dimidtwo[0] = dimid_time;
   dimidtwo[1] = dimid_Nc;

   dimidthree[0] = dimid_time;
   dimidthree[1] = dimid_Nk;
   dimidthree[2] = dimid_Nc;

   // eta
   if ((retval = nc_def_var(ncid,"eta",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Sea surface elevation");
   nc_addattr(ncid, varid,"units","m");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time yv xv");

   // cfl_count
   if ((retval = nc_def_var(ncid,"limiting_cell",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Cell frequency of CFL limit");
   nc_addattr(ncid, varid,"units","-");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time yv xv");

   // minimum time step per cell
   if ((retval = nc_def_var(ncid,"min_time_step",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Minimum time step limitation");
   nc_addattr(ncid, varid,"units","s");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time yv xv");

   //u
   if ((retval = nc_def_var(ncid,"uc",NC_DOUBLE,3,dimidthree,&varid)))
      ERR(retval);
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Eastward water velocity component");
   nc_addattr(ncid, varid,"units","m s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time z_r yv xv");

   //v
   if ((retval = nc_def_var(ncid,"vc",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval);
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Northward water velocity component");
   nc_addattr(ncid, varid,"units","m s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time z_r yv xv");

   //w
   dimidthree[1] = dimid_Nkw;
   if ((retval = nc_def_var(ncid,"w",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval); 
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Vertical water velocity component");
   nc_addattr(ncid, varid,"units","m s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time z_w yv xv");

   dimidthree[1] = dimid_Nk;

   //nu_v
   if ((retval = nc_def_var(ncid,"nu_v",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval); 
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Vertical eddy viscosity");
   nc_addattr(ncid, varid,"units","m2 s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time z_r yv xv");

   // turb_q
   if ((retval = nc_def_var(ncid,"turb_q",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval); 
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","turbulent velocity scale"); 
   nc_addattr(ncid, varid,"units","m s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time z_r yv xv");

   //turb_l
   if ((retval = nc_def_var(ncid,"turb_l",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval); 
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","turbulent length scale");
   nc_addattr(ncid, varid,"units","m");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time z_r yv xv");

   //salinity
   if(prop->beta>0){
     if ((retval = nc_def_var(ncid,"salt",NC_DOUBLE,3,dimidthree,&varid)))
      ERR(retval);
     if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
     if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Salinity");
    nc_addattr(ncid, varid,"units","ppt");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time z_r yv xv");
   }

   //temperature
   if(prop->gamma>0){
     if ((retval = nc_def_var(ncid,"temp",NC_DOUBLE,3,dimidthree,&varid)))
      ERR(retval); 
     if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
     if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Water temperature");
    nc_addattr(ncid, varid,"units","degrees C");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time z_r yv xv");
   }

   //rho
   if( (prop->gamma>0) || (prop->beta>0) ){
     if ((retval = nc_def_var(ncid,"rho",NC_DOUBLE,3,dimidthree,&varid)))
      ERR(retval);
     if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
     if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Water density");
    nc_addattr(ncid, varid,"units","kg m-3");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time z_r yv xv");
   }

   //age
   if(prop->calcage>0){
     if ((retval = nc_def_var(ncid,"agec",NC_DOUBLE,3,dimidthree,&varid)))
      ERR(retval);
     if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
     if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Age concentration");
    nc_addattr(ncid, varid,"units","");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time z_r yv xv");

    if ((retval = nc_def_var(ncid,"agealpha",NC_DOUBLE,3,dimidthree,&varid)))
      ERR(retval);
    if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
    if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Age alpha parameter");
    nc_addattr(ncid, varid,"units","seconds");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time z_r yv xv");

    /*
    //Age source term
    dimidtwo[0] = dimid_Nk;
    dimidtwo[1] = dimid_Nc;
    if ((retval = nc_def_var(ncid,"agesource",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval);
    if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
    if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Age source term (>0 =source");
    nc_addattr(ncid, varid,"units","");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","z_r yv xv");
    // Set back to time for the other variables
    dimidtwo[0] = dimid_time;
    */
   }

   //sediment
   if(prop->computeSediments){
     for(k=0;k<sediments->Nsize;k++){
       sprintf(varname,"sed%d",k+1);
     
       if ((retval = nc_def_var(ncid,varname,NC_DOUBLE,3,dimidthree,&varid)))
         ERRM(retval,"Define sed var");
       if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
         ERR(retval);
       if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
         ERR(retval);
       nc_addattr(ncid, varid,"long_name","Sediment concentration");
       nc_addattr(ncid, varid,"units","mg l-1"); // is that true?
       nc_addattr(ncid, varid,"mesh","suntans_mesh");
       nc_addattr(ncid, varid,"location","face");
       nc_addattr(ncid, varid,"coordinates","time z_r yv xv");
     }
   }

   //U
   dimidthree[2] = dimid_Ne;
   if ((retval = nc_def_var(ncid,"U",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval); 
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Edge normal velocity");
   nc_addattr(ncid, varid,"units","m s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","edge");
   nc_addattr(ncid, varid,"coordinates","time z_r ye xe");  

   // Meteorological variables (2-D) //
   
   if(prop->metmodel>0){
    // Uwind
    if ((retval = nc_def_var(ncid,"Uwind",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Eastward wind velocity component");
    nc_addattr(ncid, varid,"units","m s-1");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");

    // Vwind
    if ((retval = nc_def_var(ncid,"Vwind",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Eastward wind velocity component");
    nc_addattr(ncid, varid,"units","m s-1");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");

    // Tair
    if ((retval = nc_def_var(ncid,"Tair",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval);
    nc_addattr(ncid, varid,"long_name","Air temperature");
    nc_addattr(ncid, varid,"units","degrees C");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");

    // Pair
    if ((retval = nc_def_var(ncid,"Pair",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Air pressure");
    nc_addattr(ncid, varid,"units","millibar");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");

    // rain
    if ((retval = nc_def_var(ncid,"rain",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Rain fall rate");
    nc_addattr(ncid, varid,"units","kg m2 s-1");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");

    //RH
    if ((retval = nc_def_var(ncid,"RH",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Relative humidity");
    nc_addattr(ncid, varid,"units","Percent (%)");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");

    //cloud
    if ((retval = nc_def_var(ncid,"cloud",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Cloud cover fraction");
    nc_addattr(ncid, varid,"units","");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");
    
    // Surface flux variables //
    // Hs
    if ((retval = nc_def_var(ncid,"Hs",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Sensible heat flux");
    nc_addattr(ncid, varid,"units","W m-2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");   
    nc_addattr(ncid, varid,"positive","down");   

    // Hl
    if ((retval = nc_def_var(ncid,"Hl",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
     nc_addattr(ncid, varid,"long_name","Latent heat flux");
    nc_addattr(ncid, varid,"units","W m-2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");   
    nc_addattr(ncid, varid,"positive","down");   

    // Hlw
    if ((retval = nc_def_var(ncid,"Hlw",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
     nc_addattr(ncid, varid,"long_name","Net longwave radiation flux");
    nc_addattr(ncid, varid,"units","W m-2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");      
    nc_addattr(ncid, varid,"positive","down");   

    // Hsw
    if ((retval = nc_def_var(ncid,"Hsw",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Net shortwave radiation flux");
    nc_addattr(ncid, varid,"units","W m-2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");   
    nc_addattr(ncid, varid,"positive","down");   

    // tau_x
    if ((retval = nc_def_var(ncid,"tau_x",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Eastward component surface wind stress");
    nc_addattr(ncid, varid,"units","N m-2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");   

    // tau_y
    if ((retval = nc_def_var(ncid,"tau_y",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
     nc_addattr(ncid, varid,"long_name","Northward component surface wind stress");
    nc_addattr(ncid, varid,"units","N m-2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");   

    //EP
    if(prop->beta > 0.0){
	if ((retval = nc_def_var(ncid,"EP",NC_DOUBLE,2,dimidtwo,&varid)))
	    ERR(retval); 
	nc_addattr(ncid, varid,"long_name","Evaporation minus precipiaton");
	nc_addattr(ncid, varid,"units","m s-1");
	nc_addattr(ncid, varid,"mesh","suntans_mesh");
	nc_addattr(ncid, varid,"location","face");
	nc_addattr(ncid, varid,"coordinates","time yv xv");   
    }

   }

   //End file definition mode
   if ((retval = nc_enddef(ncid)))
     ERR(retval);

   
   /**********************************************************
   *
   * Write data (needs to be done out of definition mode for classic model)
   *
   ****************************************************************/
   nc_write_intvar(ncid,"cells",grid,grid->cells,myproc);
   nc_write_intvar(ncid,"face",grid,grid->face,myproc);
   nc_write_int(ncid,"nfaces",grid->nfaces,myproc);
   //nc_write_int(ncid,"edges",grid->edges,myproc);
   nc_write_intvar(ncid,"neigh",grid,grid->neigh,myproc);
   nc_write_int(ncid,"grad",grid->grad,myproc);
   nc_write_int(ncid,"gradf",grid->gradf,myproc);
   nc_write_int(ncid,"mark",grid->mark,myproc);
   nc_write_int(ncid,"mnptr",grid->mnptr,myproc);
   nc_write_int(ncid,"eptr",grid->eptr,myproc);
   
   nc_write_double(ncid,"xv",grid->xv,myproc);
   nc_write_double(ncid,"yv",grid->yv,myproc);
   nc_write_double(ncid,"xe",grid->xe,myproc);
   nc_write_double(ncid,"ye",grid->ye,myproc);
   nc_write_double(ncid,"xp",grid->xp,myproc);
   nc_write_double(ncid,"yp",grid->yp,myproc);

   nc_write_intvar(ncid,"normal",grid,grid->normal,myproc);
   nc_write_double(ncid,"n1",grid->n1,myproc);
   nc_write_double(ncid,"n2",grid->n2,myproc);
   nc_write_double(ncid,"df",grid->df,myproc);
   nc_write_double(ncid,"dg",grid->dg,myproc);
   nc_write_doublevar(ncid,"def",grid,grid->def,myproc);
   nc_write_double(ncid,"Ac",grid->Ac,myproc);

   nc_write_double(ncid,"dz",grid->dz,myproc);
   nc_write_double(ncid,"z_r",z_r,myproc);
   nc_write_double(ncid,"z_w",z_w,myproc);
   nc_write_int(ncid,"Nk",grid->Nk,myproc);
   nc_write_int(ncid,"Nke",grid->Nke,myproc);
   nc_write_double(ncid,"dv",grid->dv,myproc);


   // Need to convert the edge array that is stored is NUMEDGECOLUMN*Ne where
   // NUMEDGECOLUMN=3
   for(n=0;n<grid->Ne;n++){
     for(j=0;j<NUMEDGECOLUMNS-1;j++){
       edges[2*n+j] = grid->edges[NUMEDGECOLUMNS*n+j];
     }
   }
   nc_write_int(ncid,"edges",edges,myproc);

   // Free the temporary vectors
   //SunFree(z_r,grid->Nkmax*sizeof(REAL),"InitialiseOutputNCugrid");
   //SunFree(z_w,(grid->Nkmax+1)*sizeof(REAL),"InitialiseOutputNCugrid");
   //SunFree(tmpvar,grid->Nc*grid->Nkmax,"InitialiseOutputNC");
   if(VERBOSE>1 && myproc==0) printf("Done.\n");

   // RH -- this might make data more available during the simulation
   nc_sync(ncid);
}// End function

/*
* Function: InitialiseAverageNCugridMerge()
* ------------------------------------
*
* Initialises the average netcdf file/s
* Files conform to the UGRID-0.9 CF conventions
* 
* One file per processor
* The pointer to each file is stored in prop->outputNetcdfFileID
* 
*/
void InitialiseAverageNCugridMerge(propT *prop, gridT *grid, averageT *average, int myproc){
   int ncid = prop->averageNetcdfFileID;
   int retval, k;
   int varid;
   int dimid_Nc, dimid_Ne , dimid_Np, dimid_time, dimid_numsides, dimid_Two, dimid_Nkw, dimid_Nk; 
   int dimidone[1];
   int dimidtwo[2];
   int dimidthree[3];
   int nofill=0;
   const size_t starttwo[] = {0,0};
   const size_t counttwo[] = {mergedGrid->Nkmax,mergedGrid->Nc};
   REAL *z_r;
   REAL *z_w;
   const int DEFLATE=1;
   const int DEFLATELEVEL=2;
   const REAL FILLVALUE = (REAL)EMPTY;


   //REAL *tmpvar;
   // Need to write the 3-D arrays as vectors
   //tmpvar = (REAL *)SunMalloc(grid->Nc*grid->Nkmax*sizeof(REAL),"InitialiseOutputNC");  
   
   /* Initialise the depth arrays */
   z_r = (REAL *)SunMalloc((grid->Nkmax)*sizeof(REAL),"InitialiseOutputNCugrid");
   z_w = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"InitialiseOutputNCugrid");
   
   /**************
    *
    * Start writing...
    *
    **************/
   if(VERBOSE>1 && myproc==0) printf("Initialising output netcdf files...\n");
   
   // Set the netcdf time ctr to 0
   prop->avgtimectr=0;
   
   /* Define the global attributes - this should be expanded to include model input parameters*/
   nc_addattr(ncid, NC_GLOBAL,"title","SUNTANS NetCDF time-averaged file");

   nc_addattr_int(ncid,NC_GLOBAL,"ntaverage",&prop->ntaverage);
   nc_addattr_real(ncid,NC_GLOBAL,"dt",&prop->dt);
   
   /********************************************************************** 
    *
    * Define the dimensions
    *
    **********************************************************************/
   if ((retval = nc_def_dim(ncid, "Nc", mergedGrid->Nc, &dimid_Nc)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "Np", grid->Np, &dimid_Np)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "Ne", mergedGrid->Ne, &dimid_Ne)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "Nk", grid->Nkmax, &dimid_Nk)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "Nkw", grid->Nkmax+1, &dimid_Nkw)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "numsides", grid->maxfaces, &dimid_numsides)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "Two", 2, &dimid_Two)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "time", NC_UNLIMITED, &dimid_time)))
	ERR(retval);
   
    /********************************************************************** 
    *
    * Define the grid topology variables and attributes
    *
    **********************************************************************/

    //suntans_mesh
    if ((retval = nc_def_var(ncid,"suntans_mesh",NC_INT,0,0,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"cf_role","mesh_topology");
    nc_addattr(ncid, varid,"long_name","Topology data of 2D unstructured mesh");
    nc_addattr(ncid, varid,"topology_dimension","2");
    nc_addattr(ncid, varid,"node_coordinates","xp yp");
    nc_addattr(ncid, varid,"face_node_connectivity","cells");
    nc_addattr(ncid, varid,"edge_node_connectivity","edges");
    nc_addattr(ncid, varid,"face_coordinates","xv yv");
    nc_addattr(ncid, varid,"edge_coordinates","xe ye");
    nc_addattr(ncid, varid,"face_edge_connectivity","face");
    nc_addattr(ncid, varid,"edge_face_connectivity","grad");

    // cells
    dimidtwo[0] = dimid_Nc;
    dimidtwo[1] = dimid_numsides;
    if ((retval = nc_def_var(ncid,"cells",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"cf_role","face_node_connectivity");
    nc_addattr(ncid, varid,"long_name","Maps every face to its corner nodes");
    //if ((retval = nc_put_var_int(ncid,varid, grid->cells)))
    //  ERR(retval);

    //nfaces
    dimidone[0] = dimid_Nc;
    if ((retval = nc_def_var(ncid,"nfaces",NC_INT,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Number of cell faces");

    //face
    dimidtwo[0] = dimid_Nc;
    dimidtwo[1] = dimid_numsides;
    if ((retval = nc_def_var(ncid,"face",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"cf_role","face_edge_connectivity");
    nc_addattr(ncid, varid,"long_name","Maps every face to its edges");
    //if ((retval = nc_put_var_int(ncid,varid, grid->face)))
    //  ERR(retval);

    //edges
    dimidtwo[0] = dimid_Ne;
    dimidtwo[1] = dimid_Two;
    if ((retval = nc_def_var(ncid,"edges",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"cf_role","edge_node_connectivity");
    nc_addattr(ncid, varid,"long_name","Maps every edge to the two nodes it connects");
    //if ((retval = nc_put_var_int(ncid,varid, grid->edges)))
    //  ERR(retval);

    //neigh
    dimidtwo[0] = dimid_Nc;
    dimidtwo[1] = dimid_numsides;
    if ((retval = nc_def_var(ncid,"neigh",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"cf_role","face_face_connectivity");
    nc_addattr(ncid, varid,"long_name","Maps every face to its neighbouring faces");
    //if ((retval = nc_put_var_int(ncid,varid, grid->neigh)))
    //  ERR(retval);

    //grad
    dimidtwo[0] = dimid_Ne;
    dimidtwo[1] = dimid_Two;
    if ((retval = nc_def_var(ncid,"grad",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"cf_role","edge_face_connectivity");
    nc_addattr(ncid, varid,"long_name","Maps every edge to the two faces it connects ");
    //if ((retval = nc_put_var_int(ncid,varid, grid->grad)))
    //  ERR(retval);
   /********************************************************************** 
    *
    * Define the grid coordinate variables and attributes 
    *
    **********************************************************************/

    //xv
    dimidone[0] = dimid_Nc;
    if ((retval = nc_def_var(ncid,"xv",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"standard_name","Easting");
    nc_addattr(ncid, varid,"long_name","Easting of 2D mesh face");
    //if ((retval = nc_put_var_double(ncid,varid, grid->xv)))
    //   ERR(retval);
      
    //yv
    dimidone[0] = dimid_Nc;
    if ((retval = nc_def_var(ncid,"yv",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"standard_name","Northing");
    nc_addattr(ncid, varid,"long_name","Northing of 2D mesh face");
    //if ((retval = nc_put_var_double(ncid,varid, grid->yv)))
    //   ERR(retval);
       
    //xp
    dimidone[0] = dimid_Np;
    if ((retval = nc_def_var(ncid,"xp",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"standard_name","Easting");
    nc_addattr(ncid, varid,"long_name","Easting of 2D mesh node");
    //if ((retval = nc_put_var_double(ncid,varid, grid->xp)))
    //   ERR(retval);
        
    //yp
    dimidone[0] = dimid_Np;
    if ((retval = nc_def_var(ncid,"yp",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"standard_name","Northing");
    nc_addattr(ncid, varid,"long_name","Northing of 2D mesh node");
    //if ((retval = nc_put_var_double(ncid,varid, grid->yp)))
    //   ERR(retval);
         
    //xe
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"xe",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"standard_name","Easting");
    nc_addattr(ncid, varid,"long_name","Easting of 2D mesh edge");
    //if ((retval = nc_put_var_double(ncid,varid, grid->xe)))
    //   ERR(retval);
          
    //ye
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"ye",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"standard_name","Northing");
    nc_addattr(ncid, varid,"long_name","Northing of 2D mesh edge");
    //if ((retval = nc_put_var_double(ncid,varid, grid->ye)))
    //   ERR(retval);
        
    /********************************************************************** 
    *
    * Define the grid metric variables 
    *
    **********************************************************************/

    //normal
    dimidtwo[0] = dimid_Nc;
    dimidtwo[1] = dimid_numsides;
    if ((retval = nc_def_var(ncid,"normal",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Dot product of unique normal with outward normal of each edge");
    //if ((retval = nc_put_var_int(ncid,varid, grid->normal)))
    //  ERR(retval);

    //n1
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"n1",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","x-component of the edge normal");
    nc_addattr(ncid, varid,"coordinates","xe ye");

    //n2
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"n2",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","y-component of the edge normal");
    nc_addattr(ncid, varid,"coordinates","xe ye");


    //df
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"df",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","edge length");
    nc_addattr(ncid, varid,"units","m");
    nc_addattr(ncid, varid,"coordinates","xe ye");
    //if ((retval = nc_put_var_double(ncid,varid, grid->df)))
    //  ERR(retval);

    //dg
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"dg",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","distance between faces on either side of edge");
    nc_addattr(ncid, varid,"units","m");
    nc_addattr(ncid, varid,"coordinates","xe ye");
    //if ((retval = nc_put_var_double(ncid,varid, grid->dg)))
    //  ERR(retval);

    //def
    dimidtwo[0] = dimid_Nc;
    dimidtwo[1] = dimid_numsides;
    if ((retval = nc_def_var(ncid,"def",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Distance between faces and edges");
    nc_addattr(ncid, varid,"units","m");
    nc_addattr(ncid, varid,"coordinates","xe ye");

    //mark
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"mark",NC_INT,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Edge marker type");
    nc_addattr(ncid, varid,"units","0 - computational; 1 - closed; 2 flux BC; 3 - stage BC; 4 - other BC; 5 - interproc; 6 - ghost cell.");
    nc_addattr(ncid, varid,"coordinates","xe ye");

    //Ac
    dimidone[0] = dimid_Nc;
    if ((retval = nc_def_var(ncid,"Ac",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Horizontal area of 2D mesh");
    nc_addattr(ncid, varid,"units","m2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","xv yv");
    //if ((retval = nc_put_var_double(ncid,varid, grid->Ac)))
    //  ERR(retval);

    /********************************************************************** 
    *
    * Define the vertical grid variables and attributes 
    *
    **********************************************************************/
   //dz
   dimidone[0] = dimid_Nk;   
   if ((retval = nc_def_var(ncid,"dz",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","z layer spacing");
   nc_addattr(ncid, varid,"units","m");
   //if ((retval = nc_put_var_double(ncid,varid, grid->dz)))
   //  ERR(retval);

   // Calculate and write the vertical coordinate z levels 
   z_w[0]=0.0;
   for(k=0;k<grid->Nkmax;k++){
      z_w[k+1] = z_w[k] + grid->dz[k];
      if(k==0){
	 z_r[k] = grid->dz[k]*0.5;
      }else{
	 z_r[k] = z_r[k-1]+grid->dz[k];
      }
   }
   //z_r
   dimidone[0] = dimid_Nk;   
   if ((retval = nc_def_var(ncid,"z_r",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"standard_name","ocean_z_coordinate");
   nc_addattr(ncid, varid,"long_name","depth at layer mid points");
   nc_addattr(ncid, varid,"units","m");  
   nc_addattr(ncid, varid,"positive","up");  
   //if ((retval = nc_put_var_double(ncid,varid, z_r)))
   //  ERR(retval);

   //z_w
   dimidone[0] = dimid_Nkw;   
   if ((retval = nc_def_var(ncid,"z_w",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"standard_name","ocean_z_coordinate");
   nc_addattr(ncid, varid,"long_name","depth at layer edges");
   nc_addattr(ncid, varid,"units","m");  
   nc_addattr(ncid, varid,"positive","up");  
   //if ((retval = nc_put_var_double(ncid,varid, z_w)))
   //  ERR(retval);

   //Nk
   dimidone[0] = dimid_Nc;   
   if ((retval = nc_def_var(ncid,"Nk",NC_INT,1,dimidone,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Number of layers at face"); 
   //if ((retval = nc_put_var_int(ncid,varid, grid->Nk)))
   //  ERR(retval);

   //Nke
   dimidone[0] = dimid_Ne;   
   if ((retval = nc_def_var(ncid,"Nke",NC_INT,1,dimidone,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Number of layers at edge");
   //if ((retval = nc_put_var_int(ncid,varid, grid->Nke)))
   //  ERR(retval);

    //dv
    dimidone[0] = dimid_Nc;
    if ((retval = nc_def_var(ncid,"dv",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"standard_name","sea_floor_depth_below_geoid");
    nc_addattr(ncid, varid,"long_name","seafloor depth");
    nc_addattr(ncid, varid,"units","m");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","xv yv");
    //if ((retval = nc_put_var_double(ncid,varid, grid->dv)))
    //  ERR(retval);

    //time
    dimidone[0] = dimid_time;
    if ((retval = nc_def_var(ncid,"time",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"long_name","time");
    nc_addattr(ncid, varid,"units","seconds since 1990-01-01 00:00:00");  

  //  //average time
  //  dimidone[0] = 0;
  //  if ((retval = nc_def_var(ncid,"average_time",NC_DOUBLE,1,dimidone,&varid)))
  //     ERR(retval);
  //  nc_addattr(ncid, varid,"long_name","Averaging time interval");
  //  nc_addattr(ncid, varid,"units","seconds");  


   /********************************************************************** 
    * 
    * Define the physical variables and attributes 
    *
    **********************************************************************/
   
   dimidtwo[0] = dimid_time;
   dimidtwo[1] = dimid_Nc;
   
   dimidthree[0] = dimid_time;
   dimidthree[1] = dimid_Nk;
   dimidthree[2] = dimid_Nc;
   
   // eta average
   if ((retval = nc_def_var(ncid,"eta_avg",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Time-averaged sea surface elevation");
   nc_addattr(ncid, varid,"units","m");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time yv xv");
    
   // eta
   if ((retval = nc_def_var(ncid,"eta",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Instantaneous sea surface elevation");
   nc_addattr(ncid, varid,"units","m");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time yv xv");
 
   //u
   if ((retval = nc_def_var(ncid,"uc",NC_DOUBLE,3,dimidthree,&varid)))
      ERR(retval);
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Time-averaged Eastward water velocity component");
   nc_addattr(ncid, varid,"units","m s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time z_r yv xv");

   //v
   if ((retval = nc_def_var(ncid,"vc",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval);   
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Time-averaged Northward water velocity component");
   nc_addattr(ncid, varid,"units","m s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time z_r yv xv");
   
   //w
   dimidthree[1] = dimid_Nkw;
   if ((retval = nc_def_var(ncid,"w",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval); 
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Time-averaged Vertical water velocity component");
   nc_addattr(ncid, varid,"units","m s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time z_w yv xv");  

   dimidthree[1] = dimid_Nk;
   
   //nu_v
   if ((retval = nc_def_var(ncid,"nu_v",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval); 
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Time-averaged Vertical eddy viscosity");
   nc_addattr(ncid, varid,"units","m2 s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time z_r yv xv");
   
   //kappa_tv
   if ((retval = nc_def_var(ncid,"kappa_tv",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval); 
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Time-averaged vertical tracer diffusivity");
   nc_addattr(ncid, varid,"units","m2 s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time z_r yv xv");

   //salinity
   if(prop->beta>0){
     if ((retval = nc_def_var(ncid,"salt",NC_DOUBLE,3,dimidthree,&varid)))
       ERR(retval);
     if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
       ERR(retval);
     if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
       ERR(retval);
     nc_addattr(ncid, varid,"long_name","Time-averaged Salinity");
     nc_addattr(ncid, varid,"units","ppt");
     nc_addattr(ncid, varid,"mesh","suntans_mesh");
     nc_addattr(ncid, varid,"location","face");
     nc_addattr(ncid, varid,"coordinates","time z_r yv xv");
     
     // Depth-integrated salinity
     if ((retval = nc_def_var(ncid,"s_dz",NC_DOUBLE,2,dimidtwo,&varid)))
       ERR(retval);
     nc_addattr(ncid, varid,"long_name","Instantaneous depth-integrated salinity");
     nc_addattr(ncid, varid,"units","psu m");
     nc_addattr(ncid, varid,"mesh","suntans_mesh");
     nc_addattr(ncid, varid,"location","face");
     nc_addattr(ncid, varid,"coordinates","time yv xv");
   }
   
   //temperature
   if(prop->gamma>0){
     if ((retval = nc_def_var(ncid,"temp",NC_DOUBLE,3,dimidthree,&varid)))
       ERR(retval); 
     if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
       ERR(retval);
     if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
       ERR(retval);
     nc_addattr(ncid, varid,"long_name","Time-averaged Water temperature");
     nc_addattr(ncid, varid,"units","degrees C");
     nc_addattr(ncid, varid,"mesh","suntans_mesh");
     nc_addattr(ncid, varid,"location","face");
     nc_addattr(ncid, varid,"coordinates","time z_r yv xv");
     
     // Depth-integrated temperature
     if ((retval = nc_def_var(ncid,"T_dz",NC_DOUBLE,2,dimidtwo,&varid)))
       ERR(retval);
     nc_addattr(ncid, varid,"long_name","Instantaneous depth-integrated temperature");
     nc_addattr(ncid, varid,"units","degrees C m");
     nc_addattr(ncid, varid,"mesh","suntans_mesh");
     nc_addattr(ncid, varid,"location","face");
     nc_addattr(ncid, varid,"coordinates","time yv xv");
   }
   
   //rho
   if( (prop->gamma>0) || (prop->beta>0) ){
     if ((retval = nc_def_var(ncid,"rho",NC_DOUBLE,3,dimidthree,&varid)))
       ERR(retval);
     if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
       ERR(retval);
     if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
       ERR(retval);
     nc_addattr(ncid, varid,"long_name","Time-averaged Water density");
     nc_addattr(ncid, varid,"units","kg m-3");
     nc_addattr(ncid, varid,"mesh","suntans_mesh");
     nc_addattr(ncid, varid,"location","face");
     nc_addattr(ncid, varid,"coordinates","time z_r yv xv");
   }
   
   //age
   if(prop->calcage>0){
     if ((retval = nc_def_var(ncid,"agec",NC_DOUBLE,3,dimidthree,&varid)))
       ERR(retval);
     if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
       ERR(retval);
     if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
       ERR(retval);
     nc_addattr(ncid, varid,"long_name","Age concentration");
     nc_addattr(ncid, varid,"units","");
     nc_addattr(ncid, varid,"mesh","suntans_mesh");
     nc_addattr(ncid, varid,"location","face");
     nc_addattr(ncid, varid,"coordinates","time z_r yv xv");
     
     if ((retval = nc_def_var(ncid,"agealpha",NC_DOUBLE,3,dimidthree,&varid)))
       ERR(retval);
     if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
       ERR(retval);
     if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
       ERR(retval);
     nc_addattr(ncid, varid,"long_name","Age alpha parameter");
     nc_addattr(ncid, varid,"units","seconds");
     nc_addattr(ncid, varid,"mesh","suntans_mesh");
     nc_addattr(ncid, varid,"location","face");
     nc_addattr(ncid, varid,"coordinates","time z_r yv xv");
   }

   //U_F
   dimidthree[2] = dimid_Ne;
   if ((retval = nc_def_var(ncid,"U_F",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval); 
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
     ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
     ERR(retval);
   nc_addattr(ncid, varid,"long_name","Time-averaged edge flux rate");
   nc_addattr(ncid, varid,"units","m3 s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","edge");
   nc_addattr(ncid, varid,"coordinates","time z_r ye xe");  

    //s_F
   dimidthree[2] = dimid_Ne;
   if ((retval = nc_def_var(ncid,"s_F",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval); 
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Time-averaged edge salt flux rate");
   nc_addattr(ncid, varid,"units","psu m3 s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","edge");
   nc_addattr(ncid, varid,"coordinates","time z_r ye xe");  

    //T_F
   dimidthree[2] = dimid_Ne;
   if ((retval = nc_def_var(ncid,"T_F",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval); 
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Time-averaged edge temperature flux rate");
   nc_addattr(ncid, varid,"units","degreesC m3 s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","edge");
   nc_addattr(ncid, varid,"coordinates","time z_r ye xe");  

   // Meteorological variables (2-D) //
   
   if(prop->metmodel>0){
    // Uwind
    if ((retval = nc_def_var(ncid,"Uwind",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Time-averaged Eastward wind velocity component");
    nc_addattr(ncid, varid,"units","m s-1");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");

    // Vwind
    if ((retval = nc_def_var(ncid,"Vwind",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Time-averaged Eastward wind velocity component");
    nc_addattr(ncid, varid,"units","m s-1");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");

    // Tair
    if ((retval = nc_def_var(ncid,"Tair",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval);
    nc_addattr(ncid, varid,"long_name","Time-averaged Air temperature");
    nc_addattr(ncid, varid,"units","degrees C");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");

    // Pair
    if ((retval = nc_def_var(ncid,"Pair",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Time-averaged Air pressure");
    nc_addattr(ncid, varid,"units","millibar");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");

    // rain
    if ((retval = nc_def_var(ncid,"rain",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Time-averaged Rain fall rate");
    nc_addattr(ncid, varid,"units","kg m2 s-1");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");

    //RH
    if ((retval = nc_def_var(ncid,"RH",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Time-averaged Relative humidity");
    nc_addattr(ncid, varid,"units","Percent (%)");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");

    //cloud
    if ((retval = nc_def_var(ncid,"cloud",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Time-averaged Cloud cover fraction");
    nc_addattr(ncid, varid,"units","");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");
    
    // Surface flux variables //
    // Hs
    if ((retval = nc_def_var(ncid,"Hs",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Time-averaged Sensible heat flux");
    nc_addattr(ncid, varid,"units","W m-2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");   
    nc_addattr(ncid, varid,"positive","down");   

    // Hl
    if ((retval = nc_def_var(ncid,"Hl",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
     nc_addattr(ncid, varid,"long_name","Time-averaged Latent heat flux");
    nc_addattr(ncid, varid,"units","W m-2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");   
    nc_addattr(ncid, varid,"positive","down");   

    // Hlw
    if ((retval = nc_def_var(ncid,"Hlw",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
     nc_addattr(ncid, varid,"long_name","Time-averaged Net longwave radiation flux");
    nc_addattr(ncid, varid,"units","W m-2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");      
    nc_addattr(ncid, varid,"positive","down");   

    // Hsw
    if ((retval = nc_def_var(ncid,"Hsw",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Time-averaged Net shortwave radiation flux");
    nc_addattr(ncid, varid,"units","W m-2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");   
    nc_addattr(ncid, varid,"positive","down");   

    // tau_x
    if ((retval = nc_def_var(ncid,"tau_x",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Time-averaged Eastward component surface wind stress");
    nc_addattr(ncid, varid,"units","N m-2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");   

    // tau_y
    if ((retval = nc_def_var(ncid,"tau_y",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
     nc_addattr(ncid, varid,"long_name","Time-averaged Northward component surface wind stress");
    nc_addattr(ncid, varid,"units","N m-2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");   

    //EP
    if(prop->beta > 0.0){
	if ((retval = nc_def_var(ncid,"EP",NC_DOUBLE,2,dimidtwo,&varid)))
	    ERR(retval); 
	nc_addattr(ncid, varid,"long_name","Time-averaged surface salt flux (S0*EP)");
	nc_addattr(ncid, varid,"units","psu m s-1");
	nc_addattr(ncid, varid,"mesh","suntans_mesh");
	nc_addattr(ncid, varid,"location","face");
	nc_addattr(ncid, varid,"coordinates","time yv xv");   
    }

   }

   //End file definition mode
   if ((retval = nc_enddef(ncid)))
	ERR(retval);

   
   /**********************************************************
   *
   * Write data (needs to be done out of definition mode for classic model)
   *
   ****************************************************************/
   nc_write_intvar(ncid,"cells",mergedGrid,mergedGrid->cells,myproc);
   nc_write_intvar(ncid,"face",mergedGrid,mergedGrid->face,myproc);
   nc_write_int(ncid,"nfaces",mergedGrid->nfaces,myproc);
   nc_write_int(ncid,"edges",mergedGrid->edges,myproc);
   //nc_write_int(ncid,"neigh",mergedGrid->neigh,myproc);
   nc_write_int(ncid,"grad",mergedGrid->grad,myproc);
   nc_write_int(ncid,"mark",mergedGrid->mark,myproc);
   //nc_write_int(ncid,"mnptr",grid->mnptr,myproc);
   //nc_write_int(ncid,"eptr",grid->eptr,myproc);
   
   nc_write_double(ncid,"xv",mergedGrid->xv,myproc);
   nc_write_double(ncid,"yv",mergedGrid->yv,myproc);
   nc_write_double(ncid,"xe",mergedGrid->xe,myproc);
   nc_write_double(ncid,"ye",mergedGrid->ye,myproc);
   nc_write_double(ncid,"xp",grid->xp,myproc);
   nc_write_double(ncid,"yp",grid->yp,myproc);

   nc_write_intvar(ncid,"normal",mergedGrid,mergedGrid->normal,myproc);
   nc_write_double(ncid,"n1",mergedGrid->n1,myproc);
   nc_write_double(ncid,"n2",mergedGrid->n2,myproc);
   nc_write_double(ncid,"df",mergedGrid->df,myproc);
   nc_write_double(ncid,"dg",mergedGrid->dg,myproc);
   nc_write_double(ncid,"def",mergedGrid->def,myproc);
   nc_write_double(ncid,"Ac",mergedGrid->Ac,myproc);

   nc_write_double(ncid,"dz",grid->dz,myproc);
   nc_write_double(ncid,"z_r",z_r,myproc);
   nc_write_double(ncid,"z_w",z_w,myproc);
   nc_write_int(ncid,"Nk",mergedGrid->Nk,myproc);
   nc_write_int(ncid,"Nke",mergedGrid->Nke,myproc);
   nc_write_double(ncid,"dv",mergedGrid->dv,myproc);

  // nc_write_double(ncid,"average_time",(float)prop->ntaverage*prop->dt,myproc);

}// End function


/*
* Function: InitialiseAverageNCugrid()
* ------------------------------------
*
* Initialises the average netcdf file/s
* Files conform to the UGRID-0.9 CF conventions
* 
* One file per processor
* The pointer to each file is stored in prop->outputNetcdfFileID
* 
*/
void InitialiseAverageNCugrid(propT *prop, gridT *grid, averageT *average, int myproc){
   int ncid = prop->averageNetcdfFileID;
   int retval, k, n, j;
   int varid;
   int dimid_Nc, dimid_Ne , dimid_Np, dimid_time, dimid_numsides, dimid_Two, dimid_Nkw, dimid_Nk; 
   int dimidone[1];
   int dimidtwo[2];
   int dimidthree[3];
   int nofill=0;
   const size_t starttwo[] = {0,0};
   const size_t counttwo[] = {grid->Nkmax,grid->Nc};
   REAL *z_r;
   REAL *z_w;
   int *edges;
   const int DEFLATE=1;
   const int DEFLATELEVEL=2;
   const REAL FILLVALUE = (REAL)EMPTY;

   /* Initialise the depth arrays */
   z_r = (REAL *)SunMalloc((grid->Nkmax)*sizeof(REAL),"InitialiseOutputNCugrid");
   z_w = (REAL *)SunMalloc((grid->Nkmax+1)*sizeof(REAL),"InitialiseOutputNCugrid");
   
   /**************
    *
    * Start writing...
    *
    **************/
   if(VERBOSE>1 && myproc==0) printf("Initialising output netcdf files...\n");
   
   // Set the netcdf time ctr to 0
   prop->avgtimectr=0;
   prop->avgctr=0;
   
   /* Define the global attributes - this should be expanded to include model input parameters*/
   nc_addattr(ncid, NC_GLOBAL,"title","SUNTANS NetCDF time-averaged file");

   nc_addattr_int(ncid,NC_GLOBAL,"ntaverage",&prop->ntaverage);
   nc_addattr_real(ncid,NC_GLOBAL,"dt",&prop->dt);
   
   /********************************************************************** 
    *
    * Define the dimensions
    *
    **********************************************************************/
   if ((retval = nc_def_dim(ncid, "Nc", grid->Nc, &dimid_Nc)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "Np", grid->Np, &dimid_Np)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "Ne", grid->Ne, &dimid_Ne)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "Nk", grid->Nkmax, &dimid_Nk)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "Nkw", grid->Nkmax+1, &dimid_Nkw)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "numsides", grid->maxfaces, &dimid_numsides)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "Two", 2, &dimid_Two)))
	ERR(retval);
   if ((retval = nc_def_dim(ncid, "time", NC_UNLIMITED, &dimid_time)))
	ERR(retval);
   
    /********************************************************************** 
    *
    * Define the grid topology variables and attributes
    *
    **********************************************************************/

    //suntans_mesh
    if ((retval = nc_def_var(ncid,"suntans_mesh",NC_INT,0,0,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"cf_role","mesh_topology");
    nc_addattr(ncid, varid,"long_name","Topology data of 2D unstructured mesh");
    nc_addattr(ncid, varid,"topology_dimension","2");
    nc_addattr(ncid, varid,"node_coordinates","xp yp");
    nc_addattr(ncid, varid,"face_node_connectivity","cells");
    nc_addattr(ncid, varid,"edge_node_connectivity","edges");
    nc_addattr(ncid, varid,"face_coordinates","xv yv");
    nc_addattr(ncid, varid,"edge_coordinates","xe ye");
    nc_addattr(ncid, varid,"face_edge_connectivity","face");
    nc_addattr(ncid, varid,"edge_face_connectivity","grad");

    // cells
    dimidtwo[0] = dimid_Nc;
    dimidtwo[1] = dimid_numsides;
    if ((retval = nc_def_var(ncid,"cells",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"cf_role","face_node_connectivity");
    nc_addattr(ncid, varid,"long_name","Maps every face to its corner nodes");
    //if ((retval = nc_put_var_int(ncid,varid, grid->cells)))
    //  ERR(retval);

    //face
    dimidtwo[0] = dimid_Nc;
    dimidtwo[1] = dimid_numsides;
    if ((retval = nc_def_var(ncid,"face",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"cf_role","face_edge_connectivity");
    nc_addattr(ncid, varid,"long_name","Maps every face to its edges");
    //if ((retval = nc_put_var_int(ncid,varid, grid->face)))
    //  ERR(retval);

    //nfaces
    dimidone[0] = dimid_Nc;
    if ((retval = nc_def_var(ncid,"nfaces",NC_INT,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Number of cell faces");

    //edges
    dimidtwo[0] = dimid_Ne;
    dimidtwo[1] = dimid_Two;
    if ((retval = nc_def_var(ncid,"edges",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"cf_role","edge_node_connectivity");
    nc_addattr(ncid, varid,"long_name","Maps every edge to the two nodes it connects");

    //neigh
    dimidtwo[0] = dimid_Nc;
    dimidtwo[1] = dimid_numsides;
    if ((retval = nc_def_var(ncid,"neigh",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"cf_role","face_face_connectivity");
    nc_addattr(ncid, varid,"long_name","Maps every face to its neighbouring faces");
    //if ((retval = nc_put_var_int(ncid,varid, grid->neigh)))
    //  ERR(retval);

    //grad
    dimidtwo[0] = dimid_Ne;
    dimidtwo[1] = dimid_Two;
    if ((retval = nc_def_var(ncid,"grad",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"cf_role","edge_face_connectivity");
    nc_addattr(ncid, varid,"long_name","Maps every edge to the two faces it connects ");
    //if ((retval = nc_put_var_int(ncid,varid, grid->grad)))
    //  ERR(retval);

    //mnptr
    dimidone[0] = dimid_Nc;
    if ((retval = nc_def_var(ncid,"mnptr",NC_INT,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Maps face indices between partitioned and unpartioned grid");
    //if ((retval = nc_put_var_int(ncid,varid, grid->mnptr)))
    //  ERR(retval);

    //eptr
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"eptr",NC_INT,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Maps edge indices between partitioned and unpartioned grid");
    //if ((retval = nc_put_var_int(ncid,varid, grid->eptr)))
    //  ERR(retval);

   /********************************************************************** 
    *
    * Define the grid coordinate variables and attributes 
    *
    **********************************************************************/

    //xv
    dimidone[0] = dimid_Nc;
    if ((retval = nc_def_var(ncid,"xv",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"standard_name","Easting");
    nc_addattr(ncid, varid,"long_name","Easting of 2D mesh face");
    //if ((retval = nc_put_var_double(ncid,varid, grid->xv)))
    //   ERR(retval);
      
    //yv
    dimidone[0] = dimid_Nc;
    if ((retval = nc_def_var(ncid,"yv",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"standard_name","Northing");
    nc_addattr(ncid, varid,"long_name","Northing of 2D mesh face");
    //if ((retval = nc_put_var_double(ncid,varid, grid->yv)))
    //   ERR(retval);
       
    //xp
    dimidone[0] = dimid_Np;
    if ((retval = nc_def_var(ncid,"xp",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"standard_name","Easting");
    nc_addattr(ncid, varid,"long_name","Easting of 2D mesh node");
    //if ((retval = nc_put_var_double(ncid,varid, grid->xp)))
    //   ERR(retval);
        
    //yp
    dimidone[0] = dimid_Np;
    if ((retval = nc_def_var(ncid,"yp",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"standard_name","Northing");
    nc_addattr(ncid, varid,"long_name","Northing of 2D mesh node");
    //if ((retval = nc_put_var_double(ncid,varid, grid->yp)))
    //   ERR(retval);
         
    //xe
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"xe",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"standard_name","Easting");
    nc_addattr(ncid, varid,"long_name","Easting of 2D mesh edge");
    //if ((retval = nc_put_var_double(ncid,varid, grid->xe)))
    //   ERR(retval);
          
    //ye
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"ye",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"standard_name","Northing");
    nc_addattr(ncid, varid,"long_name","Northing of 2D mesh edge");
    //if ((retval = nc_put_var_double(ncid,varid, grid->ye)))
    //   ERR(retval);
        
    /********************************************************************** 
    *
    * Define the grid metric variables 
    *
    **********************************************************************/

    //normal
    dimidtwo[0] = dimid_Nc;
    dimidtwo[1] = dimid_numsides;
    if ((retval = nc_def_var(ncid,"normal",NC_INT,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Dot product of unique normal with outward normal of each edge");
    //if ((retval = nc_put_var_int(ncid,varid, grid->normal)))
    //  ERR(retval);

    //n1
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"n1",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","x-component of the edge normal");
    nc_addattr(ncid, varid,"coordinates","xe ye");
    //if ((retval = nc_put_var_double(ncid,varid, grid->n1)))
    //  ERR(retval);

    //n2
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"n2",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","y-component of the edge normal");
    nc_addattr(ncid, varid,"coordinates","xe ye");
    //if ((retval = nc_put_var_double(ncid,varid, grid->n2)))
    //  ERR(retval);

    //df
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"df",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","edge length");
    nc_addattr(ncid, varid,"units","m");
    nc_addattr(ncid, varid,"coordinates","xe ye");
    //if ((retval = nc_put_var_double(ncid,varid, grid->df)))
    //  ERR(retval);

    //dg
    dimidone[0] = dimid_Ne;
    if ((retval = nc_def_var(ncid,"dg",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","distance between faces on either side of edge");
    nc_addattr(ncid, varid,"units","m");
    nc_addattr(ncid, varid,"coordinates","xe ye");
    //if ((retval = nc_put_var_double(ncid,varid, grid->dg)))
    //  ERR(retval);

    //def
    dimidtwo[0] = dimid_Nc;
    dimidtwo[1] = dimid_numsides;
    if ((retval = nc_def_var(ncid,"def",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Distance between faces and edges");
    nc_addattr(ncid, varid,"units","m");
    nc_addattr(ncid, varid,"coordinates","xe ye");
    //if ((retval = nc_put_var_double(ncid,varid, grid->def)))
    //  ERR(retval);

    //Ac
    dimidone[0] = dimid_Nc;
    if ((retval = nc_def_var(ncid,"Ac",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Horizontal area of 2D mesh");
    nc_addattr(ncid, varid,"units","m2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","xv yv");
    //if ((retval = nc_put_var_double(ncid,varid, grid->Ac)))
    //  ERR(retval);

    /********************************************************************** 
    *
    * Define the vertical grid variables and attributes 
    *
    **********************************************************************/
   //dz
   dimidone[0] = dimid_Nk;   
   if ((retval = nc_def_var(ncid,"dz",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","z layer spacing");
   nc_addattr(ncid, varid,"units","m");
   //if ((retval = nc_put_var_double(ncid,varid, grid->dz)))
   //  ERR(retval);

   // Calculate and write the vertical coordinate z levels 
   z_w[0]=0.0;
   for(k=0;k<grid->Nkmax;k++){
      z_w[k+1] = z_w[k] + grid->dz[k];
      if(k==0){
	 z_r[k] = grid->dz[k]*0.5;
      }else{
	 z_r[k] = z_r[k-1]+grid->dz[k];
      }
   }
   //z_r
   dimidone[0] = dimid_Nk;   
   if ((retval = nc_def_var(ncid,"z_r",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"standard_name","ocean_z_coordinate");
   nc_addattr(ncid, varid,"long_name","depth at layer mid points");
   nc_addattr(ncid, varid,"units","m");  
   nc_addattr(ncid, varid,"positive","up");  
   //if ((retval = nc_put_var_double(ncid,varid, z_r)))
   //  ERR(retval);

   //z_w
   dimidone[0] = dimid_Nkw;   
   if ((retval = nc_def_var(ncid,"z_w",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"standard_name","ocean_z_coordinate");
   nc_addattr(ncid, varid,"long_name","depth at layer edges");
   nc_addattr(ncid, varid,"units","m");  
   nc_addattr(ncid, varid,"positive","up");  
   //if ((retval = nc_put_var_double(ncid,varid, z_w)))
   //  ERR(retval);

   //Nk
   dimidone[0] = dimid_Nc;   
   if ((retval = nc_def_var(ncid,"Nk",NC_INT,1,dimidone,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Number of layers at face"); 
   //if ((retval = nc_put_var_int(ncid,varid, grid->Nk)))
   //  ERR(retval);

   //Nke
   dimidone[0] = dimid_Ne;   
   if ((retval = nc_def_var(ncid,"Nke",NC_INT,1,dimidone,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Number of layers at edge");
   //if ((retval = nc_put_var_int(ncid,varid, grid->Nke)))
   //  ERR(retval);

    //dv
    dimidone[0] = dimid_Nc;
    if ((retval = nc_def_var(ncid,"dv",NC_DOUBLE,1,dimidone,&varid)))
      ERR(retval);
    nc_addattr(ncid, varid,"standard_name","sea_floor_depth_below_geoid");
    nc_addattr(ncid, varid,"long_name","seafloor depth");
    nc_addattr(ncid, varid,"units","m");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","xv yv");
    //if ((retval = nc_put_var_double(ncid,varid, grid->dv)))
    //  ERR(retval);

    //time
    dimidone[0] = dimid_time;
    if ((retval = nc_def_var(ncid,"time",NC_DOUBLE,1,dimidone,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"long_name","time");
    nc_addattr(ncid, varid,"units","seconds since 1990-01-01 00:00:00");  

  //  //average time
  //  dimidone[0] = 0;
  //  if ((retval = nc_def_var(ncid,"average_time",NC_DOUBLE,1,dimidone,&varid)))
  //     ERR(retval);
  //  nc_addattr(ncid, varid,"long_name","Averaging time interval");
  //  nc_addattr(ncid, varid,"units","seconds");  


   /********************************************************************** 
    * 
    * Define the physical variables and attributes 
    *
    **********************************************************************/
   
   dimidtwo[0] = dimid_time;
   dimidtwo[1] = dimid_Nc;
   
   dimidthree[0] = dimid_time;
   dimidthree[1] = dimid_Nk;
   dimidthree[2] = dimid_Nc;
   
   // eta
   if ((retval = nc_def_var(ncid,"eta",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Time-averaged Sea surface elevation");
   nc_addattr(ncid, varid,"units","m");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time yv xv");
    
   //u
   if ((retval = nc_def_var(ncid,"uc",NC_DOUBLE,3,dimidthree,&varid)))
      ERR(retval);
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Time-averaged Eastward water velocity component");
   nc_addattr(ncid, varid,"units","m s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time z_r yv xv");

   //v
   if ((retval = nc_def_var(ncid,"vc",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval);   
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Time-averaged Northward water velocity component");
   nc_addattr(ncid, varid,"units","m s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time z_r yv xv");
   
   //w
   dimidthree[1] = dimid_Nkw;
   if ((retval = nc_def_var(ncid,"w",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval); 
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Time-averaged Vertical water velocity component");
   nc_addattr(ncid, varid,"units","m s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time z_w yv xv");  

   dimidthree[1] = dimid_Nk;
   
   //nu_v
   if ((retval = nc_def_var(ncid,"nu_v",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval); 
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Time-averaged Vertical eddy viscosity");
   nc_addattr(ncid, varid,"units","m2 s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time z_r yv xv");
   
   //kappa_tv
   if ((retval = nc_def_var(ncid,"kappa_tv",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval); 
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Time-averaged vertical tracer diffusivity");
   nc_addattr(ncid, varid,"units","m2 s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","face");
   nc_addattr(ncid, varid,"coordinates","time z_r yv xv");

   //salinity
   if(prop->beta>0){
     if ((retval = nc_def_var(ncid,"salt",NC_DOUBLE,3,dimidthree,&varid)))
      ERR(retval);
     if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
     if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Time-averaged Salinity");
    nc_addattr(ncid, varid,"units","ppt");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time z_r yv xv");

    // Depth-integrated salinity
    if ((retval = nc_def_var(ncid,"s_dz",NC_DOUBLE,2,dimidtwo,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"long_name","Instantaneous depth-integrated salinity");
    nc_addattr(ncid, varid,"units","psu m");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");
    
   }
   
   //temperature
   if(prop->gamma>0){
     if ((retval = nc_def_var(ncid,"temp",NC_DOUBLE,3,dimidthree,&varid)))
      ERR(retval); 
     if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
     if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Time-averaged Water temperature");
    nc_addattr(ncid, varid,"units","degrees C");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time z_r yv xv");

     // Depth-integrated temperature
    if ((retval = nc_def_var(ncid,"T_dz",NC_DOUBLE,2,dimidtwo,&varid)))
       ERR(retval);
    nc_addattr(ncid, varid,"long_name","Instantaneous depth-integrated temperature");
    nc_addattr(ncid, varid,"units","degrees C m");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");
    
   }
   
   //rho
   if( (prop->gamma>0) || (prop->beta>0) ){
     if ((retval = nc_def_var(ncid,"rho",NC_DOUBLE,3,dimidthree,&varid)))
      ERR(retval);
     if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
     if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Time-averaged Water density");
    nc_addattr(ncid, varid,"units","kg m-3");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time z_r yv xv");
   }
   
   //age
   if(prop->calcage>0){
     if ((retval = nc_def_var(ncid,"agec",NC_DOUBLE,3,dimidthree,&varid)))
      ERR(retval);
     if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
     if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Age concentration");
    nc_addattr(ncid, varid,"units","");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time z_r yv xv");
    
    if ((retval = nc_def_var(ncid,"agealpha",NC_DOUBLE,3,dimidthree,&varid)))
      ERR(retval);
    if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
    if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
    nc_addattr(ncid, varid,"long_name","Age alpha parameter");
    nc_addattr(ncid, varid,"units","seconds");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time z_r yv xv");

   }
   //U_F
   dimidthree[2] = dimid_Ne;
   if ((retval = nc_def_var(ncid,"U_F",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval); 
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Time-averaged edge flux rate");
   nc_addattr(ncid, varid,"units","m3 s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","edge");
   nc_addattr(ncid, varid,"coordinates","time z_r ye xe");  

    //s_F
   dimidthree[2] = dimid_Ne;
   if ((retval = nc_def_var(ncid,"s_F",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval); 
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Time-averaged edge salt flux rate");
   nc_addattr(ncid, varid,"units","psu m3 s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","edge");
   nc_addattr(ncid, varid,"coordinates","time z_r ye xe");  

    //T_F
   dimidthree[2] = dimid_Ne;
   if ((retval = nc_def_var(ncid,"T_F",NC_DOUBLE,3,dimidthree,&varid)))
     ERR(retval); 
   if ((retval = nc_def_var_fill(ncid,varid,nofill,&FILLVALUE))) // Sets a _FillValue attribute
      ERR(retval);
   if ((retval = nc_def_var_deflate(ncid,varid,0,DEFLATE,DEFLATELEVEL))) // Compresses the variable
      ERR(retval);
   nc_addattr(ncid, varid,"long_name","Time-averaged edge temperature flux rate");
   nc_addattr(ncid, varid,"units","degreesC m3 s-1");
   nc_addattr(ncid, varid,"mesh","suntans_mesh");
   nc_addattr(ncid, varid,"location","edge");
   nc_addattr(ncid, varid,"coordinates","time z_r ye xe");  


   // Meteorological variables (2-D) //
   
   if(prop->metmodel>0){
    // Uwind
    if ((retval = nc_def_var(ncid,"Uwind",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Time-averaged Eastward wind velocity component");
    nc_addattr(ncid, varid,"units","m s-1");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");

    // Vwind
    if ((retval = nc_def_var(ncid,"Vwind",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Time-averaged Eastward wind velocity component");
    nc_addattr(ncid, varid,"units","m s-1");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");

    // Tair
    if ((retval = nc_def_var(ncid,"Tair",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval);
    nc_addattr(ncid, varid,"long_name","Time-averaged Air temperature");
    nc_addattr(ncid, varid,"units","degrees C");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");

    // Pair
    if ((retval = nc_def_var(ncid,"Pair",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Time-averaged Air pressure");
    nc_addattr(ncid, varid,"units","millibar");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");

    // rain
    if ((retval = nc_def_var(ncid,"rain",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Time-averaged Rain fall rate");
    nc_addattr(ncid, varid,"units","kg m2 s-1");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");

    //RH
    if ((retval = nc_def_var(ncid,"RH",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Time-averaged Relative humidity");
    nc_addattr(ncid, varid,"units","Percent (%)");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");

    //cloud
    if ((retval = nc_def_var(ncid,"cloud",NC_DOUBLE,2,dimidtwo,&varid)))
	ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Time-averaged Cloud cover fraction");
    nc_addattr(ncid, varid,"units","");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");
    
    // Surface flux variables //
    // Hs
    if ((retval = nc_def_var(ncid,"Hs",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Time-averaged Sensible heat flux");
    nc_addattr(ncid, varid,"units","W m-2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");   
    nc_addattr(ncid, varid,"positive","down");   

    // Hl
    if ((retval = nc_def_var(ncid,"Hl",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
     nc_addattr(ncid, varid,"long_name","Time-averaged Latent heat flux");
    nc_addattr(ncid, varid,"units","W m-2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");   
    nc_addattr(ncid, varid,"positive","down");   

    // Hlw
    if ((retval = nc_def_var(ncid,"Hlw",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
     nc_addattr(ncid, varid,"long_name","Time-averaged Net longwave radiation flux");
    nc_addattr(ncid, varid,"units","W m-2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");      
    nc_addattr(ncid, varid,"positive","down");   

    // Hsw
    if ((retval = nc_def_var(ncid,"Hsw",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Time-averaged Net shortwave radiation flux");
    nc_addattr(ncid, varid,"units","W m-2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");   
    nc_addattr(ncid, varid,"positive","down");   

    // tau_x
    if ((retval = nc_def_var(ncid,"tau_x",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
    nc_addattr(ncid, varid,"long_name","Time-averaged Eastward component surface wind stress");
    nc_addattr(ncid, varid,"units","N m-2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");   

    // tau_y
    if ((retval = nc_def_var(ncid,"tau_y",NC_DOUBLE,2,dimidtwo,&varid)))
      ERR(retval); 
     nc_addattr(ncid, varid,"long_name","Time-averaged Northward component surface wind stress");
    nc_addattr(ncid, varid,"units","N m-2");
    nc_addattr(ncid, varid,"mesh","suntans_mesh");
    nc_addattr(ncid, varid,"location","face");
    nc_addattr(ncid, varid,"coordinates","time yv xv");   

    //EP
    if(prop->beta > 0.0){
	if ((retval = nc_def_var(ncid,"EP",NC_DOUBLE,2,dimidtwo,&varid)))
	    ERR(retval); 
	nc_addattr(ncid, varid,"long_name","Time-averaged surface salt flux (S0*EP)");
	nc_addattr(ncid, varid,"units","psu m s-1");
	nc_addattr(ncid, varid,"mesh","suntans_mesh");
	nc_addattr(ncid, varid,"location","face");
	nc_addattr(ncid, varid,"coordinates","time yv xv");   
    }

   }

   //End file definition mode
   if ((retval = nc_enddef(ncid)))
	ERR(retval);

   /**********************************************************
   *
   * Write data (needs to be done out of definition mode for classic model)
   *
   ****************************************************************/
   nc_write_int(ncid,"cells",grid->cells,myproc);
   nc_write_int(ncid,"face",grid->face,myproc);
   nc_write_int(ncid,"nfaces",grid->nfaces,myproc);
   if (0) {
     nc_write_int(ncid,"edges",grid->edges,myproc);
   } else {
     // InitializeNCUgrid does this instead for edges:
     // Need to convert the edge array that is stored is NUMEDGECOLUMN*Ne where
     // NUMEDGECOLUMN=3
     /* Initialize an edge array */
     edges = (int *)SunMalloc(2*(grid->Ne)*sizeof(int),"InitialiseOutputNCugrid");
     
     for(n=0;n<grid->Ne;n++){
       for(j=0;j<NUMEDGECOLUMNS-1;j++){
         edges[2*n+j] = grid->edges[NUMEDGECOLUMNS*n+j];
       }
     }
     nc_write_int(ncid,"edges",edges,myproc);
     SunFree(edges,2*grid->Ne*sizeof(int),"InitialiseOutputNCugrid");
   }

   nc_write_int(ncid,"neigh",grid->neigh,myproc);
   nc_write_int(ncid,"grad",grid->grad,myproc);
   nc_write_int(ncid,"mnptr",grid->mnptr,myproc);
   nc_write_int(ncid,"eptr",grid->eptr,myproc);
   
   nc_write_double(ncid,"xv",grid->xv,myproc);
   nc_write_double(ncid,"yv",grid->yv,myproc);
   nc_write_double(ncid,"xe",grid->xe,myproc);
   nc_write_double(ncid,"ye",grid->ye,myproc);
   nc_write_double(ncid,"xp",grid->xp,myproc);
   nc_write_double(ncid,"yp",grid->yp,myproc);

   nc_write_int(ncid,"normal",grid->normal,myproc);
   nc_write_double(ncid,"n1",grid->n1,myproc);
   nc_write_double(ncid,"n2",grid->n2,myproc);
   nc_write_double(ncid,"df",grid->df,myproc);
   nc_write_double(ncid,"dg",grid->dg,myproc);
   nc_write_double(ncid,"def",grid->def,myproc);
   nc_write_double(ncid,"Ac",grid->Ac,myproc);

   nc_write_double(ncid,"dz",grid->dz,myproc);
   nc_write_double(ncid,"z_r",z_r,myproc);
   nc_write_double(ncid,"z_w",z_w,myproc);
   nc_write_int(ncid,"Nk",grid->Nk,myproc);
   nc_write_int(ncid,"Nke",grid->Nke,myproc);
   nc_write_double(ncid,"dv",grid->dv,myproc);

  // nc_write_double(ncid,"average_time",(float)prop->ntaverage*prop->dt,myproc);

}// End function

/*
* Function: WriteAverageNCmerge()
* -----------------------------
* Main function for writing SUNTANS output to netcdf file/s
* 
*/
void WriteAverageNCmerge(propT *prop, gridT *grid, averageT *average, physT *phys, metT *met, int blowup, int numprocs, MPI_Comm comm, int myproc){
   int ncid;// = prop->averageNetcdfFileID;
   int varid, retval, k;
   // Start and count vectors for one, two and three dimensional arrays
   size_t startone[] = {prop->avgtimectr};
   size_t countone[] = {1};
   size_t starttwo[] = {prop->avgtimectr,0};
   size_t counttwo[] = {1,grid->Nc};
   size_t startthree[] = {prop->avgtimectr,0,0};
   size_t countthree[] = {1,grid->Nkmax,grid->Nc};
   const size_t countthreew[] = {1,grid->Nkmax+1,grid->Nc};
   const REAL time[] = {prop->nctime};
   char str[2*BUFFERLENGTH], filename[BUFFERLENGTH];

   nc_set_log_level(3); // This helps with debugging errors
   
   // RH moved to UpdateAverageVariables
   //prop->avgctr+=1;

   // if (prop->n==0) {
   //   printf("We do get prop->n==0\n");
   // }
   
   // Output the first time step but don't compute the average
   // note that the first step of output from a restart will not
   // be accurate, since we don't have the integrated values from
   // the previous run.  Still, this should keep the map and
   // average outputs the same length when ntout==ntaverage, and
   // eta should be valid here.
   if(!(prop->n%prop->ntaverage) || prop->n==prop->nstart) {

     // Work out if we need to open a new averages file or not
     if( !(prop->avgtimectr%prop->nstepsperncfile) || (prop->n==1+prop->nstart) ){
       if(prop->avgfilectr>average->initialavgfilectr){
         // Close the old netcdf file
         if(myproc==0){
           printf("Closing opened output netcdf file...\n");
           MPI_NCClose(prop->averageNetcdfFileID);
         }
       }
       
       // Open the new netcdf file
       MPI_GetFile(filename,DATAFILE,"averageNetcdfFile","OpenFiles",myproc);
       sprintf(str,"%s_%04d.nc",filename,prop->avgfilectr);
       if(myproc==0){
         prop->averageNetcdfFileID = MPI_NCOpen(str,NC_CLASSIC_MODEL|NC_NETCDF4,"OpenFiles",myproc);
       }else{
         prop->averageNetcdfFileID=-1;
       }
       
       // Initialise a new output file
       if(myproc==0)
         InitialiseAverageNCugridMerge(prop, grid, average, myproc);
       
       prop->avgfilectr += 1;
       
	// Reset the time counter
       prop->avgtimectr = 0;
       startone[0] = prop->avgtimectr;
       starttwo[0] = prop->avgtimectr;
       startthree[0] = prop->avgtimectr;
     }
     ncid = prop->averageNetcdfFileID;
     
     //Compute the averages 
     ComputeAverageVariables(grid,average,phys,met,prop->avgctr,prop);
     
     //Communicate the values
     SendRecvAverages(prop,grid,average,comm,myproc); 
     
     //Reset the counter
     // RH: moved to ZeroVariables.
     // prop->avgctr=0;
     
     if(myproc==0 && VERBOSE>1){ 
       if(!blowup) 
         printf("Outputting average data to netcdf at step %d of %d\n",prop->n,prop->nsteps+prop->nstart);
       else
         printf("Outputting blowup averagedata to netcdf at step %d of %d\n",prop->n,prop->nsteps+prop->nstart);
     }
     
     /* Write the time data*/
     if(myproc==0){
       if ((retval = nc_inq_varid(ncid, "time", &varid)))
         ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, startone, countone, time )))
         ERR(retval);
       countthree[2] = mergedGrid->Nc;
       counttwo[1] = mergedGrid->Nc;
     }
     /* Write to the physical variables*/

    // 2D cell-centered variables
    nc_write_2D_merge(ncid,prop->avgtimectr,  average->h, prop, grid, "eta", numprocs, myproc, comm);
    nc_write_2D_merge(ncid,prop->avgtimectr,  average->h_avg, prop, grid, "eta_avg", numprocs, myproc, comm);
    if(prop->beta>0)
	nc_write_2D_merge(ncid,prop->avgtimectr,  average->s_dz, prop, grid, "s_dz", numprocs, myproc, comm);
    if(prop->gamma>0)
	nc_write_2D_merge(ncid,prop->avgtimectr,  average->T_dz, prop, grid, "T_dz", numprocs, myproc, comm);

    if(prop->metmodel>0){
    	// Atmospheric flux variables
	nc_write_2D_merge(ncid,prop->avgtimectr,  average->Uwind, prop, grid, "Uwind", numprocs, myproc, comm);
	nc_write_2D_merge(ncid,prop->avgtimectr,  average->Vwind, prop, grid, "Vwind", numprocs, myproc, comm);
	nc_write_2D_merge(ncid,prop->avgtimectr,  average->Tair, prop, grid, "Tair", numprocs, myproc, comm);
	nc_write_2D_merge(ncid,prop->avgtimectr,  average->Pair, prop, grid, "Pair", numprocs, myproc, comm);
	nc_write_2D_merge(ncid,prop->avgtimectr,  average->rain, prop, grid, "rain", numprocs, myproc, comm);
	nc_write_2D_merge(ncid,prop->avgtimectr,  average->RH, prop, grid, "RH", numprocs, myproc, comm);
	nc_write_2D_merge(ncid,prop->avgtimectr,  average->cloud, prop, grid, "cloud", numprocs, myproc, comm);
	nc_write_2D_merge(ncid,prop->avgtimectr,  average->Hs, prop, grid, "Hs", numprocs, myproc, comm);
	nc_write_2D_merge(ncid,prop->avgtimectr,  average->Hl, prop, grid, "Hl", numprocs, myproc, comm);
	nc_write_2D_merge(ncid,prop->avgtimectr,  average->Hlw, prop, grid, "Hlw", numprocs, myproc, comm);
	nc_write_2D_merge(ncid,prop->avgtimectr,  average->Hsw, prop, grid, "Hsw", numprocs, myproc, comm);
	nc_write_2D_merge(ncid,prop->avgtimectr,  average->tau_x, prop, grid, "tau_x", numprocs, myproc, comm);
	nc_write_2D_merge(ncid,prop->avgtimectr,  average->tau_y, prop, grid, "tau_y", numprocs, myproc, comm);
	if(prop->beta > 0)
	    nc_write_2D_merge(ncid,prop->avgtimectr,  average->EP, prop, grid, "EP", numprocs, myproc, comm);

    }

    // 3D cell-centered variables
    nc_write_3D_merge(ncid,prop->avgtimectr,  average->uc, prop, grid, "uc",0, numprocs, myproc, comm);
    nc_write_3D_merge(ncid,prop->avgtimectr,  average->vc, prop, grid, "vc",0, numprocs, myproc, comm);
    nc_write_3D_merge(ncid,prop->avgtimectr,  average->nu_v, prop, grid, "nu_v",0, numprocs, myproc, comm);
    nc_write_3D_merge(ncid,prop->avgtimectr,  average->kappa_tv, prop, grid, "kappa_tv",0, numprocs, myproc, comm);


    if(prop->beta>0)
	nc_write_3D_merge(ncid,prop->avgtimectr,  average->s, prop, grid, "salt",0, numprocs, myproc, comm);

    if(prop->gamma>0)
	nc_write_3D_merge(ncid,prop->avgtimectr,  average->T, prop, grid, "temp",0, numprocs, myproc, comm);

    if( (prop->gamma>0) || (prop->beta>0) ) 
	nc_write_3D_merge(ncid,prop->avgtimectr,  average->rho, prop, grid, "rho",0, numprocs, myproc, comm);

    if(prop->calcage){
	nc_write_3D_merge(ncid,prop->avgtimectr,  average->agec, prop, grid, "agec",0, numprocs, myproc, comm);
	nc_write_3D_merge(ncid,prop->avgtimectr,  average->agealpha, prop, grid, "agealpha",0, numprocs, myproc, comm);
    }
  
    // Vertical velocity 
    nc_write_3D_merge(ncid,prop->avgtimectr,  average->w, prop, grid, "w",1, numprocs, myproc, comm);
    
    // 3D edge-based variables 
    nc_write_3Dedge_merge(ncid,prop->avgtimectr,  average->U_F, prop, grid, "U_F", 0, numprocs, myproc, comm);
    if(prop->beta>0)
	nc_write_3Dedge_merge(ncid,prop->avgtimectr,  average->s_F, prop, grid, "s_F", 0, numprocs, myproc, comm);
    if(prop->gamma>0)
	nc_write_3Dedge_merge(ncid,prop->avgtimectr,  average->T_F, prop, grid, "T_F", 0, numprocs, myproc, comm);
     
    // Zero the arrays after they have been written(don't do it for the initial step)
    if(prop->n>1+prop->nstart) {
      if(myproc==0) {
        printf("Zeroing average variables at prop->n=%d\n",prop->n);
      }
      ZeroAverageVariables(grid,average,prop);
    } else {
      if(myproc==0) {
        printf("Skipping zeroing average variables at prop->n=%d\n",prop->n);
      }
    }

    /* Update the netcdf time index */
    prop->avgtimectr += 1;  
   }
  
} // End of function



/*
* Function: WriteAverageNC()
* -----------------------------
* Main function for writing SUNTANS output to netcdf file/s
* 
*/
void WriteAverageNC(propT *prop, gridT *grid, averageT *average, physT *phys, metT *met, int blowup, MPI_Comm comm, int myproc){
   int ncid = prop->averageNetcdfFileID;
   int varid, retval, k;
   // Start and count vectors for one, two and three dimensional arrays
   const size_t startone[] = {prop->avgtimectr};
   const size_t countone[] = {1};
   const size_t starttwo[] = {prop->avgtimectr,0};
   const size_t counttwo[] = {1,grid->Nc};
   size_t startthree[] = {prop->avgtimectr,0,0};
   size_t countthree[] = {1,grid->Nkmax,grid->Nc};
   const size_t countthreew[] = {1,grid->Nkmax+1,grid->Nc};
   const REAL time[] = {prop->nctime};

   nc_set_log_level(3); // This helps with debugging errors

   // RH moved to UpdateAverageVariables
   //prop->avgctr+=1;
   
   // Output the first time step but don't compute the average 
   if( (prop->n==1+prop->nstart) || !(prop->n%prop->ntaverage) ) {
     //Compute the averages 
     ComputeAverageVariables(grid,average,phys,met,prop->avgctr,prop);
     
     //Communicate the values
     SendRecvAverages(prop,grid,average,comm,myproc); 

     //Reset the counter
     // RH: moved to ZeroVariables.
     // prop->avgctr=0;

     if(myproc==0 && VERBOSE>1){ 
       if(!blowup) 
         printf("Outputting average data to netcdf at step %d of %d\n",prop->n,prop->nsteps+prop->nstart);
       else
         printf("Outputting blowup averagedata to netcdf at step %d of %d\n",prop->n,prop->nsteps+prop->nstart);
     }
     
     /* Write the time data*/
     if ((retval = nc_inq_varid(ncid, "time", &varid)))
       ERR(retval);
     if ((retval = nc_put_vara_double(ncid, varid, startone, countone, time )))
       ERR(retval);
     
    /* Write to the physical variables*/
    if ((retval = nc_inq_varid(ncid, "eta", &varid)))
	ERR(retval);
    if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, average->h )))
 	ERR(retval);
    
    if ((retval = nc_inq_varid(ncid, "uc", &varid)))
	ERR(retval);
    ravel(average->uc, average->tmpvar, grid);
    if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, average->tmpvar )))
	ERR(retval);
    
    if ((retval = nc_inq_varid(ncid, "vc", &varid)))
	ERR(retval);
    ravel(average->vc, average->tmpvar, grid);
    if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, average->tmpvar )))
	ERR(retval);
      
    // write w at cell top and bottom
    if ((retval = nc_inq_varid(ncid, "w", &varid)))
	ERR(retval);
    ravelW(average->w, average->tmpvarW, grid);
    if ((retval = nc_put_vara_double(ncid, varid, startthree, countthreew, average->tmpvarW )))
	ERR(retval);

    if ((retval = nc_inq_varid(ncid, "nu_v", &varid)))
	ERR(retval);
    ravel(average->nu_v, average->tmpvar, grid);
    if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, average->tmpvar )))
	ERR(retval);

    if ((retval = nc_inq_varid(ncid, "kappa_tv", &varid)))
	ERR(retval);
    ravel(average->kappa_tv, average->tmpvar, grid);
    if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, average->tmpvar )))
	ERR(retval);
    
    // Tracers
     if(prop->beta>0){
       if ((retval = nc_inq_varid(ncid, "salt", &varid)))
	  ERR(retval);
      ravel(average->s, average->tmpvar, grid);
      if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, average->tmpvar )))
	  ERR(retval);

	if ((retval = nc_inq_varid(ncid, "s_dz", &varid)))
	    ERR(retval);
	if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, average->s_dz )))
	    ERR(retval);
     }
     
     if(prop->gamma>0){
	if ((retval = nc_inq_varid(ncid, "temp", &varid)))
	  ERR(retval);
	ravel(average->T, average->tmpvar, grid);
	if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, average->tmpvar )))
	  ERR(retval);

	if ((retval = nc_inq_varid(ncid, "T_dz", &varid)))
	    ERR(retval);
	if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, average->T_dz )))
	    ERR(retval);
     }
      
     if( (prop->gamma>0) || (prop->beta>0) ){ 
	if ((retval = nc_inq_varid(ncid, "rho", &varid)))
	  ERR(retval);
	ravel(average->rho, average->tmpvar, grid);
	if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, average->tmpvar )))
	  ERR(retval);
     }

     if(prop->calcage>0){ 
	if ((retval = nc_inq_varid(ncid, "agec", &varid)))
	  ERR(retval);
	ravel(average->agec, average->tmpvar, grid);
	if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, average->tmpvar )))
	  ERR(retval);

	if ((retval = nc_inq_varid(ncid, "agealpha", &varid)))
	  ERR(retval);
	ravel(average->agealpha, average->tmpvar, grid);
	if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, average->tmpvar )))
	  ERR(retval);

    }

     // Edge fluxes
     countthree[2] = grid->Ne;
     if ((retval = nc_inq_varid(ncid, "U_F", &varid)))
	ERR(retval);
     ravelEdge(average->U_F, average->tmpvarE, grid);
     if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, average->tmpvarE )))
        ERR(retval);

     if ((retval = nc_inq_varid(ncid, "s_F", &varid)))
	ERR(retval);
     ravelEdge(average->s_F, average->tmpvarE, grid);
     if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, average->tmpvarE )))
        ERR(retval);

     if ((retval = nc_inq_varid(ncid, "T_F", &varid)))
	ERR(retval);
     ravelEdge(average->T_F, average->tmpvarE, grid);
     if ((retval = nc_put_vara_double(ncid, varid, startthree, countthree, average->tmpvarE )))
        ERR(retval);

     // Wind variables
     if(prop->metmodel>0){
       if ((retval = nc_inq_varid(ncid, "Uwind", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, average->Uwind )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "Vwind", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, average->Vwind )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "Tair", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, average->Tair )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "Pair", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, average->Pair )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "rain", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, average->rain )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "RH", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, average->RH )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "cloud", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, average->cloud )))
	 ERR(retval);
       
       // Heat flux variables
       if ((retval = nc_inq_varid(ncid, "Hs", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, average->Hs )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "Hl", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, average->Hl )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "Hlw", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, average->Hlw )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "Hsw", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, average->Hsw )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "tau_x", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, average->tau_x )))
	 ERR(retval);
       
       if ((retval = nc_inq_varid(ncid, "tau_y", &varid)))
	 ERR(retval);
       if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, average->tau_y )))
	 ERR(retval);
       
       if(prop->beta > 0.0){
	  if ((retval = nc_inq_varid(ncid, "EP", &varid)))
	     ERR(retval);
	  if ((retval = nc_put_vara_double(ncid, varid, starttwo, counttwo, average->EP )))
	     ERR(retval);
       }
     }
     
     // Zero the arrays after they have been written(don't do it for the initial step)
     if(prop->n>1+prop->nstart) {
       // printf("Zeroing average variables at prop->n=%d\n",prop->n);
       ZeroAverageVariables(grid,average,prop);
     } else {
       // printf("Skipping zeroing average variables at prop->n=%d\n",prop->n);
     }
     
     /* Update the netcdf time index */
     prop->avgtimectr += 1;  
   }
  
} // End of function



/*
 * Function: nc_addattr()
 * -------------------------
 *
 * Wrapper function to add a text attribute into a netcdf file
 *
 */
static void nc_addattr(int ncid, int varid, char *attname, char *attvalue){
    int retval;
    if ((retval = nc_put_att_text(ncid, varid, attname ,strlen(attvalue),attvalue)))
	ERR(retval);
	
}//End function

/*
 * Function: nc_addattr_int()
 * -------------------------
 *
 * Wrapper function to add a int attribute into a netcdf file
 *
 */
static void nc_addattr_int(int ncid, int varid, char *attname, int *attvalue){
    int retval;
    if ((retval = nc_put_att_int(ncid, varid, attname ,NC_INT,1,attvalue)))
	ERR(retval);
	
}

/*
 * Function: nc_addattr_real()
 * -------------------------
 *
 * Wrapper function to add a double attribute into a netcdf file
 *
 */
static void nc_addattr_real(int ncid, int varid, char *attname, REAL *attvalue){
    int retval;
    if ((retval = nc_put_att_double(ncid, varid, attname ,NC_DOUBLE,1,attvalue)))
	ERR(retval);
	
}

/* 
* Function: ravel()
* -----------------
* Unravel a 2-D SUNTANS array [Nc, Nk] into a vector 
* This is necessary for writing a 2-D array to netcdf as the missing cells need to be filled
*
*/
static void ravel(REAL **tmparray, REAL *tmpvec,gridT *grid){
  int j,k;
  int nk=grid->Nkmax, nc=grid->Nc;
  
  for(j=0;j<nc;j++){
    for(k=0;k<nk;k++){
      if(k<grid->Nk[j]){
        tmpvec[k*nc+j] = tmparray[j][k];
      }else{
	tmpvec[k*nc+j] = (REAL)EMPTY;
      }
    }
  }
}//End of function


/* 
* Function: ravelW()
* -----------------
* Unravel a 2-D SUNTANS array [Nc, Nk+1] into a vector 
* This is necessary for writing a 2-D array to netcdf as the missing cells need to be filled
*
*/
static void ravelW(REAL **tmparray, REAL *tmpvec,gridT *grid){
  int j,k;
  int nk=grid->Nkmax+1, nc=grid->Nc;
  
  for(j=0;j<nc;j++){
    for(k=0;k<nk;k++){
      if(k<grid->Nk[j]+1){
        tmpvec[k*nc+j] = tmparray[j][k];
      }else{
	tmpvec[k*nc+j] = (REAL)EMPTY;
      }
    }
  }
}//End of function

/* 
* Function: ravelEdge()
* -----------------
* Unravel a 2-D SUNTANS array [Ne, Nk] into a vector 
* This is necessary for writing a 2-D array to netcdf as the missing edges need to be filled
*
*/
static void ravelEdge(REAL **tmparray, REAL *tmpvec,gridT *grid){
  int i,k;
  int nk=grid->Nkmax, ne=grid->Ne;
  
  for(i=0;i<ne;i++){
    for(k=0;k<nk;k++){
      if(k<grid->Nke[i]){
        tmpvec[k*ne+i] = tmparray[i][k];
      }else{
	tmpvec[k*ne+i] = (REAL)EMPTY;
      }
    }
  }
}//End of function
const void* FillValue(int empty){
  /* Converts the EMPTY value expression type to match the type expected by nc_def_var_fill*/
  empty = (REAL)empty;
}
 
/*###############################################################
*
* Meteorological input NetCDF functions
*
#################################################################*/

/*
* Function: ReadMetNC()
* ---------------------
* Main function for reading in the meteorological data from the netcdf file
*
*/

void ReadMetNC(propT *prop, gridT *grid, metinT *metin,int myproc){
    int retval, j,k;
    int t0;
    int varid;
    char *vname;
    size_t start[2];
    size_t count[]={1,1};
    int ncid = prop->metncid;

    if(metin->t0==-1){
	metin->t1 = getTimeRec(prop->nctime,metin->time,(int)metin->nt);
	metin->t0 = metin->t1-1;
	metin->t2 = metin->t1+1;
    }
    t0 = metin->t0;
    
    //printf("Model time(0) = %f, time index = %d of %d\n",prop->nctime,t0,metin->nt);
    start[0] = t0;
    start[1] = 0;
    count[0] = NTmet;

    vname = "Uwind";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from netcdf file...\n",vname);
    count[1] = metin->NUwind;
    nc_read_2D(ncid,vname,start,count, metin->Uwind, myproc);

    vname = "Vwind";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from netcdf file...\n",vname);
    count[1] = metin->NVwind;
    nc_read_2D(ncid,vname,start,count, metin->Vwind, myproc);

    vname = "Tair";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from netcdf file...\n",vname);
    count[1] = metin->NTair;
    nc_read_2D(ncid,vname,start,count, metin->Tair, myproc); 

    vname = "Pair";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from netcdf file...\n",vname);
    count[1] = metin->NPair;
    nc_read_2D(ncid,vname,start,count, metin->Pair, myproc);

    vname = "rain";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from netcdf file...\n",vname);
    count[1] = metin->Nrain;
    nc_read_2D(ncid,vname,start,count, metin->rain, myproc);

    vname = "RH";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from netcdf file...\n",vname);
    count[1] = metin->NRH;
    nc_read_2D(ncid,vname,start,count, metin->RH, myproc);

    vname = "cloud";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from netcdf file...\n",vname);
    count[1] = metin->Ncloud;
    nc_read_2D(ncid,vname,start,count, metin->cloud, myproc);
} //End function

/*
* Function: ReadMetNCcoord()
* --------------------------
* Read the coordinate information from the netcdf file 
*
*/
void ReadMetNCcoord(propT *prop, gridT *grid, metinT *metin, int myproc){
    /* Read the data from the meteorological netcdf file into the metin structure */
    int retval, j;
    int varid;
    char *vname;
    int ncid = prop->metncid;

    printf("Reading MET coordinates\n"); // RH

    /* Get the horizontal coordintates*/
    vname = "x_Uwind";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
      ERRM(retval,vname);
    if ((retval = nc_get_var_double(ncid, varid,metin->x_Uwind))) 
      ERRM(retval,vname); 
    vname = "y_Uwind";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
      ERRM(retval,vname);
    if ((retval = nc_get_var_double(ncid, varid,metin->y_Uwind))) 
      ERRM(retval,vname); 
    vname = "x_Vwind";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
      ERRM(retval,vname);
    if ((retval = nc_get_var_double(ncid, varid,metin->x_Vwind)))
      ERRM(retval,vname);
    vname = "y_Vwind";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
      ERRM(retval,vname);
    if ((retval = nc_get_var_double(ncid, varid,metin->y_Vwind))) 
      ERRM(retval,vname);
    vname = "x_Tair";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
      ERRM(retval,vname);
    if ((retval = nc_get_var_double(ncid, varid,metin->x_Tair))) 
      ERRM(retval,vname);
    vname = "y_Tair";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
      ERRM(retval,vname);
    if ((retval = nc_get_var_double(ncid, varid,metin->y_Tair)))
      ERRM(retval,vname);
    vname = "x_Pair";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid,metin->x_Pair))) 
      ERR(retval); 
    vname = "y_Pair";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid,metin->y_Pair))) 
      ERR(retval); 
    vname = "x_rain";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid,metin->x_rain))) 
      ERR(retval); 
    vname = "y_rain";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid,metin->y_rain))) 
      ERR(retval); 
    vname = "x_RH";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid,metin->x_RH))) 
      ERR(retval); 
    vname = "y_RH";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid,metin->y_RH))) 
      ERR(retval); 
    vname = "x_cloud";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid,metin->x_cloud))) 
      ERR(retval); 
    vname = "y_cloud";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid,metin->y_cloud))) 
      ERR(retval); 
    
    /* Vertical coordinates */
    vname = "z_Uwind";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid,metin->z_Uwind))) 
      ERR(retval); 
    vname = "z_Vwind";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid,metin->z_Vwind))) 
      ERR(retval); 
    vname = "z_Tair";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid,metin->z_Tair))) 
      ERR(retval); 
    vname = "z_RH";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid,metin->z_RH))) 
      ERR(retval); 
    
    /* Time */
    vname = "Time";
    if(VERBOSE>2 && myproc==0) printf("Reading variable: %s...\n",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
      ERRM(retval,vname);
    if ((retval = nc_get_var_double(ncid, varid,metin->time))) 
      ERRM(retval,vname); 
    
    if(VERBOSE>2 && myproc==0) printf("Finished Reading met netcdf coordinates.\n");
}// End function

/*
* Function: returndimlen()
* ------------------------
* Returns the length of a netcdf dimension 
*/

size_t returndimlen(int ncid, char *dimname){
 int retval;
 int dimid;
 size_t dimlen;

 if ((retval =nc_inq_dimid(ncid,dimname,&dimid)))
   ERRM(retval,dimname);

 if ((retval = nc_inq_dimlen(ncid,dimid, &dimlen)))
   ERRM(retval,dimname);
 return dimlen;
} //End function

/*###############################################################
*
* Boundary contition input NetCDF functions
*
#################################################################*/

/*
 * Function: ReadBdyNC()
 * -----------------------------
 * Reads in boundary netcdf data into the forward and back time steps 
 *
 */
void ReadBdyNC(propT *prop, gridT *grid, int myproc, MPI_Comm comm){
    int retval, j, k, n, p, sendSize;
    int scal_idx;
    int t0, t1;
    int varid;
    char *vname;
    char namebuff[300];
    
    size_t start[]={0,0,0};
    size_t start2[]={0,0};
    size_t count[]={0,0,0};
    size_t count2[]={0,0};
    int ncid = prop->netcdfBdyFileID;  
    size_t Nk = bound->Nk;
    size_t Ntype3 = bound->Ntype3;
    size_t Ntype2 = bound->Ntype2;
    size_t Nseg = bound->Nseg;

    //Find the time index of the middle time step (t1) 
    if(bound->t0==-1){
       bound->t1 = getTimeRecBnd(prop->nctime,bound->time,(int)bound->Nt); //this is in met.c
       bound->t0=bound->t1-1;
       bound->t2=bound->t1+1;
       printf("myproc: %d, bound->t0: %d, nctime: %f\n",myproc,bound->t0, prop->nctime);
    }
    t0 = bound->t0;
    t1 = bound->t1;

    count[0]=NT;
    count[1]=Nk;

    count2[0]=NT;

    start[0]=t0;
    start[1]=0;
    start[2]=0;

    start2[0]=t0;
    start2[1]=0;

    //if(myproc==0) printf("t0 = %d [Nt = %d]\n",t0,bound->Nt);    
    if(bound->hasType2){
      
      count[2]=Ntype2;

      vname = "boundary_u";
      if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from boundry netcdf file...\n",vname);
      nc_read_3D(ncid, vname, start, count, bound->boundary_u_t );
      
      vname = "boundary_v";
      if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from boundry netcdf file...\n",vname);
      nc_read_3D(ncid, vname, start, count, bound->boundary_v_t );
      
      vname = "boundary_w";
      if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from boundry netcdf file...\n",vname);
      nc_read_3D(ncid, vname, start, count, bound->boundary_w_t );

      for(scal_idx=0;scal_idx<bound->num_scalars;scal_idx++) {
        sprintf(namebuff,"boundary_%s",bound->scalars[scal_idx].varname); // e.g. "boundary_T"
        if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from boundry netcdf file...\n",namebuff);
        nc_read_3D(ncid, namebuff, start, count, bound->scalars[scal_idx].boundary_scal_t );
      }
    }

    if(bound->hasType3){
      count[0]=NT;
      count[1]=Nk;
      count[2]=Ntype3;
      count2[1]=Ntype3;
      
      sendSize = count[0]*count[1]*count[2];
      
      vname = "uc";
      if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from boundry netcdf file...\n",vname);
      nc_read_3D(ncid, vname, start, count, bound->uc_t );
      
      vname = "vc";
      if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from boundry netcdf file...\n",vname);
      nc_read_3D(ncid, vname, start, count, bound->vc_t );
      
      vname = "wc";
      if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from boundry netcdf file...\n",vname);
      nc_read_3D(ncid, vname, start, count, bound->wc_t );

      for(scal_idx=0;scal_idx<bound->num_scalars;scal_idx++) {
        vname=bound->scalars[scal_idx].varname; // e.g. "T"
        if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from boundry netcdf file...\n",vname);
        nc_read_3D(ncid, vname, start, count, bound->scalars[scal_idx].scal_t );
      }
      
      vname = "h";//2D array
      if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from boundry netcdf file...\n",vname);
      nc_read_2D(ncid, vname, start2, count2, bound->h_t, myproc );
     }// End read type-3

     //Flux boundary data
    if(bound->hasType2 && bound->hasSeg){
      
      count2[1]=Nseg;
      //if(myproc==0){
      vname = "boundary_Q";//2D array
      if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from boundry netcdf file...\n",vname);
      nc_read_2D(ncid, vname, start2, count2, bound->boundary_Q_t, myproc);
      //}
      //sendSize = count2[0]*count2[1];
      //MPI_Bcast(&(bound->boundary_Q_t[0][0]),sendSize,MPI_DOUBLE,0,comm);
    }//End flux read
    
    if(bound->Npoint_source) {
       count2[1]=bound->Npoint_source;
       // these differ from the above 
       count2[0]=1; // no quadratic interpolation
       start2[0]=t1; // just read the middle time step

       if(VERBOSE>2) {
         printf("[p=%d] Reading %ld point source at time index %d\n",myproc,bound->Npoint_source,t1);
       }
       
       vname="point_Q";
       if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from boundary netcdf file...\n",vname);
       if ((retval = nc_inq_varid(ncid, vname, &varid)))
         ERRM(retval,"looking up point_Q");
       if ((retval = nc_get_vara_double(ncid, varid, start2, count2, bound->point_Q)))
         ERRM(retval,"reading point_Q"); 
       // nc_read_2D(ncid, vname, start2, count2, (REAL**)bound->point_Q, myproc);

       for(scal_idx=0;scal_idx<bound->num_scalars;scal_idx++) {
         sprintf(namebuff,"point_%s",bound->scalars[scal_idx].varname); // e.g. "point_T"
         vname=namebuff;
         
         if(VERBOSE>2 && myproc==0) printf("Reading variable: %s from boundary netcdf file...\n",vname);
         if ((retval = nc_inq_varid(ncid, vname, &varid)))
           ERRM(retval,vname);
         if ((retval = nc_get_vara_double(ncid, varid, start2, count2,
                                          bound->scalars[scal_idx].point_scal)))
           ERRM(retval,vname); 
       }       
     }

   // Wait for all processors
   //MPI_Barrier(comm);

 }//End function

/*
 * Function: ReadBndNCcoord()
 * --------------------------
 * Reads the coordinate information from the netcdf file into the boundary structure
 *
 */
void ReadBndNCcoord(int ncid, propT *prop, gridT *grid, int myproc, MPI_Comm comm){

    int retval, j;
    int varid;
    char *vname;
    //int ncid = prop->netcdfBdyFileID; 

    //if(myproc==0){
    vname = "time";
    if(VERBOSE>2 && myproc==0) printf("Reading boundary variable: %s...",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid,bound->time))) 
      ERR(retval); 
    if(VERBOSE>2 && myproc==0) printf("done.\n");

    vname = "z";
    if(VERBOSE>2 && myproc==0) printf("Reading boundary variable: %s...",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid,bound->z))) 
      ERR(retval); 
    if(VERBOSE>2 && myproc==0) printf("done.\n");

    if(bound->hasType3>0){

	vname = "xv";
	if(VERBOSE>2 && myproc==0) printf("Reading boundary variable: %s...",vname);
	if ((retval = nc_inq_varid(ncid, vname, &varid)))
	    ERR(retval);
	if ((retval = nc_get_var_double(ncid, varid,bound->xv))) 
	  ERR(retval); 
	if(VERBOSE>2 && myproc==0) printf("done.\n");

	vname = "yv";
	if(VERBOSE>2 && myproc==0) printf("Reading boundary variable: %s...",vname);
	if ((retval = nc_inq_varid(ncid, vname, &varid)))
	    ERR(retval);
	if ((retval = nc_get_var_double(ncid, varid,bound->yv))) 
	      ERR(retval); 
	if(VERBOSE>2 && myproc==0) printf("done.\n");

	vname = "cellp";
	if(VERBOSE>2 && myproc==0) printf("Reading boundary variable: %s...",vname);
	if ((retval = nc_inq_varid(ncid, vname, &varid)))
	    ERR(retval);
	if ((retval = nc_get_var_int(ncid, varid,bound->cellp))) 
	  ERR(retval); 
	if(VERBOSE>2 && myproc==0) printf("done.\n");
    }//end if

    if(bound->Npoint_source){
      vname="point_cell";
      if(VERBOSE>2 && myproc==0) printf("Reading boundary variable: %s...",vname);
      if ((retval = nc_inq_varid(ncid, vname, &varid)))
        ERR(retval);
      if ((retval = nc_get_var_int(ncid, varid,bound->point_cell))) 
        ERR(retval); 
      if(VERBOSE>2 && myproc==0) printf("done.\n");

      vname="point_layer";
      if(VERBOSE>2 && myproc==0) printf("Reading boundary variable: %s...",vname);
      if ((retval = nc_inq_varid(ncid, vname, &varid)))
        ERR(retval);
      if ((retval = nc_get_var_int(ncid, varid,bound->point_layer))) 
        ERR(retval); 
      if(VERBOSE>2 && myproc==0) printf("done.\n");
    }
    
    if(bound->hasType2>0){

	vname = "xe";
	if(VERBOSE>2 && myproc==0) printf("Reading boundary variable: %s...\n",vname);
	if ((retval = nc_inq_varid(ncid, vname, &varid)))
	    ERR(retval);
	if ((retval = nc_get_var_double(ncid, varid,bound->xe))) 
	  ERR(retval); 

	vname = "ye";
	if(VERBOSE>2 && myproc==0) printf("Reading boundary variable: %s...\n",vname);
	if ((retval = nc_inq_varid(ncid, vname, &varid)))
	    ERR(retval);
	if ((retval = nc_get_var_double(ncid, varid,bound->ye))) 
	  ERR(retval); 

	vname = "edgep";
	if(VERBOSE>2 && myproc==0) printf("Reading boundary variable: %s...\n",vname);
	if ((retval = nc_inq_varid(ncid, vname, &varid)))
	    ERR(retval);
	if ((retval = nc_get_var_int(ncid, varid,bound->edgep))) 
	  ERR(retval); 
    }//end if

    if(bound->hasSeg>0){
	vname = "segedgep";
	if(VERBOSE>2 && myproc==0) printf("Reading boundary variable: %s...\n",vname);
	if ((retval = nc_inq_varid(ncid, vname, &varid)))
	    ERR(retval);
	if ((retval = nc_get_var_int(ncid, varid,bound->segedgep))) 
	  ERR(retval); 

	vname = "segp";
	if(VERBOSE>2 && myproc==0) printf("Reading boundary variable: %s...\n",vname);
	if ((retval = nc_inq_varid(ncid, vname, &varid)))
	    ERR(retval);
	if ((retval = nc_get_var_int(ncid, varid,bound->segp))) 
	  ERR(retval); 
    }
    //} // End processor 0 read

    /*
    //Distribute the arrays
    MPI_Bcast(&(bound->time[0]),bound->Nt,MPI_DOUBLE,0,comm);
    MPI_Bcast(&(bound->z[0]),bound->Nk,MPI_DOUBLE,0,comm);
    if(bound->hasType3>0){
	MPI_Bcast(&(bound->xv[0]),bound->Ntype3,MPI_DOUBLE,0,comm);
	MPI_Bcast(&(bound->yv[0]),bound->Ntype3,MPI_DOUBLE,0,comm);
	MPI_Bcast(&(bound->cellp[0]),bound->Ntype3,MPI_INT,0,comm);
    }
    if(bound->hasType2>0){
	MPI_Bcast(&(bound->xe[0]),bound->Ntype2,MPI_DOUBLE,0,comm);
	MPI_Bcast(&(bound->ye[0]),bound->Ntype2,MPI_DOUBLE,0,comm);
	MPI_Bcast(&(bound->edgep[0]),bound->Ntype2,MPI_INT,0,comm);
    }
    if(bound->hasSeg>0){
	MPI_Bcast(&(bound->segedgep[0]),bound->Ntype2,MPI_INT,0,comm);
	MPI_Bcast(&(bound->segp[0]),bound->Nseg,MPI_INT,0,comm);
    }
    */

   // Wait for all processors
   MPI_Barrier(comm);

}//End function

/*
* Function: returndimlenBC()
* --------------------------
* Returns the length of a dimension 
* Returns a zero if the dimension is not found and does not raise an error
*/
size_t returndimlenBC(int ncid, char *dimname){
  int retval;
  int dimid;
  size_t dimlen;
  
  if ((retval =nc_inq_dimid(ncid,dimname,&dimid)))
    return 0;

  if ((retval = nc_inq_dimlen(ncid,dimid, &dimlen)))
    ERR(retval);

  return dimlen;
} // End function

/*###############################################################
*
* Initial contition input NetCDF functions
*/

/*
 * Function: ReadInitialNCcoord()
 * -----------------------------
 * Reads the dimensions from the initial condition netcdf file
 *
 */
void ReadInitialNCcoord(propT *prop, gridT *grid, int *Nci, int *Nei, int *Nki, int *T0, int myproc){

  int Nt;

   // Read the spatial dimension sizes
    *Nci = (int)returndimlenBC(prop->initialNCfileID,"Nc");
    *Nei = (int)returndimlenBC(prop->initialNCfileID,"Ne"); // will be 0 if not present.
    *Nki = (int)returndimlenBC(prop->initialNCfileID,"Nk");

    // Check the dimension with the grid
    if(*Nki != grid->Nkmax){
	printf("Error! Number of layers in initial condition file (%d) not equal to Nkmax (%d).\n",*Nki,grid->Nkmax); 
	MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    // Find the index of the closest time point, T0
    Nt = (int)returndimlenBC(prop->initialNCfileID,"time");
    *T0 = getICtime(prop,Nt, myproc);
    if (*T0>=Nt) *T0 = Nt-1;
    //*T0=0;
    return;

 } // End function

/*
 * Function: getICtime()
 * -----------------------------
 * Return the closest time index from the initial condition file
 *
 */
int getICtime(propT *prop, int Nt, int myproc){

   int retval, varid;
   int ncid = prop->initialNCfileID;
   //REAL time[Nt]; 
   REAL *ictime;
   char *vname;

   ictime = (REAL *)SunMalloc(Nt*sizeof(REAL),"getICtime");

   vname = "time";
    if(VERBOSE>2 && myproc==0) printf("Reading initial condition %s...",vname);
    if ((retval = nc_inq_varid(ncid, vname, &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid, &ictime[0] ))) 
      ERR(retval); 
    if(VERBOSE>2 && myproc==0) printf("done.\n");

    return getTimeRecBnd(prop->nctime, ictime, (int)Nt);

} // End function

/*
 * Function: ReturnFreeSurfaceNC()
 * -------------------------------
 * Reads the free surface from the initial condition netcdf array
 *
 */
void ReturnFreeSurfaceNC(propT *prop, physT *phys, gridT *grid, REAL *htmp, int Nci, int T0, int myproc){
   int i;
   size_t start[] = {T0, 0};
   size_t count[] = {1,Nci};
   //REAL htmp[Nci];

   int varid, retval;
   int ncid = prop->initialNCfileID;

   if(VERBOSE>1 && myproc==0) printf("Reading free-surface initial condition from netcdf file...\n");
   //printf("Initial condition file: T0 = %d, Nci = %d\n",T0,Nci);
   //nc_read_2D(prop->initialNCfileID, "eta", start, count, htmp , myproc);
    if ((retval = nc_inq_varid(ncid, "eta", &varid)))
	ERR(retval);
    if ((retval = nc_get_vara_double(ncid, varid, start, count, &htmp[0]))) 
	ERR(retval); 

   for(i=0;i<grid->Nc;i++) {

     phys->h[i]=htmp[grid->mnptr[i]];
     //phys->h[i]=0;
     // RH consolidate this in phys.c where this is called, to reduce potential
     // for differences between netcdf and non-netcdf
     // if(phys->h[i]<-grid->dv[i] + DRYCELLHEIGHT) 
     //   phys->h[i]=-grid->dv[i] + DRYCELLHEIGHT;
  }
} // End function


/*
 * Function: ReturnSalinityNC()
 * -------------------------------
 * Reads the salinity from the initial condition netcdf array
 *
 */
void ReturnSalinityNC(propT *prop, physT *phys, gridT *grid, REAL *htmp, int Nci, int Nki, int T0, int myproc){
   int i,k,ind;
   size_t start[] = {T0, 0, 0};
   size_t count[] = {1, Nki, Nci};
   //REAL htmp[Nki][Nci];

   int varid, retval;
   int ncid = prop->initialNCfileID;

   if(VERBOSE>1 && myproc==0) printf("Reading salinity initial condition from netcdf file...\n");
    if ((retval = nc_inq_varid(ncid, "salt", &varid)))
	ERR(retval);
    if ((retval = nc_get_vara_double(ncid, varid, start, count, &htmp[0]))) 
	ERR(retval); 

   for(i=0;i<grid->Nc;i++) {
     k=grid->ctop[i];
     if(k==grid->Nk[i]) k--; // be sure to get the bed layer even if 'dry' from ctop.
     for(;k<grid->Nk[i];k++) {
       ind = k*Nci + grid->mnptr[i]; 
       phys->s[i][k]=htmp[ind];
       phys->s0[i][k]=htmp[ind];
     }
  }
} // End function

// Pull sediment initial condition for size class sizeno from netcdf
void ReturnSedimentNC(int sizeno, propT *prop, physT *phys, gridT *grid, REAL *htmp, int Nci, int Nki, int T0, int myproc){
  int i,k,ind;
  size_t start[] = {T0, 0, 0};
  size_t count[] = {1, Nki, Nci};
  
  char varname[256];
   
  int varid, retval;
  int ncid = prop->initialNCfileID;

  // 1-based to match the names in suntans.dat
  sprintf(varname,"sed%d",sizeno+1);
  
  if(VERBOSE>1 && myproc==0) printf("Reading %s initial condition from netcdf file...\n",varname);
  if ((retval = nc_inq_varid(ncid, varname, &varid))) {
    printf("Could not find sediment initial condition %s, will default to 0.0\n",varname);
    for(i=0;i<grid->Nc;i++) {
      k=grid->ctop[i];
      if(k==grid->Nk[i]) k--; // be sure to get the bed layer even if 'dry' from ctop.
      for(;k<grid->Nk[i];k++) {
        sediments->SediC[sizeno][i][k]=0.0;
      }
    }
    return;
  }
  
  if ((retval = nc_get_vara_double(ncid, varid, start, count, &htmp[0]))) 
    ERRM(retval," could not read sediment initial condition"); 
  
  for(i=0;i<grid->Nc;i++) {
    k=grid->ctop[i];
    if(k==grid->Nk[i]) k--; // be sure to get the bed layer even if 'dry' from ctop.
    for(;k<grid->Nk[i];k++) {
      ind = k*Nci + grid->mnptr[i]; 
      sediments->SediC[sizeno][i][k]=htmp[ind];
    }
  }
} // End function


/*
 * Function: ReturnTemperatureNC()
 * -------------------------------
 * Reads the temperature from the initial condition netcdf array
 *
 */
void ReturnTemperatureNC(propT *prop, physT *phys, gridT *grid, REAL *htmp, int Nci, int Nki, int T0, int myproc){
   int i,k,ind;
   size_t start[] = {T0, 0, 0};
   size_t count[] = {1, Nki, Nci};
   //REAL htmp[Nki][Nci];

   int varid, retval;
   int ncid = prop->initialNCfileID;

   if(VERBOSE>1 && myproc==0) printf("Reading temperature initial condition from netcdf file...\n");
    if ((retval = nc_inq_varid(ncid, "temp", &varid)))
	ERR(retval);
    if ((retval = nc_get_vara_double(ncid, varid, start, count, &htmp[0]))) 
	ERR(retval); 

   for(i=0;i<grid->Nc;i++) {
     k=grid->ctop[i];
     if(k==grid->Nk[i]) k--; // be sure to get the bed layer even if 'dry' from ctop.
     for(;k<grid->Nk[i];k++) {
       ind = k*Nci + grid->mnptr[i]; 
       phys->T[i][k]=htmp[ind];
     }
  }
} // End function

/*
 * Function: ReturnAgeNC()
 * -------------------------------
 * Reads the age variables (agec & agealpha) from the initial condition netcdf array
 *
 */
void ReturnAgeNC(propT *prop, gridT *grid, REAL *htmp, int Nci, int Nki, int T0, int myproc){
   int i,k,ind;
   size_t start[] = {T0, 0, 0};
   size_t count[] = {1, Nki, Nci};
   //REAL htmp[Nki][Nci];

   int varid, retval;
   int ncid = prop->initialNCfileID;

   if(VERBOSE>1 && myproc==0) printf("Reading agec initial condition from netcdf file...\n");
    if ((retval = nc_inq_varid(ncid, "agec", &varid)))
	ERR(retval);
    if ((retval = nc_get_vara_double(ncid, varid, start, count, &htmp[0]))) 
	ERR(retval); 

   for(i=0;i<grid->Nc;i++) {
     k=grid->ctop[i];
     if (k==grid->Nk[i]) k--; // be sure to get the bed layer even if 'dry' from ctop.
     for(;k<grid->Nk[i];k++) {
       ind = k*Nci + grid->mnptr[i]; 
       age->agec[i][k]=htmp[ind];
     }
   }
   
   if(VERBOSE>1 && myproc==0) printf("Reading agealpha initial condition from netcdf file...\n");
   if ((retval = nc_inq_varid(ncid, "agealpha", &varid)))
     ERR(retval);
   if ((retval = nc_get_vara_double(ncid, varid, start, count, &htmp[0]))) 
     ERR(retval); 

   for(i=0;i<grid->Nc;i++) {
     k=grid->ctop[i];
     if (k==grid->Nk[i]) k--; // be sure to get the bed layer even if 'dry' from ctop.
     for(;k<grid->Nk[i];k++) {
       //for(k=0;k<grid->Nk[i];k++) {
       ind = k*Nci + grid->mnptr[i]; 
       age->agealpha[i][k]=htmp[ind];
     }
   }

   if(VERBOSE>1 && myproc==0) printf("Reading agesource term from netcdf file...\n");
    if ((retval = nc_inq_varid(ncid, "agesource", &varid)))
	ERR(retval);
    if ((retval = nc_get_var_double(ncid, varid, &htmp[0]))) 
	ERR(retval); 

   for(i=0;i<grid->Nc;i++) {
     k=grid->ctop[i];
     if (k==grid->Nk[i]) k--; // be sure to get the bed layer even if 'dry' from ctop.
     for(;k<grid->Nk[i];k++) {
       //for(k=0;k<grid->Nk[i];k++) {
       ind = k*Nci + grid->mnptr[i]; 
       age->agesource[i][k]=htmp[ind];
     }
   }
} // End function

/*
 * Function: ReturnZ0BNC()
 * -------------------------------
 * Reads the edge-centered roughness from the initial condition netcdf array
 *
 */
void ReturnZ0BNC(propT *prop, physT *phys, gridT *grid, REAL *htmp, int Nei, int T0, int myproc){
  int j,ind;
  size_t start[] = {T0, 0};
  size_t count[] = {1, Nei};

  int varid, retval;
  int ncid = prop->initialNCfileID;

  if ((retval = nc_inq_varid(ncid, "z0B", &varid))) {
    if(VERBOSE>1 && myproc==0) 
      printf("No roughness in netcdf file.  Will use %.3e from suntans.dat...\n",prop->z0B);
    return;
  }

  // most of IC has a time dimension with a single entry.
  // make sure that is true here to avoid an ugly bug.
  if ((retval=nc_inq_varndims(ncid, varid,&ind)))
    ERR(retval); 

  if(ind!=2) {
    printf("While reading roughness from netcdf variable 'z0B'...\n");
    printf("Expected dimensions (time,edge), but got %d dimensions\n",
           ind);
    exit(1);
  }

  if(VERBOSE>1 && myproc==0) 
    printf("Reading roughness from netcdf variable 'z0B'...\n");

  if ((retval = nc_get_vara_double(ncid, varid, start, count, &htmp[0]))) 
    ERR(retval); 

  for(j=0;j<grid->Ne;j++) {
    ind = grid->eptr[j]; // I think eptr is correct, but not 100% sure
    phys->z0B_spec[j]=htmp[ind];
  }
} // End function


/*
* Function: GetTimeRecBnd()
* ------------------
* Retuns the index of the first preceding time step in the vector time
*/
int getTimeRecBnd(REAL nctime, REAL *time, int nt){
    int j;

    for(j=0;j<nt;j++){
       if (time[j]>=nctime)
	 //return j-1;
	 return j;
    }
    return nt;
}

