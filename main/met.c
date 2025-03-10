/*
 * Routines for parsing meteorological input data into SUNTANS 
 * 
 */

#include "met.h"
//#include "phys.h"
#include "mynetcdf.h"
#include "sendrecv.h"
#include "util.h"

/* Private functions */
void calcInterpWeights(gridT *grid, propT *prop, REAL *xo, REAL *yo, int Ns, int **index, REAL **klambda,int myproc);
static REAL semivariogram(int varmodel, REAL nugget, REAL sill, REAL range, REAL D);
void FindNearestMetStations(propT *prop, gridT *grid, metinT **metin, int myproc);
void weightInterpArray(REAL **D, REAL **klambda, gridT *grid, int Ns, int **index, int nt, REAL **Dout);
void weightInterpField(REAL *D, REAL **klambda, gridT *grid, int Ns, int **index, REAL *Dout);
static REAL specifichumidity(REAL RH, REAL Ta, REAL Pair);
static REAL qsat(REAL Tw, REAL Pair);
static REAL satvap(REAL Ta, REAL Pair);
static REAL vappres(REAL Ta, REAL RH, REAL Pair);
static REAL longwave_berliand(REAL Ta, REAL Tw, REAL C_cloud, REAL RH, REAL Pair);
static REAL longwave(REAL Ta, REAL Tw, REAL C_cloud);
static void cor30a(REAL *y);
static REAL psiu_30(REAL zet);
static REAL psit_30(REAL zet);

/* Start of functions */

/*
* Function: initialiseMetFields()
* -----------------------------
* Driver function to initialise all of the meterological inputs 
*
*/
void InitialiseMetFields(propT *prop, gridT *grid, metinT *metin, metT *met, int myproc){
 
  int retval;
  int i,j;
  

 /*  Read in the coordinate data*/
 if(VERBOSE>3 && myproc==0) printf("Reading netcdf coordinate data...\n");
 ReadMetNCcoord(prop,grid,metin, myproc);

 /* Find the nearest met points to each grid point */
 FindNearestMetStations(prop, grid, &metin, myproc);
 
 /* Calculating the interpolation weights for each variable*/
 calcInterpWeights(grid,prop, metin->x_Uwind, metin->y_Uwind, 
    metin->max_nearest_Uwind, metin->nearest_Uwind, metin->WUwind, myproc);
 calcInterpWeights(grid,prop, metin->x_Vwind, metin->y_Vwind,
     metin->max_nearest_Vwind, metin->nearest_Vwind, metin->WVwind, myproc);
 calcInterpWeights(grid,prop, metin->x_Tair, metin->y_Tair,
    metin->max_nearest_Tair, metin->nearest_Tair, metin->WTair, myproc);
 calcInterpWeights(grid,prop, metin->x_Pair, metin->y_Pair,
    metin->max_nearest_Pair, metin->nearest_Pair, metin->WPair, myproc);
 calcInterpWeights(grid,prop, metin->x_rain, metin->y_rain,
    metin->max_nearest_rain, metin->nearest_rain, metin->Wrain, myproc);
 calcInterpWeights(grid,prop, metin->x_RH, metin->y_RH,
    metin->max_nearest_RH, metin->nearest_RH, metin->WRH, myproc);
 calcInterpWeights(grid,prop, metin->x_cloud, metin->y_cloud,
    metin->max_nearest_cloud, metin->nearest_cloud, metin->Wcloud, myproc);
 
 if(VERBOSE>3 && myproc==0){
    printf("Uwind weights:\n");
    for(i=0;i<grid->Nc;i++){
	    printf("xv=%f, yv=%f, Weights: ",grid->xv[i],grid->yv[i]);
	    for(j=0;j<metin->NUwind;j++){
	      printf("%f, ",metin->WUwind[i][j]);
	    }
	    printf("\n");
    }
 }
 
 /*  Interpolate the heights of some variables */
 if(VERBOSE>1 && myproc==0) printf("Interpolating height coordinates onto grid...");
 weightInterpField(metin->z_Uwind, metin->WUwind, grid,
    metin->max_nearest_Uwind, metin->nearest_Uwind, met->z_Uwind);
 weightInterpField(metin->z_Vwind, metin->WVwind, grid,
    metin->max_nearest_Vwind, metin->nearest_Vwind, met->z_Vwind);
 weightInterpField(metin->z_Tair, metin->WTair, grid,
    metin->max_nearest_Tair, metin->nearest_Tair, met->z_Tair);
 weightInterpField(metin->z_RH, metin->WRH, grid,
    metin->max_nearest_RH, metin->nearest_RH, met->z_RH);
 if(VERBOSE>1 && myproc==0) printf("Done.\n");
  
} // End of InitialiseMetFields

/*
* Function: updateMetData()
* -----------------------------
* Main function for updating the met structure and interpolating onto the model time step 
*
*/
void updateMetData(propT *prop, gridT *grid, metinT *metin, metT *met, int myproc, MPI_Comm comm){
  
  int j,i,iptr, t0, t1, t2; 
  REAL dt, r1, r2;
   
  t1 = getTimeRec(prop->nctime,metin->time,metin->nt);

  if ( t1<0 ) {
    if(myproc==0) printf("NetCDF met data starts too late\n");
    t1=0;
  }
  if ( t1>=metin->nt ) {
    if(myproc==0) printf("NetCDF met data ends too early\n");
    t1=metin->nt-1;
  }

    /* Only interpolate the data onto the grid if need to*/
    if (metin->t1!=t1){
      if(VERBOSE>3 && myproc==0) printf("Updating netcdf variable at nc timestep: %d\n",t1);

      metin->t1=t1;
      // rather than fail with nan when the data is not long enough, this will
      // keep t0,t1,t2 valid, though it must be tested for in QuadInterp.
      if(t1>0) {
        metin->t0=t1-1;
      } else {
        metin->t0=t1;
      }
      if(t1+1<metin->nt) {
        metin->t2=t1+1;
      } else {
        metin->t2=t1;
      }

      /* Read in the data two time steps*/
      /* RH: This had been before the update above */
      ReadMetNC(prop, grid, metin, myproc);

      /* Interpolate the two time steps onto the grid*/
      weightInterpArray(metin->Uwind, metin->WUwind, grid,
	  metin->max_nearest_Uwind, metin->nearest_Uwind, NTmet, met->Uwind_t);
      weightInterpArray(metin->Vwind, metin->WVwind, grid,
	  metin->max_nearest_Vwind, metin->nearest_Vwind, NTmet, met->Vwind_t);
      weightInterpArray(metin->Tair, metin->WTair, grid,
	  metin->max_nearest_Tair, metin->nearest_Tair, NTmet, met->Tair_t);
      weightInterpArray(metin->Pair, metin->WPair, grid,
      	metin->max_nearest_Pair, metin->nearest_Pair, NTmet, met->Pair_t);
      weightInterpArray(metin->rain, metin->Wrain, grid,
	  metin->max_nearest_rain, metin->nearest_rain, NTmet, met->rain_t);
      weightInterpArray(metin->RH, metin->WRH, grid,
	  metin->max_nearest_RH, metin->nearest_RH, NTmet, met->RH_t);
      weightInterpArray(metin->cloud, metin->Wcloud, grid,
	  metin->max_nearest_cloud, metin->nearest_cloud, NTmet, met->cloud_t);
    }
    
    /* Do a linear temporal interpolation */
    //dt = metin->time[metin->t1]-metin->time[metin->t0];
    //r2 = (prop->nctime - metin->time[metin->t0])/dt;
    //r1 = 1.0-r2;
    
    t0=metin->t0;
    t1=metin->t1;
    t2=metin->t2;
     //printf("tmod = %f, tlow = %f (r1=%f), thigh = %f (r2=%f)\n",prop->nctime, metin->time[metin->t0],r1,metin->time[metin->t1],r2);
    
//    for (j=0;j<grid->Nc;j++){
 for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];  
    /*
      met->Uwind[j] = met->Uwind_t[0][j]*r1 + met->Uwind_t[1][j]*r2;
      met->Vwind[j] = met->Vwind_t[0][j]*r1 + met->Vwind_t[1][j]*r2;
      met->Tair[j] = met->Tair_t[0][j]*r1 + met->Tair_t[1][j]*r2;
      met->Pair[j] = met->Pair_t[0][j]*r1 + met->Pair_t[1][j]*r2;
      met->rain[j] = met->rain_t[0][j]*r1 + met->rain_t[1][j]*r2;
      met->RH[j] = met->RH_t[0][j]*r1 + met->RH_t[1][j]*r2;
      met->cloud[j] = met->cloud_t[0][j]*r1 + met->cloud_t[1][j]*r2;
    */
       //Quadratic temporal interpolation
       met->Uwind[i] = QuadInterp(prop->nctime,metin->time[t0],metin->time[t1],metin->time[t2],met->Uwind_t[0][i],met->Uwind_t[1][i],met->Uwind_t[2][i]);
       met->Vwind[i] = QuadInterp(prop->nctime,metin->time[t0],metin->time[t1],metin->time[t2],met->Vwind_t[0][i],met->Vwind_t[1][i],met->Vwind_t[2][i]);
       met->Tair[i] = QuadInterp(prop->nctime,metin->time[t0],metin->time[t1],metin->time[t2],met->Tair_t[0][i],met->Tair_t[1][i],met->Tair_t[2][i]);
       met->Pair[i] = QuadInterp(prop->nctime,metin->time[t0],metin->time[t1],metin->time[t2],met->Pair_t[0][i],met->Pair_t[1][i],met->Pair_t[2][i]);
       met->rain[i] = QuadInterp(prop->nctime,metin->time[t0],metin->time[t1],metin->time[t2],met->rain_t[0][i],met->rain_t[1][i],met->rain_t[2][i]);
       met->RH[i] = QuadInterp(prop->nctime,metin->time[t0],metin->time[t1],metin->time[t2],met->RH_t[0][i],met->RH_t[1][i],met->RH_t[2][i]);
       met->cloud[i] = QuadInterp(prop->nctime,metin->time[t0],metin->time[t1],metin->time[t2],met->cloud_t[0][i],met->cloud_t[1][i],met->cloud_t[2][i]);

      /* Place bounds on rain, humidity and cloud variables */
       if (met->cloud[i]<0.0) 
	 met->cloud[i]=0.0;
       if (met->cloud[i]>1.0)
	 met->cloud[i]=1.0;
       if ( (prop->metmodel!=5) && (met->rain[i]<0.0) ){
         // metmodel 5 allows for evaporation as negative rain.
	 met->rain[i]=0.0;
       }
       if (met->RH[i]<0.0)
	 met->RH[i]=0.0;
       if (met->RH[i]>100.0)
	 met->RH[i]=100.0;
    }

    //Communicate the arrays
    ISendRecvCellData2D(met->Uwind,grid,myproc,comm);
    ISendRecvCellData2D(met->Vwind,grid,myproc,comm);
    ISendRecvCellData2D(met->Tair,grid,myproc,comm);
    ISendRecvCellData2D(met->Pair,grid,myproc,comm);
    ISendRecvCellData2D(met->rain,grid,myproc,comm);
    ISendRecvCellData2D(met->RH,grid,myproc,comm);
    ISendRecvCellData2D(met->cloud,grid,myproc,comm);
} // End of updateMetData

/*
* Function: AllocateMet()
* -----------------------------
* Allocates memory to the meteorological structure array on the SUNTANS grid points
*
*/
void AllocateMet(propT *prop, gridT *grid, metT **met , int myproc){
  int j, k, n;
  int Nc = grid->Nc;
  
  if(VERBOSE>3 && myproc==0) printf("Allocating met structure...\n");
  *met = (metT *)SunMalloc(sizeof(metT),"AllocateMet");
  
  (*met)->z_Uwind = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->z_Vwind = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->z_Tair = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->z_RH = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  
  (*met)->Uwind = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->Vwind = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->Tair = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->Pair = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->rain = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->RH = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->cloud = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  
  (*met)->Hs = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->Hl = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->Hlw = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->Hsw = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->tau_x = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->tau_y = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->ustar = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->Tstar = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->qstar = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->EP = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->Htmp = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  (*met)->xtmp = (REAL *)SunMalloc(19*sizeof(REAL),"AllocateMet"); 

  (*met)->Uwind_t = (REAL **)SunMalloc(NTmet*sizeof(REAL *),"AllocateMet");
  (*met)->Vwind_t = (REAL **)SunMalloc(NTmet*sizeof(REAL *),"AllocateMet");
  (*met)->Tair_t = (REAL **)SunMalloc(NTmet*sizeof(REAL *),"AllocateMet");
  (*met)->Pair_t = (REAL **)SunMalloc(NTmet*sizeof(REAL *),"AllocateMet");
  (*met)->rain_t = (REAL **)SunMalloc(NTmet*sizeof(REAL *),"AllocateMet");
  (*met)->RH_t = (REAL **)SunMalloc(NTmet*sizeof(REAL *),"AllocateMet");
  (*met)->cloud_t = (REAL **)SunMalloc(NTmet*sizeof(REAL *),"AllocateMet");
  for(n=0;n<NTmet;n++){
      (*met)->Uwind_t[n] = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");   
      (*met)->Vwind_t[n] = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
      (*met)->Tair_t[n] = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
      (*met)->Pair_t[n] = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
      (*met)->rain_t[n] = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
      (*met)->RH_t[n] = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
      (*met)->cloud_t[n] = (REAL *)SunMalloc(Nc*sizeof(REAL),"AllocateMet");
  }
  
  for(k=0;k<Nc;k++){
      (*met)->z_Uwind[k] = 0.0;
      (*met)->z_Vwind[k] = 0.0;
      (*met)->z_Tair[k] = 0.0;
      (*met)->z_RH[k] = 0.0;
      
      (*met)->Uwind[k] = 0.0;
      (*met)->Vwind[k] = 0.0;
      (*met)->Tair[k] = 0.0;
      (*met)->Pair[k] = 0.0;
      (*met)->rain[k] = 0.0;
      (*met)->RH[k] = 0.0;
      (*met)->cloud[k] = 0.0;
      
      (*met)->Hs[k] = 0.0;
      (*met)->Hl[k] = 0.0;
      (*met)->Hlw[k] = 0.0;
      (*met)->Hsw[k] = 0.0;
      (*met)->tau_x[k] = 0.0;
      (*met)->tau_y[k] = 0.0;
      (*met)->ustar[k] = 0.0;
      (*met)->Tstar[k] = 0.0;
      (*met)->qstar[k] = 0.0;
      (*met)->EP[k] = 0.0;
      (*met)->Htmp[k] = 0.0;
      
      for(n=0;n<NTmet;n++){
	  (*met)->Uwind_t[n][k] = 0.0;
	  (*met)->Vwind_t[n][k] = 0.0;
	  (*met)->Tair_t[n][k] = 0.0;
	  (*met)->Pair_t[n][k] = 0.0;
	  (*met)->rain_t[n][k] = 0.0;
	  (*met)->RH_t[n][k] = 0.0;
	  (*met)->cloud_t[n][k] = 0.0;
      }  
  }
  for(j=0;j<19;j++){
    (*met)->xtmp[j]=0.0;
  }
} // End of function
  
/*
* Function: AllocateMetIn()
* -----------------------------
* Allocates memory to the meteorological input structure array
*
*/
void AllocateMetIn(propT *prop, gridT *grid, metinT **metin, int myproc){
  int j, k, n, retval;
  size_t NUwind;
  size_t NVwind;
  size_t NTair;
  size_t NPair;
  size_t Nrain;
  size_t NRH;
  size_t Ncloud;
  size_t nt;
  int Nc = grid->Nc;

  if(VERBOSE>3 && myproc==0) printf("Allocating metin structure...\n");
  *metin = (metinT *)SunMalloc(sizeof(metinT),"AllocateMetIn");

  /* Scalars */
  (*metin)->NUwind = returndimlen(prop->metncid,"NUwind");
  (*metin)->NVwind = returndimlen(prop->metncid,"NVwind");
  (*metin)->NTair = returndimlen(prop->metncid,"NTair");
  (*metin)->NPair = returndimlen(prop->metncid,"NPair");
  (*metin)->Nrain = returndimlen(prop->metncid,"Nrain");
  (*metin)->NRH = returndimlen(prop->metncid,"NRH");
  (*metin)->Ncloud = returndimlen(prop->metncid,"Ncloud");
  (*metin)->nt = returndimlen(prop->metncid,"nt");

  (*metin)->t0 = -1;
  (*metin)->t1 = -1;
  (*metin)->t2 = -1;

  NUwind = (*metin)->NUwind;
  NVwind = (*metin)->NVwind;
  NTair = (*metin)->NTair;
  NPair = (*metin)->NPair;
  Nrain = (*metin)->Nrain;
  NRH = (*metin)->NRH;
  Ncloud = (*metin)->Ncloud;
  nt = (*metin)->nt;

  if (VERBOSE>3){
      printf("NUwind = %d\n",(int)NUwind);
      printf("NVwind = %d\n",(int)NVwind); 
      printf("NTair = %d\n",(int)NTair);
      printf("NPair = %d\n",(int)NPair);
      printf("Nrain = %d\n",(int)Nrain);
      printf("NRH = %d\n",(int)NRH);
      printf("Ncloud = %d\n",(int)Ncloud);
      printf("nt = %d\n",(int)nt);
      printf("Nc = %d\n",(int)Nc);
  }
  /* Allocate the coordinate vectors*/
  //printf("Allocating coordinates...\n");
  (*metin)->x_Uwind = (REAL *)SunMalloc(NUwind*sizeof(REAL),"AllocateMetIn");
  (*metin)->x_Vwind = (REAL *)SunMalloc(NVwind*sizeof(REAL),"AllocateMetIn");
  (*metin)->x_Tair = (REAL *)SunMalloc(NTair*sizeof(REAL),"AllocateMetIn");
  (*metin)->x_Pair = (REAL *)SunMalloc(NPair*sizeof(REAL),"AllocateMetIn");
  (*metin)->x_rain = (REAL *)SunMalloc(Nrain*sizeof(REAL),"AllocateMetIn");
  (*metin)->x_RH = (REAL *)SunMalloc(NRH*sizeof(REAL),"AllocateMetIn");
  (*metin)->x_cloud = (REAL *)SunMalloc(Ncloud*sizeof(REAL),"AllocateMetIn");
  
  (*metin)->y_Uwind = (REAL *)SunMalloc(NUwind*sizeof(REAL),"AllocateMetIn");
  (*metin)->y_Vwind = (REAL *)SunMalloc(NVwind*sizeof(REAL),"AllocateMetIn");
  (*metin)->y_Tair = (REAL *)SunMalloc(NTair*sizeof(REAL),"AllocateMetIn");
  (*metin)->y_Pair = (REAL *)SunMalloc(NPair*sizeof(REAL),"AllocateMetIn");
  (*metin)->y_rain = (REAL *)SunMalloc(Nrain*sizeof(REAL),"AllocateMetIn");
  (*metin)->y_RH = (REAL *)SunMalloc(NRH*sizeof(REAL),"AllocateMetIn");
  (*metin)->y_cloud = (REAL *)SunMalloc(Ncloud*sizeof(REAL),"AllocateMetIn");

  (*metin)->z_Uwind = (REAL *)SunMalloc(NUwind*sizeof(REAL),"AllocateMetIn");
  (*metin)->z_Vwind = (REAL *)SunMalloc(NVwind*sizeof(REAL),"AllocateMetIn");
  (*metin)->z_Tair = (REAL *)SunMalloc(NTair*sizeof(REAL),"AllocateMetIn");
  (*metin)->z_RH = (REAL *)SunMalloc(NRH*sizeof(REAL),"AllocateMetIn"); 
  
  (*metin)->time = (REAL *)SunMalloc(nt*sizeof(REAL),"AllocateMetIn");
  
  /* Allocate the 2-D grid based weights */
  (*metin)->WUwind = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocateMetIn");
  (*metin)->WVwind = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocateMetIn");
  (*metin)->WTair = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocateMetIn");
  (*metin)->WPair = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocateMetIn");
  (*metin)->Wrain = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocateMetIn");
  (*metin)->WRH = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocateMetIn");
  (*metin)->Wcloud = (REAL **)SunMalloc(Nc*sizeof(REAL *),"AllocateMetIn");
  for(j=0;j<Nc;j++){
      (*metin)->WUwind[j] = (REAL *)SunMalloc(NUwind*sizeof(REAL),"AllocateMetIn");
      (*metin)->WVwind[j] = (REAL *)SunMalloc(NVwind*sizeof(REAL),"AllocateMetIn");
      (*metin)->WTair[j] = (REAL *)SunMalloc(NTair*sizeof(REAL),"AllocateMetIn");
      (*metin)->WPair[j] = (REAL *)SunMalloc(NPair*sizeof(REAL),"AllocateMetIn");
      (*metin)->Wrain[j] = (REAL *)SunMalloc(Nrain*sizeof(REAL),"AllocateMetIn");
      (*metin)->WRH[j] = (REAL *)SunMalloc(NRH*sizeof(REAL),"AllocateMetIn");
      (*metin)->Wcloud[j] = (REAL *)SunMalloc(Ncloud*sizeof(REAL),"AllocateMetIn");
  }
  
  /* Allocate the 2-D variable data (NTmet time steps)*/
  (*metin)->Uwind = (REAL **)SunMalloc(NTmet*sizeof(REAL *),"AllocateMetIn");
  (*metin)->Vwind = (REAL **)SunMalloc(NTmet*sizeof(REAL *),"AllocateMetIn");
  (*metin)->Tair = (REAL **)SunMalloc(NTmet*sizeof(REAL *),"AllocateMetIn");
  (*metin)->Pair = (REAL **)SunMalloc(NTmet*sizeof(REAL *),"AllocateMetIn");
  (*metin)->rain = (REAL **)SunMalloc(NTmet*sizeof(REAL *),"AllocateMetIn");
  (*metin)->RH = (REAL **)SunMalloc(NTmet*sizeof(REAL *),"AllocateMetIn");
  (*metin)->cloud = (REAL **)SunMalloc(NTmet*sizeof(REAL *),"AllocateMetIn");
  for(n=0;n<NTmet;n++){
      (*metin)->Uwind[n] = (REAL *)SunMalloc(NUwind*sizeof(REAL),"AllocateMetIn");   
      (*metin)->Vwind[n] = (REAL *)SunMalloc(NVwind*sizeof(REAL),"AllocateMetIn");
      (*metin)->Tair[n] = (REAL *)SunMalloc(NTair*sizeof(REAL),"AllocateMetIn");
      (*metin)->Pair[n] = (REAL *)SunMalloc(NPair*sizeof(REAL),"AllocateMetIn");
      (*metin)->rain[n] = (REAL *)SunMalloc(Nrain*sizeof(REAL),"AllocateMetIn");
      (*metin)->RH[n] = (REAL *)SunMalloc(NRH*sizeof(REAL),"AllocateMetIn");
      (*metin)->cloud[n] = (REAL *)SunMalloc(Ncloud*sizeof(REAL),"AllocateMetIn");
  }
  
  /* Initialises all of the input meteorological arrays with zeros*/ 
  // Need to allocate variable by variable as the lengths are different
  if(VERBOSE>2 && myproc==0) printf("Uwind, nj = %d, Nc = %d...\n",(int)(*metin)->NUwind,grid->Nc);
  for(j=0;j<(*metin)->NUwind;j++){
      (*metin)->x_Uwind[j]=0.0;
      (*metin)->y_Uwind[j]=0.0;
      (*metin)->z_Uwind[j]=0.0;
      for(k=0;k<Nc;k++){
	 (*metin)->WUwind[k][j]=0.0;
      }
       for(n=0;n<NTmet;n++){
	 (*metin)->Uwind[n][j]=0.0;
      }
  }  
  if(VERBOSE>2 && myproc==0) printf("Vwind, nj = %d, Nc = %d...\n",(int)(*metin)->NVwind,grid->Nc);
  for(j=0;j<(*metin)->NVwind;j++){
      (*metin)->x_Vwind[j]=0.0;
      (*metin)->y_Vwind[j]=0.0;
      (*metin)->z_Vwind[j]=0.0;
      for(k=0;k<Nc;k++){
	 (*metin)->WVwind[k][j]=0.0;
      }
       for(n=0;n<NTmet;n++){
	 (*metin)->Vwind[n][j]=0.0;
      }
  }  
  if(VERBOSE>2 && myproc==0) printf("Tair, nj = %d, Nc = %d...\n",(int)(*metin)->NTair,grid->Nc);
  for(j=0;j<(*metin)->NTair;j++){
      (*metin)->x_Tair[j]=0.0;
      (*metin)->y_Tair[j]=0.0;
      (*metin)->z_Tair[j]=0.0;
      for(k=0;k<Nc;k++){
	 (*metin)->WTair[k][j]=0.0;
      }
       for(n=0;n<NTmet;n++){
	 (*metin)->Tair[n][j]=0.0;
      }
  }  
  if(VERBOSE>2 && myproc==0) printf("Pair, nj = %d, Nc = %d...\n",(int)(*metin)->NPair,grid->Nc);
  for(j=0;j<(*metin)->NPair;j++){
      (*metin)->x_Pair[j]=0.0;
      (*metin)->y_Pair[j]=0.0;
      for(k=0;k<Nc;k++){
	 (*metin)->WPair[k][j]=0.0;
      }
       for(n=0;n<NTmet;n++){
	 (*metin)->Pair[n][j]=0.0;
      }
  }  
  if(VERBOSE>2 && myproc==0) printf("rain, nj = %d, Nc = %d...\n",(int)(*metin)->Nrain,grid->Nc);
  for(j=0;j<(*metin)->Nrain;j++){
      (*metin)->x_rain[j]=0.0;
      (*metin)->y_rain[j]=0.0;
      for(k=0;k<Nc;k++){
	 (*metin)->Wrain[k][j]=0.0;
      }
       for(n=0;n<NTmet;n++){
	 (*metin)->rain[n][j]=0.0;
      }
  }  
  if(VERBOSE>2 && myproc==0) printf("RH, nj = %d, Nc = %d...\n",(int)(*metin)->NRH,grid->Nc);
  for(j=0;j<(*metin)->NRH;j++){
      (*metin)->x_RH[j]=0.0;
      (*metin)->y_RH[j]=0.0;
      (*metin)->z_RH[j]=0.0;
      for(k=0;k<Nc;k++){
	 (*metin)->WRH[k][j]=0.0;
      }
       for(n=0;n<NTmet;n++){
	 (*metin)->RH[n][j]=0.0;
      }
  }  
  if(VERBOSE>2 && myproc==0) printf("cloud, nj = %d, Nc = %d...\n",(int)(*metin)->Ncloud,grid->Nc);
  for(j=0;j<(*metin)->Ncloud;j++){
      (*metin)->x_cloud[j]=0.0;
      (*metin)->y_cloud[j]=0.0;
      for(k=0;k<Nc;k++){
	 (*metin)->Wcloud[k][j]=0.0;
      }
       for(n=0;n<NTmet;n++){
	 (*metin)->cloud[n][j]=0.0;
      }
  }
  if(VERBOSE>2 && myproc==0) printf("time, nt = %d ...\n",(int)(*metin)->nt);
  for(j=0;j<nt;j++){
      (*metin)->time[j]=0.0;
  }
}


/*
* Function calcInterpWeights()
* ----------------------------
* Calculates the interpolation weights for all grid points based on "Ns" interpolants
* at cooridinates (xo, yo)
*/
void calcInterpWeights(gridT *grid, propT *prop, REAL *xo, REAL *yo, int Ns, int **index, REAL **klambda, int myproc){
    
    int j, i, jj, ii, iptr, jptr;
    int Nc = grid->Nc;
    const REAL inversepower = 1.2; // RH trying something smoother than original 2.2;
    REAL sumgamma, dist, tmp;
    REAL *gamma;
    REAL **C, **Ctmp;
    // Allocate the arrays
    if(prop->varmodel==0){
      gamma = (REAL *)SunMalloc(Ns*sizeof(REAL),"CalcInterpWeights");
    }else{
	gamma = (REAL *)SunMalloc((Ns+1)*sizeof(REAL),"CalcInterpWeights");
	C = (REAL **)SunMalloc((Ns+1)*sizeof(REAL),"CalcInterpWeights");
	Ctmp = (REAL **)SunMalloc((Ns+1)*sizeof(REAL),"CalcInterpWeights");
	for (j=0;j<Ns+1;j++){
	  C[j] = (REAL *)SunMalloc((Ns+1)*sizeof(REAL),"CalcInterpWeights");
	  Ctmp[j] = (REAL *)SunMalloc((Ns+1)*sizeof(REAL),"CalcInterpWeights");
	}
    }
    
    if(prop->varmodel==0){ // Inverse distance weighting
      if(VERBOSE>1 && myproc==0) printf("Calculating interpolation weights using inverse distance weighting...\n");     
      for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
	i = grid->cellp[iptr];
        // approximate min length scale for inverse distance
        tmp=sqrt(grid->Ac[i]);
	sumgamma=0.0;
	for(j=0;j<Ns;j++){
          dist = pow(grid->xv[i]-xo[index[i][j]],2) + pow(grid->yv[i]-yo[index[i][j]],2);
          gamma[j] = 1.0/pow(dist+tmp,inversepower);
          sumgamma += gamma[j];
	}
	for(j=0;j<Ns;j++){
          klambda[i][j] = gamma[j]/sumgamma;
	}
      }
      
      //SunFree(gamma,Ns,"CalcInterpWeights");
    
    }else{ // kriging
	if(VERBOSE>1 && myproc==0)  printf("Calculating interpolation weights using kriging...\n");
	
	// Loop through each model grid point and
	for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
	  i = grid->cellp[iptr];  

	  // Construct the LHS Matrix C
	  for(ii=0;ii<Ns+1;ii++){
	    for(j=0;j<Ns+1;j++){
	      C[ii][j]=1.0;
	    }
	  }
	  for(ii=0;ii<Ns;ii++){
	    C[ii][ii]=semivariogram(prop->varmodel, prop->nugget, prop->sill, prop->range, 0.0);
	    for(j=ii+1;j<Ns;j++){
	      dist = sqrt(  pow(xo[index[i][ii]]-xo[index[i][j]],2) + pow(yo[index[i][ii]]-yo[index[i][j]],2) );
	      C[ii][j] = semivariogram(prop->varmodel, prop->nugget, prop->sill, prop->range, dist);
	      C[j][ii]=C[ii][j];
	    }
	  }
	  C[Ns][Ns]=0.0;


	  // calculate the  weights
	  for(jj=0;jj<Ns;jj++){
	    dist = sqrt( pow(grid->xv[i]-xo[index[i][jj]],2) + pow(grid->yv[i]-yo[index[i][jj]],2) );
	    gamma[jj] = semivariogram(prop->varmodel, prop->nugget, prop->sill, prop->range, dist);
	  }
	  gamma[Ns]=1.0;
	  
	  // Solve the linear system
	  for(ii=0;ii<Ns+1;ii++){
	    for(jj=0;jj<Ns+1;jj++){
	      Ctmp[ii][jj]=C[ii][jj];
	    }
	  }
	  linsolve(Ctmp,gamma,Ns+1);
	  
	  // Write to the weights array
	  for(jj=0;jj<Ns;jj++){
	    klambda[i][jj] = gamma[jj];
	  }
	  
	  // Check the weights
          sumgamma=0.0;
          // printf("W[j]:\n");
          for(jj=0;jj<Ns;jj++){ // don't include the last point
            sumgamma+=gamma[jj];
            //printf("%3.6f ",gamma[jj]);
          }
          if ( fabs(sumgamma - 1.0) > 0.001 ) {
            printf("\nSumW = %f (should equal 1.000)\n",sumgamma);
	  }
          if (sumgamma !=sumgamma ) {
            printf("\nSumW = %f (should equal 1.000)\n",sumgamma);
          }
	} 
	
      // Free up the arrays
      /*
      SunFree(gamma,Ns+1,"CalcInterpWeights");
      for (j=0;j<Ns+1;j++){
	SunFree(C[j],Ns+1,"CalcInterpWeights");
	SunFree(Ctmp[j],Ns+1,"CalcInterpWeights");
      }
      */
    }// end of kriging

    if (VERBOSE) {
      printf("[p=%d] End of calcInterpWeights, Ns=%d total memory used: %lu\n",
	     myproc,Ns,TotSpace);
    }
} // End of calcInterpWeights

/*
* Function: semivariogram()
* -------------------------
* Calculates the semivariogram function used by the kriging interpolation scheme
*
*/
static REAL semivariogram(int varmodel, REAL nugget, REAL sill, REAL range, REAL D){
  
  REAL tmp;
 if (varmodel==1){
   // Spherical model
   if(D > range){
     return sill;
   }else{
      tmp = D/range;
      return (nugget + (sill - nugget) * (1.5*tmp - 0.5 * pow(tmp,3))); 
   }
 }  
 return 0.0;
}

/*
* Function: weightInterpArray()
* -----------------------------
* Perform weighted interpolation on a field D [2d-array]
*
*/
void weightInterpArray(REAL **D, REAL **klambda, gridT *grid, int Ns, int **index, int nt, REAL **Dout){
  int i,j, k, iptr;
  for(k=0;k<nt;k++){
    //for(i=0;i<Nc;i++){
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
        i = grid->cellp[iptr];  
	Dout[k][i] = 0.0;
	for(j=0;j<Ns;j++){
	    Dout[k][i] += klambda[i][j] * D[k][index[i][j]];
	}
    }
  }
}

/* 
* Function: weightInterpField()
* -----------------------------
* Perform weighted interpolation on a field D [vector]
*
*/
void weightInterpField(REAL *D, REAL **klambda, gridT *grid, int Ns, int **index, REAL *Dout){
  int i,j,iptr;
  
//  for(i=0;i<Nc;i++){
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];  
      Dout[i] = 0.0;
      for(j=0;j<Ns;j++){
	  Dout[i] += klambda[i][j] * D[index[i][j]];
      }
  }
}

/* 
* Function: FindNearestMetStations()
* -----------------------------
* Find the N Nearest met stations to each grid cell
*
*/
void FindNearestMetStations(propT *prop, gridT *grid, metinT **metin, int myproc){

    int Nc = grid->Nc;
    int i,iptr;

    // Determine the maximum points for each variables
    (*metin)->max_nearest_Uwind =(int)Min((REAL)MAXNEAR, (REAL)(*metin)->NUwind);
    (*metin)->max_nearest_Vwind =(int)Min((REAL)MAXNEAR, (REAL)(*metin)->NVwind);
    (*metin)->max_nearest_Tair = (int)Min((REAL)MAXNEAR, (REAL)(*metin)->NTair);
    (*metin)->max_nearest_Pair = (int)Min((REAL)MAXNEAR, (REAL)(*metin)->NPair);
    (*metin)->max_nearest_RH = (int)Min((REAL)MAXNEAR, (REAL)(*metin)->NRH);
    (*metin)->max_nearest_rain = (int)Min((REAL)MAXNEAR, (REAL)(*metin)->Nrain);
    (*metin)->max_nearest_cloud =(int)Min((REAL)MAXNEAR, (REAL)(*metin)->Ncloud);

    // Allocate the indexing arrays
    (*metin)->nearest_Uwind = (int **)SunMalloc(Nc*sizeof(int *),"AllocateMetIn");
    (*metin)->nearest_Vwind = (int **)SunMalloc(Nc*sizeof(int *),"AllocateMetIn");
    (*metin)->nearest_Tair = (int **)SunMalloc(Nc*sizeof(int *),"AllocateMetIn");
    (*metin)->nearest_Pair = (int **)SunMalloc(Nc*sizeof(int *),"AllocateMetIn");
    (*metin)->nearest_rain = (int **)SunMalloc(Nc*sizeof(int *),"AllocateMetIn");
    (*metin)->nearest_RH = (int **)SunMalloc(Nc*sizeof(int *),"AllocateMetIn");
    (*metin)->nearest_cloud = (int **)SunMalloc(Nc*sizeof(int *),"AllocateMetIn");
    for(i=0;i<Nc;i++){
        (*metin)->nearest_Uwind[i] = (int *)SunMalloc((*metin)->max_nearest_Uwind*sizeof(int),"AllocateMetIn");
        (*metin)->nearest_Vwind[i] = (int *)SunMalloc((*metin)->max_nearest_Vwind*sizeof(int),"AllocateMetIn");
        (*metin)->nearest_Tair[i] = (int *)SunMalloc((*metin)->max_nearest_Tair*sizeof(int),"AllocateMetIn");
        (*metin)->nearest_Pair[i] = (int *)SunMalloc((*metin)->max_nearest_Pair*sizeof(int),"AllocateMetIn");
        (*metin)->nearest_rain[i] = (int *)SunMalloc((*metin)->max_nearest_rain*sizeof(int),"AllocateMetIn");
        (*metin)->nearest_RH[i] = (int *)SunMalloc((*metin)->max_nearest_RH*sizeof(int),"AllocateMetIn");
        (*metin)->nearest_cloud[i] = (int *)SunMalloc((*metin)->max_nearest_cloud*sizeof(int),"AllocateMetIn");
    }
    
    // Go through and find the N nearest points for each variable
    //FindNearest(int *points, REAL *x, REAL *y, int N, int np, REAL xi, REAL yi)
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
	i = grid->cellp[iptr];  

	FindNearest((*metin)->nearest_Uwind[i], (*metin)->x_Uwind, (*metin)->y_Uwind,
		(*metin)->NUwind, (*metin)->max_nearest_Uwind, grid->xv[i], grid->yv[i]);
	FindNearest((*metin)->nearest_Vwind[i], (*metin)->x_Vwind, (*metin)->y_Vwind,
		(*metin)->NVwind, (*metin)->max_nearest_Vwind, grid->xv[i], grid->yv[i]);
	FindNearest((*metin)->nearest_Tair[i], (*metin)->x_Tair, (*metin)->y_Tair,
		(*metin)->NTair, (*metin)->max_nearest_Tair, grid->xv[i], grid->yv[i]);
	FindNearest((*metin)->nearest_Pair[i], (*metin)->x_Pair, (*metin)->y_Pair,
		(*metin)->NPair, (*metin)->max_nearest_Pair, grid->xv[i], grid->yv[i]);
	FindNearest((*metin)->nearest_RH[i], (*metin)->x_RH, (*metin)->y_RH,
		(*metin)->NRH, (*metin)->max_nearest_RH, grid->xv[i], grid->yv[i]);
	FindNearest((*metin)->nearest_rain[i], (*metin)->x_rain, (*metin)->y_rain,
		(*metin)->Nrain, (*metin)->max_nearest_rain, grid->xv[i], grid->yv[i]);
	FindNearest((*metin)->nearest_cloud[i], (*metin)->x_cloud, (*metin)->y_cloud,
		(*metin)->Ncloud, (*metin)->max_nearest_cloud, grid->xv[i], grid->yv[i]);
    
    }


}// End function

/*
* Function: updateAirSeaFluxes()
* ------------------------------
* Main routine for calculating the air-sea heat and salt fluxes
*
* Computed terms are stored in the met structure array
*
* These routines are activated when "metmodel" >= 2 in suntans.dat
* (note this gets called from HeatSources *and* from phys.c)
*/ 
void updateAirSeaFluxes(propT *prop, gridT *grid, physT *phys, metT *met,REAL **T){
  
  int i, ktop, iptr, n;
  int Nc = grid->Nc;
  REAL *x=met->xtmp; // pointer to vector passed to cor30
  REAL Umag; // Wind Speed magnitude
  // Constant flux coefficients
  REAL Cd = prop->Cda; // Drag coefficient
  REAL Ch = prop->Ch; // Stanton Number
  REAL Ce = prop->Ce; //Dalton Number
  REAL cp = 4186.0; // Specific heat of water
  REAL cpa = 1004.67; // Specific heat of air
  REAL Lv = 2.50e6; // Latent heat of vaporization
  REAL rhoa = 1.20; // Density of air	
  

  /* Terms in the COARE3.0 input vector
   REAL u=x[0]; //wind speed (m/s]  at height zu [m]
   REAL us=x[1]; //surface current speed in the wind direction [m/s]
   REAL ts=x[2]; //bulk water temperature [C] if jcool=1, interface water T if jcool=0  
   REAL t=x[3]; //bulk air temperature [C], height zt
   REAL Qs=x[4]; //bulk water spec hum [kg/kg] if jcool=1, ...
   REAL Q=x[5]; //bulk air spec hum [kg/kg], height zq
   REAL Rs=x[6]; //downward solar flux [W/m^2]
   REAL Rl=x[7]; //downard IR flux [W/m^2]
   REAL rain=x[8]; //rain rate [mm/hr]
   REAL zi=x[9]; //PBL depth [m]
   REAL P=x[10]; //Atmos surface pressure [mb]
   REAL zu=x[11]; //wind speed measurement height [m]
   REAL zt=x[12]; //air T measurement height [m]
   REAL zq=x[13]; //air q measurement height [m]
   REAL lat=x[14]; //latitude [deg, N=+]
   REAL jcool=x[15]; //implement cool calculation skin switch, 0=no, 1=yes
   REAL jwave=x[16]; //implement wave dependent roughness model
   REAL twave=x[17]; //wave period [s]
   REAL hwave=x[18]; //wave height [m]
  */
  
  /* The goal is fill up vector x   */ 
  //x = (REAL *)SunMalloc(19*sizeof(REAL),"updateAirSeaFluxes");

  
  //if(myproc==0) printf(" j, Hs, Hl, tau, Hlw, Hsw\n"); 
  // for(j=0;j<Nc;j++){
 for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
    i = grid->cellp[iptr];  
    ktop = grid->ctop[i];
    // Wind speed
    //Umag = sqrt( pow( (met->Uwind[i]-phys->uc[i][ktop]) ,2) + pow( (met->Vwind[i]-phys->vc[i][ktop]),2) );
    Umag = sqrt( pow(met->Uwind[i],2) + pow(met->Vwind[i],2) );
    x[0] = Umag;

    // Surface current speed in wind direction
    // This is the projection of the water velocity vector onto the wind velocity vector
    //x[1] = 0.0;
    if ( x[0] > 1e-10 ) {
      x[1] = fabs(phys->uc[i][ktop]*met->Uwind[i]/Umag + phys->vc[i][ktop]*met->Vwind[i]/Umag);
    } else {
      // if Umag is too small, just use wind magnitude
      x[1] = sqrt(pow(phys->uc[i][ktop],2)+pow(phys->vc[i][ktop],2));
    }

    // Water temperature
    x[2] = T[i][ktop];
    
    // Air temperature
    x[3] = met->Tair[i];
    
    // Water specific humidty
    x[4] = qsat(T[i][ktop], met->Pair[i]);
    
    // Air specific humidity
    x[5] = specifichumidity(met->RH[i],met->Tair[i],met->Pair[i]);
    
    // Longwave radiation
    met->Hlw[i] = longwave(met->Tair[i],T[i][ktop],met->cloud[i]);
    //met->Hlw[i] = longwave_berliand(met->Tair[i],T[i][ktop],met->cloud[i],met->RH[i],met->Pair[i]);
    x[6] = met->Hlw[i];
    
    // Shortwave radiation
    met->Hsw[i] = shortwave(prop->nctime/86400.0,prop->latitude,met->cloud[i],prop->gmtoffset/24.0);
    x[7] = met->Hsw[i];
    
    //rain [mm/hr] (rain heat flux is not included at the moment)
    x[8] = met->rain[i]*3600;
    
    //Air pressure [mb]
    x[10] = met->Pair[i];
    
    //wind speed height
    x[11] = met->z_Uwind[i];
    
    //air temp height
    x[12] = met->z_Tair[i];
    
    //humidity measurement height
    x[13] = met->z_RH[i];
    
     //Set some constant values
     x[9] = 600.0; 
     x[14] = prop->latitude;
     x[15] = 0.0;
     x[16] = 0.0;
     x[17] = 0.0;
     x[18] = 0.0;
    //printf(" %d, x[0]:%6.6f, x[1]:%6.6f, x[2]:%6.6f, x[3]:%6.6f, x[4]:%6.6f, x[5]:%6.6f\n",j,x[0],x[1],x[2],x[3],x[4],x[5]);
    if (prop->metmodel==2){
      /* Calculate the actual fluxes */
      cor30a(x);
      /* Output variable back into the input array 
      x[0] = hsb;
      x[1] = hlb;
      x[2] = tau;
      x[3] = zo;
      x[4] = zot;
      x[5] = zoq;
      x[6] = rhoa;
      x[7] = usr;
      x[8] = tsr;
      x[9] = qsr;
      x[10] = dter;
      x[11] = dqer;
      x[12] = tkt;
      x[13] = RF;
      x[14] = ut; // Speed including gustiness
      x[15] = Cd;
      x[16] = Ch;
      x[17] = Ce;
      x[18] = ug;
      */
      // Output the fluxes to the met structure array
      met->Hs[i] = -x[0];//Note the change of sign
      met->Hl[i] = -x[1];
      met->ustar[i] = x[7];//These are not used at present
      met->Tstar[i] = x[8];
      met->qstar[i] = x[9];
    
      /* Calculate the wind stress components
      * tau_x = rhoa * Cd * S * (Ucurrent - Uwind) : Fairall et al, 1996
      */
      // *** I think this is double counting if the surface currents have already been accounted for above***
      //met->tau_x[j] = x[6] * x[15] * x[14] * (phys->uc[j][ktop] - met->Uwind[j]);
      //met->tau_y[j] = x[6] * x[15] * x[14] * (phys->vc[j][ktop] - met->Vwind[j]);
      //met->tau_x[j] = x[6] * x[15] * x[14] * (met->Uwind[j] - phys->uc[j][ktop]);
      //met->tau_y[j] = x[6] * x[15] * x[14] * (met->Vwind[j] - phys->vc[j][ktop]);
      // No gust speed in stress term
      met->tau_x[i] = x[6] * x[10] * Umag * (met->Uwind[i] - phys->uc[i][ktop]);
      met->tau_y[i] = x[6] * x[10] * Umag * (met->Vwind[i] - phys->vc[i][ktop]);

      //No surface current dependence 
      //met->tau_x[j] = x[6] * x[15] * x[14] * met->Uwind[j];
      //met->tau_y[j] = x[6] * x[15] * x[14] * met->Vwind[j];
      //met->tau_x[j] = 1.2 * 0.0011 * x[14] * met->Uwind[j];
      //met->tau_y[j] = 1.2 * 0.0011 * x[14] * met->Vwind[j];
      //printf("%10.6f, %10.6f, %10.6f, %10.6f\n",x[6],x[15],x[14],met->Vwind[j]);

    } else if(prop->metmodel>=3){// Compute fluxes with constant parameters
      met->Hs[i] = - rhoa * cpa * Ch * Umag * (x[2] - x[3]);//T_w > T_a -> Hs is negative
      met->Hl[i] = - rhoa * Lv * Ce * Umag * (x[4] - x[5]);

      // Large and Pond, 1981
      if(Umag<11.0) Cd=1.2e-3;
      else Cd=0.49e-3 + 0.065e-3*Umag;
      met->tau_x[i] = rhoa * Cd * Umag * (met->Uwind[i] - phys->uc[i][ktop]); 
      met->tau_y[i] = rhoa * Cd * Umag * (met->Vwind[i] - phys->vc[i][ktop]);
    }

    // Check for nans and dump the inputs
    for(n=0;n<8;n++){
    	if(x[n]!=x[n]){
	   printf("Error in COARE3.0 Algorithm at i = %d, x[%d] = nan.\n",i,n);
	   printf("Uwind[%d] = %6.10f, z_Uwind = %6.10f m\n",i,met->Uwind[i],met->z_Uwind[i]);
	   printf("Vwind[%d] = %6.10f, z_Vwind = %6.10f m\n",i,met->Vwind[i],met->z_Vwind[i]);
	   printf("Tair[%d] = %6.10f, z_Tair = %6.10f m\n",i,met->Tair[i],met->z_Tair[i]);
	   printf("Pair[%d] = %6.10f\n",i,met->Pair[i]);
	   printf("rain[%d] = %6.10f\n",i,met->rain[i]);
	   printf("RH[%d] = %6.10f, z_RH = %6.10f m\n",i,met->RH[i],met->z_RH[i]);
	   printf("cloud[%d] = %6.10f\n",i,met->cloud[i]);
	   printf("T[%d][%d] = %6.10f\n",i,ktop,T[i][ktop]);
           printf(" cell location: %10.4f %10.4f\n",grid->xv[i], grid->yv[i]);
           
	   MPI_Finalize();
	   exit(EXIT_FAILURE);
	}
    }
  }
  //SunFree(x,19,"updateAirSeaFluxes");
} // End updateAirFluxes


/* 
* Function: specifichumidity()
* ----------------------------
* Convert relative humidity (%) to specific humidity (kg/kg 
*/
static REAL specifichumidity(REAL RH, REAL Ta, REAL Pair){
 
 //REAL cff;
 ////REAL eps=1e-8;
 ///*
 // * Compute air saturation vapor pressure (mb), using Teten formula.
 // */
 // cff=(1.0007+3.46E-6*Pair)*6.1121*exp(17.502*Ta/(240.97+Ta+SMALL));
 // 
 // /*
 // *  Compute specific humidity at Saturation, Qair (kg/kg).
 // */
 // cff=cff*RH/100;                    // Vapor pres (mb)
 // return (0.62197*(cff/(Pair-0.378*cff+SMALL))); // Spec hum (kg/kg)
 return RH*0.01*qsat(Ta,Pair);
} // End specifichumidity


/*
* Function: qsat()
 -----------------
* Compute water saturation vapor pressure (mb), using Teten formula.
*/
static REAL qsat(REAL Tw, REAL Pair){
 REAL cff;
 //REAL eps=1e-10;

 cff=(1.0007+3.46E-6*Pair)*6.1121* exp(17.502*Tw/(240.97+Tw+SMALL));

  //  Vapor Pressure reduced for salinity (Kraus & Businger, 1994, pp 42).
  cff=cff*0.98;
 
  return (0.62197*(cff/(Pair-0.378*cff+SMALL)));
} // End qsat

/*
 * Function: satvap()
 * --------------------
 *  Calculate the saturation vapor pressure
 *
 *  Taken from the matlab airsea toolbox
 */
static REAL satvap(REAL Ta, REAL Pair){
    REAL ew, fw;

    ew = pow( (0.7859+0.03477*Ta) / (1.+0.00412*Ta),10.);
    fw = 1.0 + 1e-6*Pair*(4.5+6e-4*pow(Ta,2));

    return fw*ew;
}
/*
 * Function: vappres()
 * --------------------
 *  Calculate the  vapor pressure from relative humidity
 *
 *  Taken from the matlab airsea toolbox
 */
static REAL vappres(REAL Ta, REAL RH, REAL Pair){
    REAL ew, rw, r;
    const REAL eps_air = 0.62197; // molecular weight ratio (water/air)

    ew = satvap(Ta,Pair);
    rw = eps_air*ew / (Pair-ew);
    r = (RH*0.01) * rw;

    return r*Pair / (eps_air+r);

}
/*
* Function: longwave_berliand()()
* --------------------
* Calculate net longwave radiation into water
*
* Uses the Berliand (1952) bulk formula
*
* Taken from the matlab air sea toolbox
* 
*/
static REAL longwave_berliand(REAL Ta, REAL Tw, REAL C_cloud, REAL RH, REAL Pair){
  
  // Constants
  const REAL T_ref = 273.16;             // conversion from C to K
  const REAL sigma = 5.6697e-8;         // Boltzmann constant (W m^{-2} K^{-4})
  const REAL epsilon_w = 0.985;           // emissivity of water
  
  REAL e_a;

  e_a = vappres(Ta,RH,Pair);
  
  return -epsilon_w*sigma*pow(Ta+T_ref,4.) * (0.39 - 0.05*sqrt(e_a)) * C_cloud -
  	4*epsilon_w*sigma*pow(Ta+T_ref,3.)*(Tw-Ta);
  
}

/*
* Function: longwave()
* --------------------
* Calculate net longwave radiation into water
*
* Ref: Martin and McCutcheon, "Hydrodynamics and Transport for Water
* Quality Modeling", 1999
*/
static REAL longwave(REAL Ta, REAL Tw, REAL C_cloud){
  
  // Constants
  const REAL T_ref = 273.16;             // conversion from C to K
  const REAL sigma = 5.67e-8;         // Boltzmann constant (W m^{-2} K^{-4})
  const REAL alpha_0 = 0.937e-5;         // Proportionality constants ([alpha_0]= K^{-2}, all others dimensionless)
  const REAL alpha_LW = 0.17;            // LW cloud cover fraction coefficient
  const REAL r_LW = 0.03;                // Fraction of longwave radiation reflected at surface
  const REAL epsilon_w = 0.97;           // emissivity of water
  
  REAL H_LE, epsilon_a, H_LW;
  
  // Emitted Long Wave Radiation
  H_LE = epsilon_w*sigma*pow(Tw + T_ref,4);
  
  //  Incoming Long Wave Radiation
  
  //emissivity of air
  epsilon_a = alpha_0*(1+alpha_LW*C_cloud)*pow(Ta + T_ref,2);
  
  H_LW = epsilon_a*sigma*(1-r_LW)*pow(Ta+T_ref,4);
  
  return  (H_LW - H_LE);

} // End longwave

/*
* Function: shortwave()
* ---------------------
* Compute solar radiation flux using the Gill, 1982 formulae 
*
*/
REAL shortwave(REAL time, REAL Lat,REAL C_cloud, REAL toffset){
  
  const REAL S=1368.0;  //[W m-2]
  const REAL albedo = 0.06;
  const REAL R = 0.76;
  REAL omega1, omega0;
  REAL delta, singamma, Qsc, Hsw;
  REAL gmttime;

  omega1 = 2*PI/1.0; // Diurnal cycle
  omega0 = 2*PI/365.25; // Annual cycle

  // Remove the gmtoffset from the time (assume units are the same)
  gmttime = time + toffset;

  delta = 23.5*PI/180 * cos(omega0*gmttime - 2.95);
  singamma = sin(delta)*sin(PI*Lat/180) - cos(delta)*cos(PI*Lat/180)*cos(omega1*gmttime);

  // Clear sky radiation
  if(singamma >= 0.0){
      Qsc = R*S*singamma;
  }else{
      Qsc = 0;
  }
  // Corrected radiation
  Hsw = (1.0-albedo) * (1.0 - 0.65 * pow(C_cloud,2)) * Qsc;

  return Hsw;
  
} // End shortwave

/* Function: bulkflux()
 * ------------------------
 * Calculates air-sea fluxes using bulk flux formulation
 *
 * Adapted from the matlab air-sea toolbox
 *
 */

//static void bulkflux(REAL *y){
//    // Input variables
//    REAL ur = y[0]; //wind speed [m/s] measured at height zr [m] 
//    REAL zr = y[1];
//    REAL Ta = y[2]; //   = air temperature [C] measured at height zt [m]
//    REAL zt = y[3];
//    REAL rh = y[4]; //   = relative humidity [%] measured at height zq [m]
//    REAL zq = y[5];
//    REAL Pa = y[6]; //   = air pressure [mb]
//    REAL Ts = y[7]; //   = sea surface temperature [C]
//
//
//}

/* 
* Function:  cor30a()
* -------------------
* Calculates air-sea fluxes using bulk flux formulation
*
* References:
*	Fairall et al, 1996, JGR
*  	Fairall et al, 2003, Journal of Climate
*
* See:
*	http://coaps.fsu.edu/COARE/flux_algor/
*
* This code is adapted from cor30.m matlab function
*/
static void cor30a(REAL *y){
  
  /* Get variables from the input vector */
  REAL u=y[0]; //wind speed (m/s]  at height zu [m]
  REAL us=y[1]; //surface current speed in the wind direction [m/s]
  REAL ts=y[2]; //bulk water temperature [C] if jcool=1, interface water T if jcool=0  
  REAL t=y[3]; //bulk air temperature [C], height zt
  REAL Qs=y[4]; //bulk water spec hum [kg/kg] if jcool=1, ...
  REAL Q=y[5]; //bulk air spec hum [kg/kg], height zq
  REAL Rs=y[6]; //downward solar flux [W/m^2]
  REAL Rl=y[7]; //downard IR flux [W/m^2]
  REAL rain=y[8]; //rain rate [mm/hr]
  REAL zi=y[9]; //PBL depth [m]
  REAL P=y[10]; //Atmos surface pressure [mb]
  REAL zu=y[11]; //wind speed measurement height [m]
  REAL zt=y[12]; //air T measurement height [m]
  REAL zq=y[13]; //air q measurement height [m]
  REAL lat=y[14]; //latitude [deg, N=+]
  REAL jcool=y[15]; //implement cool calculation skin switch, 0=no, 1=yes
  REAL jwave=y[16]; //implement wave dependent roughness model
  REAL twave=y[17]; //wave period [s]
  REAL hwave=y[18]; //wave height [m]
  
  int nits =3;
  int i;
  
  REAL Beta, von, fdg, tdk, grav;
  REAL Rgas, Le, cpa, cpv, rhoa, visa;
  REAL Al, be, cpw, rhow, visw, tcw, bigc, wetc;
  REAL lwave, cwave, Rns, Rnl;
  REAL dt, du, dq, ta, ug, dter, dqer;
  REAL u10, usr, zo10, Cd10, Ch10, Ct10, zot10;
  REAL Cd, Ct, CC, Ribcu, Ribu, zetu, L10, tsr, qsr, tkt, charn, rr, Bf;
  REAL hsb, hlb, qout, dels, qcol, alq;
  REAL xlamx, dwat, dtmp, alfac, RF, wbar, hl_webb, zet;
  REAL zo, zot, zoq, L, ut, tau, Ch, Ce, Chn_10, Cdn_10, Cen_10;
  
  //printf("y[13]:%6.6f, y[14]:%6.6f, y[15]:%6.6f, y[16]:%6.6f, y[17]:%6.6f, y[18]:%6.6f\n",y[13],y[14],y[15],y[16],y[17],y[18]);
  
/**********   set constants *************/
  Beta=1.2;
  von=0.4;
  fdg=1.00;
  tdk=273.16;
  grav=9.81;
  /*************  air constants ************/
  Rgas=287.1;
  Le=(2.501-.00237*ts)*1e6;
  cpa=1004.67;
  cpv=cpa*(1+0.84*Q);
  rhoa=P*100/(Rgas*(t+tdk)*(1+0.61*Q));
  visa=1.326e-5*(1+6.542e-3*t+8.301e-6*t*t-4.84e-9*t*t*t);
   /************  cool skin constants  *******/
   Al=2.1e-5*pow(ts+3.2,0.79);
   be=0.026;
   cpw=4000;
   rhow=1022;
   visw=1e-6;
   tcw=0.6;
   bigc=16*grav*cpw*pow(rhow*visw,3.0)/(tcw*tcw*rhoa*rhoa);
   wetc=0.622*Le*Qs/(Rgas*pow(ts+tdk,2.0));
   
   zo = 0.0;
   zot=0.0;
   zoq=0.0;
   L=0.0;
   /***************   wave parameters  *********/
   //lwave=grav/2/PI*pow(twave,2);
   //cwave=grav/2/PI*twave;
   
   /**************  compute aux stuff *******/
   //Rns=Rs*.945;
   //Rnl=0.97*(5.67e-8*pow(ts-0.3*jcool+tdk,4)-Rl);
   
   /***************   Begin bulk loop *******/
   
   /***************  first guess ************/
   du=u-us;
   dt=ts-t-.0098*zt;
   dq=Qs-Q;
   ta=t+tdk;
   ug=.5;
   dter=0.3; 
   dqer=wetc*dter;
   ut=sqrt(du*du+ug*ug);
   u10=ut*log(10/1e-4)/log(zu/1e-4);
   usr=.035*u10;
   zo10=0.011*usr*usr/grav+0.11*visa/usr;
   Cd10=pow(von/log(10/zo10),2);
   Ch10=0.00115;
   Ct10=Ch10/sqrt(Cd10);
   zot10=10/exp(von/Ct10);
   Cd=pow(von/log(zu/zo10),2);
   Ct=von/log(zt/zot10);
   CC=von*Ct/Cd;
   Ribcu=-zu/zi/.004/pow(Beta,3);
   Ribu=-grav*zu/ta*((dt-dter*jcool)+.61*ta*dq)/pow(ut,2);

   if(Ribu<0){
     zetu=CC*Ribu/(1+Ribu/Ribcu);
   }else{
     zetu=CC*Ribu*(1+27/9*Ribu/CC);
   }
   
   L10=zu/zetu;
   if (zetu>50){
    nits=1;
   }
   usr=ut*von/(log(zu/zo10) - psiu_30(zu/L10));
   tsr=-dt*von*fdg/(log(zt/zot10)-psit_30(zt/L10));
   qsr=-dq*von*fdg/(log(zq/zot10)-psit_30(zq/L10));
   //tsr=-(dt-dter*jcool)*von*fdg/(log(zt/zot10)-psit_30(zt/L10));
   //qsr=-(dq-wetc*dter*jcool)*von*fdg/(log(zq/zot10)-psit_30(zq/L10));
   
   tkt=.001;
   
   charn=0.011;
   if (ut>10.0){
     charn=0.011+(ut-10)/(18-10)*(0.018-0.011);
   }
   if(ut>18.0){
     charn=0.018;
   }
   
   /***************  bulk loop ************/
  for(i=0;i<nits;i++){
     zet=von*grav*zu/ta*(tsr*(1+0.61*Q)+.61*ta*qsr)/(usr*usr)/(1+0.61*Q);
     //Hard-wire it without waves
     zo=charn*usr*usr/grav+0.11*visa/usr; // Eq. 6

     //if((int)jwave==0){
     //    zo=charn*usr*usr/grav+0.11*visa/usr;
     //}
     //if((int)jwave==1){
     //    zo=50.0/2.0/PI*lwave*pow(usr/cwave,4.5)+0.11*visa/usr;
     //} // Oost et al
     //if((int)jwave==2){
     //    zo=1200.0*hwave*pow(hwave/lwave,4.5)+0.11*visa/usr;
     //} // Taylor and Yelland
      
      rr=zo*usr/visa; // Roughness reynolds number (see Eq. 7)
      L=zu/zet;
      zoq=Min(1.15e-4,5.5e-5/pow(rr,-0.6));// Eq. 28
      zot=zoq;
      usr=ut*von/(log(zu/zo)-psiu_30(zu/L));
      tsr=-dt*von*fdg/(log(zt/zot10)-psit_30(zt/L10));
      qsr=-dq*von*fdg/(log(zq/zot10)-psit_30(zq/L10));
      //tsr=-(dt-dter*jcool)*von*fdg/(log(zt/zot)-psit_30(zt/L));
      //qsr=-(dq-wetc*dter*jcool)*von*fdg/(log(zq/zoq)-psit_30(zq/L));
      Bf=-grav/ta*usr*(tsr+.61*ta*qsr);
      
      if (Bf>0){
	ug=Beta*pow(Bf*zi,0.333);
      }else{
	ug=.2;
      }
      
      ut=sqrt(du*du+ug*ug);
      //Rnl=0.97*(5.67e-8*pow(ts-dter*jcool+tdk,4)-Rl);
      hsb=-rhoa*cpa*usr*tsr;
      hlb=-rhoa*Le*usr*qsr;
      /*
      qout=Rnl+hsb+hlb;
      dels=Rns*(.065+11*tkt-6.6e-5/tkt*(1-exp(-tkt/8.0e-4))); 	// Eq.16 Shortwave
      qcol=qout-dels;
      alq=Al*qcol+be*hlb*cpw/Le;					// Eq. 7 Buoy flux water

     if (alq>0){
     	xlamx=6/pow(1+pow(bigc*alq/pow(usr,4),0.75),0.333);				// Eq 13 Saunders
        tkt=xlamx*visw/(sqrt(rhoa/rhow)*usr);			//Eq.11 Sub. thk

     }else{
       xlamx=6.0;
       tkt=Min(.01,xlamx*visw/(sqrt(rhoa/rhow)*usr));			// Eq.11 Sub. thk
     }
     
      dter=qcol*tkt/tcw; //  Eq.12 Cool skin
      dqer=wetc*dter;
      */
     
  }//bulk iter loop
  
  tau=rhoa*usr*usr*du/ut;  //stress
  hsb=-rhoa*cpa*usr*tsr;
  hlb=-rhoa*Le*usr*qsr;
  
  
  ///****************   rain heat flux ********/
  //dwat=2.11e-5*pow((t+tdk)/tdk,1.94); // water vapour diffusivity
  //dtmp=(1.+3.309e-3*t-1.44e-6*t*t)*0.02411/(rhoa*cpa); 	//heat diffusivity
  //alfac= 1/(1+(wetc*Le*dwat)/(cpa*dtmp));      	// wet bulb factor
  //RF= rain*alfac*cpw*((ts-t-dter*jcool)+(Qs-Q-dqer*jcool)*Le/cpa)/3600;
  //
  ///****************   Webb et al. correection  ************/
  //wbar=1.61*hlb/Le/(1+1.61*Q)/rhoa+hsb/rhoa/cpa/ta;//formulation in hlb already includes webb
  //hl_webb=rhoa*wbar*Q*Le;
  //
  ///**************   compute transfer coeffs relative to ut @meas. ht **********/
  Cd=tau/rhoa/ut/Max(.1,du);
  Ch=-usr*tsr/ut/(dt);
  Ce=-usr*qsr/(dq)/ut;
  //Ch=-usr*tsr/ut/(dt-dter*jcool);
  //Ce=-usr*qsr/(dq-dqer*jcool)/ut;
  //
  ///************  10-m neutral coeff realtive to ut ********/
  //Cdn_10=von*von/log(10/zo)/log(10/zo);
  //Chn_10=von*von*fdg/log(10/zo)/log(10/zot);
  //Cen_10=von*von*fdg/log(10/zo)/log(10/zoq);
  
  /* Output variable back into the input array */
  y[0] = hsb;
  y[1] = hlb;
  y[2] = tau;
  y[3] = zo;
  y[4] = zot;
  y[5] = zoq;
  y[6] = rhoa;
  y[7] = usr;
  y[8] = tsr;
  y[9] = qsr;
  y[10] = Cd;
  //y[11] = Ch;
  //y[12] = Ce;
  //y[10] = dter;
  //y[11] = dqer;
  //y[12] = tkt;
  //y[13] = RF;
  //y[14] = ut;
  //y[15] = Cd;
  //y[16] = Ch;
  //y[17] = Ce;
  //y[18] = ug;
  
  //printf("y[13]:%6.6f, y[14]:%6.6f, y[15]:%6.6f, y[16]:%6.6f, y[17]:%6.6f, y[18]:%6.6f\n",y[13],y[14],y[15],y[16],y[17],y[18]);
  //y=[hsb hlb tau zo zot zoq L usr tsr qsr dter dqer tkt RF wbar Cd Ch Ce Cdn_10 Chn_10 Cen_10 ug ];
  //   1   2   3   4  5   6  7  8   9  10   11   12  13  14  15  16 17 18    19      20    21  22
} // End cor30a

/* Velocity stability function */
static REAL psiu_30(REAL zet){
  
  REAL x, psik, psic, f, c,  psi;

  x=pow(1.0-15.0*zet,0.25);
  psik=2.0*log((1.0+x)/2.0)+log((1+x*x)/2.0)-2*atan(x)+2.0*atan(1.0);
  x=pow(1-10.15*zet,0.3333);
  psic=1.5*log((1.0+x+x*x)/3.0)-sqrt(3.0)*atan((1.0+2.0*x)/sqrt(3.0))+4.0*atan(1.0)/sqrt(3.0);
  f=zet*zet/(1.0+zet*zet);
  psi=(1.0-f)*psik+f*psic;                                               
  
  if(zet>0){
    c=Min(50.0,0.35*zet);
    psi=-(pow(1.0+1.0*zet,1.0)+.667*(zet-14.28)/exp(c)+8.525);
  }
  return psi;
  
} // End psiu_30

/* Temperature stability function */
static REAL psit_30(REAL zet){
 
  REAL x, psik, psic, f, c,  psi;
  
  x=pow(1.0-15.0*zet,0.5);
  psik=2.0*log((1+x)/2.0);
  x=pow(1-34.15*zet,0.3333);
  psic=1.5*log((1.0+x+x*x)/3.0)-sqrt(3.0)*atan((1.0+2.0*x)/sqrt(3.0))+4.0*atan(1.0)/sqrt(3.0);
  f=zet*zet/(1+zet*zet);
  psi=(1.0-f)*psik+f*psic;  
  
  if(zet>0){
    c=Min(50,0.35*zet);
    psi=-(pow(1+2/3*zet,1.5)+0.6667*(zet-14.28)/exp(c)+8.525);
  }
  return psi;
} // End psit_30


/* Alter the volume of the water column to account for rain and evaporation.
   rain is taken as a signed quantity in m/s.  Evaporation is limited to
   a minimum depth to avoid creating hypersaline cells which can create salt mass
   during wetting/drying.

   Only active for metmodel==5.
 */ 
void RainEvapSources(gridT *grid, physT *phys, metT *met, propT *prop, int myproc, MPI_Comm comm)
{
  int i,iptr;
  REAL dz_avail;
  
  if( prop->metmodel==5 ) {
    // remove evaporation/add rain.
    for(iptr=grid->celldist[0];iptr<grid->celldist[1];iptr++) {
      i = grid->cellp[iptr];
      // store the actual 'velocity' in EP,
      // and flip sign to match other usages of EP, where it is evaporation - precip.
      met->EP[i]=-met->rain[i];
      // adjusted to avoid removing too much water.
      if( met->EP[i]>0 ) {
        // at most evaporate enough to get down to DZMIN_EVAP
        dz_avail=phys->h[i]+grid->dv[i] - DZMIN_EVAP;
        met->EP[i]=Min( met->EP[i], dz_avail/prop->dt );
      }
      phys->h[i] += -met->EP[i]*prop->dt;
    }
  }
}
                     
