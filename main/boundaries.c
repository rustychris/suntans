/* -*- tab-width: 8 -*-
 * File: boundaries.c
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * This file contains functions to impose the boundary conditions on u.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#include "sendrecv.h"
#include "boundaries.h"
#include "mynetcdf.h"
#include "sediments.h"

/* boundary-related globals */
boundT *bound;
FILE *windFID;


// Local functions
//static void GetBoundaryVelocity(REAL *ub, int *forced, REAL x, REAL y, REAL t, REAL h, REAL d, REAL omega, REAL amp);
//static void SetUVWH(gridT *grid, physT *phys, propT *prop, int ib, int j, int boundary_index, REAL boundary_flag);

static void MatchBndPoints(propT *prop, gridT *grid, int myproc);
static void FluxtoUV(propT *prop, gridT *grid, int myproc,MPI_Comm comm);
static void SegmentArea(propT *prop, gridT *grid, int myproc, MPI_Comm comm);
int isGhostEdge(int j, gridT *grid, int myproc);
int getTimeRecBnd(REAL nctime, REAL *time, int nt);

/*
 * Function: OpenBoundaryFluxes
 * Usage: OpenBoundaryFluxes(q,ubnew,ubn,grid,phys,prop);
 * ----------------------------------------------------
 * This will update the boundary flux at the edgedist[2] to edgedist[3] edges.
 *
 * Note that phys->uold,vold contain the velocity at time step n-1 and
 * phys->uc,vc contain it at time step n.
 *
 * The radiative open boundary condition does not work yet!!!  For this reason c[k] is
 * set to 0
 *
 */
void OpenBoundaryFluxes(REAL **q, REAL **ub, REAL **ubn, gridT *grid, physT *phys, propT *prop,
                        int myproc) {
  int j, jptr, ib, k, forced;
  REAL **uc = phys->uc, **vc = phys->vc, **ucold = phys->uold, **vcold = phys->vold;
  REAL z, c0, c1, C0, C1, dt=prop->dt, u0, u0new, uc0, vc0, uc0old, vc0old, ub0;

  for(jptr=grid->edgedist[2];jptr<grid->edgedist[5];jptr++) {
    j = grid->edgep[jptr];

    ib = grid->grad[2*j];

    // RH: This used to just zero out the entries that are immediately set below.
    for(k=0;k<grid->etop[j];k++) 
      ub[j][k]=0;

    for(k=grid->etop[j];k<grid->Nke[j];k++) {
      ub[j][k]=phys->boundary_u[jptr-grid->edgedist[2]][k]*grid->n1[j]+
	phys->boundary_v[jptr-grid->edgedist[2]][k]*grid->n2[j];
#if defined(DBG_PROC) && defined(DBG_EDGE)
      if(myproc==DBG_PROC && j==DBG_EDGE) {
        printf("[p=%d j=%d k=% 3d] OpenBoundaryFluxes ub(aka utmp)=%.5f\n",
               myproc, j, k, ub[j][k]);
      }
#endif
    }
  }
}

/* Function: PointSources
 * Usage: PointSources(grid,phys,prop,myproc, comm);
 * -------------------------------------------------------------
 * Updates cell volumes with point-source flows
 */
void PointSources(gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm)
{
  int i,i_src;
  REAL Q, htmp;
  
  for(i=0;i<bound->Npoint_source;i++) {
    i_src=bound->ind_point[i];
    if(i_src<0) continue;
    
    Q=bound->point_Q[i];
    phys->h[i_src] += prop->dt*Q / grid->Ac[i_src];
  }
}

void PointSourcesContinuity(REAL **w, gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm)
{
  int i,k;
  int i_src;
  int k_src;
  REAL Q;
  REAL theta=prop->theta;

  for(i=0;i<bound->Npoint_source;i++) {
    i_src=bound->ind_point[i];
    if(i_src<0) continue;
    
    
    Q=bound->point_Q[i];

    k_src=bound->point_layer[i];
    // Clamp to available layers
    if(k_src>grid->Nk[i_src]-1)
      k_src=grid->Nk[i_src]-1;
    else if(k_src<grid->ctop[i_src])
      k_src=grid->ctop[i_src];

    // Theta here is tricky -
    // Continuity is basically this:
    // theta* (w[k]-w[k+1]) + (1-theta) (wtmp2[k]-wtmp2[k+1]) = 0
    // plus horizontal fluxes.
    // a mass balance over the time step with semi-implicit w.
    // with point source, it should be
    // theta* (w[k]-w[k+1]) + (1-theta) (wtmp2[k]-wtmp2[k+1]) - Q[k]/A= 0
    // with Q centered on the [t,t+\Delta t] interval, while wtmp2 is at
    // t and w is at t+\Delta t.
    // un-rearranging:
    // w[k]=w[k+1] - (1-theta)/theta * (wtmp2[k]-wtmp2[k+1]) + Q[k]/A/theta
    // thus here we use 1/theta.
    // it's a counter-intuitive, but it’s the same 1/theta factor that the
    // horizontal fluxes get.

    for(k=k_src;k>=grid->ctop[i_src];k--) {
      w[i_src][k] += Q/grid->Ac[i_src]/theta;
    }
  }
}

void PointSourceScalar(REAL *scalar,REAL **A, REAL **B, gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm)
{
  int i;
  int i_src;
  int k_src;
  REAL Q_src;
  
  // what is the difference between A and B?
  // UpdateScalars takes src1 (B) and src2 (A)
  //  - src1 is added to the diagonal in the tridiag system,
  //    and has units of 1/time.
  //  - src2 is added to the rhs. is has units of concentration/time
  for(i=0;i<bound->Npoint_source;i++) {
    i_src=bound->ind_point[i];
    if(i_src<0) continue;
    
    k_src=bound->point_layer[i];
    // Clamp to available layers
    if(k_src>grid->Nk[i_src]-1)
      k_src=grid->Nk[i_src]-1;
    else if(k_src<grid->ctop[i_src])
      k_src=grid->ctop[i_src];
    
    Q_src=bound->point_Q[i];
    //printf("Updating scalar: A[i=%d][k=%d] += scal=%.2f\n",
    //       i_src,k_src,scalar[i]);
    A[i_src][k_src] += scalar[i]*Q_src/(grid->Ac[i_src]*grid->dzzold[i_src][k_src]);
  }
}

void PointSourceTemp(REAL **A, REAL **B, gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm)
{
  PointSourceScalar(bound->point_T,A,B,grid,phys,prop,myproc,comm);
}
void PointSourceSalt(REAL **A, REAL **B, gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm)
{
  PointSourceScalar(bound->point_S,A,B,grid,phys,prop,myproc,comm);
}

/*
 * Function: BoundaryScalars
 * Usage: BoundaryScalars(boundary_s,boundary_T,grid,phys,prop);
 * -------------------------------------------------------------
 * This will set the values of the scalars at the open boundaries.
 * 
 */
void BoundaryScalars(gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm) {
  int jptr, j, ib, k, jind;
  int iptr, i, ii;
  int nf,ne,neigh;
  REAL z;
  
  //Type-2 zero gradient (Neumann) boundary condition
  /* 
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
      j=grid->edgep[jptr];
      ib=grid->grad[2*j];

      for(k=grid->ctop[ib];k<grid->Nk[ib];k++) {
	phys->boundary_T[jptr-grid->edgedist[2]][k]=phys->T[ib][k];
	phys->boundary_s[jptr-grid->edgedist[2]][k]=phys->s[ib][k];
      }
  }
  */

  if(prop->netcdfBdy) {
    BoundaryScalarGeneric(phys->boundary_T, // [jptr-edgedist[2], k]
                          phys->T, 
                          bound->T_scal,
                          grid,phys,prop,myproc,comm);
    BoundaryScalarGeneric(phys->boundary_s, // [jptr-edgedist[2], k]
                          phys->s,
                          bound->S_scal,
                          grid,phys,prop,myproc,comm);
  } else {
    //No NetCDF
    // Type-2
    for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
      jind = jptr-grid->edgedist[2];
      j = grid->edgep[jptr];
      ib=grid->grad[2*j];
      
      for(k=grid->ctop[ib];k<grid->Nk[ib];k++) {
        phys->boundary_T[jind][k]=0;
        phys->boundary_s[jind][k]=0;
      }
    }

    for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
      i = grid->cellp[iptr];
      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
        phys->T[i][k] = 0;
        phys->s[i][k] = 0;
      }
    }
  }
  // Need to communicate the cell data for type 3 boundaries
  ISendRecvCellData3D(phys->T,grid,myproc,comm);
  ISendRecvCellData3D(phys->s,grid,myproc,comm);

  // Set the edge array to the value in the boundary array
  /*
    RH: note that ind3edge is now longer than ind3, and can have more entries
    than ind3.  Not sure how that affects this code.

    The code that this note preceded has been dropped, but leaving this
    note in place out of paranoia.
  */ 
} // End funciton

/*
 * Function: BoundaryVelocities
 * Usage: BoundaryVelocities(grid,phys,prop);
 * ------------------------------------------
 * This will set the values of u,v,w, and h at the boundaries.
 *
 */
void BoundaryVelocities(gridT *grid, physT *phys, propT *prop, int myproc, MPI_Comm comm) {
  int i, ii, j, jj, jind, iptr, jptr, n, k;
  REAL u,v,w,h,rampfac;

   //printf("#####\nUpdating boundary velocities on processor: %d\n#####\n",myproc);
  if(prop->thetaramptime>0){
      rampfac = 1-exp(-prop->rtime/prop->thetaramptime);//Tidal rampup factor 
  }else{
      rampfac = 1.0;
  }

  // Test
  // REAL amp = 0.25;
  // REAL omega = 7.27e-5;

   // Update the netcdf boundary data
   if(prop->netcdfBdy==1){
       UpdateBdyNC(prop,grid,myproc,comm);
   }

  // Type-2
  if(prop->netcdfBdy){
  ii=-1;
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
    jind = jptr-grid->edgedist[2];
    j = grid->edgep[jptr];
    ii+=1;

    //phys->boundary_h[jind]=0.0; // Not used??
    
    for(k=0;k<grid->etop[j];k++) {
      phys->boundary_u[jind][k]=0.0;
      phys->boundary_v[jind][k]=0.0;
      phys->boundary_w[jind][k]=0.0;
    }
    
    for(k=grid->etop[j];k<grid->Nke[j];k++) {
     phys->boundary_u[jind][k]=bound->boundary_u[k][bound->ind2[ii]]*rampfac;
     phys->boundary_v[jind][k]=bound->boundary_v[k][bound->ind2[ii]]*rampfac;
     phys->boundary_w[jind][k]=bound->boundary_w[k][bound->ind2[ii]]*rampfac;
    }
  }
  }else{ // No NetCDF
      for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
	jind = jptr-grid->edgedist[2];
	j = grid->edgep[jptr];
	for(k=grid->etop[j];k<grid->Nke[j];k++) {
	    phys->boundary_u[jind][k]=0;
	    phys->boundary_v[jind][k]=0;
	    phys->boundary_w[jind][k]=0;
	}
      }
  }

  // Type-3
  if(prop->netcdfBdy){
    ii=-1;
    for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
      i = grid->cellp[iptr];
      ii+=1;
      // RH: 20180724: used to have rampfac here, but since h may have a constant
      // offset, scaling toward 0 is not a good solution.
      phys->h[i]=bound->h[bound->ind3[ii]];
    }
  }else{ // No NetCDF
    for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
      i = grid->cellp[iptr];
      phys->h[i]=0;
    }
  }

  // RH - seems that type marker 3 edges are underconstrained --
  // set them to zero.  this probably isn't necessary as these
  // should just be zero for all time.
  for(jptr=grid->edgedist[3];jptr<grid->edgedist[4];jptr++)
    {
      jind = jptr-grid->edgedist[2];
      j = grid->edgep[jptr];
      for(k=grid->etop[j];k<grid->Nke[j];k++) {
        phys->boundary_u[jind][k]=0;
        phys->boundary_v[jind][k]=0;
        phys->boundary_w[jind][k]=0;
      }
    }

  // Set velocities in type-3 cells
  // Try recalculating the cell tops here
  ISendRecvCellData2D(phys->h,grid,myproc,comm);
  UpdateDZ(grid,phys,prop,0,myproc);
  if(prop->netcdfBdy){
    ii=-1;
    for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
      i = grid->cellp[iptr];
      ii+=1;
      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
        phys->uc[i][k]=bound->uc[k][bound->ind3[ii]]*rampfac;
        phys->vc[i][k]=bound->vc[k][bound->ind3[ii]]*rampfac;
        //phys->wc[i][k]=bound->wc[k][bound->ind3[ii]]*rampfac;
      }
    }
  }else{//No NetCDF
    for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
      i = grid->cellp[iptr];
      for(k=grid->ctop[i];k<grid->Nk[i];k++) {
        phys->uc[i][k]=0;
        phys->vc[i][k]=0;
      }
    }
  }
  // Need to communicate the cell data for type 3 boundaries
  ISendRecvCellData3D(phys->uc,grid,myproc,comm);
  ISendRecvCellData3D(phys->vc,grid,myproc,comm);
  //ISendRecvCellData3D(phys->wc,grid,myproc,comm);
}

/*
 * Function: WindStress
 * Usage: WindStress(grid,phys,prop,myproc);
 * -----------------------------------------
 * Set the wind stress as well as the bottom stress.
 * tau_B is not currently in use (4/1/05).
 *
 * !!!!!!!!!!! This should be moved to met.c !!!!!!!!!!!!!!!
 */
void WindStress(gridT *grid, physT *phys, propT *prop, metT *met, int myproc) {
  int j, jptr;
  int Nc=grid->Nc; 
  int i,iptr, nf, ne, nc1, nc2, neigh;
  REAL dgf, def1, def2, rampfac;

  // I don't think that wind stress needs to be ramped up, and it makes
  // short testing runs more error-prone.
  // if(prop->thetaramptime>0){
  //   rampfac = 1-exp(-prop->rtime/prop->thetaramptime);//Rampup factor 
  // }else{
  rampfac = 1.0;
  //}
  
  if(prop->metmodel>0){
    // Interpolate the spatially variable wind stress onto the edges
    //Loop through edges and find the cell neighbours
    // computational edges
    for(jptr=grid->edgedist[0];jptr<grid->edgedist[1];jptr++) {
      j = grid->edgep[jptr]; 
      
      nc1 = grid->grad[2*j];
      nc2 = grid->grad[2*j+1];
      if(nc1==-1) nc1=nc2;
      if(nc2==-1) nc2=nc1;
      
      // Note that dgf==dg only when the cells are orthogonal!
      def1 = grid->def[nc1*grid->maxfaces+grid->gradf[2*j]];
      def2 = grid->def[nc2*grid->maxfaces+grid->gradf[2*j+1]];
      dgf = def1+def2;
      
      //This assume cells are orthogonal
      //   phys->tau_T[ne] = (met->tau_x[nc1]*def1/grid->dg[ne] +
      //	   met->tau_x[nc2]*def2/grid->dg[ne])*grid->n1[ne] + 
      //       (met->tau_y[nc1]*def1/grid->dg[ne] + 
      //	met->tau_y[nc2]*def2/grid->dg[ne])*grid->n2[ne];  
      
      def1 /= dgf;
      def2 /= dgf;
      phys->tau_T[j] = (met->tau_x[nc2]*def1 + met->tau_x[nc1]*def2)*grid->n1[j] + 
        (met->tau_y[nc2]*def1 + met->tau_y[nc1]*def2)*grid->n2[j];  
      
      phys->tau_T[j] /= RHO0; 
      phys->tau_T[j] *= rampfac;
    }
  }else{
    // Set stress to constant
    for(jptr=grid->edgedist[0];jptr<grid->edgedist[5];jptr++) {
      j = grid->edgep[jptr];

      phys->tau_T[j]=grid->n2[j]*prop->tau_T;
      phys->tau_B[j]=0;
    }
  }
}

/*
* Function: InitBoundaryData()
* -----------------------------
* Initialise the boundary condition data for the model
*
* This is called from phys.c
*/
void InitBoundaryData(propT *prop, gridT *grid, int myproc, MPI_Comm comm){

    // Step 1) Allocate the structure array data
    // Moved to phys.c

    // Step 2) Read in the coordinate info
    if(VERBOSE>1 && myproc==0) printf("Reading netcdf boundary coordinate data...\n");
    ReadBndNCcoord(prop->netcdfBdyFileID, prop, grid, myproc,comm);

    // Step 3) Match each boundary point with its local grid point
    if(VERBOSE>1 && myproc==0) printf("Matching boundary points...\n");
    MatchBndPoints(prop, grid, myproc);

    // Step 4) Read in the forward and backward time steps into the boundary arrays
    if(VERBOSE>1 && myproc==0) printf("Reading netcdf boundary initial data...\n");
    ReadBdyNC(prop, grid, myproc,comm);
    
}//end function

 /*
  * Function: UpdateBdyNC()
  * -----------------------------
  * Update the boundary netcdf data and temporally interpolate onto the model time step
  *
  */     
void UpdateBdyNC(propT *prop, gridT *grid, int myproc, MPI_Comm comm){
  int n, j, k, t0, t1, t2;
  int scal_idx;
  REAL dt, r1, r2, mu;
  if(VERBOSE>2) {
    printf("[p=%d] top of UpdateBdyNC\n",myproc);
  }
  
  t1 = getTimeRecBnd(prop->nctime,bound->time,bound->Nt);
  t0=t1-1;
  t2=t1+1;
  
  /* Only update the boundary data if need to*/
  if (bound->t1!=t1){
    if(VERBOSE>3 && myproc==0) printf("Updating netcdf boundary data at nc timestep: %d\n",t1);
    /* Read in the data two time steps*/
    bound->t1=t1;
    bound->t2=t2;
    bound->t0=t0;
    ReadBdyNC(prop, grid, myproc,comm);
    //Try this to avoid parallel read errors
    /*
      for(n=0;n<64;n++){
      MPI_Barrier(comm);
      if(n==myproc)
      ReadBdyNC(prop, grid, myproc);
      }
    */
  }

    
  if(VERBOSE>2) {
    printf("[p=%d] UpdateBdyNC post ReadBdyNC\n",myproc);
  }
  
   /*Linear temporal interpolation coefficients*/
    //dt = bound->time[bound->t1]-bound->time[bound->t0];
    //r2 = (prop->nctime - bound->time[bound->t0])/dt;
    //r1 = 1.0-r2;

    /*Cosine temporal interpolation coefficients*/
    //dt = bound->time[bound->t1]-bound->time[bound->t0];
    //mu = (prop->nctime - bound->time[bound->t0])/dt;
    //r2 = (1.0 - cos(PI*mu))/2.0;
    //r1 = 1.0-r2;

  if(bound->hasType2>0){
    //printf("Updating type-2 boundaries on proc %d\n",myproc);
    for (j=0;j<bound->Ntype2;j++){
      for (k=0;k<bound->Nk;k++){
        //Quadratic temporal interpolation
        bound->boundary_u[k][j] = QuadInterp(prop->nctime,bound->time[t0],bound->time[t1],bound->time[t2],
                                             bound->boundary_u_t[0][k][j],bound->boundary_u_t[1][k][j],bound->boundary_u_t[2][k][j] );
        bound->boundary_v[k][j] = QuadInterp(prop->nctime,bound->time[t0],bound->time[t1],bound->time[t2],
                                             bound->boundary_v_t[0][k][j],bound->boundary_v_t[1][k][j],bound->boundary_v_t[2][k][j] );
        bound->boundary_w[k][j] = QuadInterp(prop->nctime,bound->time[t0],bound->time[t1],bound->time[t2],
                                             bound->boundary_w_t[0][k][j],bound->boundary_w_t[1][k][j],bound->boundary_w_t[2][k][j] );
        for(scal_idx=0;scal_idx<bound->num_scalars;scal_idx++) {
          bound->scalars[scal_idx].boundary_scal[k][j] =
            QuadInterp(prop->nctime,
                       bound->time[t0],bound->time[t1],bound->time[t2],
                       bound->scalars[scal_idx].boundary_scal_t[0][k][j],
                       bound->scalars[scal_idx].boundary_scal_t[1][k][j],
                       bound->scalars[scal_idx].boundary_scal_t[2][k][j] );
        }
      }
    }
  }
  if(bound->hasType3>0){
    //printf("Updating type-3 boundaries on proc %d\n",myproc);
    for (j=0;j<bound->Ntype3;j++){
      for (k=0;k<bound->Nk;k++){
        // Quadratic temporal interpolation
        //printf("bound->S[0][0] = %f\n",bound->S[0][0]);
        //printf("bound->S_t[1][0][0] = %f, bound->S[0][0] = %f\n",bound->S_t[1][0][0],bound->S[0][0]);
        bound->uc[k][j] = bound->uc_t[1][k][j];
        bound->uc[k][j] = QuadInterp(prop->nctime,bound->time[t0],bound->time[t1],bound->time[t2],
                                     bound->uc_t[0][k][j],bound->uc_t[1][k][j],bound->uc_t[2][k][j]);
        bound->vc[k][j] = QuadInterp(prop->nctime,bound->time[t0],bound->time[t1],bound->time[t2],
                                     bound->vc_t[0][k][j],bound->vc_t[1][k][j],bound->vc_t[2][k][j]);
        bound->wc[k][j] = QuadInterp(prop->nctime,bound->time[t0],bound->time[t1],bound->time[t2],
                                     bound->wc_t[0][k][j],bound->wc_t[1][k][j],bound->wc_t[2][k][j]);
        for(scal_idx=0;scal_idx<bound->num_scalars;scal_idx++) {
          bound->scalars[scal_idx].scal[k][j] =
            QuadInterp(prop->nctime,bound->time[t0],bound->time[t1],bound->time[t2],
                       bound->scalars[scal_idx].scal_t[0][k][j],
                       bound->scalars[scal_idx].scal_t[1][k][j],
                       bound->scalars[scal_idx].scal_t[2][k][j]);
        }
      }
      bound->h[j] = QuadInterp(prop->nctime,bound->time[t0],bound->time[t1],bound->time[t2],
                               bound->h_t[0][j],bound->h_t[1][j],bound->h_t[2][j]);
    }
  }
  
  if(VERBOSE>2) {
    printf("[p=%d] UpdateBdyNC post hasType3\n",myproc);
  }
    
  // printf("Updating flux (Q) boundaries on proc %d\n",myproc);
  // Interpolate Q and find the velocities based on the (dynamic) segment area
  if(bound->hasType2 && bound->hasSeg>0){
    for (j=0;j<bound->Nseg;j++){
      bound->boundary_Q[j] = QuadInterp(prop->nctime,bound->time[t0],bound->time[t1],bound->time[t2],bound->boundary_Q_t[0][j],bound->boundary_Q_t[1][j],bound->boundary_Q_t[2][j]);
      // printf("Interpolated Q[seg=%d]=%.3f\n",j,bound->boundary_Q[j]);
    }
    // This function does the actual conversion
    FluxtoUV(prop,grid,myproc,comm);
  }
    
  if(VERBOSE>2) {
    printf("[p=%d] Finished UpdateBdyNC\n",myproc);
  }

}//End function

 /*
  * Function: FluxtoUV()
  * -----------------------------
  * Converts boundary flux information into u and v information 
  *
  */     
static void FluxtoUV(propT *prop, gridT *grid, int myproc, MPI_Comm comm){
    int ii, j, k, jptr, jind;
    int ss,n;
    REAL dz;

    // Step 1) Find the total area of each segment
    SegmentArea(prop,grid,myproc,comm);
   //for (n=0;n<bound->Nseg;n++){
   //	 printf("Processor: %d,  segment #: %d, segment ID: %d, segment area: %10.6f [m2]\n",myproc,n,bound->segp[n],bound->segarea[n]); 
   // }
    //Step 2) Loop through again and calculate the velocity for each type-2 cell 
    for(j=0;j<bound->Ntype2;j++){
    	jind = bound->localedgep[j];

	//Only calculate if the segment flag >0 
	if(bound->segedgep[j]>0 && jind!=-1 && !isGhostEdge(jind,grid,myproc) ){
	    //Find the segment index
	    ss = -1;
	    for (n=0;n<bound->Nseg;n++){
		if( bound->segp[n] == bound->segedgep[j]) 
		    ss=n;
	    }

	    //printf("Processor: %d, j: %d, jind: %d, mark: %d, segment #: %d, segment ID: %d, segment area: %10.6f [m2]\n",myproc,j,jind,grid->mark[jind],ss,bound->segp[ss],bound->segarea[ss]); 
	    //isGhostEdge(jind,grid,myproc);

	    for(k=grid->etop[jind];k<grid->Nke[jind];k++) {
		bound->boundary_u[k][j] = grid->n1[jind] * bound->boundary_Q[ss] / bound->segarea[ss]; 
		bound->boundary_v[k][j] = grid->n2[jind] * bound->boundary_Q[ss] / bound->segarea[ss]; 
	    }
	}//End if
    }
}//End function FluxtoUV

/*
 * Function: SegmentArea()
 * ----------------------------
 * Calculates the area of each boundary segment
 * Note that this needs to be done across processors...
 */
static void SegmentArea(propT *prop, gridT *grid, int myproc, MPI_Comm comm){
    int ii, j, k, jptr, jind;
    int ss,n;
    REAL dz;
    
    //Zero the area first
    for (n=0;n<bound->Nseg;n++){
	bound->localsegarea[n]=0.0;
	//bound->segarea[n]=0.0;
    }
    for(j=0;j<bound->Ntype2;j++){
        jind = bound->localedgep[j]; // the local edge pointer will = -1 if the boundary cell is not on the processor's grid
	//Only calculate if the segment flag >0 
	if(bound->segedgep[j]>0 && jind != -1 && !isGhostEdge(jind,grid,myproc) ){
	    //Find the segment index
	    ss = -1;
	    for (n=0;n<bound->Nseg;n++){
		if(bound->segp[n] ==bound->segedgep[j]) 
		    ss=n;
	    }

	    //Loop through all vertical cells and sum the area
	    for(k=grid->etop[jind];k<grid->Nke[jind];k++) {
	    	//dz = grid->dzf[jind][k];
		dz = grid->dzz[grid->grad[2*jind]][k];//Updwind cell height
		bound->localsegarea[ss]+= grid->df[jind] * dz ;
	//	 printf("Processor: %d, k: %d, j: %d,  segment #: %d, segment ID: %d, df = %f, dzf = %f, segment area: %10.6f [m2]\n",myproc,k, j, ss,bound->segp[ss],grid->df[j],dz,bound->segarea[ss]); 
	    }
	}//end if
    }

    //Sum the area across processors
    MPI_Reduce(bound->localsegarea, bound->segarea,(int)bound->Nseg, MPI_DOUBLE, MPI_SUM, 0,comm);
    MPI_Bcast(bound->segarea,(int)bound->Nseg, MPI_DOUBLE,0,comm);

}//End function SegmentArea()

/*
* Function: isGhostEdge()
* ----------------------------
* Find whether or not an edge is part of a ghost cell on a particular processor
* I think the quickest test is to find whether one of the other edges has mark==6.
*/

int isGhostEdge(int j, gridT *grid, int myproc){
    int ib, isGhost, nf, ne;
    
    isGhost=0;
    // Cell index
    ib = grid->grad[2*j];
    // Loop through the cell edges
    //printf("Processor: %d, j: %d, ",myproc,j);
    //for(nf=0;nf<NFACES;nf++){ 
	// Edge index
    //ne = grid->face[ib*NFACES+nf];
    for(nf=0;nf<grid->nfaces[ib];nf++) {
	ne = grid->face[ib*grid->maxfaces+nf];
	//printf("ne: %d, mark: %d, ",ne,grid->mark[ne]);
	if (grid->mark[ne]==6)
	    isGhost=1;
    }
    //printf("\n");

    return isGhost;

}//End function isGhostEdge


/*
 * Function: MatchBndPoints()
 * -------------------------
 * Checks that the boundary arrays match the grid sizes
 */
static void MatchBndPoints(propT *prop, gridT *grid, int myproc){
  int iptr, jptr, jj, ii, ne, nc1, nc2, nc, j, ib;

  if(myproc==0) printf("Boundary NetCDF file grid # type 2 points = %d\n",(int)bound->Ntype2);
  if(myproc==0) printf("Boundary NetCDF file grid # type 3 points = %d\n",(int)bound->Ntype3);

  if(myproc==0) printf("Matching type-2 points...\n");
  //Type-2
  ii=-1;
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
    ii+=1;

    // Match suntans edge cell with the type-2 boundary file point
    for(jj=0;jj<bound->Ntype2;jj++){
      if(grid->eptr[grid->edgep[jptr]]==bound->edgep[jj]){
        bound->ind2[ii]=jj;
        bound->localedgep[jj]=grid->edgep[jptr]; 
        //printf("grid->eptr:%d, bound->edgep[jj]: %d, jj: %d, jptr: %d, grid->edgep[jptr]: %d, localedgep[jj]: %d\n",grid->eptr[grid->edgep[jptr]],bound->edgep[jj],jj,jptr,grid->edgep[jptr],bound->localedgep[jj]); 
      }
    }
  }
  if(myproc==0) printf("Matching type-3 points...\n");
  // Type-3
  ii=-1;
  for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
    ii+=1;
    // Match suntans grid cell with the type-3 boundary file point
    for(jj=0;jj<bound->Ntype3;jj++){
      if(grid->mnptr[grid->cellp[iptr]]==bound->cellp[jj]){
        bound->ind3[ii]=jj;
      }
    }
    //printf("Type 3 : Processor = %d, jptr = %d, cellp[jptr]=%d, bound->ind3[ii]=%d\n",myproc,jptr,grid->cellp[jptr],bound->ind3[ii]);
  }

  if(myproc==0) printf("Matching point source cells...\n");
  for(jj=0;jj<bound->Npoint_source;jj++){
    bound->ind_point[jj]=-1; // if this point is on a different proc, ind_point remains -1
    for(ii=0;ii<grid->Nc;ii++) {
      if(grid->mnptr[ii]==bound->point_cell[jj]){
        bound->ind_point[jj]=ii;
        printf("Matched point source %d (cell=%d) to proc=%d cell=%d\n",
               jj,bound->point_cell[jj],myproc,ii);
      }
    }
  }
  
  if(myproc==0) printf("Matching type-2 edges...\n");

  //Type-3 edges
  ii=-1;
  for(jptr=grid->edgedist[3];jptr<grid->edgedist[4];jptr++) {
    ii+=1;
    j = grid->edgep[jptr];
    ib=grid->mnptr[grid->grad[2*j]];
    // Match suntans edge cell with the type-3 boundary file point
    for(jj=0;jj<bound->Ntype3;jj++){
      if(ib==bound->cellp[jj]){
        bound->ind3edge[ii]=jj;
      }
    }
  }
  // Check that ind2 and ind3 do not contain any -1 (non-matching points)

  // Check that the number of vertical grid points match
  if(bound->Nk != grid->Nkmax){
    printf("Error! Number of layers in open boundary file (%d) not equal to Nkmax (%d).\n",(int)bound->Nk,grid->Nkmax);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }
} // End function


// allocate the per-scalar arrays for a single scalar species, using the next unused
// entry in bound->scalars
scalar_boundT *AllocateBoundaryScalar(const char *name, boundT *bound, int myproc){
  // Nt, Npoint_source;
  int j, k, i, n;

  int Nk=bound->Nk;
  int Ntype2=bound->Ntype2;
  int Ntype3=bound->Ntype3;
  
  if(bound->num_scalars==MAXSCALARS) {
    printf("No more boundary scalar elements available. Increase MAXSCALARS (%d)\n",MAXSCALARS);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  printf("scalar %d: %s\n",bound->num_scalars, name);
  
  scalar_boundT *scalar=&(bound->scalars[bound->num_scalars]);
  bound->num_scalars++;

  strcpy(scalar->varname,name);

  if(bound->Npoint_source) {
    scalar->point_scal=(REAL*)SunMalloc( bound->Npoint_source*sizeof(REAL),"AllocateBoundaryScalar");
    for(j=0;j<bound->Npoint_source;j++) {
      scalar->point_scal[j]=0.;
    }
  } else {
    scalar->point_scal=NULL;
  }

  if(bound->hasType2==1){
    scalar->boundary_scal = (REAL **)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryScalar");
    scalar->boundary_scal_t = (REAL ***)SunMalloc(NT*sizeof(REAL),"AllocateBoundaryScalar");
    for(k=0;k<Nk;k++){
      scalar->boundary_scal[k] = (REAL *)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryScalar");
      for(j=0;j<Ntype2;j++) {
        scalar->boundary_scal[k][j]=0.0;
      }
    }      

    for(n=0;n<NT;n++){
      scalar->boundary_scal_t[n] = (REAL **)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryScalar");
      for(k=0;k<Nk;k++){
        scalar->boundary_scal_t[n][k] = (REAL *)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryScalar");
        for(j=0;j<Ntype2;j++) {
          scalar->boundary_scal_t[n][k][j]=0.0;
        }
      }
    }
  } else {
    scalar->boundary_scal=NULL;
    scalar->boundary_scal_t=NULL;
  }

  if(bound->hasType3==1){
    scalar->scal = (REAL **)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
    scalar->scal_t = (REAL ***)SunMalloc(NT*sizeof(REAL),"AllocateBoundaryData");
    for(k=0;k<Nk;k++){
      scalar->scal[k] = (REAL *)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
      for(i=0;i<Ntype3;i++) {
        scalar->scal[k][i]=0.0;
      }
    }
    for(n=0;n<NT;n++){
      scalar->scal_t[n] = (REAL **)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
      for(k=0;k<Nk;k++){
        scalar->scal_t[n][k] = (REAL *)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
        for(i=0;i<Ntype3;i++) {
          scalar->scal_t[n][k][i]=0.0;
        }
      }
    }
  } else {
    scalar->scal=NULL;
    scalar->scal_t=NULL;
  }
  
  return scalar;
}
/*
 * Function: AllocateBoundaryData()
 * -------------------------------
 * Allocate boundary structure arrays.
 */
void AllocateBoundaryData(propT *prop, gridT *grid, boundT **bound, int myproc, MPI_Comm comm){
  int Ntype2, Ntype3, Nseg, Nt, Nk, Npoint_source;
  int j, k, i, n;
  int n3, n3e, n2; // Number of type 2 and 3 on a processor's grid
  
  //Allocate memory
  if(VERBOSE>1 && myproc==0) printf("Allocating boundary data structure...\n");
  
  *bound = (boundT *)SunMalloc(sizeof(boundT),"AllocateBoundaryData");
  
  // Read in the dimensions
  printf("Reading boundary netcdf dimensions...\n");
  Ntype3 = returndimlenBC(prop->netcdfBdyFileID,"Ntype3");
  Ntype2 = returndimlenBC(prop->netcdfBdyFileID,"Ntype2");
  Nseg = returndimlenBC(prop->netcdfBdyFileID,"Nseg");
  Nt = returndimlenBC(prop->netcdfBdyFileID,"Nt");
  Nk = returndimlenBC(prop->netcdfBdyFileID,"Nk");
  Npoint_source = returndimlenBC(prop->netcdfBdyFileID,"Npoint");
  
  (*bound)->Ntype3=Ntype3;
  (*bound)->Nseg=Nseg;
  (*bound)->Ntype2=Ntype2;
  (*bound)->Nt=Nt;
  (*bound)->Nk=Nk;
  (*bound)->Npoint_source=Npoint_source;
  (*bound)->num_scalars=0;

  // Check if boundary types are in the file
  if ( (*bound)->Ntype3==0 ) {
    (*bound)->hasType3=0;
  } else {
    (*bound)->hasType3=1;
  }
  
  if ((*bound)->Ntype2==0){
    (*bound)->hasType2=0;
  } else {
    (*bound)->hasType2=1;
  }
    
  if ((*bound)->Nseg==0){
    (*bound)->hasSeg=0;
  } else {
    (*bound)->hasSeg=1;
  }

  // Print the array sizes
  if(VERBOSE>1 && myproc==0){
    printf("Ntype 3 = %d\n",Ntype3);
    printf("Ntype 2 = %d\n",Ntype2);
    printf("Nseg= %d\n",Nseg);
    printf("Nt = %d\n",Nt);
    printf("Nk = %d\n",Nk);
    printf("Npoint_source = %d\n",Npoint_source);
    printf("hasType2 = %d\n",(*bound)->hasType2);
    printf("hasType3 = %d\n",(*bound)->hasType3);
    printf("hasSeg = %d\n",(*bound)->hasSeg);
  }

  scalar_boundT *T_scal,*S_scal;
  T_scal=AllocateBoundaryScalar("T",*bound,myproc);
  S_scal=AllocateBoundaryScalar("S",*bound,myproc);
  (*bound)->T_scal=T_scal;
  (*bound)->S_scal=S_scal;
  
  if((*bound)->Npoint_source) {
    (*bound)->point_cell=(int*)SunMalloc( Npoint_source*sizeof(int),"AllocateBoundaryData");
    (*bound)->point_layer=(int*)SunMalloc( Npoint_source*sizeof(int),"AllocateBoundaryData");
    (*bound)->ind_point=(int*)SunMalloc( Npoint_source*sizeof(int),"AllocateBoundaryData");
    
    (*bound)->point_Q=(REAL*)SunMalloc( Npoint_source*sizeof(REAL),"AllocateBoundaryData");

    for(i=0;i<(*bound)->Npoint_source;i++){
      (*bound)->ind_point[i]=-1;
      (*bound)->point_cell[i]=-1;
      (*bound)->point_layer[i]=-1;
    }
  } else {
    (*bound)->point_cell=NULL;
    (*bound)->point_layer=NULL;
    (*bound)->ind_point=NULL;
    
    (*bound)->point_Q=NULL;
  }
  
  (*bound)->point_T=T_scal->point_scal;
  (*bound)->point_S=S_scal->point_scal;
  
  // Allocate the type2 arrays
  if((*bound)->hasType2==1){
    (*bound)->edgep = (int *)SunMalloc(Ntype2*sizeof(int),"AllocateBoundaryData");
    (*bound)->localedgep = (int *)SunMalloc(Ntype2*sizeof(int),"AllocateBoundaryData");
    (*bound)->xe = (REAL *)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryData");
    (*bound)->ye = (REAL *)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryData");
    
    (*bound)->boundary_u = (REAL **)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
    (*bound)->boundary_v = (REAL **)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
    (*bound)->boundary_w = (REAL **)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
    
    (*bound)->boundary_u_t = (REAL ***)SunMalloc(NT*sizeof(REAL),"AllocateBoundaryData");
    (*bound)->boundary_v_t = (REAL ***)SunMalloc(NT*sizeof(REAL),"AllocateBoundaryData");
    (*bound)->boundary_w_t = (REAL ***)SunMalloc(NT*sizeof(REAL),"AllocateBoundaryData");
    
    for(k=0;k<Nk;k++){
      (*bound)->boundary_u[k] = (REAL *)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryData");
      (*bound)->boundary_v[k] = (REAL *)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryData");
      (*bound)->boundary_w[k] = (REAL *)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryData");
    }
    for(n=0;n<NT;n++){
      (*bound)->boundary_u_t[n] = (REAL **)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
      (*bound)->boundary_v_t[n] = (REAL **)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
      (*bound)->boundary_w_t[n] = (REAL **)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
      for(k=0;k<Nk;k++){
        (*bound)->boundary_u_t[n][k] = (REAL *)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryData");
        (*bound)->boundary_v_t[n][k] = (REAL *)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryData");
        (*bound)->boundary_w_t[n][k] = (REAL *)SunMalloc(Ntype2*sizeof(REAL),"AllocateBoundaryData");
      }
    }
    // Allocate the pointer arrays for each processor's grid
    n2 = grid->edgedist[3]-grid->edgedist[2];
    (*bound)->ind2 = (int *)SunMalloc(n2*sizeof(int),"AllocateBoundaryData");
  }//endif

  (*bound)->boundary_T_t=T_scal->boundary_scal_t;
  (*bound)->boundary_T  =T_scal->boundary_scal;
  (*bound)->boundary_S_t=S_scal->boundary_scal_t;
  (*bound)->boundary_S  =S_scal->boundary_scal;
  
  // Allocate the segment arrays (subset of type2)
  if ((*bound)->hasSeg==1){
    (*bound)->segp = (int *)SunMalloc(Nseg*sizeof(int),"AllocateBoundaryData");
    (*bound)->segedgep = (int *)SunMalloc(Ntype2*sizeof(int),"AllocateBoundaryData");
    
    (*bound)->segarea = (REAL *)SunMalloc(Nseg*sizeof(REAL),"AllocateBoundaryData");
    (*bound)->localsegarea = (REAL *)SunMalloc(Nseg*sizeof(REAL),"AllocateBoundaryData");
    (*bound)->boundary_Q  = (REAL *)SunMalloc(Nseg*sizeof(REAL),"AllocateBoundaryData");
    (*bound)->boundary_Q_t  = (REAL **)SunMalloc(NT*sizeof(REAL),"AllocateBoundaryData");
    for(n=0;n<NT;n++){
      (*bound)->boundary_Q_t[n] = (REAL *)SunMalloc(Nseg*sizeof(REAL),"AllocateBoundaryData");
    }
  }
  // Allocate the type3 arrays   
  if((*bound)->hasType3==1){
    (*bound)->cellp = (int *)SunMalloc(Ntype3*sizeof(int),"AllocateBoundaryData");
    (*bound)->xv = (REAL *)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
    (*bound)->yv = (REAL *)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
    
    (*bound)->h = (REAL *)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
    (*bound)->h_t = (REAL **)SunMalloc(NT*sizeof(REAL),"AllocateBoundaryData");
    
    (*bound)->uc = (REAL **)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
    (*bound)->vc = (REAL **)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
    (*bound)->wc = (REAL **)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
    
    (*bound)->uc_t = (REAL ***)SunMalloc(NT*sizeof(REAL),"AllocateBoundaryData");
    (*bound)->vc_t = (REAL ***)SunMalloc(NT*sizeof(REAL),"AllocateBoundaryData");
    (*bound)->wc_t = (REAL ***)SunMalloc(NT*sizeof(REAL),"AllocateBoundaryData");
    
    for(k=0;k<Nk;k++){
      (*bound)->uc[k] = (REAL *)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
      (*bound)->vc[k] = (REAL *)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
      (*bound)->wc[k] = (REAL *)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
    }
    for(n=0;n<NT;n++){
      (*bound)->uc_t[n] = (REAL **)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
      (*bound)->vc_t[n] = (REAL **)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
      (*bound)->wc_t[n] = (REAL **)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");
      for(k=0;k<Nk;k++){
        (*bound)->uc_t[n][k] = (REAL *)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
        (*bound)->vc_t[n][k] = (REAL *)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
        (*bound)->wc_t[n][k] = (REAL *)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
      }
      (*bound)->h_t[n] = (REAL *)SunMalloc(Ntype3*sizeof(REAL),"AllocateBoundaryData");
    }
    // Allocate the pointer arrays for each processor's grid
    n3 = grid->celldist[2]-grid->celldist[1];
    // RH: very conservative.  Either we have to allow for more
    n3e= n3*grid->maxfaces;
    // it is possible, in particular in mpi, for the later code to
    // find more than n3 edges to put into ind3edge.  At the moment there is no code
    // which uses ind3edge.  allocate it conservatively large and press on.
    (*bound)->ind3 = (int *)SunMalloc(n3*sizeof(int),"AllocateBoundaryData");
    (*bound)->ind3edge = (int *)SunMalloc(n3e*sizeof(int),"AllocateBoundaryData");
  }//endif

  (*bound)->T = T_scal->scal;
  (*bound)->S = S_scal->scal;
  (*bound)->T_t = T_scal->scal_t;
  (*bound)->S_t = S_scal->scal_t;
  
  (*bound)->time = (REAL *)SunMalloc(Nt*sizeof(REAL),"AllocateBoundaryData");
  (*bound)->z = (REAL *)SunMalloc(Nk*sizeof(REAL),"AllocateBoundaryData");

  // Zero the boundary arrays
  if(myproc==0) printf("Zeroing boundary arrays...\n");

  for(j=0;j<(*bound)->Nt;j++){
    (*bound)->time[j]=0.0;
  }
  (*bound)->t0=-1;
  (*bound)->t1=-1;
  (*bound)->t2=-1;
  for(k=0;k<(*bound)->Nk;k++){
    (*bound)->z[k]=0.0;
  }
  
  if((*bound)->hasType2>0){
    for(j=0;j<Ntype2;j++){
      (*bound)->edgep[j]=0;
      (*bound)->localedgep[j]=-1;
      (*bound)->xe[j]=0.0;
      (*bound)->ye[j]=0.0;
      for(k=0;k<Nk;k++){
        (*bound)->boundary_u[k][j]=0.0;
        (*bound)->boundary_v[k][j]=0.0;
        (*bound)->boundary_w[k][j]=0.0;
        for(n=0;n>NT;n++){
          (*bound)->boundary_u_t[n][k][j]=0.0;
          (*bound)->boundary_v_t[n][k][j]=0.0;
          (*bound)->boundary_w_t[n][k][j]=0.0;
        }
      }
    }
    for(i=0;i<n2;i++){
      (*bound)->ind2[i]=-1;
    }
  }
  
  if((*bound)->hasSeg>0){
    for(j=0;j<Ntype2;j++){
      (*bound)->segedgep[j]=0;
    }
    for(j=0;j<Nseg;j++){
      (*bound)->segp[j]=0;
      (*bound)->segarea[j]=0.0;
      (*bound)->localsegarea[j]=0.0;
      (*bound)->boundary_Q[j]=0.0;
      for(n=0;n<NT;n++){
        (*bound)->boundary_Q_t[n][j]=0.0;
      }
    }
  }
  if(myproc==0) printf("Finished Zeroing Type2 boundary arrays...\n");
  
  if((*bound)->hasType3>0){
    for(j=0;j<Ntype3;j++){
      (*bound)->cellp[j]=0;
      (*bound)->xv[j]=0.0;
      (*bound)->yv[j]=0.0;
      
      (*bound)->h[j]=0.0;
      for(n=0;n<NT;n++){
        (*bound)->h_t[n][j]=0.0;
      }
      for(k=0;k<Nk;k++){
        (*bound)->uc[k][j]=0.0;
        (*bound)->vc[k][j]=0.0;
        (*bound)->wc[k][j]=0.0;
        for(n=0;n<NT;n++){
          (*bound)->uc_t[n][k][j]=0.0;
          (*bound)->vc_t[n][k][j]=0.0;
          (*bound)->wc_t[n][k][j]=0.0;
        }
      }
    }
    for(i=0;i<n3;i++){
      (*bound)->ind3[i]=-1;
    }
    for(i=0;i<n3e;i++){
      (*bound)->ind3edge[i]=-1;
    }
  }
  if(myproc==0) printf("Finished Zeroing Type 3 boundary arrays...\n");
} //End function


void BoundaryScalarGeneric(REAL **boundary_scal, // [jptr-edgedist[2], k]
                           REAL **scal, // [cellp,k]
                           scalar_boundT *scalar,
                           gridT *grid, physT *phys, propT *prop,int myproc, MPI_Comm comm) {
  int k,i,ii,jptr,jind,j,ib,iptr;

  ii=-1;
  for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
    jind = jptr-grid->edgedist[2];
    j = grid->edgep[jptr];
    ib=grid->grad[2*j];
    ii+=1;
    
    for(k=grid->ctop[ib];k<grid->Nk[ib];k++) {
      boundary_scal[jind][k]=scalar->boundary_scal[k][bound->ind2[ii]];
    }
  }

  ii=-1;
  for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
    i = grid->cellp[iptr];
    ii+=1;
    for(k=grid->ctop[i];k<grid->Nk[i];k++) {
      scal[i][k] = scalar->scal[k][bound->ind3[ii]];
    }
  }
}

/*
 * Function: BoundarySediment
 * Usage: BoundarySediment(boundary_s,boundary_T,grid,phys,prop);
 * -------------------------------------------------------------
 * This will set the values of the suspended sediment concentration
 * at the open boundaries.
 * 
 */
void BoundarySediment(gridT *grid, physT *phys, propT *prop,int myproc, MPI_Comm comm) {
  int jptr, j, ib, k,nosize,i,iptr;
  REAL z;

  if ( prop->netcdfBdy ) {
    for(nosize=0;nosize<sediments->Nsize;nosize++){
      BoundaryScalarGeneric(sediments->boundary_sediC[nosize],
                            sediments->SediC[nosize],
                            sediments->sed_bounds[nosize],
                            grid,phys,prop,myproc,comm);
    }
  } else {
    // At the upstream boundary
    for(jptr=grid->edgedist[2];jptr<grid->edgedist[3];jptr++) {
      j=grid->edgep[jptr];
      ib=grid->grad[2*j];
      for(nosize=0;nosize<sediments->Nsize;nosize++){
        for(k=grid->ctop[ib];k<grid->Nk[ib];k++) {
          sediments->boundary_sediC[nosize][jptr-grid->edgedist[2]][k]=0;
        }
      }
    }

    // At the ocean boundary
    for(iptr=grid->celldist[1];iptr<grid->celldist[2];iptr++) {
      i = grid->cellp[iptr];
      for(nosize=0;nosize<sediments->Nsize;nosize++){
        for(k=0;k<grid->ctop[i];k++) {
          sediments->SediC[nosize][i][k]=0;
        } 
        for(k=grid->ctop[i];k<grid->Nk[i];k++) {
          sediments->SediC[nosize][i][k]=0;
        }
      }
    }
  }
}
