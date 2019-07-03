/*
 * Function: CalculateSubsteps
 * ----------------------------------
 * Find the global timestep limitation for scalar advection, and set
 * prop->nsubsteps for later use in UpdateScalars.
 * Only horizontal advection is considered (ignores horizontal diffusive
 * limitations), and currently only valid for upwind.
 *
 */

#include "suntans.h"
#include "phys.h"
#include "grid.h"
#include "sendrecv.h"
#include "util.h"
#include "timer.h"

#include "substepping.h"

typedef struct {
  int proc;
  int i,k;
  REAL x,y;
  REAL dt_min;
  int nsubsteps;
} substep_info;

void CalculateSubsteps(gridT *grid, physT *phys, propT *prop, int myproc,
                       int numprocs, MPI_Comm comm)
{
  int i,j,k,kf,ktop,nf,proc,normal;
  REAL Qout,u_semi;
  REAL dt_min_i;
  REAL theta = prop->theta;
  substep_info mine,theirs;
  MPI_Status status;

  mine.dt_min = 1e10;
  mine.i = mine.k = -1;
  mine.proc = myproc;

  for(i=0;i<grid->Nc;i++) {
    // it's possible that this should only consider celldist[0] to celldist[1] cells.

    // In UpdateScalars, the fluxes are through face heights dzfold, and into layers
    // k = ctopold,..., Nk-1, with starting volumes dzzold*Ac

    // etopold can be higher than ctopold, though, and even worse, these edges can
    // be active, but have semi-implicit flow *out* of the cell.  basically we have to
    // recreate the ktop indexing scheme of UpdateScalars.

    // so ktop is the top layer of the cell for which fluxes will be considered.
    if ( grid->ctopold[i] == grid->Nk[i] ) {
      // even in dry cells, we still consider the bed cell.
      ktop = grid->Nk[i] - 1;
    } else {
      ktop = grid->ctopold[i];
    }

    for(k=ktop;k<grid->Nk[i];k++) {
      Qout = 0;
      for(nf=0;nf<grid->nfaces[i];nf++) {
        j = grid->face[i*grid->maxfaces+nf];

        // if normal is +1, then a positive velocity
        // is out of this cell.
        normal = grid->normal[i*grid->maxfaces+nf];

        // for the top layer, we include any higher-up flux faces, too
        // careful that k maybe less than etopold (since etop's come from upwind
        // cells, not necessarily taller cells).  Shouldn't be a problem, though,
        // since if we wind up with a kf that is above etopold, it will just have 
        // zero velocity anyway.

        for(kf = (k==ktop)? grid->etopold[j]:k ;
            kf <= k && kf < grid->Nke[j] ;
            kf ++ ) {
          
          // Use the semi-implicit velocity, the one that is used for
          // the actual flux in UpdateScalars
          u_semi = (theta*phys->u[j][kf]+(1-theta)*phys->utmp2[j][kf]) * normal;
          if ( u_semi > 0.0 ) {
            // outward flow - include it
            Qout += u_semi * grid->dzfold[j][kf] * grid->df[j];
          }
          
#ifdef DBG_CELL
          if( DBG_CELL==i && DBG_PROC==myproc ) {
            printf("p=%d i=%d, for cell layer k=%d  outward u_semiimplicit[j=%d][kf=%d] = %.10f through dzfold=%.10f\n",
                   myproc,i,k,j,kf,u_semi,grid->dzfold[j][kf]);
          }
#endif
        }
      }

      // RH 2019-01-17: include vertical fluxes, too.
      if (k>ktop ) {
        // vertical velocity on top face, positive out
        u_semi=theta*phys->wnew[i][k] + (1-theta)*phys->wtmp2[i][k];
        if(u_semi>0.0) Qout += u_semi*grid->Ac[i];
#ifdef DBG_CELL
        if( DBG_CELL==i && DBG_PROC==myproc ) {
          printf("p=%d i=%d, for cell layer k=%d   upward w_semiimplicit[i=%d][k=%d] = %.10f\n",
                 myproc,i,k,i,k,u_semi);
        }
#endif // DBG_CELL
      }
      // vertical velocity on lower face.  here we can assume that
      // the boundary is fixed, and use w without testing k<Nk[i]-1
      u_semi=-( theta*phys->wnew[i][k+1] + (1-theta)*phys->wtmp2[i][k+1] );
      if(u_semi>0.0) Qout += u_semi*grid->Ac[i];
#ifdef DBG_CELL
      if( DBG_CELL==i && DBG_PROC==myproc ) {
        printf("p=%d i=%d, for cell layer k=%d downward w_semiimplicit[i=%d][k=%d] = %.10f\n",
               myproc,i,k,i,k+1,u_semi);
      }
#endif // DBG_CELL
      
#ifdef DBG_CELL
      if( DBG_CELL==i && DBG_PROC==myproc ) {
        printf("p=%d i=%d  k=%d Qout=%.10f dzzold*Ac=%.10f\n",
               myproc,i,k,Qout,grid->dzzold[i][k]*grid->Ac[i]);
      }
#endif
      
      if ( Qout > 0.0 ) {
        // careful for dzzold of a bed cell when that watercolumn has been dried out.
        // UpdateDZ() should never zero out the bed dzz, so it should be okay (as long as
        // dzmin is nonzero)
        if ( grid->dzzold[i][k] == 0.0 ) {
          printf("DANGER: while checking number of substeps, p=%d,i=%d,k=%d has outflow, but dzzold=%f\n",
                 myproc,i,k,grid->dzzold[i][k] );
        } else {
          // note that while it seems like multiplying Qout by 2 is unnecessary, in the case
          // of TVD this is the constraint published by Casulli & Zanolli (OceMod 2005)
          dt_min_i = grid->dzzold[i][k] * grid->Ac[i] / (2*Qout);
          if ( dt_min_i < mine.dt_min ) {
            mine.i = i;
            mine.k = k;
            mine.dt_min = dt_min_i;
            mine.x = grid->xv[i];
            mine.y = grid->yv[i];
          }
          if (dt_min_i < phys->min_time_step[i] ) {
            phys->min_time_step[i]=dt_min_i;
          }
        }
      } else if ( Qout < 0.0 ) {
        printf("p=%d: Qout[i=%d][k=%d] = %f  NEGATIVE!\n",myproc,i,k,Qout);
      }
    }
  }

  mine.nsubsteps = (int) ceil( prop->dt / mine.dt_min );

  // now we have to communicate this with everyone else -

  if (myproc == 0) {
    for ( proc=1 ; proc < numprocs ; proc++ ) {
      MPI_Recv((void *)(&theirs),sizeof(theirs),MPI_BYTE,proc,37,comm,&status);

      if ( theirs.dt_min < mine.dt_min ) {
        mine.dt_min = theirs.dt_min;
        mine.proc = theirs.proc;
        mine.i = theirs.i;
        mine.k = theirs.k;
        mine.nsubsteps = theirs.nsubsteps;
        mine.x = theirs.x;
        mine.y = theirs.y;
      }

      if ( mine.nsubsteps > 500 ) {
        printf(" WARNING: proc %d i=%d k=%d yielded %d substeps. Changing that to 250\n",
               mine.proc,mine.i,mine.k, mine.nsubsteps);
        mine.nsubsteps = 250;
        mine.dt_min = prop->dt / mine.nsubsteps;
      }
    }
  } else {
    MPI_Send((void *)(&mine),sizeof(mine),MPI_BYTE,0,37,comm);
  }

  // And redistribute the result:
  MPI_Bcast((void*)(&mine),sizeof(mine),MPI_BYTE,0,comm);

  if (myproc == 0 && VERBOSE) {
    if( VERBOSE > 1 || mine.nsubsteps > 1 ) {
      printf("Substep limitation: dt_min=%f, %d substeps, from p=%d,i=%d,k=%d  %.2f %.2f\n",
             mine.dt_min,mine.nsubsteps,mine.proc,mine.i,mine.k,mine.x,mine.y);
      fflush(stdout);
    }
  }
  if ( myproc == 0 ) {
    if( prop->n == prop->nstart+1 ) {
      fprintf(prop->SubstepFID,"base_dt,limit_dt,nsubsteps,limit_proc,limit_cell,limit_layer\n");
    }
    fprintf(prop->SubstepFID,"%f,%f,%d,%d,%d,%d\n",prop->dt,mine.dt_min,mine.nsubsteps,mine.proc,mine.i,mine.k);
    // DBG
    fflush(prop->SubstepFID);
  }

  if ( mine.proc==myproc ) {
    phys->limiting_cell[mine.i]++;
  }
  
  if ( mine.nsubsteps > 0 )
    prop->nsubsteps = mine.nsubsteps;
  else
    prop->nsubsteps = 1;

  MPI_Barrier(comm);
}
