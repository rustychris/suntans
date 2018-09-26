/*
 * Rusty Holleman
 * calculate how many substeps would be required for upwind
 * advection. This does not yet perform any substepping,
 * but this calculation is useful both in anticipation of
 * substepping and also for calculating a finite-volume
 * Courant number.
 */

void CalculateSubsteps(gridT *grid, physT *phys, propT *prop, int myproc,
                       int numprocs, MPI_Comm comm);
