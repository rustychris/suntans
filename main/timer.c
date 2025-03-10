/*
 * File: timer.c
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Contains functions used for wallclock timing.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#include "mympi.h"
#include "timer.h"

/* global state for timer */
REAL t_start, t_source, t_predictor, t_nonhydro, t_turb, t_transport, t_io, t_comm,
  t_check, t_tictoc;


/*
 * Function: Timer
 * Usage: printf("Time = %f\n",Timer()-t0);
 * ----------------------------------------
 * Returns the time in seconds and uses the timer
 * defined in timer.h.
 *
 */
extern REAL Timer(void) {
  return (REAL)MPI_Wtime();
}

/*
 * Function: Tic
 * ----------------------------------------
 * Initializes the timer (tic-toc timer)
 * 
 */
void Tic(void) {
  t_tictoc = (REAL)MPI_Wtime();
}

/*
 * Function: Toc
 * ----------------------------------------
 * Returns the time span from the Tic()
 * initialization.
 * 
 */
REAL Toc(void) {
  return ((REAL)MPI_Wtime() - t_tictoc);
}
