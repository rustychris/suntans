/*
 * File: sendrecv.h
 * Author: Oliver B. Fringer
 * Institution: Stanford University
 * --------------------------------
 * Header file for sendrecv.c.
 *
 * Copyright (C) 2005-2006 The Board of Trustees of the Leland Stanford Junior 
 * University. All Rights Reserved.
 *
 */
#ifndef _sendrecv_h
#define _sendrecv_h

#include "grid.h"
#include "mympi.h"

void AllocateTransferArrays(gridT **grid, int myproc, int numprocs, MPI_Comm comm);
void FreeTransferArrays(gridT *grid, int myproc, int numprocs, MPI_Comm comm);

void ISendRecvCellData2DTag(REAL *celldata, gridT *grid, int myproc, MPI_Comm comm,int tag);
// void ISendRecvCellData2D(REAL *celldata, gridT *grid, int myproc, MPI_Comm comm);
#define ISendRecvCellData2D(celldata,grid,myproc,comm) \
  ISendRecvCellData2DTag(celldata,grid,myproc,comm,__LINE__)

void ISendRecvCellData3DTag(REAL **celldata, gridT *grid, int myproc, MPI_Comm comm,int tag);
// void ISendRecvCellData3D(REAL **celldata, gridT *grid, int myproc, MPI_Comm comm);
#define ISendRecvCellData3D(celldata, grid,myproc,comm) \
  ISendRecvCellData3DTag(celldata,grid,myproc,comm,__LINE__)

void ISendRecvWDataTag(REAL **celldata, gridT *grid, int myproc, MPI_Comm comm,int tag);
#define ISendRecvWData(celldata,grid,myproc,comm) \
  ISendRecvWDataTag(celldata,grid,myproc,comm,__LINE__)

void ISendRecvEdgeData3DTag(REAL **edgedata, gridT *grid, int myproc, MPI_Comm comm, int tag);
#define ISendRecvEdgeData3D(edgedata,grid,myproc,comm) \
  ISendRecvEdgeData3DTag(edgedata,grid,myproc,comm,__LINE__)

void ISendRecvEdgeData2DTag(REAL *edgedata, gridT *grid, int myproc, MPI_Comm comm, int tag);
#define ISendRecvEdgeData2D(edgedata,grid,myproc,comm) \
  ISendRecvEdgeData2DTag(edgedata,grid,myproc,comm,__LINE__)

void CheckCommunicateCells(gridT *maingrid, gridT *localgrid, int myproc, MPI_Comm comm);
void CheckCommunicateEdges(gridT *maingrid, gridT *localgrid, int myproc, MPI_Comm comm);

void SyncBarrier(int tag,int myproc, MPI_Comm comm);


#endif

