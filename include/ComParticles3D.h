/***************************************************************************
  ComParticles.h  -  Library to manage communication of particles among processors
  -------------------
begin                : May 2008
copyright            : (C) 2008 KUL
developers           : Stefano Markidis, Giovanni Lapenta
 ***************************************************************************/

#ifndef ComParticles3D_H
#define ComParticles3D_H

#include "MPIdata.h"
#include "ComBasic3D.h"

/** comunicate particles and receive particles to and from 6 processors */
void communicateParticles(long long buffer_size_x, long long buffer_size_y, long long buffer_size_z, double *b_Xleft, double *b_Xright, double *b_Yleft, double *b_Yright, double *b_Zleft, double *b_Zright, VirtualTopology3D * vct);

/** communicate the global sum */
inline int globalSum(int value) {
  int sum = 0;
//  MPI_Barrier(MPI_COMM_WORLD); // This is time-consuming and should be debug code, only!
  MPI_Allreduce(&value, &sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  return (sum);
}

/** communicate the global maximum */
inline int globalMaximum(int value) {
  int max = 0;
//  MPI_Barrier(MPI_COMM_WORLD); // This is time-consuming and should be debug code, only!
  MPI_Allreduce(&value, &max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  return (max);
}


#endif
