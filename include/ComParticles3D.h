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
void communicateParticles(int buffer_size, double *b_Xleft, double *b_Xright, double *b_Yleft, double *b_Yright, double *b_Zleft, double *b_Zright, VirtualTopology3D * vct);

/** communicate the number of particles are not in the right domain*/
inline int reduceNumberParticles(int rightDomain) {
  int result = 0;
//  MPI_Barrier(MPI_COMM_WORLD); // This is time-consuming and should be debug code, only!
  MPI_Allreduce(&rightDomain, &result, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  return (result);
}

/** communicate the maximum number of particles from a domain */
inline int reduceMaxNpExiting(int npExitingMax) {
  int result = 0;
//  MPI_Barrier(MPI_COMM_WORLD); // This is time-consuming and should be debug code, only!
  MPI_Allreduce(&npExitingMax, &result, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  return (result);
}


#endif
