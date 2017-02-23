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

/*! mlmd: i need the communicator also */
/** communicate the number of particles are not in the right domain*/
//inline int reduceNumberParticles(int rightDomain) {
inline int reduceNumberParticles(int rightDomain, MPI_Comm Comm) {
  int result = 0;
  MPI_Barrier(Comm);
  MPI_Allreduce(&rightDomain, &result, 1, MPI_INT, MPI_SUM, Comm);
  return (result);
}

/** communicate the maximum number of particles from a domain */
/*! mlmd: i need the communicator also */
//inline int reduceMaxNpExiting(int npExitingMax) {
inline int reduceMaxNpExiting(int npExitingMax, MPI_Comm Comm) {
  int result = 0;
  MPI_Barrier(Comm);
  MPI_Allreduce(&npExitingMax, &result, 1, MPI_INT, MPI_MAX, Comm);
  return (result);
}


#endif
