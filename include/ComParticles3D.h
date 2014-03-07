/***************************************************************************
  ComParticles.h  -  Library to manage communication of particles among processors
  -------------------
begin                : May 2008
copyright            : (C) 2008 KUL
developers           : Stefano Markidis, Giovanni Lapenta
 ***************************************************************************/

#ifndef ComParticles3D_H
#define ComParticles3D_H

#include "ComBasic3D.h"

/** comunicate particles and receive particles to and from 6 processors */
void communicateParticles(int buffer_size, double *b_Xleft, double *b_Xright, double *b_Yleft, double *b_Yright, double *b_Zleft, double *b_Zright, VirtualTopology3D * vct);

/** communicate the number of particles are not in the right domain*/
int reduceNumberParticles(int rightDomain);

/** communicate the maximum number of particles from a domain */
int reduceMaxNpExiting(int npExitingMax);

#endif
